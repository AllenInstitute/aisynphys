"""
Low-level relational database / sqlalchemy interaction.

The actual schemas for database tables are implemented in other files in this subpackage.
"""
from __future__ import division, print_function

import os, sys, io, time, json, threading, gc, re
from datetime import datetime
from collections import OrderedDict
import numpy as np
try:
    import queue
except ImportError:
    import Queue as queue

import sqlalchemy
from distutils.version import LooseVersion
if LooseVersion(sqlalchemy.__version__) < '1.2':
    raise Exception('requires at least sqlalchemy 1.2')

from sqlalchemy import create_engine, Column, Integer, String, Boolean, Float, Date, DateTime, LargeBinary, ForeignKey, or_, and_
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship, deferred, sessionmaker, aliased, reconstructor
from sqlalchemy.types import TypeDecorator
from sqlalchemy.sql.expression import func

from .. import config

# database version should be incremented whenever the schema has changed
db_version = 12
db_name = '{database}_{version}'.format(database=config.synphys_db, version=db_version)
default_app_name = ('mp_a:' + ' '.join(sys.argv))[:60]

def db_address(host, db_name=None, app_name=None):
    if app_name is None:
        app_name = default_app_name
    if host.startswith('postgres'):
        return "{host}/{db_name}?application_name={app_name}".format(host=host, db_name=db_name, app_name=app_name)
    else:
        return host


db_address_ro = db_address(config.synphys_db_host, db_name)
if config.synphys_db_host_rw is None:
    db_address_rw = None
else:
    db_address_rw = db_address(config.synphys_db_host_rw, db_name)

# make a passwordless version of the address for printing
db_address_rw_clean = db_address_rw
if db_address_rw is not None:
    m = re.match(r"(\w+\:/+)(([^\:]+)(\:\S+)?(\@))?([^\?]+)", db_address_rw)
    if m is not None:
        g = m.groups()
        db_address_rw_clean = g[0] + (g[2] or '') + (g[4] or '') + g[5]


default_sample_rate = 20000


_sample_rate_str = '%dkHz' % (default_sample_rate // 1000)


class TableGroup(object):
    """Class used to manage a group of tables that act as a single unit--tables in a group
    are always created and deleted together.
    """    
    def __init__(self, tables):
        self.tables = OrderedDict([(t.__table__.name,t) for t in tables])

    def __getitem__(self, item):
        return self.tables[item]

    def drop_tables(self):
        global engine_rw
        drops = []
        for k in self.tables:
            if k in engine_rw.table_names():
                drops.append(k)
        if len(drops) == 0:
            return
            
        if engine_rw.name == 'sqlite':
            engine_rw.execute('drop table %s' % (','.join(drops)))
        else:
            engine_rw.execute('drop table %s cascade' % (','.join(drops)))

    def create_tables(self):
        global engine_rw, engine_ro
        
        # if we do not have write access, just verify that the tables exist
        if engine_rw is None:
            for k in self.tables:
                if k not in engine_ro.table_names():
                    raise Exception("Table %s not found in database %s" % (k, db_address_ro))
            return

        create_tables(tables=[self[k].__table__ for k in self.tables])
        
    def enable_triggers(self, enable):
        """Enable or disable triggers for all tables in this group.
        
        This can be used to temporarily disable constraint checking on tables that are under development,
        allowing the rest of the pipeline to continue operating (for example, if removing an object from 
        the pipeline would violate a foreign key constraint, disabling triggers will allow this constraint
        to go unchecked).
        """
        s = Session(readonly=False)
        enable = 'enable' if enable else 'disable'
        for table in self.tables.keys():
            s.execute("alter table %s %s trigger all;" % (table, enable))
        s.commit()


#----------- define ORM classes -------------

class NDArray(TypeDecorator):
    """For marshalling arrays in/out of binary DB fields.
    """
    impl = LargeBinary
    
    def process_bind_param(self, value, dialect):
        if value is None:
            return b'' 
        buf = io.BytesIO()
        np.save(buf, value, allow_pickle=False)
        return buf.getvalue()
        
    def process_result_value(self, value, dialect):
        if value == b'':
            return None
        buf = io.BytesIO(value)
        return np.load(buf, allow_pickle=False)


class JSONObject(TypeDecorator):
    """For marshalling objects in/out of json-encoded text.
    """
    impl = String
    
    def process_bind_param(self, value, dialect):
        return json.dumps(value)
        
    def process_result_value(self, value, dialect):
        return json.loads(value)


class FloatType(TypeDecorator):
    """For marshalling float types (including numpy).
    """
    impl = Float
    
    def process_bind_param(self, value, dialect):
        if value is None:
            return None
        return float(value)
        
    #def process_result_value(self, value, dialect):
        #buf = io.BytesIO(value)
        #return np.load(buf, allow_pickle=False)


column_data_types = {
    'int': Integer,
    'float': FloatType,
    'bool': Boolean,
    'str': String,
    'date': Date,
    'datetime': DateTime,
    'array': NDArray,
#    'object': JSONB,  # provides support for postges jsonb, but conflicts with sqlite
    'object': JSONObject,
}


ORMBase = declarative_base()

all_tables = OrderedDict()
def make_table(name, columns, base=None, **table_args):
    """Generate an ORM mapping class from a simplified schema format.

    Columns named 'id' (int) and 'meta' (object) are added automatically.

    Parameters
    ----------
    name : str
        Name of the table, used to set __tablename__ in the new class
    base : class or None
        Base class on which to build the new table class
    table_args : keyword arguments
        Extra keyword arguments are used to set __table_args__ in the new class
    columns : list of tuple
        List of column specifications. Each column is given as a tuple:
        ``(col_name, data_type, comment, {options})``. Where *col_name* and *comment* 
        are strings, *data_type* is a key in the column_data_types global, and
        *options* is a dict providing extra initialization arguments to the sqlalchemy
        Column (for example: 'index', 'unique'). Optionally, *data_type* may be a 'tablename.id'
        string indicating that this column is a foreign key referencing another table.
    """
    global all_tables
    props = {
        '__tablename__': name,
        '__table_args__': table_args,
        'id': Column(Integer, primary_key=True),
    }

    for column in columns:
        colname, coltype = column[:2]
        kwds = {} if len(column) < 4 else column[3]
        kwds['comment'] = None if len(column) < 3 else column[2]
        defer_col = kwds.pop('deferred', False)
        ondelete = kwds.pop('ondelete', None)

        if coltype not in column_data_types:
            if not coltype.endswith('.id'):
                raise ValueError("Unrecognized column type %s" % coltype)
            props[colname] = Column(Integer, ForeignKey(coltype, ondelete=ondelete), **kwds)
        else:
            ctyp = column_data_types[coltype]
            props[colname] = Column(ctyp, **kwds)

        if defer_col:
            props[colname] = deferred(props[colname])

    # props['time_created'] = Column(DateTime, default=func.now())
    # props['time_modified'] = Column(DateTime, onupdate=func.current_timestamp())
    props['meta'] = Column(column_data_types['object'])

    if base is None:
        new_table = type(name, (ORMBase,), props)
    else:
        # need to jump through a hoop to allow __init__ on table classes;
        # see: https://docs.sqlalchemy.org/en/latest/orm/constructors.html
        if hasattr(base, '_init_on_load'):
            @reconstructor
            def _init_on_load(self, *args, **kwds):
                base._init_on_load(self)
            props['_init_on_load'] = _init_on_load
        new_table = type(name, (base,ORMBase), props)

    all_tables[name] = new_table
    return new_table


#-------------- initial DB access ----------------
engine_ro = None
engine_rw = None
engine_pid = None  # pid of process that created this engine. 
def init_engines():
    """Initialize ro and rw (if possible) database engines.
    
    If any engines are currently defined, they will be disposed first.
    """
    global engine_ro, engine_rw, engine_pid
    dispose_engines()
    
    if db_address_ro.startswith('postgres'):
        opts_ro = {'pool_size': 10, 'max_overflow': 40, 'isolation_level': 'AUTOCOMMIT'}
    else:
        opts_ro = {}        
    engine_ro = create_engine(db_address_ro, **opts_ro)

    if db_address_rw is not None:
        if db_address_rw.startswith('postgres'):
            opts_rw = {'pool_size': 10, 'max_overflow': 40}
        else:
            opts_rw = {}        
        engine_rw = create_engine(db_address_rw, **opts_rw)

    engine_pid = os.getpid()


def dispose_engines():
    global engine_ro, engine_rw, engine_pid, _sessionmaker_ro, _sessionmaker_rw
    if engine_ro is not None:
        engine_ro.dispose()
        engine_ro = None
        _sessionmaker_ro = None
    if engine_rw is not None:
        engine_rw.dispose()
        engine_rw = None
        _sessionmaker_rw = None
        
    # collect now or else we might try to collect engine-related garbage in forked processes,
    # which can lead to "OperationalError: server closed the connection unexpectedly"
    # Note: if this turns out to be flaky as well, we can just disable connection pooling.
    gc.collect()
    
    engine_pid = None    


def get_engines():
    """Return ro and rw (if possible) database engines.
    """
    global engine_ro, engine_rw, engine_pid
    if os.getpid() != engine_pid:
        # In forked processes, we need to re-initialize the engine before
        # creating a new session, otherwise child processes will
        # inherit and muck with the same connections. See:
        # https://docs.sqlalchemy.org/en/latest/faq/connections.html#how-do-i-use-engines-connections-sessions-with-python-multiprocessing-or-os-fork
        if engine_pid is not None:
            print("Making new session for subprocess %d != %d" % (os.getpid(), engine_pid))
        dispose_engines()
    
    if engine_ro is None:
        init_engines()
    
    return engine_ro, engine_rw
    

init_engines()


_sessionmaker_ro = None
_sessionmaker_rw = None
# external users should create sessions from here.
def Session(readonly=True):
    """Create and return a new database Session instance.
    
    If readonly is True, then the session is created using read-only credentials and has autocommit enabled.
    This prevents idle-in-transaction timeouts that occur when GUI analysis tools would otherwise leave transactions
    open after each request.
    """
    global _sessionmaker_ro, _sessionmaker_rw, engine_ro, engine_rw, engine_pid
    engine_ro, engine_rw = get_engines()
            
    if readonly:
        if _sessionmaker_ro is None:
            _sessionmaker_ro = sessionmaker(bind=engine_ro)
        return _sessionmaker_ro()
    else:
        if engine_rw is None:
            raise RuntimeError("Cannot start read-write DB session; no write access engine is defined (see config.synphys_db_host_rw)")
        if _sessionmaker_rw is None:
            _sessionmaker_rw = sessionmaker(bind=engine_rw)
        return _sessionmaker_rw()


def reset_db():
    """Drop the existing synphys database and initialize a new one.
    """
    global engine_rw, db_name
    
    dispose_engines()
    
    drop_database(db_name)
    create_database(db_name)

    # reconnect to DB
    init_engines()

    create_tables()
    grant_readonly_permission(db_name)


def list_databases(host=None):
    if host is None:
        host = config.synphys_db_host
    pg_engine = create_engine(host + '/postgres')
    with pg_engine.begin() as conn:
        conn.connection.set_isolation_level(0)
        return [rec[0] for rec in conn.execute('SELECT datname FROM pg_catalog.pg_database;')]


def drop_database(db_name, host=None):
    if host is None:
        host = config.synphys_db_host_rw
    
    if host.startswith('sqlite:///'):
        dbfile = host[10:]
        if os.path.isfile(dbfile):
            os.remove(dbfile)
    elif host.startswith('postgres'):
        pg_engine = create_engine(host + '/postgres')
        with pg_engine.begin() as conn:
            conn.connection.set_isolation_level(0)
            try:
                conn.execute('drop database %s' % db_name)
            except sqlalchemy.exc.ProgrammingError as err:
                if 'does not exist' not in err.message:
                    raise
    else:
        raise TypeError("Unsupported database backend %s" % host.split(':')[0])


def create_database(db_name, host=None):
    if host is None:
        host = config.synphys_db_host_rw
    if host.startswith('sqlite'):
        return
    elif host.startswith('postgres'):
        ro_user = config.synphys_db_readonly_user

        # connect to postgres db just so we can create the new DB
        pg_engine = create_engine(db_address(host, 'postgres'))
        with pg_engine.begin() as conn:
            conn.connection.set_isolation_level(0)
            conn.execute('create database %s' % db_name)
            # conn.execute('ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT SELECT ON TABLES TO %s;' % ro_user)

    else:
        raise TypeError("Unsupported database backend %s" % host.split(':')[0])


def grant_readonly_permission(db_name, host=None):
    if host is None:
        host = config.synphys_db_host_rw
    if host.startswith('sqlite'):
        return
    elif host.startswith('postgres'):
        ro_user = config.synphys_db_readonly_user

        # grant readonly permissions
        pg_engine = create_engine(db_address(host, db_name))
        with pg_engine.begin() as conn:
            conn.connection.set_isolation_level(0)
            for cmd in [
                ('GRANT CONNECT ON DATABASE %s TO %s' % (db_name, ro_user)),
                ('GRANT USAGE ON SCHEMA public TO %s' % ro_user),
                ('GRANT SELECT ON ALL TABLES IN SCHEMA public to %s' % ro_user)]:
                conn.execute(cmd)


def create_tables(tables=None, engine=None):
    """Create tables in the database.
    
    A list of *tables* may be optionally specified (see sqlalchemy.schema.MetaData.create_all) to 
    create a subset of known tables.
    """
    # Create all tables
    global ORMBase
    if engine is None:
        _, engine = get_engines()
    ORMBase.metadata.create_all(bind=engine, tables=tables)
    grant_readonly_permission(db_name)


def vacuum(tables=None):
    """Cleans up database and analyzes table statistics in order to improve query planning.
    Should be run after any significant changes to the database.
    """
    engine_ro, engine_rw = get_engines()
    with engine_rw.begin() as conn:
        conn.connection.set_isolation_level(0)
        if tables is None:
            conn.execute('vacuum analyze')
        else:
            for table in tables:
                conn.execute('vacuum analyze %s' % table)


_default_session = None
def default_session(fn):
    """Decorator used to auto-fill `session` keyword arguments
    with a global default Session instance.
    
    If the global session is used, then it will be rolled back after 
    the decorated function returns (to prevent idle-in-transaction timeouts).
    """
    
    def wrap_with_session(*args, **kwds):
        used_default_session = False
        if kwds.get('session', None) is None:
            kwds['session'] = get_default_session()
            used_default_session = True
        try:
            ret = fn(*args, **kwds)
            return ret
        finally:
            if used_default_session:
                get_default_session().rollback()
    return wrap_with_session    


def get_default_session():
    global _default_session
    if _default_session is None:
        _default_session = Session(readonly=True)
    return _default_session


def bake_sqlite(sqlite_file, **kwds):
    """Dump a copy of the database to an sqlite file.
    """
    sqlite_addr = "sqlite:///%s" % sqlite_file
    sqlite_engine = create_engine(sqlite_addr)
    create_tables(engine=sqlite_engine)
    read_engine = get_engines()[0]
    
    last_size = 0
    for table in copy_tables(read_engine, sqlite_engine, **kwds):
        size = os.stat(sqlite_file).st_size
        diff = size - last_size
        last_size = size
        print("   sqlite file size:  %0.2fGB  (+%0.2fGB for %s)" % (size*1e-9, diff*1e-9, table))


def clone_database(db_name, host=None, overwrite=False, **kwds):
    """Copy the current database to a new database of the given name.
    """
    if db_name in list_databases(host):
        if overwrite:
            drop_database(db_name, host)
        else:
            raise Exception("Destination database %s already exists." % db_name)

    if host is None:
        host = config.synphys_db_host_rw
    
    source_engine = get_engines()[0]
    dest_address = db_address(host, db_name)
    dest_engine = create_engine(dest_address)
    create_database(db_name, host)
    create_tables(engine=dest_engine)
    
    for table in copy_tables(source_engine, dest_engine, **kwds):
        pass


def copy_tables(source_engine, dest_engine, tables=None, skip_tables=(), skip_errors=False):
    """Iterator that copies all modules in a database from one engine to another.
    
    Yields each module and table name as it is completed.
    
    This function does not create tables in dest_engine; use create_tables if needed.
    """
    global all_tables
    read_session = sessionmaker(bind=source_engine)()
    write_session = sessionmaker(bind=dest_engine)()
    
    for table_name, table in all_tables.items():
        if (table_name in skip_tables) or (tables is not None and table_name not in tables):
            print("Skipping %s.." % table_name)
            continue
        print("Cloning %s.." % table_name)
        
        # read from table in background thread, write to table in main thread.
        reader = TableReadThread(table)
        for i,rec in enumerate(reader):
            try:
                write_session.execute(table.__table__.insert(rec))
            except Exception:
                if skip_errors:
                    print("Skip record %d:" % i)
                    sys.excepthook(*sys.exc_info())
                else:
                    raise
            if i%200 == 0:
                print("%d/%d   %0.2f%%\r" % (i, reader.max_id, (100.0*(i+1.0)/reader.max_id)), end="")
                sys.stdout.flush()
            
        print("   committing %d rows..                    " % i)
        write_session.commit()
        read_session.rollback()
        
        yield table_name

    print("Optimizing database..")    
    write_session.execute("analyze")
    write_session.commit()
    print("All finished!")


class TableReadThread(threading.Thread):
    """Iterator that yields records (all columns) from a table.
    
    Records are queried chunkwise and queued in a background thread to enable more efficient streaming.
    """
    def __init__(self, table, chunksize=1000):
        threading.Thread.__init__(self)
        self.daemon = True
        
        self.table = table
        self.chunksize = chunksize
        self.queue = queue.Queue(maxsize=5)
        self.max_id = Session().query(func.max(table.id)).all()[0][0]        
        self.start()
        
    def run(self):
        session = Session()
        table = self.table
        chunksize = self.chunksize
        all_columns = table.__table__.c
        for i in range(0, self.max_id, chunksize):
            query = session.query(*all_columns).filter((table.id >= i) & (table.id < i+chunksize))
            records = query.all()
            self.queue.put(records)
        self.queue.put(None)
        session.rollback()
        session.close()
    
    def __iter__(self):
        while True:
            recs = self.queue.get()
            if recs is None:
                break
            for rec in recs:
                yield rec
    
