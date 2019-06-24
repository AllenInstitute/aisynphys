import os, glob, re, pickle, time
from datetime import datetime
from collections import OrderedDict
from acq4.util.DataManager import getDirHandle
from .. import database as db
from ..database import slice_tables
from .pipeline_module import DatabasePipelineModule
from .. import config, lims, constants
from ..util import datetime_to_timestamp, timestamp_to_datetime


class SlicePipelineModule(DatabasePipelineModule):
    """Imports per-slice metadata into DB.
    """
    name = 'slice'
    dependencies = []
    table_group = slice_tables
    
    @classmethod
    def create_db_entries(cls, job_id, session):
        slices = all_slices()
        path = slices[job_id]
        dh = getDirHandle(path)
        info = dh.info()
        parent_info = dh.parent().info()
        
        # pull some metadata from LIMS
        sid = info['specimen_ID'].strip()
        limsdata = lims.specimen_info(sid)

        quality = info.get('slice quality', None)
        try:
            quality = int(quality)
        except Exception:
            quality = None

        # Interpret slice time
        slice_time = parent_info.get('time_of_dissection', None)
        if slice_time == '':
            slice_time = None
        if slice_time is not None:
            m = re.match(r'((20\d\d)-(\d{1,2})-(\d{1,2}) )?(\d+):(\d+)', slice_time.strip())
            if m is not None:
                _, year, mon, day, hh, mm = m.groups()
                if year is None:
                    date = datetime.fromtimestamp(dh.parent().info()['__timestamp__'])
                    slice_time = datetime(date.year, date.month, date.day, int(hh), int(mm))
                else:
                    slice_time = datetime(int(year), int(mon), int(day), int(hh), int(mm))

        # construct full genotype string 
        genotype = limsdata['genotype'] or ''
        for info in (parent_info, info):
            inj = info.get('injections')
            if inj in (None, ''):
                continue
            if inj not in constants.INJECTIONS:
                raise KeyError("Injection %r is unknown in constants.INJECTIONS" % inj)
            genotype = ';'.join(genotype.split(';') + [constants.INJECTIONS[inj]])

        fields = {
            'acq_timestamp': info['__timestamp__'],
            'species': limsdata['organism'],
            'date_of_birth': limsdata['date_of_birth'],
            'age': limsdata['age'],
            'sex': limsdata['sex'],
            'genotype': genotype,
            'orientation': limsdata['plane_of_section'],
            'surface': limsdata['exposed_surface'],
            'hemisphere': limsdata['hemisphere'],
            'quality': quality,
            'slice_time': slice_time,
            'slice_conditions': {},
            'lims_specimen_name': sid,
            'storage_path': dh.name(relativeTo=dh.parent().parent()),
        }

        sl = db.Slice(**fields)
        session.add(sl)

    @classmethod
    def job_records(cls, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        """
        return session.query(db.Slice).filter(db.Slice.acq_timestamp.in_(job_ids)).all()

    @classmethod
    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        slices = all_slices()
        ready = OrderedDict()
        for ts, path in slices.items():
            mtime = os.stat(os.path.join(path, '.index')).st_mtime
            # test file updates:
            # import random
            # if random.random() > 0.8:
            #     mtime *= 2
            ready[ts] = timestamp_to_datetime(mtime)
        return ready


_all_slices = None
def all_slices():
    """Return a dict mapping {slice_timestamp: path} for all known slices.
    
    This is only generated once per running process; set _all_slices = None
    to force the list to be regenerated.
    """
    global _all_slices
    if _all_slices is not None:
        return _all_slices
        
    # Speed things up by caching this list with a 4 hour timeout
    cachefile = os.path.join(config.cache_path, 'all_slices.pkl')
    if os.path.exists(cachefile):
        age = time.time() - os.stat(cachefile).st_mtime
        if age < 4 * 3600:
            print("Loaded slice timestamps from cache (%0.1f hours old)" % (age/3600.))
            return pickle.load(open(cachefile, 'r'))
    
    slice_dirs = sorted(glob.glob(os.path.join(config.synphys_data, '*', 'slice_*')))
    _all_slices = OrderedDict()
    for path in slice_dirs:
        dh = getDirHandle(path)
        ts = dh.info().get('__timestamp__')
        if ts is None:
            print("MISSING TIMESTAMP: %s" % path)
            continue
        _all_slices[ts] = path
        
    try:
        tmpfile = cachefile+'.tmp'
        pickle.dump(_all_slices, open(tmpfile, 'w'))
        os.rename(tmpfile, cachefile)
    except:
        if os.path.exists(tmpfile):
            os.remove(tmpfile)
    
    return _all_slices
