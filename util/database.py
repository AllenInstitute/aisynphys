from __future__ import division, print_function

import os, sys, user, argparse
from multipatch_analysis import database as db
from multipatch_analysis.config import synphys_db
from multipatch_analysis import config
from multipatch_analysis.pipeline import all_modules

parser = argparse.ArgumentParser()
parser.add_argument('--reset-db', action='store_true', default=False, help="Drop all tables in the database.", dest='reset_db')
parser.add_argument('--vacuum', action='store_true', default=False, help="Ask the database to clean/optimize itself.")
parser.add_argument('--bake', type=str, default=None, help="Bake current database into an sqlite file.")
parser.add_argument('--clone', type=str, default=None, help="Clone current database into a new database with the given name.")
parser.add_argument('--tables', type=str, default=None, help="Comma-separated list of tables to include while baking.")
parser.add_argument('--skip-tables', type=str, default="", help="Comma-separated list of tables to skip while baking.", dest="skip_tables")
parser.add_argument('--overwrite', action='store_true', default=False, help="Overwrite existing sqlite file.")
parser.add_argument('--update', action='store_true', default=False, help="Update existing sqlite file.")
parser.add_argument('--drop', type=str, default=None, help="Drop database with the given name.")
parser.add_argument('--dbg', action='store_true', default=False, help="Start debugging console.")


args = parser.parse_args(sys.argv[1:])

if args.dbg:
    import pyqtgraph as pg
    pg.dbg()

if args.reset_db:
    ans = raw_input('Reset database "%s"? ' % db.database.db_address_rw_clean)
    if ans == 'y':
        print("  Ok, here we go..")
        db.reset_db()
        print("    ..done.")
    else:
        print("  Oh very well. Some other time, perhaps.")

if args.vacuum:
    # cleans up DB and analyzes column statistics to improve query performance
    print("Mopping up %s.." % db.db_name)
    db.vacuum()
    print("   ..done.")


tables = None if args.tables is None else args.tables.split(',')
if args.bake is not None:
    if os.path.exists(args.bake):
        if args.overwrite:
            os.remove(args.bake)
        elif not args.update:
            print("sqlite database file %s already exists" % args.bake)
            sys.exit(0)
        
    db.bake_sqlite(args.bake, tables=tables, skip_tables=args.skip_tables.split(','))


if args.clone is not None:
    db.clone_database(args.clone, tables=tables, skip_tables=args.skip_tables.split(','))


if args.drop is not None:
    ans = raw_input("Seriously? I'm gonna dump the entire \"%s\" database? (y/n) " % args.drop)
    if ans == 'y':
        print("  Ok, don't look..")
        db.database.drop_database(args.drop)
        print("    ..done.")
    else:
        print("  Whew! Close one.")
    