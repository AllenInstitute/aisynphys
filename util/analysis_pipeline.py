from __future__ import print_function
import argparse, sys, os, logging
import pyqtgraph as pg 
from multipatch_analysis.pipeline import all_modules
import multipatch_analysis.database as db
from multipatch_analysis import config


if __name__ == '__main__':
    logging.basicConfig()
    
    all_modules = all_modules()
    
    parser = argparse.ArgumentParser(description="Process analysis pipeline jobs")
    parser.add_argument('modules', type=str, nargs='+', help="The name of the analysis module(s) to run: %s" % ', '.join(list(all_modules.keys())))
    parser.add_argument('--rebuild', action='store_true', default=False, help="Remove and rebuild tables for this analysis")
    parser.add_argument('--retry', action='store_true', default=False, help="Retry processing jobs that previously failed")
    parser.add_argument('--workers', type=int, default=None, help="Set the number of concurrent processes during update")
    parser.add_argument('--local', action='store_true', default=False, help="Disable concurrent processing to make debugging easier")
    parser.add_argument('--raise-exc', action='store_true', default=False, help="Disable catching exceptions encountered during processing", dest='raise_exc')
    parser.add_argument('--limit', type=int, default=None, help="Limit the number of experiments to process")
    parser.add_argument('--uids', type=lambda s: [float(x) for x in s.split(',')], default=None, help="Select specific IDs to analyze (or drop)", )
    parser.add_argument('--drop', action='store_true', default=False, help="Drop selected analysis results (do not run updates)", )
    parser.add_argument('--vacuum', action='store_true', default=False, help="Run VACUUM ANALYZE on the database to optimize its query planner", )
    parser.add_argument('--bake', action='store_true', default=False, help="Bake an sqlite file after the pipeline update completes", )
    
    
    args = parser.parse_args(sys.argv[1:])

    if args.local:
        pg.dbg()
    
    if 'all' in args.modules:
        modules = list(all_modules.values())
    else:
        modules = []
        for mod in args.modules:
            try:
                if mod.startswith(':'):
                    i = all_modules.keys().index(mod[1:])
                    modules.extend(all_modules.values()[:i])
                elif mod.endswith(':'):
                    i = all_modules.keys().index(mod[:-1])
                    modules.extend(all_modules.values()[i:])
                else:
                    modules.append(all_modules[mod])
            except (KeyError, ValueError):
                print('Unknown analysis module "%s"; options are: %s' % (mod, list(all_modules.keys())))
                sys.exit(-1)

    # sort topologically
    modules = [m for m in list(all_modules.values()) if m in modules]
    
    if args.rebuild:
        mod_names = ', '.join([module.name for module in modules])
        if raw_input("Rebuild modules: %s? (y/n) " % mod_names) != 'y':
            print("  Nuts.")
            sys.exit(-1)

    if args.bake and os.path.exists(config.synphys_db_sqlite):
        msg = "sqlite database file %s already exists; ok to overwrite? (y/n) " % config.synphys_db_sqlite
        ans = raw_input(msg)
        if ans == 'y':
            print("  Ok, you asked for it..")
            os.remove(config.synphys_db_sqlite)
        else:
            print("  Phooey.")
            args.bake = False

    if args.rebuild:
        for module in modules:
            print("Dropping module %s" % module.name)
            module.drop_all(reinitialize=True)
        print("  done.")

    if args.drop:
        for module in modules:
            if args.uids is None:
                print("Dropping and reinitializing module %s" % module.name)
                module.drop_all(reinitialize=True)
            else:
                print("Dropping %d jobs in module %s" % (len(args.uids), module.name))
                module.drop_jobs(job_ids=args.uids)
    else:
        report = []
        for module in modules:
            print("=============================================")
            result = module.update(job_ids=args.uids, retry_errors=args.retry, limit=args.limit, parallel=not args.local, workers=args.workers, raise_exceptions=args.raise_exc)
            report.append((module, result))
            
        if args.vacuum:
            print("Starting vacuum..")
            db.vacuum()
            print("   ..finished vacuum.")

        print("\n================== Error Report ===========================")
        for module, result in report:
            print("------ %s : %d errors -------" % (module.name, result['n_errors']))
            for job, err in result['errors'].items():
                print("    %0.3f : %s" % (job, err))
        
            
        print("\n================== Update Report ===========================")
        for module, result in report:
            print("{name:20s}  dropped: {n_dropped:6d}  updated: {n_updated:6d} ({n_retry:6d} retry)  errors: {n_errors:6d}".format(name=module.name, **result))

    if args.bake:
        print("\n================== Bake Sqlite ===========================")
        db.bake_sqlite(config.synphys_db_sqlite)
