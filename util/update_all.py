"""
Runs all pipeline analysis stages on a daily schedule. 

Why not use cron like a normal human being?
Because I like having the script run in the console where I cam monitor and
debug problems more easily.
"""

import os, sys, time, argparse
from datetime import datetime, timedelta
from collections import OrderedDict


def delay(hour=2):
    """Sleep until *hour*"""
    now = datetime.now()
    tomorrow = now + timedelta(days=1)
    next_run = datetime(tomorrow.year, tomorrow.month, tomorrow.day, 3, 0)
    delay = (next_run - now).total_seconds()

    print("Sleeping %d seconds until %s.." % (delay, next_run))
    time.sleep(delay)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run all analysis pipeline stages to import / analyze new data on a schedule.")
    parser.add_argument('--now', default=False, action='store_true', help="Run once immediately before starting scheduled updates.")
    parser.add_argument('--skip', default='', help="comma-separated list of stages to skip")
    args = parser.parse_args(sys.argv[1:])

    if not args.now:
        delay()

    stages = OrderedDict([
        ('sync',                    ('python util/sync_rigs_to_server.py', 'sync raw data to server')),
        ('pipeline',                ('python util/analysis_pipeline.py all', 'run analysis pipeline')),
        ('vacuum',                  ('python util/database.py --vacuum', 'vacuum database')),
        ('bake',                    ('python util/database.py --bake=synphys_current.sqlite --overwrite', 'bake sqlite')),
    ])

    skip = [] if args.skip == '' else args.skip.split(',')
    for name in skip:
        if name not in stages:
            print("Unknown stage %r. Options are: %r" % (name, list(stages.keys())))
            sys.exit(-1)

    while True:
        logfile = 'update_logs/' + time.strftime('%Y-%m-%d_%H-%M-%S') + '.log'
        for name, cmd in stages.items():
            cmd, msg = cmd
            msg = ("======================================================================================\n" + 
                   "    " + msg + "\n" + 
                   "======================================================================================\n")
            print(msg)
            open(logfile, 'a').write(msg)
            
            if name in skip:
                msg = "   [ skipping ]\n"
                print(msg)
                open(logfile, 'a').write(msg)
                
                skip.remove(name)  # only skip once
                continue

            os.system(cmd + " 2>&1 | tee -a " + logfile)
        delay()


