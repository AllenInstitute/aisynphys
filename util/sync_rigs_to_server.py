"""
Used to synchronize raw data from rigs to central server.


Practical note: The file synchronization implements an rsync-like
functionality in Python. In theory we would have preferred to just use 
rsync itself, but in practice this presents a few difficult issues:

* Rsync on windows requires cygwin and sshd. These work in most cases,
  but they are not as thoroughly developed as their linux counterparts.
  Bugs encountered in this process can be showstoppers.
* The actual file synchronization code here is simple and transparent,
  allowing us to customize versioning, error checking, logging,
  and file name transformations. Achieving the equivalent with rsync
  actually requires more code than we have written here, just to handle
  subprocessing, CLI flag generation, and fragile pipe communication.
"""

import os, sys, shutil, glob, traceback, pickle, time, re
from acq4.util.DataManager import getDirHandle

from multipatch_analysis import config
from multipatch_analysis.util import sync_file


def sync_experiment(site_dir):
    """Synchronize all files for an experiment to the server.

    Argument must be the path of an experiment _site_ folder. This will also cause
    synchronization for the parent (slice) and grandparent (day) folders to ensure
    that all slice images and metadata are copied. Sibling site and slice folders
    will _not_ be copied.

    Return a list of changes made.
    """
    site_dh = getDirHandle(site_dir)
    changes = []
    slice_dh = site_dh.parent()
    expt_dh = slice_dh.parent()
    
    now = time.strftime('%Y-%m-%d_%H:%M:%S')
    log("========== %s : Sync %s to server" % (now, site_dh.name()))
    skipped = 0
    
    try:
        # Decide how the top-level directory will be named on the remote server
        # (it may already be there from a previous slice/site, or the current
        # name may already be taken by another rig.)
        server_expt_path = os.path.join(config.synphys_data, get_server_path(expt_dh))
        
        log("    using server path: %s" % server_expt_path)
        skipped += _sync_paths(expt_dh.name(), server_expt_path, changes)
        
        # Copy slice files if needed
        server_slice_path = os.path.join(server_expt_path, slice_dh.shortName())
        skipped += _sync_paths(slice_dh.name(), server_slice_path, changes)

        # Copy site files if needed
        server_site_path = os.path.join(server_slice_path, site_dh.shortName())
        skipped += _sync_paths(site_dh.name(), server_site_path, changes)
        
        log("    Done; skipped %d files." % skipped)
        
    except Exception:
        err = traceback.format_exc()
        changes.append(('error', site_dh.name(), err))
        log(err)

    return changes


def log(msg):
    print(msg)
    with open(os.path.join(config.synphys_data, 'sync_log'), 'ab') as log_fh:
        log_fh.write(msg+'\n')


def _sync_paths(source, target, changes):
    """Non-recursive directory sync.

    Return the number of skipped files.
    """
    skipped = 0
    if not os.path.isdir(target):
        os.mkdir(target)
        changes.append(('mkdir', source, target))

    # Leave a note about the source of this data
    open(os.path.join(target, 'sync_source'), 'wb').write(source)

    for fname in os.listdir(source):
        src_path = os.path.join(source, fname)
        if os.path.isfile(src_path):
            dst_path = os.path.join(target, fname)

            # Skip Igor temporary files
            if fname.endswith('.pxpT0'):
                continue
            
            # Skip large files:
            #   - pxp > 20GB
            #   - others > 5GB
            src_size = os.stat(src_path).st_size
            # extension may be buried behind a backup date like "somefile.pxp_2018-11-20_02-01-20_0"
            m = re.match(r'.*(.[a-z]{3})(_2.*)?', os.path.split(src_path)[1])
            ext = '' if m is None else m.groups()[0]
            max_size = {'.pxp': 20e9, '.nwb': 7e9}.get(ext, 5e9)

            if src_size > max_size:
                log("    err! %s => %s" % (src_path, dst_path))
                changes.append(('error', src_path, 'file too large'))
                continue
            
            status = sync_file(src_path, dst_path)
            if status == 'skip':
                skipped += 1
            elif status == 'copy':
                log("    copy %s => %s" % (src_path, dst_path))
                changes.append(('copy', src_path, dst_path))
            elif status == 'update':
                log("    updt %s => %s" % (src_path, dst_path))
                changes.append(('update', src_path, dst_path))

    return skipped


def get_server_path(dh):
    """Given a directory handle to an experiment storage folder on a rig,
    return the path to the corresponding location on the server, relative
    to the server's root storage path.
    """
    # Find parent with dirType='Day'
    root = dh
    while True:
        if root.info().get('dirType') == 'Day':
            break
        parent = root.parent()
        if parent is root:
            raise Exception("Can't find the experiment root folder for %s" % dh.name())
        root = parent
    
    root_ts = root.info()['__timestamp__']
    return os.path.join('%0.3f'%root_ts, dh.name(relativeTo=root))


def find_all_sites(root):
    sites = glob.glob(os.path.join(root, '*', 'slice_*', 'site_*'))
    sites.sort(reverse=True)
    return sites


def sync_all(source='archive'):
    """Synchronize all known rig data paths to the server

    *source* should be either 'primary' or 'archive', referring to the paths
    specified in config.rig_data_paths.
    """
    log = []
    synced_paths = []
    # Loop over all rigs
    for rig_name, data_paths in config.rig_data_paths.items():
        # Each rig may have multiple paths to check
        for data_path in data_paths:
            data_path = data_path[source]

            # Get a list of all experiments stored in this path
            paths = find_all_sites(data_path)

            # synchronize files for each experiment to the server
            new_log, changed_paths = sync_experiments(paths)
            log.extend(new_log)
            synced_paths.append((rig_name, data_path, len(changed_paths), len(paths)))
    
    return log, synced_paths


def sync_experiments(paths):
    """Given a list of paths to experiment site folders, synchronize all to the server
    """
    log = []
    changed_paths = []
    for site_dir in paths:
        try:
            changes = sync_experiment(site_dir)
            if len(changes) > 0:
                log.append((site_dir, changes))
                changed_paths.append(site_dir)
        except Exception:
            exc = traceback.format_exc()
            print(exc)
            log.append((site_dir, [], exc, []))
    return log, changed_paths


if __name__ == '__main__':
    
    paths = sys.argv[1:]
    if len(paths) == 0:
        # Synchronize all known rig data paths
        log, synced_paths = sync_all(source='archive')
        print("==========================\nSynchronized files from:")
        for rig_name, data_path, n_expts_changed, n_expts_found in synced_paths:
            print("%s  :  %s  (%d/%d expts updated)" % (rig_name, data_path, n_expts_changed, n_expts_found))

    else:
        # synchronize just the specified path(s)
        log = sync_experiments(paths)
    
    errs = [change for site in log for change in site[1] if change[0] == 'error']
    print("\n----- DONE ------\n   %d errors" % len(errs))
    
    for err in errs:
        print(err[1], err[2])
