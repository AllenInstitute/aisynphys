from aisynphys.pipeline.pipeline_module import DatabasePipelineModule
from aisynphys import config
from .opto_slice import OptoSlicePipelineModule
from collections import OrderedDict
import csv, codecs, glob, os
from neuroanalysis.data.experiment import AI_Experiment
from neuroanalysis.data.loaders.opto_experiment_loader import OptoExperimentLoader
from optoanalysis import data_model
from ... import config
from neuroanalysis.util.optional_import import optional_import
getDirHandle = optional_import('acq4.util.DataManager', 'getDirHandle')


class OptoExperimentPipelineModule(DatabasePipelineModule):
    """Imports per-experiment metadata into DB.
    """
    name = 'opto_experiment'
    dependencies = [OptoSlicePipelineModule]
    table_group = ['experiment', 'electrode', 'cell', 'pair']

    @classmethod
    def create_db_entries(cls, job, session, expt=None):
        """Generate DB entries for *job_id* and add them to *session*.
        """
        job_id = job['job_id']
        db = job['database']

        try:
            if expt is None:
                expt = load_experiment(job_id)

            # look up slice record in DB
            ts = expt.info.get('slice_info', {}).get('__timestamp__', 0.0)
            slice_entry = db.slice_from_timestamp(ts, session=session)

            fields = {
                'ext_id': expt.ext_id,
                'project_name': expt.project_name,
                'date': expt.datetime,
                'target_region': expt.target_region,
                'internal': expt.expt_info.get('internal'),
                'acsf': expt.expt_info.get('solution'),
                'target_temperature': expt.target_temperature,
                'rig_name': expt.rig_name,
                'operator_name': expt.rig_operator,
                'storage_path': None if expt.path is None else os.path.relpath(expt.path, config.synphys_data), 
                'ephys_file': None if expt.loader.get_ephys_file() is None else os.path.relpath(expt.loader.get_ephys_file(), expt.path),
                'acq_timestamp': expt.info.get('site_info',{}).get('__timestamp__')     
            }

            expt_entry = db.Experiment(**fields)
            expt_entry.slice = slice_entry
            session.add(expt_entry)

            cell_entries = {}
            for name, cell in expt.cells.items():
                if cell.electrode is not None:
                    elec = cell.electrode
                    elec_entry = db.Electrode(
                        experiment=expt_entry, 
                        ext_id=elec.electrode_id, 
                        patch_status=elec.patch_status,
                        start_time=elec.start_time,
                        stop_time=elec.stop_time,
                        device_id=elec.device_id)
                    session.add(elec_entry)

                cell_entry = db.Cell(
                    experiment=expt_entry,
                    ext_id=cell.cell_id,
                    electrode=elec_entry if cell.electrode is not None else None,
                    cre_type=cell.cre_type,
                    target_layer=cell.target_layer,
                    position=cell.position,
                    depth=cell.depth,
                    cell_class=None, ## fill in cell class fields later in Morphology module (that doesn't currently exist for opto)
                    cell_class_nonsynaptic=None,
                    meta=cell.info
                )
                session.add(cell_entry)
                cell_entries[cell] = cell_entry

            pair_entries = {}
            for name, pair in expt.pairs.items():
                pair_entry = db.Pair(
                    experiment=expt_entry,
                    pre_cell=cell_entries[pair.preCell],
                    post_cell=cell_entries[pair.postCell],
                    has_synapse=pair.isSynapse(),
                    # has_polysynapse,
                    has_electrical=None,
                    # crosstalk_artifact,
                    n_ex_test_spikes=0,  # will be counted in opto_dataset pipeline module
                    n_in_test_spikes=0,  # will be counted in opto_dataset pipeline module
                    distance=pair.distance,
                    # lateral_distance,  # should these be filled in in the opto_cortical_location pipeline module?
                    # vertical_distance, 
                    # reciprocal_id, # fill this in below
                )
                pair_entries[name] = pair_entry
                session.add(pair_entry)

            ## fill in reciprocal ids
            for name, pair in expt.pairs.items():
                pair_entries[name].reciprocal = pair_entries.get((name[1],name[0]))

        except:
            session.rollback()
            raise
        
    
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        # only need to return from experiment table; other tables will be dropped automatically.
        db = self.database
        return session.query(db.Experiment).filter(db.Experiment.ext_id.in_(job_ids)).all()

    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """

        slice_module = self.pipeline.get_module('opto_slice')
        finished_slices = slice_module.finished_jobs()
        
        # cache = synphys_cache.get_cache()
        # all_expts = cache.list_experiments()
        db = self.database
        session = db.session()
        slices = session.query(db.Slice.storage_path).all()
        slice_paths = [s[0] for s in slices]
        
        #ymls = []
        #for rec in slices:
        #    path = rec[0]
        #    ymls.extend(glob.glob(os.path.join(config.synphys_data, path, 'site_*', 'pipettes.yml')))
        expts = read_expt_csvs()
        
        n_errors = {}
        n_no_slice = []
        ready = OrderedDict()
        print("checking for ready expts....")
        for i, expt in enumerate(expts['expt_list']):
            #print("Checking experiment %i/%i"%(i, len(expts['expt_list'])))
            expt['connections_dir'] = config.connections_dir
            try:
                if expt['site_path'] == '':
                    cnx_json = os.path.join(config.connections_dir, expt['experiment'])
                    if not os.path.exists(cnx_json):
                        n_errors[expt['experiment']] = "File not found: %s" % cnx_json
                        continue
                    site_path = cnx_json
                    ex = AI_Experiment(loader=OptoExperimentLoader(load_file=cnx_json, meta_info=expt), verify=True)
                else:
                    site_path = os.path.join(config.synphys_data, expt['rig_name'].lower(), 'phys', expt['site_path'])
                    slice_path = getDirHandle(os.path.split(site_path)[0]).name(relativeTo=getDirHandle(config.synphys_data))
                    if not slice_path in slice_paths:
                        n_no_slice.append(expt['experiment'])
                        continue
                    ex = AI_Experiment(loader=OptoExperimentLoader(site_path=site_path, meta_info=expt), verify=True)

                raw_data_mtime = ex.last_modification_time
                
                slice_ts = ex.info.get('slice_info', {}).get('__timestamp__')
                if slice_ts is None:
                    slice_ts = 0.0
                slice_mtime, slice_success = finished_slices.get('%.3f'%slice_ts, (None, None))
                #print('found expt for path:', site_path)
            except Exception as exc:
                n_errors[expt['experiment']] = exc
                continue
            if slice_mtime is None or slice_success is False:
            #    slice_mtime = 0
                n_no_slice.append(expt['experiment'])
                continue

            ready[ex.ext_id] = {'dep_time':max(raw_data_mtime, slice_mtime), 'meta':{'source':site_path}}
        
        print("Found %d experiments; %d are able to be processed, %d will be skipped due to errors, %d will be skipped due to missing or failed slice entries." % (len(expts['expt_list']), len(ready), len(n_errors), len(n_no_slice)))
        if len(n_errors) > 0 or len(n_no_slice) > 0:
            print("-------- skipped experiments: ----------")
            for e, exc in n_errors.items():
                print('     %s: Error - %s' %(e.split('_conn')[0],exc))
            for e in n_no_slice:
                print('     %s: skipped due to problem with slice' % e.split('_conn')[0])
        return ready

    def dependent_job_ids(self, module, job_ids):
        """Return a list of all finished job IDs in this module that depend on 
        specific jobs from another module.
        """
        if type(module) not in self.dependencies:
            raise ValueError("%s does not depend on module %s" % (cls, module))
        
        db = self.database
        session = db.session()
        dep_ids = session.query(db.Experiment.ext_id).join(db.Slice).filter(db.Slice.acq_timestamp.in_(job_ids)).all()
        session.rollback()
        return [rec.ext_id for rec in dep_ids]


_csv_data = None
def read_expt_csvs():

    global _csv_data
    if _csv_data is not None:
        return _csv_data

    _csv_data = OrderedDict()

    expt_csv = config.experiment_csv
    distance_csv = config.distance_csv

    _csv_data['expt_list'] = []
    with open(expt_csv, 'r') as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            _csv_data['expt_list'].append(row)

    _csv_data['distances']=[]
    if distance_csv is not None:
        with codecs.open(distance_csv, "r", encoding="utf-8-sig") as csv_file:
            reader = csv.DictReader(csv_file)
            for row in reader:
                _csv_data['distances'].append(row)

    return _csv_data

def load_experiment(job_id):
    all_expts = read_expt_csvs()
    indices = [i for i, e in enumerate(all_expts['expt_list']) if job_id in e['experiment']]
    if len(indices) > 1:
        raise Exception("Cannot resolve job_id: %s. Found %s" % (job_id, [all_expts['expt_list'][i]['experiment'] for i in indices]))
    elif len(indices) == 0:
        raise Exception("Could not find csv entry for %s"%job_id)

    entry = all_expts['expt_list'][indices[0]]
    entry['distances'] = [e for e in all_expts['distances'] if e['exp_id']==job_id]
    entry['connections_dir'] = config.connections_dir ## pass this in so we can look here for a connections file, even if we pass in a site path

    if entry['site_path'] != '':
        site_path = os.path.join(config.synphys_data, entry['rig_name'].lower(), 'phys', entry['site_path'])
        if not os.path.exists(site_path):
            raise Exception('%s does not exist' % site_path)
        expt = AI_Experiment(loader=OptoExperimentLoader(site_path=site_path, meta_info=entry))
    else:
        cnx_json = os.path.join(config.connections_dir, entry['experiment'])
        expt = AI_Experiment(loader=OptoExperimentLoader(load_file=cnx_json, meta_info=entry))

    return expt

