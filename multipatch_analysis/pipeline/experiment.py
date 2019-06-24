from __future__ import division, print_function
import os, sys, glob, re, time
import numpy as np
from datetime import datetime
from collections import OrderedDict
from acq4.util.DataManager import getDirHandle
from .. import database as db
from ..database import experiment_tables
from .pipeline_module import DatabasePipelineModule
from .. import config, synphys_cache
from .. import lims
from ..util import datetime_to_timestamp
from ..experiment import Experiment
from .slice import SlicePipelineModule


class ExperimentPipelineModule(DatabasePipelineModule):
    """Imports per-experiment metadata into DB.
    """
    name = 'experiment'
    dependencies = [SlicePipelineModule]
    table_group = experiment_tables
    
    
    @classmethod
    def create_db_entries(cls, job_id, session):
        cache = synphys_cache.get_cache()
        all_expts = cache.list_experiments()
        site_path = all_expts[job_id]
        expt = Experiment(site_path=site_path)
        
        # look up slice record in DB
        ts = expt.slice_timestamp
        slice_entry = db.slice_from_timestamp(ts, session=session)
        
        expt_info = expt.expt_info
        expt_lims_id = lims.expt_cluster_ids(slice_entry.lims_specimen_name, expt.timestamp)

        if len(expt_lims_id) == 1:
            expt_lims_id = expt_lims_id[0]
        elif len(expt_lims_id) == 0:
            expt_lims_id = None
        else:
            raise Exception ('Too many LIMS specimens %d' % len(expt_lims_id))
        fields = {
            'original_path': expt.original_path,
            'storage_path': expt.server_path,
            'ephys_file': None if expt.nwb_file is None else os.path.relpath(expt.nwb_file, expt.path),
            'rig_name': expt.rig_name,
            'project_name': expt.project_name,
            'acq_timestamp': expt.timestamp,
            'target_region': expt.target_region,
            'internal': expt_info.get('internal'),
            'acsf': expt_info.get('solution'),
            'target_temperature': expt.target_temperature,
            'lims_specimen_id': expt_lims_id,
        }

        # Create entry in experiment table
        expt_entry = db.Experiment(**fields)
        expt_entry.slice = slice_entry
        session.add(expt_entry)

        # create pipette and cell entries
        cell_entries = {}
            
        for e_id, elec in expt.electrodes.items():
            elec_entry = db.Electrode(experiment=expt_entry, ext_id=elec.electrode_id, device_id=elec.device_id)
            for k in ['patch_status', 'start_time', 'stop_time',  
                      'initial_resistance', 'initial_current', 'pipette_offset',
                      'final_resistance', 'final_current']:
                if hasattr(elec, k):
                    setattr(elec_entry, k, getattr(elec, k))
            session.add(elec_entry)

            if elec.cell is not None:
                cell = elec.cell
                if expt_lims_id is not None:
                    cell_specimens = lims.child_specimens(expt_lims_id)
                    if len(cell_specimens) != 0:
                        lims_cells = lims.cluster_cells(expt_lims_id)
                        cell_meta = {}
                        cell_id_map = {int(lims_cell.external_specimen_name): lims_cell.id for lims_cell in lims_cells if lims_cell is not None}
                        cell_lims_id = cell_id_map.get(cell.cell_id)
                        cell_meta['lims_specimen_id'] = cell_lims_id
                    else:
                        cell_meta['lims_specimen_id'] = None
                cell_entry = db.Cell(
                    experiment=expt_entry,
                    electrode=elec_entry,
                    ext_id=cell.cell_id,
                    cre_type=cell.cre_type,
                    target_layer=cell.target_layer,
                    is_excitatory=cell.is_excitatory,
                    depth=cell.depth,
                    position=cell.position,
                    meta=cell_meta
                )
                session.add(cell_entry)
                cell_entries[cell] = cell_entry

        # create pairs
        for i, pre_cell in expt.cells.items():
            for j, post_cell in expt.cells.items():
                if i == j:
                    continue
                # check to see if i,j is in manual connection calls
                # (do not use expt.connections, which excludes some connections based on QC)
                syn_calls = expt.connection_calls
                synapse = None if syn_calls is None else ((i, j) in syn_calls)
                gap_calls = expt.gap_calls
                electrical = None if gap_calls is None else ((i, j) in gap_calls)

                pre_cell_entry = cell_entries[pre_cell]
                post_cell_entry = cell_entries[post_cell]
                p1, p2 = pre_cell.position, post_cell.position
                if None in [p1, p2]:
                    distance = None
                else:
                    distance = np.linalg.norm(np.array(p1) - np.array(p2))
                
                pair_entry = db.Pair(
                    experiment=expt_entry,
                    pre_cell=pre_cell_entry,
                    post_cell=post_cell_entry,
                    synapse=synapse,
                    electrical=electrical,
                    n_ex_test_spikes=0,  # will be counted later
                    n_in_test_spikes=0,
                    distance=distance,
                )
                session.add(pair_entry)

                pre_id = pre_cell_entry.electrode.device_id
                post_id = post_cell_entry.electrode.device_id
        
    @classmethod
    def job_records(cls, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        # only need to return from experiment table; other tables will be dropped automatically.
        return session.query(db.Experiment).filter(db.Experiment.acq_timestamp.in_(job_ids)).all()

    @classmethod
    def dependent_job_ids(cls, module, job_ids):
        """Return a list of all finished job IDs in this module that depend on 
        specific jobs from another module.
        """
        if module not in cls.dependencies:
            raise ValueError("%s does not depend on module %s" % (cls, module))
        
        session = db.Session()
        dep_ids = session.query(db.Experiment.acq_timestamp).join(db.Slice).filter(db.Slice.acq_timestamp.in_(job_ids)).all()
        session.rollback()
        return [rec.acq_timestamp for rec in dep_ids]

    @classmethod
    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        finished_slices = SlicePipelineModule.finished_jobs()
        
        # cache = synphys_cache.get_cache()
        # all_expts = cache.list_experiments()
        
        session = db.Session()
        slices = session.query(db.Slice.storage_path).all()
        
        ymls = []
        for rec in slices:
            path = rec[0]
            ymls.extend(glob.glob(os.path.join(config.synphys_data, path, 'site_*', 'pipettes.yml')))
        
        n_errors = 0
        ready = OrderedDict()
        for i,yml_path in enumerate(ymls):
            print("  checking experiment %d/%d          \r" % (i, len(ymls)), end='')
            sys.stdout.flush()
            site_path = os.path.dirname(yml_path)
            try:
                expt = Experiment(site_path=site_path, verify=False)
                raw_data_mtime = expt.last_modification_time
                slice_ts = expt.slice_timestamp
                slice_mtime, slice_success = finished_slices.get(slice_ts, None)
            except Exception:
                n_errors += 1
                continue
            if slice_mtime is None or slice_success is False:
                continue
            ready[expt.timestamp] = max(raw_data_mtime, slice_mtime)
        
        print("Found %d experiments; %d are able to be processed, %d were skipped due to errors." % (len(ymls), len(ready), n_errors))
        return ready
