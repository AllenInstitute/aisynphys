import os, glob, re, time
import numpy as np
from datetime import datetime
from collections import OrderedDict
from acq4.util.DataManager import getDirHandle
from .. import config, synphys_cache
from .. import lims
from .. import qc
from ..util import timestamp_to_datetime
from ..experiment import Experiment
from .. import database as db
from ..database import dataset_tables
from .pipeline_module import DatabasePipelineModule
from .experiment import ExperimentPipelineModule
from ..connection_detection import PulseStimAnalyzer, MultiPatchSyncRecAnalyzer, BaselineDistributor
from neuroanalysis.baseline import float_mode
from neuroanalysis.data import PatchClampRecording
from ..data import MultiPatchExperiment, MultiPatchProbe


class DatasetPipelineModule(DatabasePipelineModule):
    """Imports NWB data per-experiment
    """
    name = 'dataset'
    dependencies = [ExperimentPipelineModule]
    table_group = dataset_tables

    # datasets are large and NWB access leaks memory
    # when running parallel, each child process may run only one job before being killed
    maxtasksperchild = 1  
    
    @classmethod
    def create_db_entries(cls, job_id, session):
        
        # Load experiment from DB
        expt_entry = db.experiment_from_timestamp(job_id, session=session)
        elecs_by_ad_channel = {elec.device_id:elec for elec in expt_entry.electrodes}
        pairs_by_device_id = {}
        for pair in expt_entry.pairs.values():
            pre_dev_id = pair.pre_cell.electrode.device_id
            post_dev_id = pair.post_cell.electrode.device_id
            pairs_by_device_id[(pre_dev_id, post_dev_id)] = pair
        
        # load NWB file
        path = os.path.join(config.synphys_data, expt_entry.storage_path)
        expt = Experiment(path)
        nwb = expt.data
        
        # Load all data from NWB into DB
        for srec in nwb.contents:
            temp = srec.meta.get('temperature', None)
            srec_entry = db.SyncRec(ext_id=srec.key, experiment=expt_entry, temperature=temp)
            session.add(srec_entry)
            
            srec_has_mp_probes = False
            
            rec_entries = {}
            all_pulse_entries = {}
            for rec in srec.recordings:
                
                # import all recordings
                electrode_entry = elecs_by_ad_channel[rec.device_id]  # should probably just skip if this causes KeyError?
                rec_entry = db.Recording(
                    sync_rec=srec_entry,
                    electrode=electrode_entry,
                    start_time=rec.start_time,
                )
                session.add(rec_entry)
                rec_entries[rec.device_id] = rec_entry
                
                # import patch clamp recording information
                if not isinstance(rec, PatchClampRecording):
                    continue
                qc_pass = qc.recording_qc_pass(rec)
                pcrec_entry = db.PatchClampRecording(
                    recording=rec_entry,
                    clamp_mode=rec.clamp_mode,
                    patch_mode=rec.patch_mode,
                    stim_name=rec.stimulus.description,
                    baseline_potential=rec.baseline_potential,
                    baseline_current=rec.baseline_current,
                    baseline_rms_noise=rec.baseline_rms_noise,
                    qc_pass=qc_pass,
                )
                session.add(pcrec_entry)

                # import test pulse information
                tp = rec.nearest_test_pulse
                if tp is not None:
                    indices = tp.indices or [None, None]
                    tp_entry = db.TestPulse(
                        electrode=electrode_entry,
                        recording=rec_entry,
                        start_index=indices[0],
                        stop_index=indices[1],
                        baseline_current=tp.baseline_current,
                        baseline_potential=tp.baseline_potential,
                        access_resistance=tp.access_resistance,
                        input_resistance=tp.input_resistance,
                        capacitance=tp.capacitance,
                        time_constant=tp.time_constant,
                    )
                    session.add(tp_entry)
                    pcrec_entry.nearest_test_pulse = tp_entry
                    
                # import information about STP protocol
                if not isinstance(rec, MultiPatchProbe):
                    continue
                srec_has_mp_probes = True
                psa = PulseStimAnalyzer.get(rec)
                ind_freq, rec_delay = psa.stim_params()
                mprec_entry = db.MultiPatchProbe(
                    patch_clamp_recording=pcrec_entry,
                    induction_frequency=ind_freq,
                    recovery_delay=rec_delay,
                )
                session.add(mprec_entry)
            
                # import presynaptic stim pulses
                pulses = psa.pulses()
                
                pulse_entries = {}
                all_pulse_entries[rec.device_id] = pulse_entries
                
                rec_tvals = rec['primary'].time_values

                for i,pulse in enumerate(pulses):
                    # Record information about all pulses, including test pulse.
                    t0 = rec_tvals[pulse[0]]
                    t1 = rec_tvals[pulse[1]]
                    data_start = max(0, t0 - 10e-3)
                    data_stop = t0 + 10e-3
                    pulse_entry = db.StimPulse(
                        recording=rec_entry,
                        pulse_number=i,
                        onset_time=t0,
                        amplitude=pulse[2],
                        duration=t1-t0,
                        data=rec['primary'].time_slice(data_start, data_stop).resample(sample_rate=20000).data,
                        data_start_time=data_start,
                    )
                    session.add(pulse_entry)
                    pulse_entries[i] = pulse_entry
                    

                # import presynaptic evoked spikes
                # For now, we only detect up to 1 spike per pulse, but eventually
                # this may be adapted for more.
                spikes = psa.evoked_spikes()
                for i,sp in enumerate(spikes):
                    pulse = pulse_entries[sp['pulse_n']]
                    if sp['spike'] is not None:
                        spinfo = sp['spike']
                        extra = {
                            'peak_time': rec_tvals[spinfo['peak_index']],
                            'max_dvdt_time': rec_tvals[spinfo['rise_index']],
                            'max_dvdt': spinfo['max_dvdt'],
                        }
                        if 'peak_diff' in spinfo:
                            extra['peak_diff'] = spinfo['peak_diff']
                        if 'peak_value' in spinfo:
                            extra['peak_value'] = spinfo['peak_value']
                        
                        pulse.n_spikes = 1
                    else:
                        extra = {}
                        pulse.n_spikes = 0
                    
                    spike_entry = db.StimSpike(
                        stim_pulse=pulse,
                        **extra
                    )
                    session.add(spike_entry)
                    pulse.first_spike = spike_entry
            
            if not srec_has_mp_probes:
                continue
            
            # import postsynaptic responses
            mpa = MultiPatchSyncRecAnalyzer(srec)
            for pre_dev in srec.devices:
                for post_dev in srec.devices:
                    if pre_dev == post_dev:
                        continue

                    # get all responses, regardless of the presence of a spike
                    responses = mpa.get_spike_responses(srec[pre_dev], srec[post_dev], align_to='pulse', require_spike=False)
                    post_tvals = srec[post_dev]['primary'].time_values
                    for resp in responses:
                        # base_entry = db.Baseline(
                        #     recording=rec_entries[post_dev],
                        #     start_index=resp['baseline_start'],
                        #     stop_index=resp['baseline_stop'],
                        #     data=resp['baseline'].resample(sample_rate=20000).data,
                        #     mode=float_mode(resp['baseline'].data),
                        # )
                        # session.add(base_entry)
                        pair_entry = pairs_by_device_id.get((pre_dev, post_dev), None)
                        if pair_entry is None:
                            continue  # no data for one or both channels
                        if resp['ex_qc_pass']:
                            pair_entry.n_ex_test_spikes += 1
                        if resp['in_qc_pass']:
                            pair_entry.n_in_test_spikes += 1
                        resp_entry = db.PulseResponse(
                            recording=rec_entries[post_dev],
                            stim_pulse=all_pulse_entries[pre_dev][resp['pulse_n']],
                            pair=pair_entry,
                            start_time=post_tvals[resp['rec_start']],
                            data=resp['response'].resample(sample_rate=20000).data,
                            ex_qc_pass=resp['ex_qc_pass'],
                            in_qc_pass=resp['in_qc_pass'],
                        )
                        session.add(resp_entry)
                        
            # generate up to 20 baseline snippets for each recording
            for dev in srec.devices:
                rec = srec[dev]
                rec_tvals = rec['primary'].time_values
                dist = BaselineDistributor.get(rec)
                for i in range(20):
                    base = dist.get_baseline_chunk(20e-3)
                    if base is None:
                        # all out!
                        break
                    start, stop = base
                    data = rec['primary'][start:stop].resample(sample_rate=20000).data

                    ex_qc_pass, in_qc_pass = qc.pulse_response_qc_pass(rec, [start, stop], None, [])

                    base_entry = db.Baseline(
                        recording=rec_entries[dev],
                        start_time=rec_tvals[start],
                        data=data,
                        mode=float_mode(data),
                        ex_qc_pass=ex_qc_pass,
                        in_qc_pass=in_qc_pass,
                    )
                    session.add(base_entry)
        
    @classmethod
    def job_records(cls, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        # only need to return from syncrec table; other tables will be dropped automatically.
        return session.query(db.SyncRec).filter(db.SyncRec.experiment_id==db.Experiment.id).filter(db.Experiment.acq_timestamp.in_(job_ids)).all()

    @classmethod
    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        # All experiments and their creation times in the DB
        expts = ExperimentPipelineModule.finished_jobs()
        
        # Look up nwb file locations for all experiments
        session = db.Session()
        expt_recs = session.query(db.Experiment.acq_timestamp, db.Experiment.storage_path, db.Experiment.ephys_file).filter(db.Experiment.ephys_file != None).all()
        expt_paths = {rec.acq_timestamp: rec for rec in expt_recs}
        session.rollback()
        
        # Return the greater of NWB mod time and experiment DB record mtime
        ready = OrderedDict()
        for expt_id, (expt_mtime, expt_success) in expts.items():
            if expt_id not in expt_paths or expt_success is False:
                # no NWB file; ignore
                continue
            rec = expt_paths[expt_id]
            ephys_file = os.path.join(config.synphys_data, rec.storage_path, rec.ephys_file)
            nwb_mtime = timestamp_to_datetime(os.stat(ephys_file).st_mtime)
            ready[rec.acq_timestamp] = max(expt_mtime, nwb_mtime)
        return ready
