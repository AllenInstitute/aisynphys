# coding: utf8
"""
For generating a DB table describing short term dynamics.

"""
from __future__ import print_function, division

import os
import numpy as np
from collections import OrderedDict
from ..util import timestamp_to_datetime
from .. import database as db
from .pipeline_module import DatabasePipelineModule
from .pulse_response import PulseResponsePipelineModule
from .connection_strength import ConnectionStrengthPipelineModule


class DynamicsPipelineModule(DatabasePipelineModule):
    """Generates dynamics analysis for each pair
    """
    name = 'dynamics'
    dependencies = [PulseResponsePipelineModule, ConnectionStrengthPipelineModule]
    table_group = db.dynamics_tables
    
    @classmethod
    def create_db_entries(cls, job_id, session):
        delays = [125, 250, 500, 1000, 2000, 4000]
        # Load experiment from DB
        expt = db.experiment_from_timestamp(job_id, session=session)
        for pair in expt.pairs.values():
            if pair.synapse is False:
                continue
            recs = session.query(db.Recording).join(db.PulseResponse).join(db.Pair).filter(db.Pair.id==pair.id).all()
            pulse_amps = {}
            for rec in recs:
                q = session.query(db.PulseResponseStrength, db.PulseResponse, db.StimPulse.pulse_number, db.MultiPatchProbe.induction_frequency, db.MultiPatchProbe.recovery_delay)
                q = q.join(db.PulseResponse, db.PulseResponseStrength.pulse_response)
                q = q.join(db.StimPulse, db.PulseResponse.stim_pulse)
                q = q.join(db.PatchClampRecording, db.PatchClampRecording.recording_id==db.PulseResponse.recording_id)
                q = q.join(db.MultiPatchProbe)
                q = q.filter(db.PulseResponse.pair_id==pair.id).filter(db.PatchClampRecording.clamp_mode=='ic').filter(db.PulseResponse.recording_id==rec.id)
                results = q.all()
                if len(results) == 0:
                    continue
                ind_freq = results[0].induction_frequency
                rec_delay = results[0].recovery_delay
                # round measured recovery delay to known delays defined in experiments, if the rounded value differs by more than 5ms ignore
                if rec_delay is None:
                    rec_delay_rounded = None
                else:
                    rec_delay_rounded = min(delays, key=lambda x:abs(x-rec_delay*1000))
                    delay_dist = abs(rec_delay_rounded - rec_delay*1000)
                    pulse_amps.setdefault(rec_delay_rounded, {})
                    if delay_dist > 5:
                        rec_delay_rounded = None
                pulse_amps.setdefault(ind_freq, {})
                sign = pair.connection_strength.synapse_type
                qc = sign+'_qc_pass'
                # make sure all 12 pulses pass qc before moving on
                qc_check = [getattr(r.pulse_response, qc) for r in results]
                if all(qc_check) is False:
                    continue
                amp_field = 'pos_dec_amp' if sign == 'ex' else 'neg_dec_amp'
                for result in results:
                    pulse_number = result.pulse_number
                    if ind_freq == 50 and rec_delay_rounded is not None:
                        pulse_amps[rec_delay_rounded].setdefault(pulse_number, [])
                        pulse_amps[rec_delay_rounded][pulse_number].append(getattr(result.pulse_response_strength, amp_field))
                    pulse_amps[ind_freq].setdefault(pulse_number, [])
                    pulse_amps[ind_freq][pulse_number].append(getattr(result.pulse_response_strength, amp_field))
            if any(pulse_amps):
                ## pulse amps is a nested dictionary with floats for keys {induction_frequency (ex 50): {pulse_number (ex 1): [list of pulse amps], 2: [pulse_amps]...}}
                pulse_ratio_8_1_50hz = np.nanmean(pulse_amps.get(50, {}).get(8, np.nan)) / np.nanmean(pulse_amps.get(50, {}).get(1, np.nan))
                pulse_ratio_9_1_125ms = np.nanmean(pulse_amps.get(125, {}).get(9, np.nan)) / np.nanmean(pulse_amps.get(125, {}).get(1, np.nan))
                pulse_ratio_9_1_250ms = np.nanmean(pulse_amps.get(250, {}).get(9, np.nan)) / np.nanmean(pulse_amps.get(250, {}).get(1, np.nan))
                pulse_ratio_9_1_500ms = np.nanmean(pulse_amps.get(500, {}).get(9, np.nan)) / np.nanmean(pulse_amps.get(500, {}).get(1, np.nan))
                pulse_ratio_9_1_1000ms = np.nanmean(pulse_amps.get(1000, {}).get(9, np.nan)) / np.nanmean(pulse_amps.get(1000, {}).get(1, np.nan))
                pulse_ratio_9_1_2000ms = np.nanmean(pulse_amps.get(2000, {}).get(9, np.nan)) / np.nanmean(pulse_amps.get(2000, {}).get(1, np.nan))
                pulse_ratio_9_1_4000ms = np.nanmean(pulse_amps.get(4000, {}).get(9, np.nan)) / np.nanmean(pulse_amps.get(4000, {}).get(1, np.nan))
                pulse_ratio_2_1_50hz = np.nanmean(pulse_amps.get(50, {}).get(2, np.nan)) / np.nanmean(pulse_amps.get(50, {}).get(1, np.nan))
                pulse_ratio_5_1_50hz = np.nanmean(pulse_amps.get(50, {}).get(5, np.nan)) / np.nanmean(pulse_amps.get(50, {}).get(1, np.nan)) 
                pulse_ratio_8_1_10hz = np.nanmean(pulse_amps.get(10, {}).get(8, np.nan)) / np.nanmean(pulse_amps.get(10, {}).get(1, np.nan))  
                pulse_ratio_8_1_20hz = np.nanmean(pulse_amps.get(20, {}).get(8, np.nan)) / np.nanmean(pulse_amps.get(20, {}).get(1, np.nan))
                pulse_ratio_8_1_100hz = np.nanmean(pulse_amps.get(100, {}).get(8, np.nan)) / np.nanmean(pulse_amps.get(100, {}).get(1, np.nan))   
                pulse_ratio_8_1_200hz = np.nanmean(pulse_amps.get(200, {}).get(8, np.nan)) / np.nanmean(pulse_amps.get(200, {}).get(1, np.nan))   

                # Write new record to DB
                dynamics = db.Dynamics(pair_id=pair.id, 
                    pulse_ratio_2_1_50hz=pulse_ratio_2_1_50hz, 
                    pulse_ratio_8_1_50hz=pulse_ratio_8_1_50hz, 
                    pulse_ratio_5_1_50hz=pulse_ratio_5_1_50hz,
                    pulse_ratio_9_1_125ms=pulse_ratio_9_1_125ms,
                    pulse_ratio_9_1_250ms=pulse_ratio_9_1_250ms,
                    pulse_ratio_9_1_500ms=pulse_ratio_9_1_500ms,
                    pulse_ratio_9_1_1000ms=pulse_ratio_9_1_1000ms,
                    pulse_ratio_9_1_2000ms=pulse_ratio_9_1_2000ms,
                    pulse_ratio_9_1_4000ms=pulse_ratio_9_1_4000ms,
                    pulse_ratio_8_1_10hz=pulse_ratio_8_1_10hz,
                    pulse_ratio_8_1_20hz=pulse_ratio_8_1_20hz,
                    pulse_ratio_8_1_100hz=pulse_ratio_8_1_100hz,
                    pulse_ratio_8_1_200hz=pulse_ratio_8_1_200hz)
                session.add(dynamics)
        
    @classmethod
    def job_records(cls, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        return session.query(db.Dynamics).filter(db.Dynamics.pair_id==db.Pair.id).filter(db.Pair.experiment_id==db.Experiment.id).filter(db.Experiment.acq_timestamp.in_(job_ids)).all()
