# coding: utf8
from __future__ import print_function, division

import os
import pyqtgraph as pg
from ... import config
from ..pipeline_module import DatabasePipelineModule
from .dataset import DatasetPipelineModule
from ...pulse_response_strength import baseline_query, response_query, measure_response, analyze_response_strength


class PulseResponsePipelineModule(DatabasePipelineModule):
    """Analyze postsynaptic responses for all presynaptic evoked spikes
    """
    name = 'pulse_response'
    dependencies = [DatasetPipelineModule]
    table_group = ['pulse_response_strength', 'baseline_response_strength']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        expt_id = job['job_id']
        _compute_strength('pulse_response', expt_id, session, db)
        _compute_strength('baseline', expt_id, session, db)
        
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        q = session.query(db.PulseResponseStrength)
        q = q.filter(db.PulseResponseStrength.pulse_response_id==db.PulseResponse.id)
        q = q.filter(db.PulseResponse.pair_id==db.Pair.id)
        q = q.filter(db.Pair.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.ext_id.in_(job_ids))
        prs = q.all()
        
        q = session.query(db.BaselineResponseStrength)
        q = q.filter(db.BaselineResponseStrength.baseline_id==db.Baseline.id)
        q = q.filter(db.Baseline.recording_id==db.Recording.id)
        q = q.filter(db.Recording.sync_rec_id==db.SyncRec.id)
        q = q.filter(db.SyncRec.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.ext_id.in_(job_ids))
        brs = q.all()
        
        return prs+brs


def _compute_strength(source, expt_id, session, db):
    """Compute per-pulse-response strength metrics
    """
    if source == 'baseline':
        q = baseline_query(session)
    elif source == 'pulse_response':
        q = response_query(session)
    else:
        raise ValueError("Invalid source %s" % source)

    # select just data for the selected experiment
    q = q.join(db.SyncRec).join(db.Experiment).filter(db.Experiment.ext_id==expt_id)

    prof = pg.debug.Profiler(delayed=False)
    
    recs = q.all()
    prof('fetch')
        
    new_recs = []

    for rec in recs:
        new_rec = {'%s_id'%source: rec.response_id}
        
        # best estimate of response amplitude using known latency for this synapse
        result = measure_response(rec, source)
        for k in ['amplitude', 'dec_fit_amp', 'dec_fit_latency', 'dec_fit_sigma']:
            new_rec[k] = result[k]
        
        # "unbiased" response analysis used to predict connectivity
        result = analyze_response_strength(rec, source)
        # copy a subset of results over to new record
        for k in ['pos_amp', 'neg_amp', 'pos_dec_amp', 'neg_dec_amp', 'pos_dec_latency', 'neg_dec_latency', 'crosstalk']:
            new_rec[k] = result[k]
        new_recs.append(new_rec)
    
    prof('process')

    # Bulk insert is not safe with parallel processes
    # if source == 'pulse_response':
    #     session.bulk_insert_mappings(PulseResponseStrength, new_recs)
    # else:
    #     session.bulk_insert_mappings(BaselineResponseStrength, new_recs)
    
    rec_type = db.PulseResponseStrength if source == 'pulse_response' else db.BaselineResponseStrength
    for rec in new_recs:
        session.add(rec_type(**rec))

    # just to collect error messages here in case we have made a mistake:
    session.flush()
    
    prof('insert')
