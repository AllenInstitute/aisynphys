import multipatch_analysis.database as db
import multipatch_analysis.connection_strength as cs 
from multipatch_analysis.database.database import TableGroup
import pandas
import multipatch_analysis.fit_average_first_pulse as afpf
import matplotlib.pyplot as plt


def join_pulse_response_to_expt(query):
    pre_rec = db.aliased(db.Recording)
    post_rec = db.aliased(db.Recording)
    joins = [
        (post_rec, db.PulseResponse.recording),
        (db.PatchClampRecording,),
        (db.MultiPatchProbe,),
        (db.SyncRec,),
        (db.Experiment,),
        (db.StimPulse, db.PulseResponse.stim_pulse),
        (pre_rec, db.StimPulse.recording),
    ]
    for join_args in joins:
        query = query.join(*join_args)

    return query, pre_rec, post_rec

def get_single_pulse_data(session, pair):
    """Select records from pulse_response_strength table
    This is a version grabbed from cs.get_amps altered to get
    additional info.  Note that if cs.get_baseline_amps is also
    called there could be some unexpected results/correspondense
    problems if cs.get_amps or cs.get_baseline is altered.
    """
    cols = [
        db.PulseResponse.id,
        db.PulseResponse.ex_qc_pass,
        db.PulseResponse.in_qc_pass,
        db.PatchClampRecording.clamp_mode,
        db.PatchClampRecording.baseline_potential,
        db.PatchClampRecording.baseline_current,
        db.StimPulse.pulse_number,
        db.StimSpike.max_dvdt_time.label('spike_time'), #note that in my fit_average_first_pulse code I filter for spikes[0] so there may be an error with multiple here
        db.PulseResponse.start_time.label('response_start_time'),
        db.MultiPatchProbe.induction_frequency.label('stim_freq'),
        db.PulseResponse.data,
    ]


    q = session.query(*cols)
    #q = q.join(db.PulseResponse)
    
    q, pre_rec, post_rec = join_pulse_response_to_expt(q)
    q = q.join(db.StimSpike)
    # note that these are aliased so needed to be joined outside the join_pulse_response_to_expt() function
    q = q.add_columns(post_rec.start_time.label('rec_start_time'))
    q = q.add_columns(post_rec.id.label('post_rec_id'))  #this will identify the train the pulse is in
        
    filters = [
        (pre_rec.electrode==pair.pre_cell.electrode,),
        (post_rec.electrode==pair.post_cell.electrode,),
        (db.PatchClampRecording.qc_pass==True,),
        (db.StimPulse.pulse_number==1,),
    ]
    for filter_args in filters:
        q = q.filter(*filter_args)
    
    # should result in chronological order
    q = q.order_by(db.PulseResponse.id).all()

    return(q)
    # df = pandas.read_sql_query(q.statement, q.session.bind)
    # recs = df.to_records()
    # return recs


#uid=1527020350.517
uid=1483748888.290
pre_cell_id=3
post_cell_id=1

# query a pair
session = db.Session()
expt = db.experiment_from_timestamp(uid, session=session)
pair = expt.pairs[(pre_cell_id, post_cell_id)]

# get pulse ids from connection strength
pulse_responses = get_single_pulse_data(session, pair)
for pr in pulse_responses:
    if pr.id != 317966:
        continue 
        
    print('pulse response id', pr.id)
    out = afpf.fit_single_first_pulse(pr, pair)
    if pr.clamp_mode == 'vc':
        plt.plot(out['vc_psp_data'], 'b')
        plt.plot(out['vc_psp_fit'], 'r')
        plt.show()
    

# # get pulse ids from pulse_responses table
# pulse_stim_response_ids=[]
# for pr in pair.pulse_responses:
#     pulse_stim_response_ids.append(pr.stim_pulse_id)
# print(pulse_stim_response_ids)
