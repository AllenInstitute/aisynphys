"""Fit average of first pulses in for voltage and current clamp and place in 
avg_first_pulse_fit table.
"""

import numpy as np
from neuroanalysis.data import Trace, TraceList
from neuroanalysis.fitting import fit_psp
from . import database as db


time_before_spike = 10.e-3 #time in seconds before spike to start trace waveforms


def measure_amp(v_array, baseline_index_window, psp_amp_index_window):
    '''measures the max of a trace within a window of the predicted psp.
    
    Input
    -----
    v_array: array-like
        voltage wave-form
    baseline_index_window: list containing two integers
        integers specify index of the beginning and end of baseline window
    psp_amp_index_window: list containing two integers
        integers specify index of the beginning and end of window to look for max/min (psp amplitude)
        
    Returns
    -------
    relative_amp: float
        measured amplitude relative to baseline
    baseline: float
        value of baseline
    '''

    # get baseline
    baseline = v_array[baseline_index_window[0]: baseline_index_window[1]].mean()  #baseline voltage
    # get the maximum value in a region around where psp amp should be
    max_psp_region = v_array[psp_amp_index_window[0]: psp_amp_index_window[1]].max()
    # get the minimum value in a region around where psp amp should be
    min_psp_region = v_array[psp_amp_index_window[0]: psp_amp_index_window[1]].min()
    # subtract the baseline value from min and max values
    bsub=np.array([max_psp_region-baseline, min_psp_region-baseline])
    # find the absolute maximum deflection from base_line
    relative_amp=bsub[np.where(abs(bsub)==max(abs(bsub)))[0][0]] 

    return relative_amp, baseline


def extract_first_pulse_info_from_Pair_object(pair, desired_clamp='ic'):
    """Extract first pulse responses and relevant information 
    from entry in the pair database. Screen out pulses that are
    not current clamp or do not pass the corresponding
    inhibitory or excitatory qc.
    
    Input
    -----
    pair: multipatch_analysis.database.database.Pair object
    desired_clamp: string
        Specifies whether current or voltage clamp sweeps are desired.
        Options are:
            'ic': current clamp
            'vc': voltage clamp
    
    Return
    ------
    pulse_responses: TraceList 
        traces where the start of each trace is 10 ms before the spike 
    pulse_ids: list of ints
        pulse ids of *pulse_responses*
    psp_amps_measured: list of floats
        amplitude of *pulse_responses* from the *pulse_response* table
    stim_freq: list of floats
        the stimulation frequency corresponding to the *pulse_responses* 
    """

    if pair.connection_strength is None:
        # print ("\t\tSKIPPING: pair_id %s, is not yielding pair.connection_strength" % pair.id)
        return [], [], [], []
    if pair.connection_strength.synapse_type is None:
        # print ("\t\tSKIPPING: pair_id %s, is not yielding pair.connection_strength.synapse_type" % pair.id)
        return [], [], [], []
    synapse_type = pair.connection_strength.synapse_type
    pulse_responses = []
    psp_amps_measured = []
    pulse_ids = []
    stim_freqs = []
    if len(pair.pulse_responses)==0:
        # print ("\t\tSKIPPING: pair_id %s, no pulse responses in pair table" % (pair.id))
        return [], [], [], []
    for pr in pair.pulse_responses:
        stim_pulse = pr.stim_pulse
        n_spikes = stim_pulse.n_spikes
        pulse_number = stim_pulse.pulse_number
        pulse_id = pr.stim_pulse_id
        ex_qc_pass = pr.ex_qc_pass
        in_qc_pass = pr.in_qc_pass
        pcr = stim_pulse.recording.patch_clamp_recording
        stim_freq = pcr.multi_patch_probe[0].induction_frequency
        clamp_mode = pcr.clamp_mode
        # current clamp
        if clamp_mode != desired_clamp:
            continue
        # ensure that there was only 1 presynaptic spike
        if n_spikes != 1:
            continue
        # we only want the first pulse of the train
        if pulse_number != 1:
            continue

        data = pr.data
        start_time = pr.start_time
        spike_time = stim_pulse.spikes[0].max_dvdt_time        
        data_trace = Trace(data=data, t0= start_time-spike_time+time_before_spike, sample_rate=db.default_sample_rate).time_slice(start=0, stop=None) #start of the data is the spike time

        # append to output lists if neurons pass qc
        if (synapse_type == 'ex' and ex_qc_pass is True) or (synapse_type == 'in' and in_qc_pass is True):
            pulse_responses.append(data_trace)
            pulse_ids.append(pulse_id)
            stim_freqs.append(stim_freq)        
        if synapse_type == 'in' and in_qc_pass is True:
            psp_amps_measured.append(pr.pulse_response_strength.neg_amp)
        if synapse_type == 'ex' and ex_qc_pass is True:
            psp_amps_measured.append(pr.pulse_response_strength.pos_amp)

    return pulse_responses, pulse_ids, psp_amps_measured, stim_freq


def get_average_pulse_response(pair, desired_clamp='ic'):
    """
    Inputs
    ------
    pair: multipatch_analysis.database.database.Pair object

    desired_clamp: string
        Specifies whether current or voltage clamp sweeps are desired.
        Options are:
            'ic': current clamp
            'vc': voltage clamp

    Returns
    -------
    Note that all returned variables are set to None if there are no acceptable (qc pasing) sweeps
    pulse_responses: TraceList 
        traces where the start of each trace is 10 ms before the spike 
    pulse_ids: list of ints
        pulse ids of *pulse_responses*
    psp_amps_measured: list of floats
        amplitude of *pulse_responses* from the *pulse_response* table
    freq: list of floats
        the stimulation frequency corresponding to the *pulse_responses* 
    avg_psp: Trace
        average of the pulse_responses
    measured_relative_amp: float
        measured amplitude relative to baseline
    measured_baseline: float
        value of baseline
    """
    # get pulses that pass qc
    pulse_responses, pulse_ids, psp_amps_measured, freq = extract_first_pulse_info_from_Pair_object(pair, desired_clamp=desired_clamp)

    # if pulses are returned take the average
    if len(pulse_responses)>0:
        avg_psp=TraceList(pulse_responses).mean()
    else:
        return None, None, None, None, None, None, None

    # get the measured baseline and amplitude of psp
    measured_relative_amp, measured_baseline=measure_amp(avg_psp.data, 
                        [0, int((time_before_spike-1.e-3)/avg_psp.dt)], 
                        [int((time_before_spike+.5e-3)/avg_psp.dt), -1])

    return pulse_responses, pulse_ids, psp_amps_measured, freq, avg_psp, measured_relative_amp, measured_baseline 


def fit_trace(waveform, excitation, clamp_mode='ic', weight=None, latency=None, latency_jitter=None):
    """
    Input
    -----
    waveform: Trace Object
        contains data to be fit
    clamp_mode: string
        'vc' denotes voltage clamp
        'ic' denotes current clamp
    excitation: str
        'ex' or 'in' specifying excitation of synapse
    latency: float or None
        Amount of time that has passed in reference to the time of
        the pre-synaptic spike.  Note that this value has to be transformed in
        reference to the start of the data waveform being fit within
        this function.
    latency_jitter: None or float
        Amount of jitter to allow before and after the latency value. If
        *latency* is None, this value must be none.
    weight: numpy array
        Relative weighting of different sections of the voltage array
        for fitting.  If specified it must be the same length as the 
        waveform
    Note there is a 'time_before_spike' variable in this code that is
    not passed in and therefore must be global.  It specifies the
    the amount of data before the spike that is being considered in 
    the here and in the rest of the code.  Note that the value is positive,
    i.e. if we start looking at the trace 10 ms before the spike we use
    10e-3 not -10e-3.
    Returns
    -------
    self.ave_psp_fit: lmfit.model.ModelResult
        fit of the average psp waveform
    weight: numpy.ndarray
        the weight assigned to each index of the input waveform for fitting
    """
    # weighting
    if weight is None:
        weight = np.ones(len(waveform.data))  #set everything to ones

    # set fitting sign to positive or negative based on excitation and clamp state
    if (excitation == 'in') and (clamp_mode == 'ic'):
        sign = '-'
    elif (excitation == 'in') and (clamp_mode == 'vc'):
        sign = '+'
    elif (excitation == 'ex') and (clamp_mode == 'ic'):
        sign = '+'
    elif (excitation == 'ex') and (clamp_mode == 'vc'):
        sign = '-'
    elif excitation == 'any':
        sign = 'any'
    else:
        raise Exception('synaptic sign not defined')
    

    if latency is None and latency_jitter:
        raise Exception('latency_jitter cannot be specified if latency is not')
    if latency_jitter:
        if latency_jitter < .1e-3:
            raise Exception('specified latency jitter less than .0e-3 may have implications for initial conditions')
        if latency < latency_jitter:
            lower_bound = time_before_spike
        else:
            lower_bound = latency + time_before_spike - latency_jitter
        xoffset=([lower_bound, latency+time_before_spike, latency+time_before_spike+latency_jitter-.1e-3], lower_bound, latency+time_before_spike+latency_jitter)  
    elif latency:
        xoffset=(latency+time_before_spike, 'fixed')
    else:
        #since these are spike aligned the psp should not happen before the spike that happens at pre_pad by definition
        xoffset=([time_before_spike+1e-3, time_before_spike+4e-3], time_before_spike, time_before_spike+5e-3)

    fit = fit_psp(waveform, 
                    clamp_mode=clamp_mode,
                    xoffset=xoffset,
                    sign=sign, 
                    weight=weight) 

    if clamp_mode == 'ic':
        scale_factor = 1.e3
        ylabel='voltage (mV)'
    if clamp_mode == 'vc':
        scale_factor = 1.e12
        ylabel='current (pA)'        

    return fit


def fit_average_first_pulses(pair):
    # get response latency from average of all pulse responses
    message = None #initialize error message 
    xoffset = pair.connection_strength.ic_fit_xoffset
    if not xoffset:
        # too much noise:
        # return {'error': 'no ic_fit_offset from connection_strength'}
        return {'error': None}
    # excitatory or inhibitory?
    excitation = pair.connection_strength.synapse_type

    # -----------fit current clamp data---------------------        
    # get pulses
    (pulse_responses_i, pulse_ids_i, psp_amps_measured_i, freq, avg_psp_i, 
        measured_relative_amp_i, measured_baseline_i) = get_average_pulse_response(pair, desired_clamp='ic')

    if pulse_responses_i:
        # weight and fit the trace
        weight_i = np.ones(len(avg_psp_i.data)) * 10.  #set everything to ten initially
        weight_i[int((time_before_spike-3e-3) / avg_psp_i.dt):int(time_before_spike / avg_psp_i.dt)] = 0.   #area around stim artifact note that since this is spike aligned there will be some blur in where the cross talk is
        weight_i[int((time_before_spike+.0001+xoffset) / avg_psp_i.dt):int((time_before_spike+.0001+xoffset+4e-3) / avg_psp_i.dt)] = 30.  #area around steep PSP rise 
        avg_fit_i = fit_trace(avg_psp_i, excitation=excitation, weight=weight_i, latency=xoffset, latency_jitter=.5e-3)
        latency_i = avg_fit_i.best_values['xoffset'] - time_before_spike
        amp_i = avg_fit_i.best_values['amp']
        rise_time_i = avg_fit_i.best_values['rise_time']
        decay_tau_i = avg_fit_i.best_values['decay_tau']
        avg_data_waveform_i = avg_psp_i.data
        avg_fit_waveform_i = avg_fit_i.best_fit
        dt_i = avg_psp_i.dt
        nrmse_i = avg_fit_i.nrmse()
    else:
        message = 'no suitable first pulses found in current clamp'
        weight_i = np.array([0])
        latency_i = None
        amp_i = None
        rise_time_i = None
        decay_tau_i = None
        avg_data_waveform_i = np.array([0])
        avg_fit_waveform_i = np.array([0])
        dt_i = None
        nrmse_i = None
    # --------------fit voltage clamp data---------------------        
    # get pulses
    (pulse_responses_v, pulse_ids_v, psp_amps_measured_v, freq_v, avg_psp_v,  
        measured_relative_amp_v, measured_baseline_v) = get_average_pulse_response(pair, desired_clamp='vc')

    if pulse_responses_v:
        # weight and fit the trace    
        weight_v = np.ones(len(avg_psp_v.data))*10.  #set everything to ten initially
        weight_v[int((time_before_spike+.0001+xoffset)/avg_psp_v.dt):int((time_before_spike+.0001+xoffset+4e-3)/avg_psp_v.dt)] = 30.  #area around steep PSP rise 
        avg_fit_v = fit_trace(avg_psp_v, excitation=excitation, clamp_mode = 'vc', weight=weight_v, latency=xoffset, latency_jitter=.5e-3)
        latency_v = avg_fit_v.best_values['xoffset'] - time_before_spike
        amp_v = avg_fit_v.best_values['amp']
        rise_time_v = avg_fit_v.best_values['rise_time']
        decay_tau_v = avg_fit_v.best_values['decay_tau']
        avg_data_waveform_v = avg_psp_v.data
        avg_fit_waveform_v = avg_fit_v.best_fit
        dt_v = avg_psp_v.dt
        nrmse_v = avg_fit_v.nrmse()

    else:
        message = 'no suitable first pulses found in voltage clamp'
        weight_v = np.array([0])
        latency_v = None
        amp_v = None
        rise_time_v = None
        decay_tau_v = None
        avg_data_waveform_v = np.array([0])
        avg_fit_waveform_v = np.array([0])
        dt_v = None
        nrmse_v = None
    #------------ done with fitting section -----------

    # dictionary for ease of translation into the outpu table
    out_dict = {
        'ic_amp': amp_i,
        'ic_latency': latency_i,
        'ic_rise_time': rise_time_i,
        'ic_decay_tau': decay_tau_i,
        'ic_avg_psp_data': avg_data_waveform_i,
        'ic_avg_psp_fit': avg_fit_waveform_i,
        'ic_dt': dt_i,
        'ic_pulse_ids': pulse_ids_i,
        'ic_nrmse': nrmse_i,
        'ic_measured_baseline': measured_baseline_i,
        'ic_measured_amp': measured_relative_amp_i,
        'ic_weight': weight_i,

        'vc_amp': amp_v,
        'vc_latency': latency_v,
        'vc_rise_time': rise_time_v,
        'vc_decay_tau': decay_tau_v,
        'vc_avg_psp_data':avg_data_waveform_v,
        'vc_avg_psp_fit': avg_fit_waveform_v,
        'vc_dt': dt_v,
        'vc_pulse_ids': pulse_ids_v,
        'vc_nrmse': nrmse_v,
        'vc_measured_baseline': measured_baseline_v,
        'vc_measured_amp': measured_relative_amp_v,
        'vc_weight': weight_v,
        'error': message
    } 
    
    return out_dict

def fit_single_first_pulse(pr, pair):
    #TODO: HAS THE APPROPRIATE QC HAPPENED?
    message = None #initialize error message for downstream processing 
    # excitatory or inhibitory?
    excitation = pair.connection_strength.synapse_type
    if not excitation:
        raise Exception('there is no synapse_type in connection_strength')

    if excitation == 'in':
        if not pr.in_qc_pass:
            return {'error': 'this pulse does not pass inhibitory qc'}    
    if excitation == 'ex':
        if not pr.ex_qc_pass:
            return {'error': 'this pulse does not pass excitatory qc'}

    # get response latency from average first pulse table
    if not pair.avg_first_pulse_fit:
        return {'error': 'no entry in avg_first_pulse_fit table for this pair'}
        

    if pr.clamp_mode == 'vc':
        weight_i = np.array([0])
        latency_i = None
        amp_i = None
        rise_time_i = None
        decay_tau_i = None
        data_waveform_i = np.array([0])
        fit_waveform_i = np.array([0])
        dt_i = None
        nrmse_i = None
        if pair.avg_first_pulse_fit.vc_latency:
            data_trace = Trace(data=pr.data, 
                t0= pr.response_start_time - pr.spike_time + time_before_spike, 
                sample_rate=db.default_sample_rate).time_slice(start=0, stop=None)
            xoffset = pair.avg_first_pulse_fit.vc_latency
            # weight and fit the trace    
            weight_v = np.ones(len(data_trace.data))*10.  #set everything to ten initially
            weight_v[int((time_before_spike+.0001+xoffset)/data_trace.dt):int((time_before_spike+.0001+xoffset+4e-3)/data_trace.dt)] = 30.  #area around steep PSP rise 
            fit_v = fit_trace(data_trace, excitation=excitation, clamp_mode = 'vc', weight=weight_v, latency=xoffset, latency_jitter=.5e-3)
            latency_v = fit_v.best_values['xoffset'] - time_before_spike
            amp_v = fit_v.best_values['amp']
            rise_time_v = fit_v.best_values['rise_time']
            decay_tau_v = fit_v.best_values['decay_tau']
            data_waveform_v = data_trace.data
            fit_waveform_v = fit_v.best_fit
            dt_v = data_trace.dt
            nrmse_v = fit_v.nrmse()

        else:
            return {'error': 'no vc_latency available from avg_first_pulse_fit table'} #no row will be made in the table because the error message is not none               

    elif pr.clamp_mode == 'ic':
        # set voltage to none since this is current clamp
        weight_v = np.array([0])
        latency_v = None
        amp_v = None
        rise_time_v = None
        decay_tau_v = None
        data_waveform_v = np.array([0])
        fit_waveform_v = np.array([0])
        dt_v = None
        nrmse_v = None
        if pair.avg_first_pulse_fit.ic_latency:
            data_trace = Trace(data=pr.data, 
                t0= pr.response_start_time - pr.spike_time + time_before_spike, 
                sample_rate=db.default_sample_rate).time_slice(start=0, stop=None)  #TODO: annoys me that this is repetitive in vc code above.
            xoffset = pair.avg_first_pulse_fit.ic_latency
            # weight and fit the trace
            weight_i = np.ones(len(data_trace.data)) * 10.  #set everything to ten initially
            weight_i[int((time_before_spike-3e-3) / data_trace.dt):int(time_before_spike / data_trace.dt)] = 0.   #area around stim artifact note that since this is spike aligned there will be some blur in where the cross talk is
            weight_i[int((time_before_spike+.0001+xoffset) / data_trace.dt):int((time_before_spike+.0001+xoffset+4e-3) / data_trace.dt)] = 30.  #area around steep PSP rise 
            fit_i = fit_trace(data_trace, excitation=excitation, weight=weight_i, latency=xoffset, latency_jitter=.5e-3)
            latency_i = fit_i.best_values['xoffset'] - time_before_spike
            amp_i = fit_i.best_values['amp']
            rise_time_i = fit_i.best_values['rise_time']
            decay_tau_i = fit_i.best_values['decay_tau']
            data_waveform_i = data_trace.data
            fit_waveform_i = fit_i.best_fit
            dt_i = data_trace.dt
            nrmse_i = fit_i.nrmse()

        else:
            return {'error': 'no ic_latency available from avg_first_pulse_fit table'} #no row will be made in the table because the error message is not none

    else:
        raise Exception('There is no clamp mode associated with this pulse')

    #------------ done with fitting section ------------------------------

    # dictionary for ease of translation into the output table
    out_dict = {
        'ic_amp': amp_i,
        'ic_latency': latency_i,
        'ic_rise_time': rise_time_i,
        'ic_decay_tau': decay_tau_i,
        'ic_psp_data': data_waveform_i,
        'ic_psp_fit': fit_waveform_i,
        'ic_dt': dt_i,
        'ic_nrmse': nrmse_i,

        'vc_amp': amp_v,
        'vc_latency': latency_v,
        'vc_rise_time': rise_time_v,
        'vc_decay_tau': decay_tau_v,
        'vc_psp_data':data_waveform_v,
        'vc_psp_fit': fit_waveform_v,
        'vc_dt': dt_v,
        'vc_nrmse': nrmse_v,

        'error': message
    } 
    
    return out_dict