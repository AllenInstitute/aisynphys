# coding: utf8
"""
Analyses that measure the strength of synaptic connections.

"""
from __future__ import print_function, division

import sys, multiprocessing, time

import numpy as np
import pyqtgraph as pg

from neuroanalysis.data import Trace
from neuroanalysis import filter
from neuroanalysis.event_detection import exp_deconvolve
from neuroanalysis.baseline import float_mode

from . import database as db


def measure_peak(trace, sign, spike_time, pulse_times, spike_delay=1e-3, response_window=4e-3):
    # Start measuring response after the pulse has finished, and no earlier than 1 ms after spike onset
    # response_start = max(spike_time + spike_delay, pulse_times[1])

    # Start measuring after spike and hope that the pulse offset doesn't get in the way
    # (if we wait for the pulse to end, then we miss too many fast rise / short latency events)
    response_start = spike_time + spike_delay
    response_stop = response_start + response_window

    # measure baseline from beginning of data until 50µs before pulse onset
    baseline_start = 0
    baseline_stop = pulse_times[0] - 50e-6

    baseline = float_mode(trace.time_slice(baseline_start, baseline_stop).data)
    response = trace.time_slice(response_start, response_stop)

    if sign == '+':
        i = np.argmax(response.data)
    else:
        i = np.argmin(response.data)
    peak = response.data[i]
    latency = response.time_values[i] - spike_time
    return peak - baseline, latency


def measure_sum(trace, sign, baseline=(0e-3, 9e-3), response=(12e-3, 17e-3)):
    baseline = trace.time_slice(*baseline).data.sum()
    peak = trace.time_slice(*response).data.sum()
    return peak - baseline

        
def deconv_filter(trace, pulse_times, tau=15e-3, lowpass=24000., lpf=True, remove_artifacts=False, bsub=True):
    if tau is not None:
        dec = exp_deconvolve(trace, tau)
    else:
        dec = trace

    if remove_artifacts:
        # after deconvolution, the pulse causes two sharp artifacts; these
        # must be removed before LPF
        cleaned = remove_crosstalk_artifacts(dec, pulse_times)
    else:
        cleaned = dec

    if bsub:
        baseline = np.median(cleaned.time_slice(cleaned.t0+5e-3, cleaned.t0+10e-3).data)
        b_subbed = cleaned - baseline
    else:
        b_subbed = cleaned

    if lpf:
        return filter.bessel_filter(b_subbed, lowpass)
    else:
        return b_subbed


def remove_crosstalk_artifacts(data, pulse_times):
    dt = data.dt
    r = [-50e-6, 250e-6]
    edges = [(int((t+r[0])/dt), int((t+r[1])/dt)) for t in pulse_times]
    # If window is too shortm then it becomes seneitive to sample noise.
    # If window is too long, then it becomes sensitive to slower signals (like the AP following pulse onset)
    return filter.remove_artifacts(data, edges, window=100e-6)


def response_query(session):
    """
    Build a query to get all pulse responses along with presynaptic pulse and spike timing
    """
    q = session.query(
        db.PulseResponse.id.label('response_id'),
        db.PulseResponse.data,
        db.PulseResponse.start_time.label('rec_start'),
        db.StimPulse.onset_time.label('pulse_start'),
        db.StimPulse.duration.label('pulse_dur'),
        db.StimSpike.max_dvdt_time.label('spike_time'),
        db.PatchClampRecording.clamp_mode,
        db.PulseResponse.ex_qc_pass,
        db.PulseResponse.in_qc_pass,
    )
    q = q.join(db.StimPulse, db.PulseResponse.stim_pulse)
    q = q.join(db.StimSpike, db.StimSpike.stim_pulse_id==db.StimPulse.id)
    q = q.join(db.Recording, db.PulseResponse.recording)
    q = q.join(db.PatchClampRecording)

    # return qc-failed records as well so we can verify qc is working
    #q = q.filter(((db.PulseResponse.ex_qc_pass==True) | (db.PulseResponse.in_qc_pass==True)))

    return q


def baseline_query(session):
    """
    Build a query to get all baseline responses
    """
    q = session.query(
        db.Baseline.id.label('response_id'),
        db.Baseline.data,
        db.PatchClampRecording.clamp_mode,
        db.Baseline.ex_qc_pass,
        db.Baseline.in_qc_pass,
    )
    q = q.join(db.Recording, db.Baseline.recording)
    q = q.join(db.PatchClampRecording, db.PatchClampRecording.recording_id==db.Recording.id)

    # return qc-failed records as well so we can verify qc is working
    # q = q.filter(((db.Baseline.ex_qc_pass==True) | (db.Baseline.in_qc_pass==True)))

    return q


def analyze_response_strength(rec, source, remove_artifacts=False, deconvolve=True, lpf=True, bsub=True, lowpass=1000):
    """Perform a standardized strength analysis on a record selected by response_query or baseline_query.

    1. Determine timing of presynaptic stimulus pulse edges and spike
    2. Measure peak deflection on raw trace
    3. Apply deconvolution / artifact removal / lpf
    4. Measure peak deflection on deconvolved trace
    """
    data = Trace(rec.data, sample_rate=db.default_sample_rate)
    if source == 'pulse_response':
        # Find stimulus pulse edges for artifact removal
        start = rec.pulse_start - rec.rec_start
        pulse_times = [start, start + rec.pulse_dur]
        if rec.spike_time is None:
            # these pulses failed QC, but we analyze them anyway to make all data visible
            spike_time = 11e-3
        else:
            spike_time = rec.spike_time - rec.rec_start
    elif source == 'baseline':
        # Fake stimulus information to ensure that background data receives
        # the same filtering / windowing treatment
        pulse_times = [10e-3, 12e-3]
        spike_time = 11e-3
    else:
        raise ValueError("Invalid source %s" % source)

    results = {}

    results['raw_trace'] = data
    results['pulse_times'] = pulse_times
    results['spike_time'] = spike_time

    # Measure crosstalk from pulse onset
    p1 = data.time_slice(pulse_times[0]-200e-6, pulse_times[0]).median()
    p2 = data.time_slice(pulse_times[0], pulse_times[0]+200e-6).median()
    results['crosstalk'] = p2 - p1

    # crosstalk artifacts in VC are removed before deconvolution
    if rec.clamp_mode == 'vc' and remove_artifacts is True:
        data = remove_crosstalk_artifacts(data, pulse_times)
        remove_artifacts = False

    # Measure deflection on raw data
    results['pos_amp'], _ = measure_peak(data, '+', spike_time, pulse_times)
    results['neg_amp'], _ = measure_peak(data, '-', spike_time, pulse_times)

    # Deconvolution / artifact removal / filtering
    if deconvolve:
        tau = 15e-3 if rec.clamp_mode == 'ic' else 5e-3
    else:
        tau = None
    dec_data = deconv_filter(data, pulse_times, tau=tau, lpf=lpf, remove_artifacts=remove_artifacts, bsub=bsub, lowpass=lowpass)
    results['dec_trace'] = dec_data

    # Measure deflection on deconvolved data
    results['pos_dec_amp'], results['pos_dec_latency'] = measure_peak(dec_data, '+', spike_time, pulse_times)
    results['neg_dec_amp'], results['neg_dec_latency'] = measure_peak(dec_data, '-', spike_time, pulse_times)
    
    return results
