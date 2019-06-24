import numpy as np

from neuroanalysis.miesnwb import MiesNwb, MiesSyncRecording, MiesRecording
from neuroanalysis.stimuli import find_square_pulses
from neuroanalysis.spike_detection import detect_evoked_spike


class MultiPatchExperiment(MiesNwb):
    """Extension of neuroanalysis data abstraction layer to include
    multipatch-specific metadata.
    """
    def create_sync_recording(self, sweep_id):
        return MultiPatchSyncRecording(self, sweep_id)

        
class MultiPatchSyncRecording(MiesSyncRecording):
    def __init__(self, nwb, sweep_id):
        MiesSyncRecording.__init__(self, nwb, sweep_id)
        self._baseline_mask = None
        try:
            self.meta['temperature'] = self.recordings[0].meta['notebook']['Async AD 1: Bath Temperature']
        except Exception:
            self.meta['temperature'] = None
    
    def create_recording(self, sweep_id, ch):
        miesrec = MiesRecording(self, sweep_id, ch)
        stim = miesrec.meta['notebook']['Stim Wave Name'].lower()
        if 'pulsetrain' in stim or 'recovery' in stim:
            return MultiPatchProbe(miesrec)
        else:
            return miesrec

    def baseline_regions(self, settling_time=100e-3):
        """Return a list of start,stop pairs indicating regions during the recording that are expected to be quiescent
        due to absence of pulses.
        """
        if self._baseline_mask is None:
            mask = np.zeros(len(self.recordings[0]['primary']), dtype=bool)
            dt = self.recordings[0]['primary'].dt
            settle_size = int(settling_time / dt)
            for rec in self.recordings:
                pa = PulseStimAnalyzer.get(rec)
                for pulse in pa.pulses():
                    mask[pulse[0]:pulse[1] + settle_size] = True
            self._baseline_mask = mask

            starts = list(np.argwhere(~mask[1:] & mask[:-1])[:,0])
            stops = list(np.argwhere(mask[1:] & ~mask[:-1])[:,0])
            if starts[0] > stops[0]:
                starts.insert(0, 0)
            if stops[-1] < starts[-1]:
                stops.append(len(mask))
            self._baseline_regions = [r for r in zip(starts, stops) if r[1] > r[0]]

        return self._baseline_regions


class MultiPatchProbe(MiesRecording):
    def __init__(self, recording):
        self._parent_rec = recording
        self._base_regions = None
        
    def __len__(self):
        return len(self._parent_rec)

    def __getattr__(self, attr):
        if '_parent_rec' not in self.__dict__:
            raise AttributeError(attr)
        return getattr(self._parent_rec, attr)

    @property
    def baseline_regions(self):
        # detect baseline regions from pulses rather than labnotebook
        # (can't always count on there being delay periods in the recording)
        if self._base_regions is None:
            self._base_regions = self._parent_rec.parent.baseline_regions()
        return self._base_regions


class Analyzer(object):
    @classmethod
    def get(cls, obj):
        """Get the analyzer attached to a recording, or create a new one.
        """
        analyzer = getattr(obj, '_' + cls.__name__, None)
        if analyzer is None:
            analyzer = cls(obj)
        return analyzer

    def _attach(self, obj):
        attr = '_' + self.__class__.__name__
        if hasattr(obj, attr):
            raise TypeError("Object %s already has attached %s" % (obj, self.__class__.__name__))
        setattr(obj, attr, self)


class PulseStimAnalyzer(Analyzer):
    """Used for analyzing a patch clamp recording with square-pulse stimuli.
    """
    def __init__(self, rec):
        self._attach(rec)
        self.rec = rec
        self._pulses = None
        self._evoked_spikes = None
        
    def pulses(self):
        """Return a list of (start, stop, amp) tuples describing square pulses
        in the stimulus.
        """
        if self._pulses is None:
            trace = self.rec['command']
            pulses = find_square_pulses(trace)
            self._pulses = []
            for p in pulses:
                start = trace.index_at(p.global_start_time)
                stop = trace.index_at(p.global_start_time + p.duration)
                self._pulses.append((start, stop, p.amplitude))
        return self._pulses

    def evoked_spikes(self):
        """Given presynaptic Recording, detect action potentials
        evoked by current injection or unclamped spikes evoked by a voltage pulse.
        """
        if self._evoked_spikes is None:
            pre_trace = self.rec['primary']

            # Detect pulse times
            pulses = self.pulses()

            # detect spike times
            spike_info = []
            for i,pulse in enumerate(pulses):
                on, off, amp = pulse
                if amp < 0:
                    # assume negative pulses do not evoke spikes
                    # (todo: should be watching for rebound spikes as well)
                    continue
                spike = detect_evoked_spike(self.rec, [on, off])
                spike_info.append({'pulse_n': i, 'pulse_ind': on, 'pulse_len': off-on, 'spike': spike})
            self._evoked_spikes = spike_info
        return self._evoked_spikes

    def stim_params(self):
        """Return induction frequency and recovery delay.
        """
        pulses = [p[0] for p in self.pulses() if p[2] > 0]
        if len(pulses) < 2:
            return None, None
        dt = self.rec['command'].dt
        ind_freq = np.round(1.0 / (dt * (pulses[1] - pulses[0])))
        rec_delay = np.round(dt*np.diff(pulses).max(), 3)
        
        return ind_freq, rec_delay


class PulseResponse(object):
    """Represents a chunk of postsynaptic recording taken during a presynaptic pulse stimulus.

    This class provides access to:
    - presynaptic stimulus waveform and metadata
    - presynaptic recording (if available)
    - postsynaptic recording
    - cell pair (if available)
    - postsynaptic cell
    """
    def __init__(self, recording=None, stim_pulse=None, pair=None, start_time=None, ex_qc_pass=None, in_qc_pass=None, post_tseries=None):
        self.recording = recording
        self.stim_pulse = stim_pulse
        self.pair = pair
        self.start_time = start_time
        self.ex_qc_pass = ex_qc_pass
        self.in_qc_pass = in_qc_pass
        self.post_tseries = post_tseries

    @property
    def pre_tseries(self):
        return self.stim_pulse.recorded_tseries

    @property
    def stim_tseries(self):
        return self.stim_pulse.stimulus_tseries
    

class StimPulse(object):
    """Represents a single stimiulus pulse intended to evoke a synaptic response.

    This class provides access to:
    - presynaptic stimulus waveform and metadata
    - presynaptic recording (if available)
    - presynaptic cell
    """
    def __init__(self, recording=None, pulse_number=None, onset_time=None, next_pulse_time=None, amplitude=None, duration=None, 
                 n_spikes=None, recorded_tseries=None):
        self.recording = recording
        self.pulse_number = pulse_number
        self.onset_time = onset_time
        self.next_pulse_time = next_pulse_time
        self.amplitude = amplitude
        self.duration = duration
        self.n_spikes = n_spikes
        self.recorded_tseries = recorded_tseries

    @property
    def recorded_tseries(self):
        if self._rec_tseries is None:
            self._rec_tseries = Trace(self.data, sample_rate=default_sample_rate, t0=self.data_start_time)
        return self._rec_tseries

    @property
    def stimulus_tseries(self):
        # can we work Stimulus objects into here, rather than generating manually??
        if self._stim_tseries is None:
            rec_ts = self.recorded_tseries  # todo: avoid loading full data for this step
            data = np.zeros(shape=rec_ts.shape)
            pstart = rec_ts.index_at(self.onset_time)
            pstop = rec_ts.index_at(self.onset_time + self.duration)
            data[pstart:pstop] = self.amplitude            
            self._stim_tseries = rec_ts.copy(data=data)
        return self._stim_tseries
