from sqlalchemy.orm import relationship
from .database import make_table, TableGroup
from .experiment import Pair

__all__ = ['avg_first_pulse_fit_table', 'AvgFirstPulseFit']


AvgFirstPulseFit = make_table(
    name='avg_first_pulse_fit',
    comment="""Contains results of psp_fit on spike aligned, average first pulse PSP for each
            connection that passed qc in current clamp. The latency is forced
            to be within +/-.5 ms to the found fitting all of the pulses in the train (available in the 
            connection_strength.ic_fit_xoffset). The heavily weighted section (meant to 
            place more importance of the wave form during the rise time) is shifted to 
            begin at the latency. Created via fit_average_first_pulse.py. 

            All units in SI.""",
    columns=[
        ('pair_id', 'pair.id', 'The ID of the entry in the pair table to which these results apply', {'index': True}),

        # current clamp
        ('ic_amp', 'float', 'fit amplitude of current clamp average first pulses'),
        ('ic_latency', 'float', 'fit time elapsed since the time of presynaptic spike (max dv/dt) of current clamp data'),
        ('ic_rise_time', 'float', 'fit rise time of psp of current clamp data'),
        ('ic_decay_tau', 'float', 'fit decay of psp of current clamp data'),
        ('ic_avg_psp_data', 'array', 'array of the data voltage waveform used in fitting; starts 10 ms before pre-synaptic spike'),
        ('ic_avg_psp_fit', 'array', 'fit array of the best fit voltage waveform starting 10 ms before pre-synaptic spike'),
        ('ic_dt', 'float', 'time step of *avg_psp* array from current clamp data'),
        ('ic_pulse_ids', 'object', 'data base pulse ids included in the current clamp fit'),
        ('ic_nrmse', 'float', 'error of fit of current clamp fit'),
        ('ic_measured_baseline', 'float', 'average voltage measured between 10 and 1 ms before a spike'),
        ('ic_measured_amp', 'float', 'voltage amplitude within a window of 0.5 ms after spike initiation (max dv/dt) until end of array specified in the pulse_response table'),
        ('ic_weight', 'array', 'weighting used during fitting of current clamp data'),

        # voltage clamp
        ('vc_amp', 'float', 'fit amplitude of voltage clamp average first pulses'),
        ('vc_latency', 'float', 'fit time elapsed since the time of presynaptic spike (max dv/dt) of voltage clamp data'),
        ('vc_rise_time', 'float', 'fit rise time of psp measured in voltage clamp'),
        ('vc_decay_tau', 'float', 'fit decay of psp measured in voltage clamp'),
        ('vc_avg_psp_data', 'array', 'array of the data current waveform used in fitting; starts 10 ms before pre-synaptic spike'),            
        ('vc_avg_psp_fit', 'array', 'fit array of the best fit current waveform starting 10 ms before pre-synaptic spike'),
        ('vc_dt', 'float', 'time step of *avg_psp* array from voltage clamp data'),
        ('vc_pulse_ids', 'object', 'data base pulse ids included in the voltage clamp fit'),
        ('vc_nrmse', 'float', 'error of fit of voltage clamp fit'),
        ('vc_measured_baseline', 'float', 'average current measured between 10 and 1 ms before a spike'),
        ('vc_measured_amp', 'float', 'current amplitude within a window of 0.5 ms after spike initiation (max dv/dt) until end of array specified in the pulse_response table'),
        ('vc_weight', 'array', 'weighting used during fitting of voltage clamp data'),
    ]
)




Pair.avg_first_pulse_fit = relationship(AvgFirstPulseFit, back_populates="pair", cascade="delete", single_parent=True, uselist=False)
AvgFirstPulseFit.pair = relationship(Pair, back_populates="avg_first_pulse_fit", single_parent=True)

avg_first_pulse_fit_table = TableGroup([AvgFirstPulseFit])
