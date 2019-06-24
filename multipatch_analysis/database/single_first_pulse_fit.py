from sqlalchemy.orm import relationship
from .database import make_table, TableGroup
from .dataset import PulseResponse


__all__ = ['single_first_pulse_fit_table', 'SingleFirstPulseFit']


SingleFirstPulseFit = make_table(
    name='single_first_pulse_fit',
    comment="""Contains results of psp_fit on individual first pulses. The latency is forced
            to the value found via fitting the average of the first pulses (available in the
            avg_first_pulse_fit table). The heavily weighted section (meant to 
            place more importance of the wave form during the rise time) is shifted to 
            begin at the latency. Created via fit_single_first_pulse.py. 
            All units in SI.""",
    columns=[
        ('pulse_response_id', 'pulse_response.id', 'The ID of the entry in the pulse_response table', {'index': True}),

        # current clamp
        ('ic_amp', 'float', 'fit amplitude of current clamp first pulses'),
        ('ic_latency', 'float', 'fit time elapsed since the time of presynaptic spike (max dv/dt) of current clamp data'),
        ('ic_rise_time', 'float', 'fit rise time of psp of current clamp data'),
        ('ic_decay_tau', 'float', 'fit decay of psp of current clamp data'),
        ('ic_psp_data', 'array', 'array of the data voltage waveform used in fitting; starts 10 ms before pre-synaptic spike'),
        ('ic_psp_fit', 'array', 'fit array of the best fit voltage waveform starting 10 ms before pre-synaptic spike'),
        ('ic_dt', 'float', 'time step of *avg_psp* array from current clamp data'),
        ('ic_nrmse', 'float', 'error of fit of current clamp fit'),

        # voltage clamp
        ('vc_amp', 'float', 'fit amplitude of voltage clamp first pulses'),
        ('vc_latency', 'float', 'fit time elapsed since the time of presynaptic spike (max dv/dt) of voltage clamp data'),
        ('vc_rise_time', 'float', 'fit rise time of psp measured in voltage clamp'),
        ('vc_decay_tau', 'float', 'fit decay of psp measured in voltage clamp'),
        ('vc_psp_data', 'array', 'array of the data current waveform used in fitting; starts 10 ms before pre-synaptic spike'),            
        ('vc_psp_fit', 'array', 'fit array of the best fit current waveform starting 10 ms before pre-synaptic spike'),
        ('vc_dt', 'float', 'time step of *avg_psp* array from voltage clamp data'),
        ('vc_nrmse', 'float', 'error of fit of voltage clamp fit'),
    ]
)

PulseResponse.single_first_pulse_fit = relationship(SingleFirstPulseFit, back_populates="pulse_response", cascade="delete", single_parent=True, uselist=False)
SingleFirstPulseFit.pulse_response = relationship(PulseResponse, back_populates="single_first_pulse_fit", single_parent=True)

single_first_pulse_fit_table = TableGroup([SingleFirstPulseFit])
