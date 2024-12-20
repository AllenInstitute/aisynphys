from sqlalchemy.orm import relationship
from . import make_table
from .experiment import Pair


__all__ = ['Dynamics']


Dynamics = make_table(
    name='dynamics',
    comment="Describes short term dynamics of synaptic connections.",
    columns=[
        ('pair_id', 'pair.id', 'The ID of the cell pair described by each record', {'index': True, 'unique': True}),
        ('qc_pass', 'bool', 'Indicates whether dynamics records pass quality control', {'index': True}),
        ('n_source_events', 'int', 'Number of qc-passed pulse response amplitudes from which dynamics metrics were generated', {'index': True}),
        ('paired_pulse_ratio_50hz', 'float', 'The median ratio of 2nd / 1st pulse amplitudes for 50Hz pulse trains.', {'index': True}),
        ('stp_initial_50hz', 'float', 'The median relative change from 1st to 2nd pulse for 50Hz pulse trains', {'index': True}),
        ('stp_initial_50hz_n', 'float', 'Number of samples represented in stp_initial_50Hz', {'index': True}),
        ('stp_initial_50hz_std', 'float', 'Standard deviation of samples represented in stp_initial_50Hz', {'index': True}),
        ('stp_induction_50hz', 'float', 'The median relative change from 1st to 6th-8th pulses for 50Hz pulse trains', {'index': True}),
        ('stp_induction_50hz_n', 'float', 'Number of samples represented in stp_induction_50Hz', {'index': True}),
        ('stp_induction_50hz_std', 'float', 'Standard deviation of samples represented in stp_induction_50Hz', {'index': True}),
        ('stp_recovery_250ms', 'float', 'The median relative change from 1st-4th to 9th-12th pulses for pulse trains with a 250 ms recovery period', {'index': True}),
        ('stp_recovery_250ms_n', 'float', 'Number of samples represented in stp_recovery_250ms', {'index': True}),
        ('stp_recovery_250ms_std', 'float', 'Standard deviation of samples represented in stp_recovery_250ms', {'index': True}),
        ('stp_recovery_single_250ms', 'float', 'The median relative change from 1st to 9th pulses for pulse trains with a 250 ms recovery period', {'index': True}),
        ('stp_recovery_single_250ms_n', 'float', 'Number of samples represented in stp_recovery_single_250ms', {'index': True}),
        ('stp_recovery_single_250ms_std', 'float', 'Standard deviation of samples represented in stp_recovery_single_250ms', {'index': True}),
        ('pulse_amp_90th_percentile', 'float', 'The 90th-percentile largest pulse amplitude, used to normalize change values in this table', {}),
        ('noise_amp_90th_percentile', 'float', 'The 90th-percentile largest amplitude measured from background noise, used for comparison to pulse_amp_90th_percentile', {}),
        
        ('pulse_amp_first_50hz', 'float', 'Median amplitude of first pulse on 50 hz trains', {}),
        ('pulse_amp_first_50hz_n', 'float', 'Number of samples represented in pulse_amp_first_50hz', {}),
        ('pulse_amp_first_50hz_std', 'float', 'Standard deviation of samples represented in pulse_amp_stp_initial_50hz', {}),

        ('pulse_amp_stp_initial_50hz', 'float', 'Median amplitude of second pulse on 50 hz trains', {}),
        ('pulse_amp_stp_initial_50hz_n', 'float', 'Number of samples represented in pulse_amp_stp_initial_50hz', {}),
        ('pulse_amp_stp_initial_50hz_std', 'float', 'Standard deviation of samples represented in pulse_amp_stp_initial_50hz', {}),

        ('pulse_amp_stp_induction_50hz', 'float', 'Median amplitude of 6th-8th pulses on 50 hz trains', {}),
        ('pulse_amp_stp_induction_50hz_n', 'float', 'Number of samples represented in pulse_amp_stp_induction_50hz', {}),
        ('pulse_amp_stp_induction_50hz_std', 'float', 'Standard deviation of samples represented in pulse_amp_stp_induction_50hz', {}),

        ('pulse_amp_stp_recovery_250ms', 'float', 'Median amplitude of 9th-12th pulses on 50 hz trains', {}),
        ('pulse_amp_stp_recovery_250ms_n', 'float', 'Number of samples represented in pulse_amp_stp_recovery_250ms', {}),
        ('pulse_amp_stp_recovery_250ms_std', 'float', 'Standard deviation of samples represented in pulse_amp_stp_recovery_250ms', {}),

        ('pulse_amp_stp_recovery_single_250ms', 'float', 'Median amplitude of 9th pulse on 50 hz trains', {}),
        ('pulse_amp_stp_recovery_single_250ms_n', 'float', 'Number of samples represented in pulse_amp_stp_recovery_single_250ms', {}),
        ('pulse_amp_stp_recovery_single_250ms_std', 'float', 'Standard deviation of samples represented in pulse_amp_stp_recovery_single_250ms', {}),

        ('noise_std', 'float', 'Standard deviation of PSP amplitudes measured from background noise', {}),
        ('variability_resting_state', 'float', 'Variability of PSP amplitudes only from events with no preceding spikes for at least 8 seconds, corrected for background noise.', {}),
        ('variability_second_pulse_50hz', 'float', 'Variability of PSP amplitudes in 2nd pulses of 50Hz trains', {}),
        ('variability_stp_induced_state_50hz', 'float', 'Variability of PSP amplitudes in 5th-8th pulses of 50Hz trains', {}),
        ('variability_change_initial_50hz', 'float', 'Difference between variability of 1st and 2nd pulses in 50Hz trains, corrected for background noise.', {}),
        ('variability_change_induction_50hz', 'float', 'Difference between variability of 1st and 5th-8th pulses in 50Hz trains, corrected for background noise.', {}),
        ('paired_event_correlation_1_2_r', 'float', 'Pearson correlation coefficient for amplitudes of 1st:2nd pulses in 50Hz trains.', {}),
        ('paired_event_correlation_1_2_p', 'float', 'Pearson correlation p-value related to paired_event_correlation_1_2_r.', {}),
        ('paired_event_correlation_2_4_r', 'float', 'Pearson correlation coefficient for amplitudes of 1st:2nd pulses in 50Hz trains.', {}),
        ('paired_event_correlation_2_4_p', 'float', 'Pearson correlation p-value related to paired_event_correlation_1_2_r.', {}),
        ('paired_event_correlation_4_8_r', 'float', 'Pearson correlation coefficient for amplitudes of 1st:2nd pulses in 50Hz trains.', {}),
        ('paired_event_correlation_4_8_p', 'float', 'Pearson correlation p-value related to paired_event_correlation_1_2_r.', {}),
        ('stp_all_stimuli', 'object', 'list of initial, induction, and recovery measurements for all stimuli presented', {'deferred': True}),
    ]
)

Pair.dynamics = relationship(Dynamics, back_populates="pair", cascade="delete", single_parent=True, uselist=False)
Dynamics.pair = relationship(Pair, back_populates="dynamics", single_parent=True)
