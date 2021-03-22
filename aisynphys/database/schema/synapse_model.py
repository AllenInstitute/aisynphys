from sqlalchemy.orm import relationship
from . import make_table
from .experiment import Pair


__all__ = ['SynapseModel']


SynapseModel = make_table(
    name='synapse_model',
    comment="Summary of stochastic model outputs per synapse",
    columns=[
        ('pair_id', 'pair.id', 'The ID of the cell pair described by each record', {'index': True, 'unique': True}),
        ('n_source_events', 'int', 'Number of qc-passed pulse response amplitudes used to fit the model', {'index': True}),

        ('max_likelihood', 'object', 'Contains maximum likelihood parameter values'),
        ('marginal_distributions', 'object', 'Contains marginal distributions for all model parameters'),
        ('confidence_intervals', 'object', 'Contains confidence intervals for all model parameters'),
        ('ml_quanta_per_spike', 'float', 'maximum likelihood value of n_release_sites * base_release_probability'),
        ('ml_sites_pr_ratio', 'float', 'maximum likelihood ratio of n_release_sites : base_release_probability'),
        

        # STP metrics generated from model simulation with max likelihood parameters
        ('paired_pulse_ratio_50hz', 'float', 'The median ratio of 2nd / 1st pulse amplitudes for 50Hz pulse trains.', {'index': True}),
        ('stp_initial_50hz', 'float', 'The median relative change from 1st to 2nd pulse for 50Hz pulse trains', {'index': True}),
        ('stp_initial_50hz_n', 'float', 'Number of samples represented in stp_initial_50Hz', {'index': True}),
        ('stp_initial_50hz_std', 'float', 'Standard deviation of samples represented in stp_initial_50Hz', {'index': True}),
        ('stp_induction_50hz', 'float', 'The median relative change from 1st to 5th-8th pulses for 50Hz pulse trains', {'index': True}),
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
        ('stp_all_stimuli', 'object', 'list of initial, induction, and recovery measurements for all stimuli presented'),
    ]
)

Pair.synapse_model = relationship(SynapseModel, back_populates="pair", cascade="delete", single_parent=True, uselist=False)
SynapseModel.paor = relationship(Pair, back_populates="synapse_model", single_parent=True)
