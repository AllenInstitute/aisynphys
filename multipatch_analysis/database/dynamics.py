from sqlalchemy.orm import relationship
from .database import make_table, TableGroup
from .experiment import Pair


__all__ = ['dynamics_tables', 'Dynamics']


Dynamics = make_table(
    name='dynamics',
    comment="Describes short term dynamics of synaptic connections.",
    columns=[
        ('pair_id', 'pair.id', 'The ID of the cell pair described by each record', {'index': True, 'unique': True}),
        ('pulse_ratio_8_1_50hz', 'float', '8:1 pulse ratio for 50Hz induction', {'index': True}),
        ('pulse_ratio_2_1_50hz', 'float', '2:1 pulse ratio for 50Hz induction', {'index': True}),
        ('pulse_ratio_5_1_50hz', 'float', '5:1 pulse ratio for 50Hz induction', {'index': True}),
        ('pulse_ratio_9_1_125ms', 'float', '9:1 pulse ratio for 50Hz induction at 125ms delay', {'index': True}),
        ('pulse_ratio_9_1_250ms', 'float', '9:1 pulse ratio for 50Hz induction at 250ms delay', {'index': True}),
        ('pulse_ratio_9_1_500ms', 'float', '9:1 pulse ratio for 50Hz induction at 500ms delay', {'index': True}),
        ('pulse_ratio_9_1_1000ms', 'float', '9:1 pulse ratio for 50Hz induction at 1000ms delay', {'index': True}),
        ('pulse_ratio_9_1_2000ms', 'float', '9:1 pulse ratio for 50Hz induction at 2000ms delay', {'index': True}),
        ('pulse_ratio_9_1_4000ms', 'float', '9:1 pulse ratio for 50Hz induction at 4000ms delay', {'index': True}),
        ('pulse_ratio_8_1_10hz', 'float', '8:1 pulse ratio for 10Hz induction', {'index': True}),
        ('pulse_ratio_8_1_20hz', 'float', '8:1 pulse ratio for 20Hz induction', {'index': True}),
        ('pulse_ratio_8_1_100hz', 'float', '8:1 pulse ratio for 100Hz induction', {'index': True}),
        ('pulse_ratio_8_1_200hz', 'float', '8:1 pulse ratio for 200Hz induction', {'index': True}),

    ]
)

Pair.dynamics = relationship(Dynamics, back_populates="pair", cascade="delete", single_parent=True, uselist=False)
Dynamics.pair = relationship(Pair, back_populates="dynamics", single_parent=True)


dynamics_tables = TableGroup([Dynamics])
