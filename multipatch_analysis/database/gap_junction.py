from sqlalchemy.orm import relationship
from .database import make_table, TableGroup
from .experiment import Pair


__all__ = ['gap_junction_tables', 'GapJunction']


GapJunction = make_table(
    name='gap_junction',
    comment= "Describes the presence of gap junction.",
    columns=[
        ('pair_id', 'pair.id', 'The ID of the entry in the pair table to which these results apply', {'index': True}),
        ('gap_junction', 'bool', 'Whether gap junction is present', {'index': True}),
        ('strength', 'float', 'Strenght of gap junction', {'index': True})
    ]
)

Pair.gap_junction = relationship(GapJunction, back_populates="pair", cascade="delete", single_parent=True, uselist=False)
GapJunction.pair = relationship(Pair, back_populates="gap_junction", single_parent=True)

gap_junction_tables = TableGroup([GapJunction])
