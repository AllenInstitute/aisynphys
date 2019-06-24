from sqlalchemy.orm import relationship
from .database import TableGroup, make_table, Column, Integer, Boolean, ForeignKey
from .experiment import Cell


__all__ = ['morphology_tables', 'Morphology']

Morphology = make_table(name='morphology', comment="Describes morphological properties of cells.", columns=[
    ('cell_id', 'cell.id', 'The ID of the cell described by each record', {'index': True, 'unique': True}),
    ('morpho_db_hash', 'int', 'hash of Morphology team database for update tracking'),
    ('pyramidal', 'bool', 'whether the experimenter labeled the cell as pyramidal', {'index': True}),
    ('cortical_layer', 'str', 'cortical layer of cell defined by layer drawing annotation', {'index': True}),
    ('qual_morpho_type', 'str', 'qualitative desctription of cell morphology', {'index': True}),
    ('dendrite_type', 'str', 'dendrite type of cell: spiny, aspiny, sparsely spiny', {'index': True}),
    ('apical_trunc_distance', 'float', 'measured distance to truncation of apical dendrite', {'index': True}),
    ('axon_trunc_distance', 'float', 'measured distance to truncation of axon', {'index': True}),
    ('apical_truncation', 'str', 'qualitative description of apical dendrite truncation', {'index': True}),
    ('axon_truncation', 'str', 'qualitative description of axon truncation', {'index': True}),
    ('axon_origin', 'str', 'origination of axon; soma, dendrite, etc', {'index': True}),

])

Cell.morphology = relationship(Morphology, back_populates="cell", cascade="delete", single_parent=True, uselist=False)
Morphology.cell = relationship(Cell, back_populates="morphology", single_parent=True)

morphology_tables = TableGroup([Morphology])
