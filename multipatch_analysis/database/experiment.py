import os
from collections import OrderedDict
from .database import TableGroup, make_table
from sqlalchemy.orm import relationship, deferred, sessionmaker, aliased
from .slice import Slice
from .. import config


__all__ = ['experiment_tables', 'Experiment', 'Electrode', 'Cell', 'Pair']


class ExperimentBase(object):
    def __getitem__(self, item):
        # Easy cell/pair getters.
        # They're inefficient, but meh.
        if isinstance(item, int):
            for cell in self.cells:
                if cell.ext_id == item:
                    return cell
        elif isinstance(item, tuple):
            for pair in self.pair_lists:
                if item == (pair.pre_cell.ext_id, pair.post_cell.ext_id):
                    return pair
    
    @property
    def cells(self):
        #return {elec.cell.ext_id: elec.cell for elec in self.electrodes if elec.cell is not None}
        return {cell.ext_id: cell for cell in self.cell_list}

    @property
    def pairs(self):
        return {(pair.pre_cell.ext_id, pair.post_cell.ext_id): pair for pair in self.pair_list}

    @property
    def nwb_file(self):
        return os.path.join(config.synphys_data, self.storage_path, self.ephys_file)

    @property
    def nwb_cache_file(self):
        from ..synphys_cache import SynPhysCache
        return SynPhysCache().get_cache(self.nwb_file)

    @property
    def data(self):
        """Data object from NWB file. 
        
        Contains all ephys recordings.
        """

        if not hasattr(self, '_data'):
            from ..data import MultiPatchExperiment
            self._data = MultiPatchExperiment(self.nwb_cache_file)
        return self._data

    @property
    def source_experiment(self):
        """Return the original Experiment object that was used to import
        data into the DB, if available.
        """
        from ..experiment_list import cached_experiments
        return cached_experiments()[self.acq_timestamp]

    def __repr__(self):
        return "<%s %0.3f>" % (self.__class__.__name__, self.acq_timestamp)


Experiment = make_table(
    name='experiment', 
    base=ExperimentBase, 
    comment= "A group of cells patched simultaneously in the same slice.", 
    columns=[
        ('original_path', 'str', 'Original location of raw data on rig.'),
        ('storage_path', 'str', 'Location of data within server or cache storage.'),
        ('ephys_file', 'str', 'Name of ephys NWB file relative to storage_path.'),
        ('rig_name', 'str', 'Identifier for the rig that generated these results.'),
        ('project_name', 'str', 'Name of the project to which this experiment belongs.', {'index': True}),
        ('acq_timestamp', 'float', 'Creation timestamp for site data acquisition folder.', {'unique': True, 'index': True}),
        ('slice_id', 'slice.id', 'ID of the slice used for this experiment', {'index': True}),
        ('target_region', 'str', 'The intended brain region for this experiment'),
        ('internal', 'str', 'The name of the internal solution used in this experiment '
                            '(or "mixed" if more than one solution was used). '
                            'The solution should be described in the pycsf database.', {'index': True}),
        ('acsf', 'str', 'The name of the ACSF solution used in this experiment. '
                        'The solution should be described in the pycsf database.', {'index': True}),
        ('target_temperature', 'float', 'The intended temperature of the experiment (but actual recording temperature is stored elsewhere)'),
        ('date', 'datetime', 'The date of this experiment'),
        ('lims_specimen_id', 'int', 'ID of LIMS "CellCluster" specimen.'),
        #('ext_id', 'str', 'Unique external identifier string for the experiment.')
    ]
)

Slice.experiments = relationship(Experiment, order_by=Experiment.id, back_populates="slice", cascade='save-update,merge,delete')
Experiment.slice = relationship(Slice, back_populates="experiments")


Electrode = make_table(
    name='electrode', 
    comment="Each electrode records a patch attempt, whether or not it resulted in a successful cell recording.",
    columns=[
        ('experiment_id', 'experiment.id', '', {'index': True}),
        ('patch_status', 'str', 'Status of the patch attempt: no seal, low seal, GOhm seal, tech fail, or no attempt'),
        ('start_time', 'datetime', 'The time when recording began for this electrode.'),
        ('stop_time', 'datetime', 'The time when recording ended for this electrode.'),
        ('device_id', 'int', 'External identifier for the device attached to this electrode (usually the MIES A/D channel)'),
        ('internal', 'str', 'The name of the internal solution used in this electrode.'),
        ('initial_resistance', 'float'),
        ('initial_current', 'float'),
        ('pipette_offset', 'float'),
        ('final_resistance', 'float'),
        ('final_current', 'float'),
        ('notes', 'str'),
        ('ext_id', 'int', 'Electrode ID (usually 1-8) referenced in external metadata records'),
    ]
)

Experiment.electrodes = relationship(Electrode, order_by=Electrode.id, back_populates="experiment", cascade='save-update,merge,delete', single_parent=True)
Electrode.experiment = relationship(Experiment, back_populates="electrodes")


Cell = make_table(
    name='cell', 
    comment="Each row represents a single cell in an experiment.",
    columns=[
        ('experiment_id', 'experiment.id', '', {'index': True}),
        ('electrode_id', 'electrode.id', 'ID of the electrode used to patch this cell, if any.', {'index': True}),
        ('cre_type', 'str', 'Comma-separated list of cre drivers apparently expressed by this cell', {'index': True}),
        ('target_layer', 'str', 'The intended cortical layer for this cell (used as a placeholder until the actual layer call is made)', {'index': True}),
        ('is_excitatory', 'bool', 'True if the cell is determined to be excitatory by synaptic current, cre type, or morphology', {'index': True}),
        ('synapse_sign', 'int', 'The sign of synaptic potentials produced by this cell: excitatory=+1, inhibitory=-1, mixed=0'),
        ('patch_start', 'float', 'Time at which this cell was first patched'),
        ('patch_stop', 'float', 'Time at which the electrode was detached from the cell'),
        ('seal_resistance', 'float', 'The seal resistance recorded for this cell immediately before membrane rupture'),
        ('has_biocytin', 'bool', 'If true, then the soma was seen to be darkly stained with biocytin (this indicates a good reseal, but does may not indicate a high-quality fill)'),
        ('has_dye_fill', 'bool', 'Indicates whether the cell was filled with fluorescent dye during the experiment'),
        ('depth', 'float', 'Depth of the cell (in m) from the cut surface of the slice.'),
        ('position', 'object', '3D location of this cell in the arbitrary coordinate system of the experiment'),
        ('ext_id', 'int', 'Cell ID (usually 1-8) referenced in external metadata records', {'index': True}),
    ]
)

Experiment.cell_list = relationship(Cell, order_by=Cell.ext_id, back_populates="experiment", cascade='save-update,merge,delete', single_parent=True)
Cell.experiment = relationship(Experiment, back_populates="cell_list")
Electrode.cell = relationship(Cell, back_populates="electrode", cascade='save-update,merge,delete', single_parent=True, uselist=False)
Cell.electrode = relationship(Electrode, back_populates="cell", single_parent=True)


class PairBase(object):
    def __repr__(self):
        uid = self.experiment.ext_id
        if uid is None or uid == '':
            uid = str('%0.3f'%self.experiment.acq_timestamp if self.experiment.acq_timestamp is not None else None)
        return "<%s %s %d %d>" % (self.__class__.__name__, uid, self.pre_cell.ext_id, self.post_cell.ext_id)


Pair = make_table(
    name='pair',
    base=PairBase,
    comment= "All possible putative synaptic connections. Each pair represents a pre- and postsynaptic cell that were recorded from simultaneously.",
    columns=[
        ('experiment_id', 'experiment.id', '', {'index': True}),
        ('pre_cell_id', 'cell.id', 'ID of the presynaptic cell', {'index': True}),
        ('post_cell_id', 'cell.id', 'ID of the postsynaptic cell', {'index': True}),
        ('synapse', 'bool', 'Whether the experimenter thinks there is a synapse', {'index': True}),
        ('electrical', 'bool', 'Whether the experimenter thinks there is a gap junction', {'index': True}),
        ('crosstalk_artifact', 'float', 'Amplitude of crosstalk artifact measured in current clamp'),
        ('n_ex_test_spikes', 'int', 'Number of QC-passed spike-responses recorded for this pair at excitatory holding potential', {'index': True}),
        ('n_in_test_spikes', 'int', 'Number of QC-passed spike-responses recorded for this pair at inhibitory holding potential', {'index': True}),
        ('synapse_sign', 'int', 'Sign of synaptic current amplitude (+1 for excitatory, -1 for inhibitory', {'index': True}),
        ('distance', 'float', 'Distance between somas (in m)'),
    ]
)

Experiment.pair_list = relationship(Pair, back_populates="experiment", cascade='save-update,merge,delete', single_parent=True)
Pair.experiment = relationship(Experiment, back_populates="pair_list")

Pair.pre_cell = relationship(Cell, foreign_keys=[Pair.pre_cell_id])
#Cell.pre_pairs = relationship(Pair, back_populates="pre_cell", single_parent=True, foreign_keys=[Pair.pre_cell])

Pair.post_cell = relationship(Cell, foreign_keys=[Pair.post_cell_id])
#Cell.post_pairs = relationship(Pair, back_populates="post_cell", single_parent=True, foreign_keys=[Pair.post_cell])        


experiment_tables = TableGroup([Experiment, Electrode, Cell, Pair])

