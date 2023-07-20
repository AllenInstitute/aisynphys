import os
from collections import OrderedDict
from sqlalchemy.orm import relationship, deferred, sessionmaker, aliased
from ... import config, constants, synphys_cache
from . import make_table
from .slice import Slice


__all__ = ['Experiment', 'Electrode', 'Cell', 'Pair']


class ExperimentBase(object):
    @property
    def cells(self):
        #return {elec.cell.ext_id: elec.cell for elec in self.electrodes if elec.cell is not None}
        return {cell.ext_id: cell for cell in self.cell_list}

    @property
    def pairs(self):
        return {(pair.pre_cell.ext_id, pair.post_cell.ext_id): pair for pair in self.pair_list}

    @property
    def nwb_file(self):
        if self.ephys_file is None:
            return None
        if config.synphys_data is not None:
            # Return path from local file repo
            return os.path.join(config.synphys_data, self.storage_path, self.ephys_file)
        else:
            # return file cached from download
            return synphys_cache.get_nwb_path(self.ext_id)

    @property
    def data(self):
        """Data object from NWB file. 
        
        Contains all ephys recordings.
        """

        if not hasattr(self, '_data'):
            from ...data import MultiPatchDataset
            if self.nwb_file is None:
                self._data = None
            else:
                self._data = MultiPatchDataset(self.nwb_file)
        return self._data

    @property
    def path(self):
        """Filesystem path to the root of this experiment.
        """
        return os.path.join(config.synphys_data, self.storage_path)
        
    @property
    def original_path(self):
        """The original path where this experiment was acquired. 
        """
        ss = os.path.join(self.path, 'sync_source')
        if os.path.isfile(ss):
            return open(ss, 'rb').read().decode('latin1')
        else:
            return self.path

    def __repr__(self):
        if self.ext_id is not None:
            return "<%s %s>" % (self.__class__.__name__, self.ext_id)
        else:
            return "<%s %0.3f>" % (self.__class__.__name__, self.acq_timestamp)


Experiment = make_table(
    name='experiment', 
    base=ExperimentBase, 
    comment= "A group of cells patched simultaneously in the same slice.", 
    columns=[
        ('ext_id', 'str', 'Unique external identifier string for the experiment.', {'unique': True, 'index': True}),
        ('slice_id', 'slice.id', 'ID of the slice used for this experiment', {'index': True}),
        ('project_name', 'str', 'Name of the project to which this experiment belongs.', {'index': True}),
        ('date', 'datetime', 'The date of this experiment'),
        ('target_region', 'str', 'The intended brain region for this experiment'),
        ('internal', 'str', 'The name of the internal solution used in this experiment '
                            '(or "mixed" if more than one solution was used). '
                            'The solution should be described in the pycsf database.', {'index': True}),
        ('acsf', 'str', 'The name of the ACSF solution used in this experiment. '
                        'The solution should be described in the pycsf database.', {'index': True}),
        ('target_temperature', 'float', 'The intended temperature of the experiment (measured temperature per-recording is stored elsewhere)'),
        ('rig_name', 'str', 'Identifier for the rig that generated these results.'),
        ('operator_name', 'str', 'Opertator that generated these results.'),
        ('storage_path', 'str', 'Location of data within server or cache storage.'),
        ('ephys_file', 'str', 'Name of ephys NWB file relative to storage_path.'),
        ('acq_timestamp', 'float', 'Creation timestamp for site data acquisition folder.', {'unique': True, 'index': True}),
    ]
)

Slice.experiments = relationship(Experiment, order_by=Experiment.id, back_populates="slice", cascade='save-update,merge,delete')
Experiment.slice = relationship(Slice, back_populates="experiments")


Electrode = make_table(
    name='electrode', 
    comment="Each electrode records a patch attempt, whether or not it resulted in a successful cell recording.",
    columns=[
        ('experiment_id', 'experiment.id', '', {'index': True}),
        ('ext_id', 'str', 'Electrode ID (usually 1-8) referenced in external metadata records'),
        ('patch_status', 'str', 'Status of the patch attempt: no seal, low seal, GOhm seal, tech fail, or no attempt'),
        ('start_time', 'datetime', 'The time when recording began for this electrode.'),
        ('stop_time', 'datetime', 'The time when recording ended for this electrode.'),
        ('device_id', 'int', 'External identifier for the device attached to this electrode (usually the MIES A/D channel)'),
        # ('internal', 'str', 'The name of the internal solution used in this electrode.'),
        # ('initial_resistance', 'float'),
        # ('initial_current', 'float'),
        # ('pipette_offset', 'float'),
        # ('final_resistance', 'float'),
        # ('final_current', 'float'),
        # ('notes', 'str'),
    ]
)

Experiment.electrodes = relationship(Electrode, order_by=Electrode.id, back_populates="experiment", cascade='save-update,merge,delete', single_parent=True)
Electrode.experiment = relationship(Experiment, back_populates="electrodes")


class CellBase(object):
    def __repr__(self):
        expt = self.experiment
        if expt is None:
            uid = "<no experiment>"
        else:
            uid = getattr(self.experiment, 'ext_id', None)
            if uid is None or uid == '':
                uid = str('%0.3f'%self.experiment.acq_timestamp if self.experiment.acq_timestamp is not None else None)
        return "<%s %s %s>" % (self.__class__.__name__, uid, self.ext_id)

    def _infer_cell_classes(self):
        """Return cell_class and cell_class_nonsynaptic based on the current contents 
        of cell.meta. 

        This is used internally by pipeline modules to update the cell_class fields whenever new evidence
        is added.
        """
        tm_classes = set()
        tms_classes = set()
        for name in ['transgenic', 'morpho', 'synaptic', 'transcriptomic']:
            cell_cls = self.meta.get(name + '_cell_class', None)
            if cell_cls is None:
                continue
            tms_classes.add(cell_cls)
            if name != 'synaptic':
                tm_classes.add(cell_cls)
            
        if len(tm_classes) == 0:
            tm_class = None
        else:
            tm_class = 'mixed' if len(tm_classes) > 1 else list(tm_classes)[0]
        
        if len(tms_classes) == 0:
            tms_class = None
        else:
            tms_class = 'mixed' if len(tms_classes) > 1 else list(tms_classes)[0]

        return tms_class, tm_class

    @property
    def pairs(self):
        cell_pairs = {}
        expt_pairs = self.experiment.pair_list
        for pair in expt_pairs:
            if self.ext_id in [pair.pre_cell.ext_id, pair.post_cell.ext_id]:
                cell_pairs[(pair.pre_cell.ext_id, pair.post_cell.ext_id)] = pair
        return cell_pairs


Cell = make_table(
    name='cell', 
    comment="Each row represents a single cell in an experiment.",
    base=CellBase,
    columns=[
        ('experiment_id', 'experiment.id', '', {'index': True}),
        ('ext_id', 'str', 'Cell ID (usually 1-8) referenced in external metadata records', {'index': True}),
        ('electrode_id', 'electrode.id', 'ID of the electrode used to patch this cell, if any.', {'index': True}),
        ('cre_type', 'str', 'Comma-separated list of cre drivers apparently expressed by this cell', {'index': True}),
        ('target_layer', 'str', 'The intended cortical layer for this cell (used as a placeholder until the actual layer call is made)', {'index': True}),
        ('position', 'object', '3D location of this cell in the arbitrary coordinate system of the experiment'),
        ('depth', 'float', 'Depth of the cell (in m) from the cut surface of the slice.'),
        ('cell_class', 'str', 'Cell class "ex", "in", or "mixed" determined by synaptic current, cre type, or morphology. '
         'This property makes use of synaptic currents to define cell class; it should _not_ be used when measuring connection probability.', {'index': True}),
        ('cell_class_nonsynaptic', 'str', 'Cell class "ex", "in", or "mixed" determined by cre type or morphology. '
         'Unlike `cell_class`, this property excludes synaptic currents as a determinant so that it can be used in measurements of connectivity.', {'index': True}),
        # ('patch_start', 'float', 'Time at which this cell was first patched'),
        # ('patch_stop', 'float', 'Time at which the electrode was detached from the cell'),
        # ('seal_resistance', 'float', 'The seal resistance recorded for this cell immediately before membrane rupture'),
        # ('has_biocytin', 'bool', 'If true, then the soma was seen to be darkly stained with biocytin (this indicates a good reseal, but does may not indicate a high-quality fill)'),
        # ('has_dye_fill', 'bool', 'Indicates whether the cell was filled with fluorescent dye during the experiment'),
    ]
)

Experiment.cell_list = relationship(Cell, order_by=Cell.ext_id, back_populates="experiment", cascade='save-update,merge,delete', single_parent=True)
Cell.experiment = relationship(Experiment, back_populates="cell_list")
Electrode.cell = relationship(Cell, back_populates="electrode", cascade='save-update,merge,delete', single_parent=True, uselist=False)
Cell.electrode = relationship(Electrode, back_populates="cell", single_parent=True)


class PairBase(object):
    @property
    def ext_id(self):
        """The external ID for this cell pair (this ID is stable across database versions)
        """
        return (self.experiment.ext_id, self.pre_cell.ext_id, self.post_cell.ext_id)
        
    def __repr__(self):
        uid = getattr(self.experiment, 'ext_id', None)
        if uid is None or uid == '':
            uid = str('%0.3f'%self.experiment.acq_timestamp if self.experiment.acq_timestamp is not None else None)
        return "<%s %s %s %s>" % (self.__class__.__name__, uid, self.pre_cell.ext_id, self.post_cell.ext_id)


Pair = make_table(
    name='pair',
    base=PairBase,
    comment= "An ordered pair of cells, possibly connected by a synapse or gap junction.",
    columns=[
        ('experiment_id', 'experiment.id', '', {'index': True}),
        ('pre_cell_id', 'cell.id', 'ID of the presynaptic cell', {'index': True}),
        ('post_cell_id', 'cell.id', 'ID of the postsynaptic cell', {'index': True}),
        ('has_synapse', 'bool', 'Whether a chemical monosynaptic connection was manually detected for this cell pair', {'index': True}),
        ('has_polysynapse', 'bool', 'Whether a polysynaptic connection was manually detected for this cell pair', {'index': True}),
        ('has_electrical', 'bool', 'Whether an electrical synapse / gap junction was manually detected for this cell pair', {'index': True}),
        ('crosstalk_artifact', 'float', 'Amplitude of crosstalk artifact measured in current clamp'),
        ('n_ex_test_spikes', 'int', 'Number of QC-passed spike-responses recorded for this pair at excitatory holding potential', {'index': True}),
        ('n_in_test_spikes', 'int', 'Number of QC-passed spike-responses recorded for this pair at inhibitory holding potential', {'index': True}),
        ('distance', 'float', 'Distance between somas (in m)'),
        ('lateral_distance', 'float', 'Distance between somas perpendicular to the pia-wm axis (in m)'),
        ('vertical_distance', 'float', 'Distance between somas along the pia-wm axis (in m)'),
        ('reciprocal_id', 'pair.id', 'ID of the reciprocal to this cell pair (the pair with pre_cell and post_cell swapped)', {'index': True}),
    ]
)

Experiment.pair_list = relationship(Pair, back_populates="experiment", cascade='save-update,merge,delete', single_parent=True)
Pair.experiment = relationship(Experiment, back_populates="pair_list")
Pair.pre_cell = relationship(Cell, foreign_keys=[Pair.pre_cell_id], uselist=False)
Pair.post_cell = relationship(Cell, foreign_keys=[Pair.post_cell_id], uselist=False)
# docs on handling mutually-dependent relationships: 
# https://docs.sqlalchemy.org/en/13/orm/relationship_persistence.html#post-update
Pair.reciprocal = relationship(Pair, foreign_keys=[Pair.reciprocal_id], uselist=False, post_update=True)
