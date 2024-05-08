import datetime
import os.path
from sqlalchemy.orm import aliased, contains_eager, selectinload
import pandas as pd
from collections import OrderedDict
from .database import Database
from .schema import schema_version, default_sample_rate
from ..synphys_cache import get_db_path, list_db_versions


class SynphysDatabase(Database):
    """Augments the Database class with convenience methods for querying the synphys database.
    """
    default_sample_rate = default_sample_rate
    schema_version = schema_version

    mouse_projects = ["mouse V1 coarse matrix", "mouse V1 pre-production"]
    human_projects = ["human coarse matrix"]

    @classmethod
    def load_sqlite(cls, sqlite_file, readonly=True):
        """Return a SynphysDatabase instance connected to an existing sqlite file.
        """
        ro_host = 'sqlite:///'
        rw_host = None if readonly else ro_host
        return cls(ro_host, rw_host, db_name=sqlite_file)
    
    @classmethod
    def load_current(cls, db_size):
        """Load the most recent version of the database that is supported by this version
        of aisynphys.

        The database file will be downloaded and cached, if an existing cache file is not found.

        Parameters
        ----------
        db_size : str
            Must be one of 'small', 'medium', or 'full'.
        """
        assert db_size in ('small', 'medium', 'full'), "db_size argument must be one of 'small', 'medium', or 'full'"
        current = cls.list_current_versions()
        if db_size not in current:
            raise Exception(f"This version of aisynphys requires database schema version {cls.schema_version}, "
                            "but no released database files were found with this schema version.")
        if current[db_size] is None:
            raise Exception(f"There are no published {db_size} database files for the current schema version ({cls.schema_version}).")
        return cls.load_version(current[db_size]['db_file'])

    @classmethod
    def list_current_versions(cls):
        """Return a dict of the most recent DB versions for each size.

        If no published DB file is available for a particular size, then the value will be set to None
        in the returned dictionary.
        """
        versions_available = list_db_versions()
        current_versions = {}
        for size in ('small', 'medium', 'full'):
            versions_with_size = [desc for desc in versions_available if desc['db_size'] == size and desc['schema_version'] == cls.schema_version]
            if len(versions_with_size) == 0:
                current_versions[size] = None
            else:
                current_versions[size] = versions_with_size[-1]

        return current_versions

    @classmethod
    def list_versions(cls, only_supported=False):
        """Return a list of all available database versions.

        Each item in the list is a dictionary with keys db_file, release_version, db_size, and schema_version.
        """
        vers = list_db_versions()
        if only_supported:
            vers = [v for v in vers if v['schema_version'] == cls.schema_version]
        return vers
    
    @classmethod
    def load_version(cls, db_version):
        """Load a named database version.
        
        Available database names can be listed using :func:`list_versions`.
        The database file will be downloaded and cached, if an existing cache file is not found.        

        Example::

            >>> from aisynphys.database import SynphysDatabase
            >>> SynphysDatabase.list_versions()
            [
                {'db_file': 'synphys_r1.0_small.sqlite'},
                {'db_file': 'synphys_r1.0_medium.sqlite'},
                {'db_file': 'synphys_r1.0_full.sqlite'},
                ...
            ]
            >>> db = SynphysDatabase.load_version('synphys_r1.0_small.sqlite')
            Downloading http://api.brain-map.org/api/v2/well_known_file_download/937779595 =>
              /home/luke/docs/aisynphys/doc/cache/database/synphys_r1.0_2019-08-29_small.sqlite
              [####################]  100.00% (73.13 MB / 73.1 MB)  4.040 MB/s  0:00:00 remaining
              done.

        """
        db_file = get_db_path(db_version)
        db = SynphysDatabase.load_sqlite(db_file, readonly=True)
        return db

    def __init__(self, ro_host, rw_host, db_name, check_schema=True):
        from .schema import ORMBase
        Database.__init__(self, ro_host, rw_host, db_name, ORMBase)
        self._project_names = None
        if check_schema:
            try:
                self._check_version()
            except Exception:
                # need to clean up carefully since this can introduce a persistent connection
                self.default_session.rollback()
                self.dispose_engines()
                raise

    @property
    def version_name(self):
        """The version name of this database, as accepted by SynphysDatabase.load_version()
        """
        return os.path.split(self.db_name)[-1]

    def initialize_database(self):
        """Optionally called after create_tables        
        """
        # initialize or verify db version
        mrec = self.metadata_record()
        if mrec is None:
            mrec = self.Metadata(meta={
                'db_version': schema_version,
                'creation_date': datetime.datetime.now().strftime('%Y-%m-%d'),
                'origin': "Allen Institute for Brain Science / Synaptic Physiology",
            })
            s = self.session(readonly=False)
            s.add(mrec)
            s.commit()
        else:
            self._check_version()

    def _check_version(self):
        mrec = self.metadata_record()
        if mrec is not None:
            ver = mrec.meta['db_version']
            assert ver == schema_version, "Database {self} has unsupported schema version {ver} (expected {schema_version})".format(locals())
    
    def metadata_record(self, session=None):
        if session is None:
            session = self.session()
            close_session = True
        else:
            close_session = False

        try:
            recs = session.query(self.Metadata).all()
            if len(recs) == 0:
                return None
            elif len(recs) > 1:
                raise Exception("Multiple metadata records found.")
            return recs[0]
        finally:
            if close_session:
                session.close()
        
    def slice_from_timestamp(self, ts, session=None):
        session = session or self.default_session
        slices = session.query(self.Slice).filter(self.Slice.acq_timestamp==ts).all()
        if len(slices) == 0:
            raise KeyError("No slice found for timestamp %0.3f" % ts)
        elif len(slices) > 1:
            raise KeyError("Multiple slices found for timestamp %0.3f" % ts)
        
        return slices[0]

    def experiment_from_timestamp(self, ts, session=None):
        session = session or self.default_session
        expts = session.query(self.Experiment).filter(self.Experiment.acq_timestamp==ts).all()
        if len(expts) == 0:
            # For backward compatibility, check for timestamp truncated to 2 decimal places
            for expt in session.query(self.Experiment).all():
                if abs((expt.acq_timestamp - ts)) < 0.01:
                    return expt
            
            raise KeyError("No experiment found for timestamp %0.3f" % ts)
        elif len(expts) > 1:
            raise RuntimeError("Multiple experiments found for timestamp %0.3f" % ts)
        
        return expts[0]

    def experiment_from_ext_id(self, ext_id, session=None):
        session = session or self.default_session
        expts = session.query(self.Experiment).filter(self.Experiment.ext_id==ext_id).all()
        if len(expts) == 0:
            raise KeyError('No experiment found for ext_id %s' %ext_id)
        elif len(expts) > 1:
            raise RuntimeError("Multiple experiments found for ext_id %s" %ext_id)

        return expts[0]

    def slice_from_ext_id(self, ext_id, session=None):
        session = session or self.default_session
        slices = session.query(self.Slice).filter(self.Slice.ext_id==ext_id).all()
        if len(slices) == 0:
            raise KeyError("No slice found for ext_id %s" % ext_id)
        elif len(slices) > 1:
            raise KeyError("Multiple slices found for ext_id %s" % ext_id)
        
        return slices[0]

    def pair_from_ext_id(self, ext_id: "str|tuple", session=None):
        if isinstance(ext_id, str):
            ext_id = ext_id.split(' ')
        session = session or self.default_session
        expt = self.experiment_from_ext_id(ext_id[0], session=session)
        return expt.pairs[ext_id[1], ext_id[2]]
    
    def cell_from_ext_id(self, ext_id: "str|tuple", session=None):
        if isinstance(ext_id, str):
            ext_id = ext_id.split(' ')
        session = session or self.default_session
        expt = self.experiment_from_ext_id(ext_id[0], session=session)
        return expt.cells[ext_id[1]]

    def list_experiments(self, session=None):
        session = session or self.default_session
        return session.query(self.Experiment).all()

    def list_project_names(self, session=None, cache=True):
        session = session or self.default_session
        if cache is False or self._project_names is None:
            self._project_names = [rec[0] for rec in session.query(self.Experiment.project_name).distinct().all()]
        return self._project_names

    def pair_query(self, pre_class=None, post_class=None, synapse=None, synapse_type=None, synapse_probed=None, electrical=None,
                   experiment_type=None, project_name=None, acsf=None, age=None, species=None, distance=None, internal=None,
                   preload=(), session=None, filter_exprs=None):
        """Generate a query for selecting pairs from the database.

        Parameters
        ----------
        pre_class : :class:`aisynphys.cell_class.CellClass` | None
            Filter for pairs where the presynaptic cell belongs to this class
        post_class : :class:`aisynphys.cell_class.CellClass` | None
            Filter for pairs where the postsynaptic cell belongs to this class
        synapse : bool | None
            Include only pairs that are (or are not) connected by a chemical synapse
        synapse_type : str | None
            Include only synapses of a particular type ('ex' or 'in')
        synapse_probed : bool | None
            If True, include only pairs that were probed for a synaptic connection (regardless
            of whether a connectin was found)
        electrical : bool | None
            Include only pairs that are (or are not) connected by an electrical synapse (gap junction)
        experiment_type : str | None
            Include only data from specific types of experiments
        project_name : str | list | None
            Value(s) to match from experiment.project_name (e.g. "mouse V1 coarse matrix" or "human coarse matrix")
        acsf : str | list | None
            Filter for ACSF recipe name(s)
        age : tuple | None
            (min, max) age ranges to filter for. Either limit may be None to disable
            that check.
        species : str | None
            Species ('mouse' or 'human') to filter for
        distance : tuple | None
            (min, max) intersomatic distance in meters
        internal : str | list | None
            Electrode internal solution recipe name(s)
        preload : list
            List of strings specifying resources to preload along with the queried pairs. 
            This can speed up performance in cases where these would otherwise be 
            individually queried later on. Options are:
            - "experiment" (includes experiment and slice)
            - "cell" (includes cell, morphology, cortical_location, and patch_seq)
            - "synapse" (includes synapse, resting_state, dynamics, conductance, and synapse_model)
            - "synapse_prediction" (includes only synapse_prediction)
        filter_exprs : list | None
            List of sqlalchemy expressions, each of which will restrict the query
            via a call to query.filter(expr)
        """

        if experiment_type == 'standard multipatch':
            assert project_name is None, "cannot specify both experiment_type and project_name"
            project_name = ['mouse V1 coarse matrix', 'mouse V1 pre-production', 'human coarse matrix']

        session = session or self.default_session
        pre_cell = aliased(self.Cell, name='pre_cell')
        post_cell = aliased(self.Cell, name='post_cell')
        pre_morphology = aliased(self.Morphology, name='pre_morphology')
        post_morphology = aliased(self.Morphology, name='post_morphology')
        pre_patch_seq = aliased(self.PatchSeq, name='pre_patch_seq')
        post_patch_seq = aliased(self.PatchSeq, name='post_patch_seq')
        pre_intrinsic = aliased(self.Intrinsic, name='pre_intrinsic')
        post_intrinsic = aliased(self.Intrinsic, name='post_intrinsic')
        pre_location = aliased(self.CorticalCellLocation, name='pre_cortical_cell_location')
        post_location = aliased(self.CorticalCellLocation, name='post_cortical_cell_location')

        query = (
            session.query(self.Pair)
            .join(pre_cell, pre_cell.id==self.Pair.pre_cell_id)
            .join(post_cell, post_cell.id==self.Pair.post_cell_id)
            .outerjoin(pre_morphology, pre_morphology.cell_id==pre_cell.id)
            .outerjoin(post_morphology, post_morphology.cell_id==post_cell.id)
            .outerjoin(pre_patch_seq, pre_patch_seq.cell_id==pre_cell.id)
            .outerjoin(post_patch_seq, post_patch_seq.cell_id==post_cell.id)
            .outerjoin(pre_intrinsic, pre_intrinsic.cell_id==pre_cell.id)
            .outerjoin(post_intrinsic, post_intrinsic.cell_id==post_cell.id)
            .outerjoin(pre_location, pre_location.cell_id==pre_cell.id)
            .outerjoin(post_location, post_location.cell_id==post_cell.id)
            .join(self.Experiment, self.Pair.experiment_id==self.Experiment.id)
            .outerjoin(self.Slice, self.Experiment.slice_id==self.Slice.id) ## don't want to drop all pairs if we don't have slice or connection strength entries
            .outerjoin(self.Synapse, self.Synapse.pair_id==self.Pair.id)
            .outerjoin(self.SynapsePrediction, self.SynapsePrediction.pair_id==self.Pair.id)
            .outerjoin(self.Dynamics, self.Dynamics.pair_id==self.Pair.id)
            .outerjoin(self.RestingStateFit, self.RestingStateFit.synapse_id==self.Synapse.id)
            .outerjoin(self.SynapseModel, self.SynapseModel.pair_id==self.Pair.id)
            .outerjoin(self.Conductance, self.Conductance.synapse_id==self.Synapse.id)
            # .outerjoin(self.PolySynapse)
            # .outerjoin(self.GapJunction)
        )

        if pre_class is not None:
            query = pre_class.filter_query(query, pre_cell, db=self)

        if post_class is not None:
            query = post_class.filter_query(query, post_cell, db=self)

        if synapse is not None:
            query = query.filter(self.Pair.has_synapse==synapse)

        if synapse_type is not None:
            assert synapse_type in ['ex', 'in', 'mixed'], "synapse_type must be 'ex', 'in', or 'mixed'"
            query = query.filter(self.Synapse.synapse_type==synapse_type)

        if synapse_probed is True:
            from ..connectivity import probed_pair_test_spike_limit
            query = query.filter(
                ((pre_cell.cell_class == 'ex') & (self.Pair.n_ex_test_spikes > probed_pair_test_spike_limit)) |
                ((pre_cell.cell_class == 'in') & (self.Pair.n_in_test_spikes > probed_pair_test_spike_limit)) |
                ((self.Pair.n_ex_test_spikes > probed_pair_test_spike_limit) & (self.Pair.n_in_test_spikes > probed_pair_test_spike_limit))
            )

        if electrical is not None:
            query = query.filter(self.Pair.has_electrical==electrical)

        if project_name is not None:
            names = [project_name] if isinstance(project_name, str) else project_name
            for name in names:
                assert name in self.list_project_names(), "project_name %r not found in database (see SynphysDatabase.list_project_names)" % name

            if isinstance(project_name, str):
                query = query.filter(self.Experiment.project_name==project_name)
            else:
                query = query.filter(self.Experiment.project_name.in_(project_name))

        if acsf is not None:
            if isinstance(acsf, str):
                query = query.filter(self.Experiment.acsf==acsf)
            else:
                query = query.filter(self.Experiment.acsf.in_(acsf))

        if age is not None:
            if age[0] is not None:
                query = query.filter(self.Slice.age>=age[0])
            if age[1] is not None:
                query = query.filter(self.Slice.age<=age[1])

        if distance is not None:
            if distance[0] is not None:
                query = query.filter(self.Pair.distance>=distance[0])
            if distance[1] is not None:
                query = query.filter(self.Pair.distance<=distance[1])

        if species is not None:
            query = query.filter(self.Slice.species==species)

        if internal is not None:
            if isinstance(internal, str):
                query = query.filter(self.Experiment.internal==internal)
            else:
                query = query.filter(self.Experiment.internal.in_(internal))

        if filter_exprs is not None:
            for expr in filter_exprs:
                query = query.filter(expr)
        
        preload_allowed = ['experiment', 'cell', 'synapse', 'synapse_prediction']
        if any([x not in preload_allowed for x in preload]):
            raise ValueError(f"Preload allowed values are: {preload_allowed}")

        if 'experiment' in preload:
            query = (
                query
                .add_entity(self.Experiment)
                .add_entity(self.Slice)
            )
            query = query.options(
                contains_eager(self.Pair.experiment),
                contains_eager(self.Experiment.slice),
            )

        if 'cell' in preload:
            query = (
                query
                .add_entity(pre_cell)
                .add_entity(post_cell)
                .add_entity(pre_morphology)
                .add_entity(post_morphology)
                .add_entity(pre_patch_seq)
                .add_entity(post_patch_seq)
                .add_entity(pre_intrinsic)
                .add_entity(post_intrinsic)
                .add_entity(pre_location)
                .add_entity(post_location)
            )
            query = query.options(
                contains_eager(self.Pair.pre_cell, alias=pre_cell), 
                contains_eager(self.Pair.post_cell, alias=post_cell), 
                contains_eager(pre_cell.morphology, alias=pre_morphology), 
                contains_eager(post_cell.morphology, alias=post_morphology), 
                contains_eager(pre_cell.patch_seq, alias=pre_patch_seq), 
                contains_eager(post_cell.patch_seq, alias=post_patch_seq), 
                contains_eager(pre_cell.intrinsic, alias=pre_intrinsic),
                contains_eager(post_cell.intrinsic, alias=post_intrinsic),
                contains_eager(pre_cell.cortical_location, alias=pre_location),
                contains_eager(post_cell.cortical_location, alias=post_location),
            )

        if 'synapse' in preload:
            query = (
                query
                .add_entity(self.Synapse)
                .add_entity(self.RestingStateFit)
                .add_entity(self.Dynamics)
                .add_entity(self.Conductance)
                .add_entity(self.SynapseModel)
            )
            query = query.options(
                contains_eager(self.Pair.synapse),
                contains_eager(self.Synapse.resting_state_fit),
                contains_eager(self.Pair.dynamics),
                contains_eager(self.Synapse.conductance),
                contains_eager(self.Pair.synapse_model),
            )

        if 'synapse_prediction' in preload:
            query = (
                query
                .add_entity(self.SynapsePrediction)
            )
            query = query.options(
                contains_eager(self.Pair.synapse_prediction),
            )

        # package the aliased cells
        query.pre_cell = pre_cell
        query.post_cell = post_cell
        query.pre_morphology = pre_morphology
        query.post_morphology = post_morphology
        query.pre_location = pre_location
        query.post_location = post_location
        query.pre_intrinsic = pre_intrinsic
        query.post_intrinsic = post_intrinsic
        query.pre_patch_seq = pre_patch_seq
        query.post_patch_seq = post_patch_seq

        return query

    def matrix_pair_query(self, pre_classes, post_classes, columns=None, pair_query_args=None):
        """Returns the concatenated result of running pair_query over every combination
        of presynaptic and postsynaptic cell class.
        """
        if pair_query_args is None:
            pair_query_args = {}

        pairs = None
        for pre_name, pre_class in pre_classes.items():
            for post_name, post_class in post_classes.items():
                pair_query = self.pair_query(
                    pre_class=pre_class,
                    post_class=post_class,
                    **pair_query_args
                )
                
                if columns is not None:
                    pair_query = pair_query.add_columns(*columns)
                
                df = pair_query.dataframe(rename_columns=False)
                df['pre_class'] = pre_name
                df['post_class'] = post_name
                if pairs is None:
                    pairs = df
                else:
                    pairs = pd.concat([pairs, df], axis=0, join='outer')
        
        return pairs

    def __getstate__(self):
        """Allows DB to be pickled and passed to subprocesses.
        """
        return {
            'ro_host': self.ro_host, 
            'rw_host': self.rw_host, 
            'db_name': self.db_name,
        }

    def __setstate__(self, state):
        self.__init__(ro_host=state['ro_host'], rw_host=state['rw_host'], db_name=state['db_name'])
