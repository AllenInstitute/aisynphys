from .database import Session, aliased, or_, and_, default_session, get_default_session, reset_db, vacuum, dispose_engines, default_sample_rate, db_name, bake_sqlite, clone_database

# Import table definitions from DB modules
from .pipeline import *
from .slice import *
from .experiment import *
from .morphology import *
from .dataset import *
from .pulse_response_strength import *
from .dynamics import *
from .connection_strength import *
from .avg_first_pulse_fit import *
from .single_first_pulse_fit import *
from .gap_junction import *


@default_session
def slice_from_timestamp(ts, session=None):
    slices = session.query(Slice).filter(Slice.acq_timestamp==ts).all()
    if len(slices) == 0:
        raise KeyError("No slice found for timestamp %0.3f" % ts)
    elif len(slices) > 1:
        raise KeyError("Multiple slices found for timestamp %0.3f" % ts)
    
    return slices[0]


@default_session
def experiment_from_timestamp(ts, session=None):
    expts = session.query(Experiment).filter(Experiment.acq_timestamp==ts).all()
    if len(expts) == 0:
        # For backward compatibility, check for timestamp truncated to 2 decimal places
        for expt in session.query(Experiment).all():
            if abs((expt.acq_timestamp - ts)) < 0.01:
                return expt
        
        raise KeyError("No experiment found for timestamp %0.3f" % ts)
    elif len(expts) > 1:
        raise RuntimeError("Multiple experiments found for timestamp %0.3f" % ts)
    
    return expts[0]

@default_session
def experiment_from_uid(uid, session=None):
    expts = session.query(Experiment).filter(Experiment.ext_id==uid).all()
    if len(expts) == 0:
        raise KeyError('No experiment found for uid %s' %uid)
    elif len(expts) > 1:
        raise RuntimeError("Multiple experiments found for uid %s" %uid)

    return expts[0]


@default_session
def list_experiments(session=None):
    return session.query(Experiment).all()


def query(*args, **kwds):
    return get_default_session().query(*args, **kwds)


def rollback():
    get_default_session().rollback()


@default_session
def pair_query(pre_class=None, post_class=None, synapse=None, electrical=None, project_name=None, acsf=None, age=None, species=None, distance=None, internal=None, session=None):
    """Generate a query for selecting pairs from the database.

    Parameters
    ----------
    pre_class : CellClass | None
        Filter for pairs where the presynaptic cell belongs to this class
    post_class : CellClass | None
        Filter for pairs where the postsynaptic cell belongs to this class
    synapse : bool | None
        Include only pairs that are (or are not) connected by a chemical synapse
    electrical : bool | None
        Include only pairs that are (or are not) connected by an electrical synapse (gap junction)
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
    
    """
    pre_cell = aliased(Cell, name='pre_cell')
    post_cell = aliased(Cell, name='post_cell')
    pre_morphology = aliased(Morphology, name='pre_morphology')
    post_morphology = aliased(Morphology, name='post_morphology')
    query = session.query(
        Pair,
        # pre_cell,
        # post_cell,
        # pre_morphology,
        # post_morphology,
        # Experiment,
        # ConnectionStrength,
    )
    query = query.join(pre_cell, pre_cell.id==Pair.pre_cell_id)
    query = query.join(post_cell, post_cell.id==Pair.post_cell_id)
    query = query.join(pre_morphology, pre_morphology.cell_id==pre_cell.id)
    query = query.join(post_morphology, post_morphology.cell_id==post_cell.id)
    query = query.join(Experiment, Pair.experiment_id==Experiment.id)
    query = query.outerjoin(Slice, Experiment.slice_id==Slice.id) ## don't want to drop all pairs if we don't have slice or connection strength entries
    query = query.outerjoin(ConnectionStrength)
    

    if pre_class is not None:
        query = pre_class.filter_query(query, pre_cell)

    if post_class is not None:
        query = post_class.filter_query(query, post_cell)

    if synapse is not None:
        query = query.filter(Pair.synapse==synapse)

    if electrical is not None:
        query = query.filter(Pair.electrical==electrical)

    if project_name is not None:
        if isinstance(project_name, str):
            query = query.filter(Experiment.project_name==project_name)
        else:
            query = query.filter(Experiment.project_name.in_(project_name))

    if acsf is not None:
        if isinstance(acsf, str):
            query = query.filter(Experiment.acsf==acsf)
        else:
            query = query.filter(Experiment.acsf.in_(acsf))

    if age is not None:
        if age[0] is not None:
            query = query.filter(Slice.age>=age[0])
        if age[1] is not None:
            query = query.filter(Slice.age<=age[1])

    if distance is not None:
        if distance[0] is not None:
            query = query.filter(Pair.distance>=distance[0])
        if distance[1] is not None:
            query = query.filter(Pair.distance<=distance[1])

    if species is not None:
        query = query.filter(Slice.species==species)

    if internal is not None:
        if isinstance(internal, str):
            query = query.filter(Experiment.internal==internal)
        else:
            query = query.filter(Experiment.internal.in_(internal))

    return query


def dataframe(query):
    import pandas
    return pandas.read_sql(query.statement, query.session.bind)

