# coding: utf8
"""
For generating a DB table describing cell morphology.

"""
from __future__ import print_function, division

import os, pyodbc, datetime
from collections import OrderedDict
from ..util import timestamp_to_datetime
from .. import database as db
from ..pipette_metadata import PipetteMetadata
from .. import config, lims
from .pipeline_module import DatabasePipelineModule
from .experiment import ExperimentPipelineModule

col_names = {
                'Qual_Morpho_Type': {'name': 'qual_morpho_type', 'type': 'str'},
                'dendrite_type': {'name': 'dendrite_type', 'type': ['spiny', 'aspiny', 'sparsely spiny', 'NEI']},
                'Apical_truncation_distance': {'name': 'apical_trunc_distance', 'type': 'float'},
                'Axon_truncation_distance': {'name': 'axon_trunc_distance', 'type': 'float'},
                'Axon origination': {'name': 'axon_origin', 'type': ['soma', 'dendrite', 'NEI']},
                'Axon_truncation': {'name': 'axon_truncation', 'type': ['truncated', 'borderline', 'intact', 'unclear', 'NEI']},
                'Apical_truncation': {'name': 'apical_truncation', 'type': ['truncated', 'borderline', 'intact','unclear', 'NEI']},

            }

class MorphologyPipelineModule(DatabasePipelineModule):
    """Imports cell morphology data for each experiment
    """
    name = 'morphology'
    dependencies = [ExperimentPipelineModule]
    table_group = db.morphology_tables
    
    @classmethod
    def create_db_entries(cls, job_id, session):
        
        # Load experiment from DB
        expt = db.experiment_from_timestamp(job_id, session=session)
        morpho_results = morpho_db()
        
        path = os.path.join(config.synphys_data, expt.storage_path)
        pip_meta = PipetteMetadata(path)

        for cell_id,cell in expt.cells.items():
            # How the experimenter described the morphology
            user_morpho = pip_meta.pipettes[cell.ext_id].get('morphology')
            cell_specimen_id = cell.meta.get('lims_specimen_id')
            cell_morpho = morpho_results.get(cell_specimen_id)
            if cell_specimen_id is not None:
                cortical_layer = lims.cell_layer(cell_specimen_id)
                if cortical_layer is not None:
                    cortical_layer = cortical_layer.lstrip('Layer')
            else:
                cortical_layer = None

            if user_morpho in (None, ''):
                pyramidal = None
            elif user_morpho == 'pyr':
                pyramidal = True
            else:
                print("Unknown morphology string: %s" % user_morpho)
                pyramidal = None

            results = {
                'pyramidal': pyramidal,
                'cortical_layer': cortical_layer,
                'morpho_db_hash': hash(None),

            }
            
            if cell_morpho is not None:
                morpho_db_hash = hash(';'.join(filter(None, cell_morpho)))
                results['morpho_db_hash'] = morpho_db_hash
                for morpho_db_name, result in col_names.items():
                    col_name = result['name']
                    col_type = result['type']
                    data = getattr(cell_morpho, morpho_db_name)
                    if data is not None:
                        if isinstance(col_type, list):
                            d = [t for t in col_type if t == data]
                            if len(d) == 1:
                                data = d[0]
                            else:
                                raise Exception ('Error parsing morphology annotation %s for cell %d from column %s' % (data, cell_specimen_id, morpho_db_name))
                        elif col_type == 'float':
                            try:
                                data = float(data)
                            except ValueError:
                                raise Exception ('Error parsing morphology annotation %s for cell %d from column %s' % (data, cell_specimen_id, morpho_db_name))
                    results[col_name] = data

            # Write new record to DB
            morphology = db.Morphology(cell_id=cell.id, **results)
            session.add(morphology)
        
    @classmethod
    def job_records(cls, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        return session.query(db.Morphology).filter(db.Morphology.cell_id==db.Cell.id).filter(db.Cell.experiment_id==db.Experiment.id).filter(db.Experiment.acq_timestamp.in_(job_ids)).all()

    @classmethod
    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        # All experiments and their creation times in the DB
        expts = ExperimentPipelineModule.finished_jobs()

        # Look up nwb file locations for all experiments
        session = db.Session()
        # expts = session.query(db.Experiment).filter(db.Experiment.acq_timestamp==1521667891.153).all()
        session.rollback()
        morpho_results = morpho_db()
        # Return the greater of NWB mod time and experiment DB record mtime
        ready = OrderedDict()
        for expt_id, (expt_mtime, success) in expts.items():
            if success is not True:
                continue

            expt = session.query(db.Experiment).filter(db.Experiment.acq_timestamp==expt_id).all()[0]
            cluster = expt.lims_specimen_id
            if cluster is None:
                continue
            cluster_cells = lims.cluster_cells(cluster)
            if cluster_cells is None:
                continue
            cell_hash_compare = []
            for cell in cluster_cells:
                if cell.id is None:
                    continue
                morpho_db_hash = hash(';'.join(filter(None, morpho_results.get(cell.id, [None]))))
                cell_morpho_rec = session.query(db.Morphology).join(db.Cell).filter(db.Cell.meta.info.get('lims_specimen_id')==str(cell.id)).all()
                if len(cell_morpho_rec) == 1:   
                    cell_rec_hash = cell_morpho_rec[0].morpho_db_hash
                    cell_hash_compare.append(morpho_db_hash == cell_rec_hash)
                else:
                    cell_hash_compare.append(False)
            if all(cell_hash_compare) is False:
                ready[expt.acq_timestamp] = datetime.datetime.now()
    
        return ready

def morpho_db():
    cnxn_str = r'DRIVER={Microsoft Access Driver (*.mdb, *.accdb)}; DBQ=%s' % config.morpho_address
    cnxn = pyodbc.connect(cnxn_str)
    cursor = cnxn.cursor()
    morpho_table = cursor.execute('select * from MPATCH_CellsofCluster')
    results = morpho_table.fetchall()
    morpho_results = {int(r.cell_specimen_id): r for r in results}

    return morpho_results