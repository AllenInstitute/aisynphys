from __future__ import print_function, division
from datetime import datetime
import pyqtgraph as pg

from multipatch_analysis import database as db


class ExperimentBrowser(pg.TreeWidget):
    """TreeWidget showing a list of experiments with cells and pairs.
    """
    # TODO: add filtering options, cell info, context menu of actions, etc.
    
    def __init__(self):
        pg.TreeWidget.__init__(self)
        self.setColumnCount(7)
        self.setHeaderLabels(['date', 'timestamp', 'rig', 'organism', 'region', 'genotype', 'acsf'])
        self.setDragDropMode(self.NoDragDrop)
        self.populate()
        self._last_expanded = None
        
    def populate(self, experiments=None):
        self.items_by_pair_id = {}
        
        self.session = db.Session()
        
        if experiments is None:
            experiments = db.list_experiments(session=self.session)
            # preload all cells,pairs so they are not queried individually later on
            pairs = self.session.query(db.Pair, db.Experiment, db.Cell, db.Slice).join(db.Experiment, db.Pair.experiment_id==db.Experiment.id).join(db.Cell, db.Cell.id==db.Pair.pre_cell_id).join(db.Slice).all()
        
        experiments.sort(key=lambda e: e.acq_timestamp)
        for expt in experiments:
            date = expt.acq_timestamp
            date_str = datetime.fromtimestamp(date).strftime('%Y-%m-%d')
            slice = expt.slice
            expt_item = pg.TreeWidgetItem(map(str, [date_str, '%0.3f'%expt.acq_timestamp, expt.rig_name, slice.species, expt.target_region, slice.genotype, expt.acsf]))
            expt_item.expt = expt
            self.addTopLevelItem(expt_item)

            for pair in expt.pair_list:
                if pair.n_ex_test_spikes == 0 and pair.n_in_test_spikes == 0:
                    continue
                cells = '%d => %d' % (pair.pre_cell.ext_id, pair.post_cell.ext_id)
                conn = {True:"syn", False:"-", None:"?"}[pair.synapse]
                types = 'L%s %s => L%s %s' % (pair.pre_cell.target_layer or "?", pair.pre_cell.cre_type, pair.post_cell.target_layer or "?", pair.post_cell.cre_type)
                pair_item = pg.TreeWidgetItem([cells, conn, types])
                expt_item.addChild(pair_item)
                pair_item.pair = pair
                pair_item.expt = expt
                self.items_by_pair_id[pair.id] = pair_item
                # also allow select by ext id
                self.items_by_pair_id[(expt.acq_timestamp, pair.pre_cell.ext_id, pair.post_cell.ext_id)] = pair_item
                
        self.verticalScrollBar().setValue(self.verticalScrollBar().maximum())
                
    def select_pair(self, pair_id):
        """Select a specific pair from the list
        """
        if self._last_expanded is not None:
            self._last_expanded.setExpanded(False)
        item = self.items_by_pair_id[pair_id]
        self.clearSelection()
        item.setSelected(True)
        parent = item.parent()
        if not parent.isExpanded():
            parent.setExpanded(True)
            self._last_expanded = parent
        else:
            self._last_expanded = None
        self.scrollToItem(item)
