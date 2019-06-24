from __future__ import print_function
import os, sys, datetime, re, glob, traceback, time, atexit, threading
try:
    import queue
except ImportError:
    import Queue as queue
from pprint import pprint
from collections import OrderedDict
import numpy as np
from .. import config, lims
from ..experiment import Experiment
from .. import database
from ..genotypes import Genotype
from .actions import ExperimentActions
from ..yaml_local import yaml

from acq4.util.DataManager import getDirHandle
import pyqtgraph as pg
import pyqtgraph.console
from pyqtgraph.Qt import QtGui, QtCore


fail_color = (255, 200, 200)
pass_color = (200, 255, 200)


class Dashboard(QtGui.QWidget):
    def __init__(self, limit=0, no_thread=False, filter_defaults=None):
        QtGui.QWidget.__init__(self)

        # fields displayed in ui
        self.visible_fields = [
            ('timestamp', float), 
            ('rig', 'U100'), 
            ('operator', 'U100'), 
            ('path', 'U100'), 
            ('project', 'U100'),
            ('description', 'U100'), 
            ('primary', 'U100'), 
            ('archive', 'U100'), 
            ('NAS', 'U100'), 
            ('backup', 'U100'), 
            ('data', 'U100'),
            ('submitted', 'U100'),
            ('connections', 'U100'),
            ('site.mosaic', 'U100'), 
            ('DB', 'U100'), 
            ('LIMS', 'U100'), 
            ('20x', 'U100'), 
            ('cell map', 'U100'), 
            ('63x', 'U100'), 
            ('morphology', 'U100')
        ]

        # data tracked but not displayed
        self.hidden_fields = [
            ('experiment', object),
            ('lims_slice_name', object),
            ('item', object),
            ('error', object),
            ('db_errors', object),
        ]

        # maps field name : index (column number)
        self.field_indices = {self.visible_fields[i][0]:i for i in range(len(self.visible_fields))}

        self.records = GrowingArray(dtype=self.visible_fields + self.hidden_fields)
        self.records_by_expt = {}  # maps expt:index

        self.selected = None

        # set up UI
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.layout.setContentsMargins(0, 0, 0, 0)

        self.splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.splitter, 0, 0)

        self.left_vbox_widget = QtGui.QWidget()
        self.splitter.addWidget(self.left_vbox_widget)
        self.left_vbox = QtGui.QVBoxLayout()
        self.left_vbox_widget.setLayout(self.left_vbox)
        self.left_vbox.setContentsMargins(0, 0, 0, 0)

        self.search_text = QtGui.QLineEdit()
        self.search_text.setPlaceholderText('search')
        # self.search_text.setClearButtonEnabled(True)   # Qt 5
        self.left_vbox.addWidget(self.search_text)
        self.search_text.textChanged.connect(self.search_text_changed)

        self.filter = pg.DataFilterWidget()
        self.filter_fields = OrderedDict([(field[0], {'mode': 'enum', 'values': {}}) for field in self.visible_fields])
        self.filter_fields['timestamp'] = {'mode': 'range'}
        self.filter_fields.pop('path')
        if filter_defaults is not None:
            for k, defs in filter_defaults.items():
                self.filter_fields[k]['values'].update(defs)

        self.filter.setFields(self.filter_fields)
        self.left_vbox.addWidget(self.filter)
        self.filter.sigFilterChanged.connect(self.filter_changed)
        self.filter.addFilter('operator')
        self.filter.addFilter('rig')
        self.filter.addFilter('project')

        self.right_splitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.splitter.addWidget(self.right_splitter)

        self.expt_tree = pg.TreeWidget()
        self.expt_tree.setSortingEnabled(True)
        self.expt_tree.setColumnCount(len(self.visible_fields))
        self.expt_tree.setHeaderLabels([f[0] for f in self.visible_fields])
        self.right_splitter.addWidget(self.expt_tree)
        self.expt_tree.itemSelectionChanged.connect(self.tree_selection_changed)

        console_text = """
        Variables:
            records : array of all known records in the dashboard
            records['experiment'] : references to all known Experiment instances
            sel : currently selected dashboard record
            expt : currently selected Experiment
            filter : DataFilterParameter
            lims : multipatch_analysis.lims module
            db : multipatch_analysis.database module
            config : multipatch_analysis.config module
        """
        self.console = pg.console.ConsoleWidget(text=console_text, namespace={
            'filter': self.filter.params,
            'sel': None,
            'records': self.records,
            'config': config,
            'lims': lims,
            'db': database,
            'dashboard': self,
            'np': np,
            'pg': pg,
        })
        self.right_splitter.addWidget(self.console)
        self.console.hide()

        self.console_btn = QtGui.QPushButton('console')
        self.console_btn.setCheckable(True)
        self.console_btn.toggled.connect(self.console_toggled)

        self.poll_btn = QtGui.QPushButton('poll')
        self.poll_btn.setCheckable(True)
        self.poll_btn.setChecked(True)
        self.poll_btn.toggled.connect(self.poll_toggled)

        self.status = QtGui.QStatusBar()
        self.status.setMaximumHeight(25)
        self.layout.addWidget(self.status, 1, 0)
        self.status.addPermanentWidget(self.poll_btn)
        self.status.addPermanentWidget(self.console_btn)

        self.resize(1000, 900)
        self.splitter.setSizes([200, 800])
        colsizes = [150, 50, 100, 230, 150, 200]
        for i in range(len(self.visible_fields)):
            size = colsizes[i] if i < len(colsizes) else 50
            self.expt_tree.setColumnWidth(i, size)

        # Add reload context menu action
        # self.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.reload_action = QtGui.QAction('Reload', self)
        self.reload_action.triggered.connect(self.reload_clicked)
        # self.addAction(self.reload_action)

        self.menu = QtGui.QMenu()
        self.menu.addAction(self.reload_action)
        
        # Add generic experiment context menu actions
        self.expt_actions = ExperimentActions()
        for act in self.expt_actions.actions.values():
            self.menu.addAction(act)
        
        # Queue of experiments to be checked
        self.expt_queue = queue.PriorityQueue()

        # collect a list of all data sources to search
        search_paths = [config.synphys_data]
        for rig_name, rig_path_sets in config.rig_data_paths.items():
            for path_set in rig_path_sets:
                search_paths.extend(list(path_set.values()))

        # Each poller searches a single data source for new experiments, adding them to the queue
        self.pollers = []
        self.poller_status = {}
        for search_path in search_paths:
            if not os.path.exists(search_path):
                print("Ignoring search path:", search_path)
                continue
            poll_thread = PollThread(self.expt_queue, search_path, limit=limit)
            poll_thread.update.connect(self.poller_update)
            if no_thread:
                poll_thread.poll()  # for local debugging
            else:
                poll_thread.start()
            self.pollers.append(poll_thread)

        # Checkers pull experiments off of the queue and check their status
        self._incoming_checker_records = []
        if no_thread:
            self.checker = ExptCheckerThread(self.expt_queue)
            self.checker.update.connect(self.checker_update)
            self.checker.run(block=False)            
        else:
            self.checkers = []
            for i in range(3):
                self.checkers.append(ExptCheckerThread(self.expt_queue))
                self.checkers[-1].update.connect(self.checker_update)
                self.checkers[-1].start()

        # shut down threads nicely on exit
        atexit.register(self.quit)

        self.status_timer = QtCore.QTimer()
        self.status_timer.timeout.connect(self.update_status)
        self.status_timer.start(1000)

        self.handle_checker_record_timer = QtCore.QTimer()
        self.handle_checker_record_timer.timeout.connect(self.handle_all_checker_records)

    def poll_toggled(self):
        if self.poll_btn.isChecked():
            self.start_polling()
        else:
            self.stop_polling()

    def start_polling(self):
        for poller in self.pollers:
            poller.enable_polling = True
            poller.waker.set()

    def stop_polling(self):
        for poller in self.pollers:
            poller.enable_polling = False

        # empty out the queue
        while self.expt_queue.qsize() > 0:
            self.expt_queue.get()

    def contextMenuEvent(self, event):
        self.menu.popup(event.globalPos())

    def tree_selection_changed(self):
        sel = self.expt_tree.selectedItems()[0]
        rec = self.records[sel.index]
        self.console.localNamespace['sel'] = rec
        self.selected = rec
        expt = rec['experiment']
        self.console.localNamespace['expt'] = expt
        self.expt_actions.experiment = expt

        msg = [
            "\n======= Selected: ======================",
            "       timestamp: %0.3f" % expt.timestamp,
            "     description: %s" % rec['description'],
        ]
        print_fields = [
            ('spec name', 'specimen_name'),
            ('genotype', 'genotype'),
            ('NAS path', 'nas_path'),
            ('primary path', 'primary_path'),
            ('archive path', 'archive_path'),
            ('backup path', 'backup_path'),
            ('biocytin URL', 'biocytin_image_url'),
            ('drawing tool', 'lims_drawing_tool_url'),
            ('cluster ID', 'cluster_id'),
            ('slice fixed', 'slice_fixed'),
            ('LIMS status', 'lims_message'), 
            ('Cell Map', 'map_message')
        ]
        for name,attr in print_fields:
            try:
                val = getattr(expt, attr)
                msg.append("%16s: %s" % (name, val))
            except Exception:
                msg.append("Error getting experiment attribute '%s':" % attr)
                msg.append(traceback.format_exc())
                break
        err = rec['error']
        if err is not None:
            msg.append("--------------------------------\nError checking experiment:")
            msg.extend([line.rstrip() for line in traceback.format_exception(*err)])

        db_errors = rec['db_errors']
        if isinstance(db_errors, dict) and len(db_errors) > 0:
            msg.append("--------------------------------\nErrors in DB pipeline:")
            for module, err in db_errors.items():
                msg.append("--- %s ---\n%s\n" % (module, err))

        msg = '\n'.join(msg)
        print(msg)
        self.console.write(msg+'\n')

    def poller_update(self, path, status):
        """Received an update on the status of a poller
        """
        self.poller_status[path] = status
        paths = sorted(self.poller_status.keys())
        msg = '    '.join(['%s: %s' % (p,self.poller_status[p]) for p in paths])
        # self.long_status.showMessage(msg)

    def update_status(self):
        self.status.showMessage("%d sites in queue" % self.expt_queue.qsize())

    def checker_update(self, rec):
        """Received an update from a worker thread describing information about an experiment
        """
        # collect incoming records and process on a timer to ensure the ui doesn't get blocked.
        self._incoming_checker_records.append(rec)
        self.handle_checker_record_timer.start(0)

    def handle_all_checker_records(self):
        recs = self._incoming_checker_records
        count = 0
        while len(recs) > 0:
            rec = recs.pop(0)
            self.handle_checker_record(rec)
            count += 1
            if count > 2:
                # yield to the event loop
                return
        self.handle_checker_record_timer.stop()

    def handle_checker_record(self, rec):
        expt = rec['experiment']
        if expt in self.records_by_expt:
            # use old record / item
            index = self.records_by_expt[expt]
            item = self.records[index]['item']
        else:
            # add new record / item
            item = pg.TreeWidgetItem()
            self.expt_tree.addTopLevelItem(item)
            rec['item'] = item
            index = self.records.add_record({})
            item.index = index

        record = self.records[index]
        record['item'] = item
        self.records_by_expt[expt] = index

        # update item/record fields
        update_filter = False
        for field, val in rec.items():
            if field in self.field_indices and isinstance(val, tuple):
                # if a tuple was given, interpret it as (text, color)
                val, color = val
            else:
                # otherwise make a guess on a good color
                color = 'w'
                if val is True:
                    color = pass_color
                elif val is False:
                    color = fail_color
                elif val in ('ERROR', 'MISSING', 'FAILED'):
                    color = fail_color
                elif val == '-':
                    # dash means item is incomplete but still acceptable
                    color = [c * 0.5 + 128 for c in pass_color]
            
            display_val = str(val)

            if field == 'timestamp' and rec.get('error') is not None:
                color = fail_color

            # update this field in the record
            record[field] = val
            
            # update this field in the tree item
            try:
                i = self.field_indices[field]
            except KeyError:
                continue

            item.setText(i, display_val)
            if color is not None:
                item.setBackgroundColor(i, pg.mkColor(color))

            # update filter fields
            filter_field = self.filter_fields.get(field)
            if filter_field is not None and filter_field['mode'] == 'enum' and display_val not in filter_field['values']:
                filter_field['values'][display_val] = True
                update_filter = True

        if update_filter:
            self.filter.setFields(self.filter_fields)

        self.filter_items(self.records[index:index+1])

    def filter_changed(self, *args):
        self.filter_items(self.records)

    def search_text_changed(self):
        self.filter_items(self.records)

    def filter_items(self, records):
        search = str(self.search_text.text()).strip()
        mask = self.filter.generateMask(records)
        for i,rec in enumerate(records):
            hidden = not mask[i]
            if search != '':
                search_hit = False
                if rec['lims_slice_name'] is not None and search in rec['lims_slice_name']:
                    search_hit = True
                elif search in str(rec['timestamp']) or search in str(np.round(rec['timestamp'], 2)):
                    search_hit = True
                elif rec['description'] is not None and search in rec['description']:
                    search_hit = True
                elif rec['path'] is not None and search in rec['path']:
                    search_hit = True
                hidden = hidden or not search_hit
            rec['item'].setHidden(hidden)

    def quit(self):
        for t in self.pollers:
            if t.isRunning():
                t.stop()
                t.wait()

        stopped = []
        for t in self.checkers:
            if t.isRunning():
                t.stop()
                stopped.append(t)
        for t in stopped:
            t.wait()  # don't start waiting until all checkers have been requested to stop!

    def closeEvent(self, ev):
        self.quit()

    def reload_clicked(self, *args):
        expt = self.selected['experiment']
        self.expt_queue.put((-expt.timestamp, expt))

    def console_toggled(self):
        self.console.setVisible(self.console_btn.isChecked())


class GrowingArray(object):
    def __init__(self, dtype, init_size=1000):
        self.size = 0
        self._data = np.empty(init_size, dtype=dtype)
        self._view = self._data[:0]

    def __len__(self):
        return self.size

    @property
    def shape(self):
        return (self.size,)

    def __getitem__(self, i):
        return self._view[i]

    def add_record(self, rec):
        index = self.size
        self._grow(self.size+1)
        self.update_record(index, rec)
        return index

    def update_record(self, index, rec):
        for k,v in rec.items():
            self._data[index][k] = v

    def _grow(self, size):
        if size > len(self._data):
            self._data = np.resize(self._data, len(self._data)*2)
        self._view = self._data[:size]
        self.size = size


class PollThread(QtCore.QThread):
    """Used to check in the background for changes to experiment status.
    """
    update = QtCore.Signal(object, object)  # search_path, status_message
    known_expts = {}
    known_expts_lock = threading.Lock()

    def __init__(self, expt_queue, search_path, limit=0, interval=3600):
        QtCore.QThread.__init__(self)
        self.expt_queue = expt_queue
        self.search_path = search_path
        self.limit = limit
        self.interval = interval
        self._stop = False
        self.waker = threading.Event()
        self.enable_polling = True   # set False to temporarily disable polling
        
    def stop(self):
        self._stop = True
        self.waker.set()

    def run(self):
        while True:
            try:
                # check for new experiments hourly
                self.poll()
                self.waker.wait(self.interval)
                if self._stop:
                    return
            except Exception:
                sys.excepthook(*sys.exc_info())
                
    def poll(self):
        # Find all available site paths across all data sources
        count = 0
        path = self.search_path

        self.update.emit(path, "Updating...")
        root_dh = getDirHandle(path)

        # iterate over all expt sites in this path
        for day_name in sorted(os.listdir(root_dh.name()), reverse=True):
            for expt_path in glob.iglob(os.path.join(root_dh.name(), day_name, 'slice_*', 'site_*')):
                if self._stop or not self.enable_polling:
                    return

                try:
                    expt = ExperimentMetadata(path=expt_path)
                    ts = expt.timestamp
                except:
                    print ('Error loading %s, ignoring and moving on...' % expt_path)    
                    continue
                # Couldn't get timestamp; show an error message
                if ts is None:
                    print("Error getting timestamp for %s" % expt)
                    continue

                with self.known_expts_lock:
                    now = time.time()
                    if ts in self.known_expts:
                        expt, last_update = self.known_expts[ts]
                        if now - last_update < self.interval:
                            # We've already seen this expt recently; skip
                            continue
                    self.known_expts[ts] = (expt, now)

                # Add this expt to the queue to be checked
                self.expt_queue.put((-ts, expt))

                count += 1
                if self.limit > 0 and count >= self.limit:
                    return
        self.update.emit(path, "Finished")


class ExptCheckerThread(QtCore.QThread):
    update = QtCore.Signal(object)

    def __init__(self, expt_queue):
        QtCore.QThread.__init__(self)
        self.expt_queue = expt_queue
        self._stop = False

    def stop(self):
        self._stop = True
        self.expt_queue.put(('stop', None))

    def run(self, block=True):
        while True:
            ts, expt = self.expt_queue.get(block=block)
            if self._stop or ts == 'stop':
                return
            rec = expt.check()
            self.update.emit(rec)


class ExperimentMetadata(Experiment):
    """Handles reading experiment metadata from several possible locations.
    """
    def __init__(self, path=None):
        Experiment.__init__(self, verify=False)
        self._site_path = path

        self.site_dh = getDirHandle(path)
        self._rig_name = None
        self._primary_path = None
        self._archive_path = None
        self._backup_path = None

        # reassign path based on order of most likely to be updated
        for path_typ in ['primary', 'archive', 'nas', 'backup']:
            path = getattr(self, path_typ + '_path')
            if path is not None and os.path.exists(path):
                self._site_path = path
                self.site_dh = getDirHandle(path)
                # reset values loaded while determining path
                self._expt_info = None
                self._slice_info = None
                self._site_info = None
                break

    def check(self):
        """Check the status of this experiment, return a dict describing
        the state of each stage in the pipeline.
        """
        try:
            rec = {'experiment': self, 'error': None}

            rec['timestamp'] = '%0.3f' % self.timestamp
            rec['rig'] = self.rig_name
            rec['operator'] = self.rig_operator
            rec['path'] = self.expt_subpath
            rec['project'] = self.slice_info.get('project', None)
            rec['primary'] = False if self.primary_path is None else (True if os.path.exists(self.primary_path) else "-")
            rec['archive'] = False if self.archive_path is None else (True if os.path.exists(self.archive_path) else "MISSING")
            rec['backup'] = False if self.backup_path is None else (True if os.path.exists(self.backup_path) else "MISSING")
            rec['NAS'] = False if self.nas_path is None else (True if os.path.exists(self.nas_path) else "MISSING")

            self.lims_message = None
            self.map_message = None
            org = self.organism
            if org is None:
                description = ("no LIMS spec info", fail_color)
            elif org == 'human':
                description = org
            else:
                try:
                    gtyp = self.genotype
                    description = ','.join(gtyp.all_drivers)
                except Exception:
                    description = (org + ' (unknown: %s)'%self.lims_record.get('genotype', None), fail_color)
            rec['description'] = description

            rec['lims_slice_name'] = self.specimen_name
            
            search_paths = [self.primary_path, self.archive_path, self.nas_path]
            submitted = False
            connections = None
            for path in search_paths:
                if path is None:
                    continue
                yml_file = os.path.join(path, 'pipettes.yml')
                if os.path.isfile(yml_file):
                    submitted = True
                    connections = (0, pass_color)
                    pips = yaml.load(open(yml_file, 'rb'))
                    for k,v in pips.items():
                        if v['got_data']:
                            if v['synapse_to'] is None:
                                connections = False
                                break
                            else:
                                connections = (connections[0] + len(v['synapse_to']), pass_color)
                    break


            rec['submitted'] = submitted
            rec['data'] = '-' if self.nwb_file is None else True
            slice_fixed = self.slice_fixed
            if slice_fixed is True:
                image_20x = self.biocytin_20x_file
                rec['20x'] = image_20x is not None
            else:
                rec['20x'] = '-'

            if rec['submitted']:
                rec['connections'] = connections
                if rec['data'] is True:
                    rec['site.mosaic'] = self.mosaic_file is not None
                    has_run, expt_success, all_success, errors = self.db_status
                    if not has_run:
                        rec['DB'] = ("MISSING", (255, 255, 100))
                    elif not expt_success:
                        rec['DB'] = ("ERROR", fail_color)
                    elif not all_success:
                        rec['DB'] = ("ERROR", (255, 255, 100))
                    else:
                        rec['DB'] = True
                    rec['db_errors'] = errors

                    cell_cluster = self.lims_cell_cluster_id
                    lims_ignore_path = os.path.join(self.archive_path, '.mpe_ignore')
                    if cell_cluster is None and os.path.isfile(lims_ignore_path):
                        lims_ignore_file = os.path.join(self.archive_path, '.mpe_ignore')
                        in_lims = "FAILED"
                        self.lims_message = open(lims_ignore_file, 'r').read()
                    else:
                        in_lims = cell_cluster is not None
                        self.lims_message = self.lims_submissions
 
                    rec['LIMS'] = in_lims
                    image_63x = self.biocytin_63x_files

                    if in_lims is True and slice_fixed is True and image_20x is not None:
                        cell_specimens = lims.child_specimens(cell_cluster)
                        if len(cell_specimens) != 0:
                            cell_info = lims.cluster_cells(cell_cluster)
                            x_coord = all([ci['x_coord'] is not None for ci in cell_info])
                            x_coord_values = all([ci['x_coord'] != 0 for ci in cell_info])
                            y_coord = all([ci['y_coord'] is not None  for ci in cell_info])
                            y_coord_values = all([ci['y_coord'] != 0 for ci in cell_info])
                            polygons = [lims.cell_polygon(cell.id) for cell in cell_info if cell is not None]
                            polygons = all(polygon is not None for polygon in polygons)
                            mapped = len(cell_info) > 0 and x_coord is True and y_coord is True
                            if mapped is True and x_coord_values is False and y_coord_values is False:
                                mapped = ('Incomplete', (255, 255, 102))
                                self.map_message = 'Cell positions contain 0-values, please re-map'
                            if mapped is True and polygons is False:
                                mapped = ('Incomplete', [c * 0.5 + 128 for c in pass_color])
                                self.map_message = 'Cell positions submitted to LIMS but polygons have not been drawn yet'
                            rec['cell map'] = mapped
                        else:
                            rec['cell map'] = False
                        image_tags = lims.specimen_tags(cell_cluster)
                        if image_tags is not None:
                            if 'cell map no go' in image_tags:
                                rec['cell map'] = '-'
                                self.map_message = 'Site marked unmappable by user'  
                            if image_63x is None:
                                if '63x no go' in image_tags:
                                    rec['63x'] = '-'
                                    rec['cell map'] = '-'
                                if '63x go' in image_tags:
                                    rec['63x'] = ('63x GO', [c * 0.5 + 128 for c in pass_color])
                                if '63x imaging started' in image_tags:
                                    rec['63x'] = ('Imaging started...', [c * 0.5 + 128 for c in pass_color])
                            else:        
                                rec['63x'] = image_63x is not None
            else:
                if self.mosaic_file is not None:
                    rec['site.mosaic'] = True
        
        except Exception as exc:
            rec['error'] = sys.exc_info()
        
        return rec

    def _get_raw_paths(self):
        expt_subpath = self.expt_subpath

        # find the local primary/archive paths that contain this experiment
        found_paths = False
        rig_name = self.rig_name
        if rig_name is None:
            return
        rig_data_paths = config.rig_data_paths.get(rig_name, [])
        for path_set in rig_data_paths:
            for root in path_set.values():
                test_path = os.path.join(root, expt_subpath)
                if not os.path.isdir(test_path):
                    continue
                dh = getDirHandle(test_path)
                if self.site_info is None:
                    raise Exception ('%s %s missing index file' % (self, self.path))
                if dh.info()['__timestamp__'] == self.site_info['__timestamp__']:
                    found_paths = True
                    # set self._primary_path, self._archive_path, etc.
                    for k,v in path_set.items():
                        setattr(self, '_'+k+'_path', os.path.join(v, expt_subpath))
                    break
            if found_paths:
                break

    @property
    def expt_subpath(self):
        path = self.site_dh.name()

        # get raw data subpath  (eg: 2018.01.20_000/slice_000/site_000)
        if os.path.abspath(path).startswith(os.path.abspath(config.synphys_data).rstrip(os.path.sep) + os.path.sep):
            # this is a server path; need to back-convert to rig path
            source_path = open(os.path.join(path, 'sync_source')).read()
            expt_subpath = os.path.join(*source_path.split('/')[-3:])
            assert os.path.abspath(path) == os.path.abspath(self.nas_path), "Expected equal paths:\n  %s\n  %s" % (os.path.abspath(path), os.path.abspath(self.nas_path))
        else:
            expt_subpath = os.path.join(*path.split(os.path.sep)[-3:])

        return expt_subpath

    @property
    def primary_path(self):
        if self._primary_path is None:
            self._get_raw_paths()
        return self._primary_path

    @property
    def archive_path(self):
        if self._archive_path is None:
            self._get_raw_paths()
        return self._archive_path

    @property
    def backup_path(self):
        if self._backup_path is None:
            self._get_raw_paths()
        return self._backup_path

    @property
    def slice_dh(self):
        return self.site_dh.parent()

    @property
    def expt_dh(self):
        return self.slice_dh.parent()

    @property
    def nas_path(self):
        return os.path.abspath(os.path.join(config.synphys_data, self.server_path))

    @property
    def organism(self):
        org = self.expt_info.get('organism')
        if org is not None:
            return org
        spec_info = self.lims_record
        if spec_info is None:
            return None
        return spec_info['organism']

    @property
    def rig_name(self):
        if self._rig_name is None:
            name = self.expt_info.get('rig_name')
            if name is not None:
                self._rig_name = name.lower()
            else:
                # infer rig name from paths
                nas_path = self.nas_path
                if nas_path is not None:
                    sync_file = os.path.join(self.nas_path, 'sync_source')
                else:
                    sync_file = None
                if sync_file is not None and os.path.isfile(sync_file):
                    path = open(sync_file, 'rb').read()
                else:
                    path = self.site_dh.name()
                m = re.search(r'(/|\\)(mp\d)(/|\\)', path)
                if m is None:
                    # last resort: see if we can find the experiment path in the config list of rig data locations
                    path = os.path.abspath(self.site_dh.name())
                    for rig,sources in config.rig_data_paths.items():
                        for data_source in sources:
                            for data_path in data_source.values():
                                if path.startswith(os.path.abspath(data_path)):
                                    self._rig_name = rig
                                    break
                        if self._rig_name is not None:
                            break
                    
                    return None
                self._rig_name = m.groups()[1].lower()
        return self._rig_name

    @property
    def db_status(self):
        jobs = database.query(database.Pipeline).filter(database.Pipeline.job_id==self.timestamp).all()
        has_run = len(jobs) > 0
        success = {j.module_name:j.success for j in jobs}
        errors = {j.module_name:j.error for j in jobs}
        expt_success = success.get('experiment', False)
        all_success = all(success.values())
        return has_run, expt_success, all_success, errors

    @property
    def lims_submissions(self):
        if self.lims_record is None:
            return None
        spec_id = self.lims_record['specimen_id']
        if spec_id is None:
            return None
        return lims.expt_submissions(spec_id, self.timestamp)

    @property
    def lims_cell_cluster_id(self):
        subs = self.lims_submissions
        if subs is None or len(subs) != 1:
            return None
        return subs[0][1]

    @property
    def slice_fixed(self):
        return self.slice_info.get('plate_well_ID') != 'not fixed'
