# -*- coding: utf-8 -*-

"""
Prototype code for analyzing connectivity and synaptic properties between cell classes.


"""

from __future__ import print_function, division

from collections import OrderedDict
import numpy as np
import pyqtgraph as pg
import multipatch_analysis.database as db
import pandas as pd
import re, cProfile, os, json, sys, copy
from analyzers import ConnectivityAnalyzer, StrengthAnalyzer, DynamicsAnalyzer, get_all_output_fields
from multipatch_analysis import constants
from multipatch_analysis.cell_class import CellClass, classify_cells, classify_pairs
from matrix_display import MatrixDisplay, MatrixWidget
from scatter_plot_display import ScatterPlotTab
from distance_plot_display import DistancePlotTab
from histogram_trace_display import HistogramTab
from pyqtgraph import parametertree as ptree
from pyqtgraph.parametertree import Parameter
from pyqtgraph.widgets.DataFilterWidget import DataFilterParameter

class SignalHandler(pg.QtCore.QObject):
        """Because we can't subclass from both QObject and QGraphicsRectItem at the same time
        """
        sigOutputChanged = pg.QtCore.Signal(object) #self

class MainWindow(pg.QtGui.QWidget):
    def __init__(self):
        pg.QtGui.QWidget.__init__(self)
        self.layout = pg.QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.h_splitter = pg.QtGui.QSplitter()
        self.h_splitter.setOrientation(pg.QtCore.Qt.Horizontal)
        self.layout.addWidget(self.h_splitter, 0, 0)
        self.control_panel_splitter = pg.QtGui.QSplitter()
        self.control_panel_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.h_splitter.addWidget(self.control_panel_splitter)
        self.update_button = pg.QtGui.QPushButton("Update Results")
        self.control_panel_splitter.addWidget(self.update_button)
        self.ptree = ptree.ParameterTree(showHeader=False)
        self.control_panel_splitter.addWidget(self.ptree)
        self.matrix_widget = MatrixWidget()
        self.h_splitter.addWidget(self.matrix_widget)
        self.tabs = Tabs()
        self.h_splitter.addWidget(self.tabs)
        self.h_splitter.setSizes([300, 600, 400])        

class Tabs(pg.QtGui.QTabWidget):
    def __init__(self, parent=None):
        pg.QtGui.QTabWidget.__init__(self)

        self.hist_tab = HistogramTab()
        self.addTab(self.hist_tab, 'Histogram and Traces')
        self.scatter_tab = ScatterPlotTab()
        self.addTab(self.scatter_tab, 'Scatter Plots')
        self.distance_tab = DistancePlotTab()
        self.addTab(self.distance_tab, 'Distance Plots')


class ExperimentFilter(object):
    def __init__(self):  
        s = db.Session()
        self._signalHandler = SignalHandler()
        self.sigOutputChanged = self._signalHandler.sigOutputChanged
        self.pairs = None
        self.acsf = None
        projects = s.query(db.Experiment.project_name).distinct().all()
        project_list = [{'name': str(record[0]), 'type': 'bool'} for record in projects]
        acsf = s.query(db.Experiment.acsf).distinct().all()
        acsf_list = [{'name': str(record[0]), 'type': 'bool'} for record in acsf]
        internal = s.query(db.Experiment.internal).distinct().all()
        internal_list = [{'name': str(record[0]), 'type': 'bool'} for record in internal]
        self.params = Parameter.create(name='Data Filters', type='group', children=[
            {'name': 'Projects', 'type': 'group', 'children':project_list},
            {'name': 'ACSF', 'type': 'group', 'children':acsf_list, 'expanded': False},
            {'name': 'Internal', 'type': 'group', 'children': internal_list, 'expanded': False},
        ])
        self.params.sigTreeStateChanged.connect(self.invalidate_output)

    def get_pair_list(self, session):
        """ Given a set of user selected experiment filters, return a list of pairs.
        Internally uses multipatch_analysis.db.pair_query.
        """
        if self.pairs is None:
            project_names = [child.name() for child in self.params.child('Projects').children() if child.value() is True]
            project_names = project_names if len(project_names) > 0 else None 
            acsf_recipes = [child.name() for child in self.params.child('ACSF').children() if child.value() is True]
            acsf_recipes = acsf_recipes if len(acsf_recipes) > 0 else None 
            internal_recipes = [child.name() for child in self.params.child('Internal').children() if child.value() is True]
            internal_recipes = internal_recipes if len(internal_recipes) > 0 else None 
            self.pairs = db.pair_query(project_name=project_names, acsf=acsf_recipes, session=session, internal=internal_recipes).all()
        return self.pairs

    def invalidate_output(self):
        self.pairs = None
        self.sigOutputChanged.emit(self)


class CellClassFilter(object):
    def __init__(self, cell_class_groups):
        self.cell_groups = None
        self.cell_classes = None
        self._signalHandler = SignalHandler()
        self.sigOutputChanged = self._signalHandler.sigOutputChanged
        self.experiment_filter = ExperimentFilter() 

        self.cell_class_groups = cell_class_groups
        combo_def = {'name': 'pre/post', 'type':'list', 'value':'both', 'values':['both', 'presynaptic', 'postsynaptic']}
        cell_group_list = [{'name': group, 'type': 'bool', 'children':[combo_def], 'expanded':False} for group in self.cell_class_groups.keys()]
        layer = [{'name': 'Define layer by:', 'type': 'list', 'values': ['target layer', 'annotated layer'], 'value': 'target layer'}]
        children = layer + cell_group_list
        self.params = Parameter.create(name="Cell Classes", type="group", children=children)
        for p in self.params.children():
            p.sigValueChanged.connect(self.expand_param)

        self.params.sigTreeStateChanged.connect(self.invalidate_output)

    def get_cell_groups(self, pairs):
        """Given a list of cell pairs, return a dict indicating which cells
        are members of each user selected cell class.
        This internally calls cell_class.classify_cells
        """
        ccg = copy.deepcopy(self.cell_class_groups)
        if self.cell_groups is None:
            self.cell_classes = []
            for group in self.params.children()[1:]:
                if group.value() is True:
                    self.cell_classes.extend(ccg[group.name()])
            self.cell_classes = self.layer_call() 
            self.cell_classes = [CellClass(**c) for c in self.cell_classes]
            self.cell_groups = classify_cells(self.cell_classes, pairs=pairs)
        return self.cell_groups, self.cell_classes

    def layer_call(self):
        layer_def = self.params['Define layer by:']
        if layer_def == 'target layer':
            for c in self.cell_classes:
                if c.get('cortical_layer') is not None:
                    del c['cortical_layer']
        elif layer_def == 'annotated layer':
            for c in self.cell_classes:
                if c.get('target_layer') is not None:
                    del c['target_layer']
        return self.cell_classes

    def get_pre_or_post_classes(self, key):
        """Return a list of postsynaptic cell_classes. This will be a subset of self.cell_classes."""
        if self.cell_classes is None:
            return []
        classes = []
        for group in self.params.children():
            if group.value() is True:
                if group['pre/post'] in ['both', key]:
                    classes.extend(self.cell_class_groups[group.name()])
        classes = [CellClass(**c) for c in classes]
        return classes

    def expand_param(self, param, value):
        if value:
            param.items.keys()[0].setExpanded(value)

    def invalidate_output(self):
        self.cell_groups = None
        self.cell_classes = None
        self.sigOutputChanged.emit(self)


class MatrixAnalyzer(object):
    sigClicked = pg.QtCore.Signal(object, object, object, object, object) # self, matrix_item, row, col
    def __init__(self, session, cell_class_groups=None, default_preset=None, preset_file=None):
        
        self.main_window = MainWindow()
        self.main_window.setGeometry(280, 130, 1500, 900)
        self.main_window.setWindowTitle('MatrixAnalyzer')
        self.main_window.show()
        self.tabs = self.main_window.tabs
        self.hist_tab = self.tabs.hist_tab
        self.scatter_tab = self.tabs.scatter_tab
        self.distance_tab = self.tabs.distance_tab
        self.distance_plot = self.distance_tab.distance_plot
        self.hist_plot = self.hist_tab.hist
        self.trace_panel = self.hist_tab.trace_plot
        self.selected = 0
        self.colors = [(0, 255, 0), (255, 0, 0), (0, 0, 255), (254, 169, 0), (170, 0, 127), (0, 230, 230)]
        
        self.analyzers = [ConnectivityAnalyzer(), StrengthAnalyzer(), DynamicsAnalyzer()]
        self.active_analyzers = []
        self.preset_file = preset_file

        self.field_map = {}
        for analyzer in self.analyzers:
            for field in analyzer.output_fields():
                if field[0] == 'None':
                    continue
                self.field_map[field[0]] = analyzer

        self.output_fields = get_all_output_fields(self.analyzers)
        self.element_scatter = self.scatter_tab.element_scatter
        self.pair_scatter = self.scatter_tab.pair_scatter
        self.experiment_filter = ExperimentFilter()
        self.cell_class_groups = cell_class_groups
        self.cell_class_filter = CellClassFilter(self.cell_class_groups)
        self.matrix_display = MatrixDisplay(self.main_window, self.output_fields, self.field_map)
        self.matrix_display_filter = self.matrix_display.matrix_display_filter
        self.element_scatter.set_fields(self.output_fields, self.field_map)
        pair_fields = [f for f in self.output_fields if f[0] not in ['connection_probability', 'gap_junction_probability', 'matrix_completeness']]
        self.pair_scatter.set_fields(pair_fields, self.field_map)
        self.visualizers = [self.matrix_display_filter, self.hist_plot, self.element_scatter, self.pair_scatter, self.distance_plot]

        self.default_preset = default_preset
        self.session = session
        self.cell_groups = None
        self.cell_classes = None

        self.presets = self.analyzer_presets()
        preset_list = [p for p in self.presets.keys()]
        self.preset_params = Parameter.create(name='Presets', type='group', children=[
            {'name': 'Analyzer Presets', 'type': 'list', 'values': preset_list, 'value': 'None'},
            {'name': 'Save as Preset', 'type': 'action', 'expanded': True, 'children': [
                {'name': 'Preset Name', 'type': 'text', 'expanded': False}],
            }])
        self.params = ptree.Parameter.create(name='params', type='group', children=[
            self.preset_params,
            self.experiment_filter.params, 
            self.cell_class_filter.params,
            self.matrix_display_filter.params,
        ])

        self.main_window.ptree.setParameters(self.params, showTop=False)
        self.params.child('Presets', 'Save as Preset').sigActivated.connect(self.save_preset)
        
        self.main_window.update_button.clicked.connect(self.update_clicked)
        self.matrix_display.matrix_widget.sigClicked.connect(self.display_matrix_element_data)
        
        self.experiment_filter.sigOutputChanged.connect(self.cell_class_filter.invalidate_output)
        self.params.child('Presets', 'Analyzer Presets').sigValueChanged.connect(self.presetChanged)

        # connect up analyzers
        for analyzer in self.analyzers:
            for visualizer in self.visualizers:
                analyzer.sigOutputChanged.connect(visualizer.invalidate_output)
            self.cell_class_filter.sigOutputChanged.connect(analyzer.invalidate_output)

    def save_preset(self):
        name = self.params['Presets', 'Save as Preset', 'Preset Name']
        if not name:
            msg = pg.QtGui.QMessageBox.warning(self.main_window, "Preset Name", "Please enter a name for your preset in the drop down under 'Preset Name'.",
            pg.QtGui.QMessageBox.Ok)
        else:
            self.presets = self.load_presets()
            cm_state = self.colormap_to_json()

            new_preset = {name: {
                'data filters': self.params.child('Data Filters').saveState(filter='user'),
                'cell classes': self.params.child('Cell Classes').saveState(filter='user'),
                'matrix colormap': cm_state,
                'text format': self.params.child('Matrix Display', 'Text format').saveState(filter='user'),
                'show confidence': self.params.child('Matrix Display', 'Show Confidence').saveState(filter='user'),
                'log scale': self.params.child('Matrix Display', 'log_scale').saveState(filter='user')
                }}
            self.presets.update(new_preset)
            self.write_presets()
            new_preset_list = [p for p in self.presets.keys()]
            self.params.child('Presets', 'Analyzer Presets').sigValueChanged.disconnect(self.presetChanged)
            self.params.child('Presets', 'Analyzer Presets').setLimits(new_preset_list)
            self.params.child('Presets', 'Analyzer Presets').sigValueChanged.connect(self.presetChanged)
            self.presets = self.analyzer_presets()

    def colormap_to_json(self):
        cm_state = self.params.child('Matrix Display', 'Color Map').saveState()
        # colormaps cannot be stored in a JSON, convert format
        cm_state.pop('fields')
        for item, state in cm_state['items'].items():
            cm = state['value']
            cm_state['items'][item]['value'] = {'pos': cm.pos.tolist(), 'color': cm.color.tolist()}
        return cm_state

    def load_presets(self):
        if os.path.exists(self.preset_file):
            try:
                with open(self.preset_file) as json_file:
                    loaded_presets = json.load(json_file, object_pairs_hook=OrderedDict)
            except:
                loaded_presets = {}
                sys.excepthook(*sys.exc_info())
                print ('Error loading analyzer presets')
        else:
            loaded_presets = {}

        return loaded_presets

    def write_presets(self):
        json.dump(self.presets, open(self.preset_file + '.new', 'wb'), indent=4)
        if os.path.exists(self.preset_file):
            os.remove(self.preset_file)
        os.rename(self.preset_file + '.new', self.preset_file)

    def analyzer_presets(self):
        self.presets = self.load_presets()
        self.json_to_colormap()

        return self.presets

    def json_to_colormap(self):
        for name, preset in self.presets.items():
            cm_state = preset['matrix colormap']
            for item, state in cm_state['items'].items():
                cm = cm_state['items'][item]['value']
                self.presets[name]['matrix colormap']['items'][item]['value'] = pg.ColorMap(np.array(cm['pos']), np.array(cm['color']))

    def presetChanged(self):
        self.clear_preset_selections()
        selected = self.params['Presets', 'Analyzer Presets']
        self.set_preset_selections(selected)
        
    def clear_preset_selections(self):
        
        for field in self.experiment_filter.params.children():
            for item in field.children():
                item.setValue(False)

        [item.setValue(False) for item in self.cell_class_filter.params.children()]

        self.matrix_display_filter.params.child('Color Map').clearChildren()
        self.matrix_display_filter.params.child('Text format').setValue('')
        self.matrix_display_filter.params.child('Show Confidence').setValue('None')
        self.matrix_display_filter.params.child('log_scale').setValue(False)

    def set_preset_selections(self, selected):
        preset_state = self.presets[selected]
        self.params.child('Data Filters').restoreState(preset_state['data filters'])
        self.params.child('Cell Classes').restoreState(preset_state['cell classes'])
        self.params.child('Matrix Display', 'Color Map').restoreState(preset_state['matrix colormap'])
        self.params.child('Matrix Display', 'Text format').restoreState(preset_state['text format'])
        self.params.child('Matrix Display', 'Show Confidence').restoreState(preset_state['show confidence'])
        self.params.child('Matrix Display', 'log_scale').restoreState(preset_state['log scale'])
    
    def analyzers_needed(self):
        ## go through all of the visualizers
        data_needed = set(['connected', 'distance'])
        for metric in self.matrix_display_filter.colorMap.children():
            data_needed.add(metric.name())
        text_fields = re.findall('\{(.*?)\}', self.matrix_display_filter.params['Text format'])
        for metric in text_fields:
            if ':' in metric:
                data_needed.add(metric.split(':')[0])
            elif '.' in metric:
                data_needed.add(metric.split('.')[0])
            else:
                data_needed.add(metric)
        data_needed.add(self.matrix_display_filter.params['Show Confidence'])
        for metric in self.element_scatter.fieldList.selectedItems():
            data_needed.add(str(metric.text()))
        for metric in self.element_scatter.filter.children():
            data_needed.add(metric.name())
        for metric in self.element_scatter.colorMap.children():
            data_needed.add(metric.name())
        for metric in self.pair_scatter.fieldList.selectedItems():
            data_needed.add(str(metric.text()))
        for metric in self.pair_scatter.filter.children():
            data_needed.add(metric.name())
        for metric in self.pair_scatter.colorMap.children():
            data_needed.add(metric.name())
        
        analyzers = set([self.field_map.get(field, None) for field in data_needed])
        self.active_analyzers = [analyzer for analyzer in analyzers if analyzer is not None]

        self.data_needed = data_needed

        print ('Active analyzers:')
        print (self.active_analyzers)

        return self.active_analyzers

    def display_matrix_element_data(self, matrix_item, event, row, col):
        field_name = self.matrix_display.matrix_display_filter.get_colormap_field()
        pre_class, post_class = [k for k, v in self.matrix_display.matrix_map.items() if v==[row, col]][0]

        #element = self.results.groupby(['pre_class', 'post_class']).get_group((pre_class, post_class))
        element = self.results.loc[self.pair_groups[(pre_class, post_class)]]
        if len(element) == 0:
            print ('%s->%s has no data, please select another element' % (pre_class, post_class))
            return
        analyzer = self.field_map[field_name]
        analyzer.print_element_info(pre_class, post_class, element, field_name)
        # from here row and col are tuples (row, pre_class) and (col, post_class) respectively
        row = (row, pre_class)
        col = (col, post_class)
        if int(event.modifiers() & pg.QtCore.Qt.ControlModifier)>0:
            self.selected += 1
            if self.selected >= len(self.colors):
                self.selected = 0
            color = self.colors[self.selected]
            self.matrix_display.color_element(row, col, color)
            self.hist_plot.plot_element_data(element, analyzer, color, self.trace_panel)
            self.distance_plot.element_distance(element, color)
        else:
            self.display_matrix_element_reset() 
            color = self.colors[self.selected]
            self.matrix_display.color_element(row, col, color)
            self.hist_plot.plot_element_data(element, analyzer, color, self.trace_panel)
            self.distance_plot.element_distance(element, color)

    def display_matrix_element_reset(self):
        self.selected = 0
        self.hist_plot.plot_element_reset()
        self.matrix_display.element_color_reset()
        self.distance_plot.element_distance_reset(self.results, color=(128, 128, 128), name='All Connections', suppress_scatter=True)


    def update_clicked(self):
        p = cProfile.Profile()
        p.enable()
        with pg.BusyCursor():
            self.analyzers_needed()
            self.update_results()
            pre_cell_classes = self.cell_class_filter.get_pre_or_post_classes('presynaptic')
            post_cell_classes = self.cell_class_filter.get_pre_or_post_classes('postsynaptic')
            self.matrix_display.update_matrix_display(self.results, self.group_results, self.cell_groups, self.field_map, pre_cell_classes=pre_cell_classes, post_cell_classes=post_cell_classes)
            self.hist_plot.matrix_histogram(self.results, self.group_results, self.matrix_display.matrix_display_filter.colorMap, self.field_map)
            self.element_scatter.set_data(self.group_results)
            self.pair_scatter.set_data(self.results)
            self.dist_plot = self.distance_plot.plot_distance(self.results, color=(128, 128, 128), name='All Connections', suppress_scatter=True)
            if self.main_window.matrix_widget.matrix is not None:
                self.display_matrix_element_reset()
        p.disable()
        # p.print_stats(sort='cumulative')

    def update_results(self):
        # Select pairs 
        self.pairs = self.experiment_filter.get_pair_list(self.session)

        # Group all cells by selected classes
        self.cell_groups, self.cell_classes = self.cell_class_filter.get_cell_groups(self.pairs)
        

        # Group pairs into (pre_class, post_class) groups
        self.pair_groups = classify_pairs(self.pairs, self.cell_groups)

        # analyze matrix elements
        for a, analysis in enumerate(self.active_analyzers):
            results = analysis.measure(self.pair_groups)
            group_results = analysis.group_result(self.pair_groups)
            if a == 0:
                self.results = results
                self.group_results = group_results
            else:
                merge_results = pd.concat([self.results, results], axis=1)
                self.results = merge_results.loc[:, ~merge_results.columns.duplicated(keep='first')]
                merge_group_results = pd.concat([self.group_results, group_results], axis=1)
                self.group_results = merge_group_results.loc[:, ~merge_group_results.columns.duplicated(keep='first')]