"""
Scatter plot displays histogram of matrix data as well as element data when element is selected
TSeries plot displays average trace responses from each pair in a selected element
"""

from __future__ import print_function, division

from collections import OrderedDict
import numpy as np
import pandas as pd
import pyqtgraph as pg
from pyqtgraph import parametertree as ptree
from pyqtgraph.parametertree import Parameter
from pyqtgraph.widgets.ColorMapWidget import ColorMapParameter
from aisynphys.ui.graphics import MatrixItem
from .analyzers import results_scatter, FormattableNumber
from ..ui.actions import PairActions


class HistogramTab(pg.QtWidgets.QWidget):
    def __init__(self):
        pg.QtWidgets.QWidget.__init__(self)
        self.layout = pg.QtWidgets.QGridLayout()
        self.setLayout(self.layout)
        self.v_splitter = pg.QtWidgets.QSplitter()
        self.v_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.layout.addWidget(self.v_splitter)
        self.hist = MatrixHistogramPlot()
        self.trace_plot = MatrixTSeriesPlot()
        self.v_splitter.addWidget(self.hist)
        self.v_splitter.addWidget(self.trace_plot)

class MatrixHistogramPlot(pg.GraphicsLayoutWidget):
    def __init__(self):
        pg.GraphicsLayoutWidget.__init__(self)
        self.setRenderHints(self.renderHints() | pg.QtGui.QPainter.Antialiasing)
        self.hist_plot = None
        self.line = None
        self.scatter = None
        self.trace_plots = None
        self.trace_plot_list = []
        self.menu = pg.QtWidgets.QMenu()
        self.pair_actions = PairActions()
        for act in self.pair_actions.actions.values():
            self.menu.addAction(act)

    def matrix_histogram(self, results, group_results, colormap, field_map):
         # plot hist or scatter of data in side panel, if there are multiple colormaps take the first enabled one
        color_map_fields = colormap.children()
        if len(color_map_fields) > 1:
            self.field_name = [field.name() for field in color_map_fields if field['Enabled'] is True][0]
        else:
            self.field_name = color_map_fields[0].name()
        field = colormap.fields[self.field_name]
        if self.hist_plot is not None:
            self.removeItem(self.hist_plot)
        self.hist_plot = self.addPlot()
        if field_map[self.field_name].name == 'connectivity':
            vals = group_results[group_results[self.field_name]['metric_summary'].notnull()][self.field_name]['metric_summary']
            y_label = 'Number of Elements' 
        else:
            vals = results[results[self.field_name].notnull()][self.field_name]
            y_label = 'Number of Pairs'
        neg = vals[vals < 0].min()
        pos = vals[vals > 0].max()
        if pos is np.nan or neg is np.nan:
            bins = np.linspace(vals.min(), vals.max(), 10)
        else:
            span = pos - neg
            pos_bins = int(round(10*(pos/span)))
            neg_bins = 10 - pos_bins
            bins = np.linspace(vals.min(), 0, neg_bins)
            bins = sorted(set(np.append(bins, np.linspace(0, vals.max(), pos_bins), axis=0)))
        y, x = np.histogram(vals, bins=bins)
        self.hist_plot.plot(x, y, stepMode=True, fillLevel=0, brush=(255,255,255,150))
        line = pg.InfiniteLine(0, pen={'color': 'w', 'width': 1, 'style': pg.QtCore.Qt.DotLine}, movable=False)
        self.hist_plot.addItem(line)
        units = field.get('units', '')
        self.hist_plot.setLabels(left=y_label, bottom=(self.field_name, units))

        return self.hist_plot

    def plot_element_data(self, element, analyzer, color, trace_panel):
        self.analyzer = analyzer
        self.trace_panel = trace_panel
        pre_class = element['pre_class'][0].name
        post_class = element['post_class'][0].name
        self.trace_pltA = self.trace_panel.addPlot()
        self.trace_pltB = self.trace_panel.addPlot()
        self.trace_plots = [self.trace_pltA, self.trace_pltB]
        
        if self.analyzer.name == 'dynamics':
            self.trace_pltB.setYLink(self.trace_pltA) 
        self.trace_plot_list.append(self.trace_plots)
        self.line, self.scatter = analyzer.plot_element_data(pre_class, post_class, element, self.field_name, color=color, trace_plt=self.trace_plots, pair_actions=self.pair_actions)
        if self.analyzer.name == 'strength' and self.field_name != 'Latency':
            self.trace_panel.removeItem(self.trace_pltA)
        else:
            self.trace_pltB.setXLink(self.trace_pltA)
            self.trace_pltB.hideAxis('left')
        if len(self.trace_plot_list) > 1:
            first_plot = self.trace_plot_list[0][0]
            for grid in self.trace_plot_list[1:]:
                for plot in grid:
                    plot.setXLink(first_plot)
        self.hist_plot.addItem(self.line)
        if self.scatter is not None:
            self.hist_plot.addItem(self.scatter)

    def plot_element_reset(self):
        self.selected = 0
        if self.hist_plot is not None:
            [self.hist_plot.removeItem(item) for item in self.hist_plot.items[2:]]
        if self.trace_plots is not None:
           self.trace_panel.clear()
           self.trace_plot = None
           self.analyzer.pair_items = {}
        self.trace_plot_list = []

    def invalidate_output(self):
        self.hist_plot.clear()
        self.trace_panel.clear()
        self.hist_plot = None
        self.line = None
        self.scatter = None
        self.trace_plot = None
        self.trace_plot_list = []

    def contextMenuEvent(self, event):
        self.menu.popup(event.globalPos())


class MatrixTSeriesPlot(pg.GraphicsLayoutWidget):
    def __init__(self):
        pg.GraphicsLayoutWidget.__init__(self)
        self.setRenderHints(self.renderHints() | pg.QtGui.QPainter.Antialiasing)

        self.menu = pg.QtWidgets.QMenu()
        self.pair_actions = PairActions()
        for act in self.pair_actions.actions.values():
            self.menu.addAction(act)

        def contextMenuEvent(self, event):
            self.menu.popup(event.globalPos())
