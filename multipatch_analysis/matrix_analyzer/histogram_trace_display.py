"""
Scatter plot displays histogram of matrix data as well as element data when element is selected
Trace plot displays average trace responses from each pair in a selected element
"""

from __future__ import print_function, division


import numpy as np
import pyqtgraph as pg
import pandas as pd
from collections import OrderedDict
from analyzers import results_scatter, FormattableNumber
from multipatch_analysis.ui.graphics import MatrixItem
from pyqtgraph import parametertree as ptree
from pyqtgraph.parametertree import Parameter
from pyqtgraph.widgets.ColorMapWidget import ColorMapParameter

class HistogramTab(pg.QtGui.QWidget):
    def __init__(self):
        pg.QtGui.QWidget.__init__(self)
        self.layout = pg.QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.v_splitter = pg.QtGui.QSplitter()
        self.v_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.layout.addWidget(self.v_splitter)
        self.hist = MatrixHistogramPlot()
        self.trace_plot = MatrixTracePlot()
        self.v_splitter.addWidget(self.hist)
        self.v_splitter.addWidget(self.trace_plot)

class MatrixHistogramPlot(pg.GraphicsLayoutWidget):
    def __init__(self):
        pg.GraphicsLayoutWidget.__init__(self)
        self.setRenderHints(self.renderHints() | pg.QtGui.QPainter.Antialiasing)
        self.hist_plot = None
        self.line = None
        self.scatter = None
        self.trace_plot = None
        self.trace_plot_list = []

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
        else:
            vals = results[results[self.field_name].notnull()][self.field_name]
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
        self.hist_plot.setLabels(left='Count', bottom=(self.field_name, units))

        return self.hist_plot

    def plot_element_data(self, element, analyzer, color, trace_panel):
        self.trace_panel = trace_panel
        pre_class = element['pre_class'][0].name
        post_class = element['post_class'][0].name
        self.trace_plot = self.trace_panel.addPlot()
        self.trace_plot_list.append(self.trace_plot)
        self.line, self.scatter = analyzer.plot_element_data(pre_class, post_class, element, self.field_name, color=color, trace_plt=self.trace_plot)
        if len(self.trace_plot_list) > 1:
            first_plot = self.trace_plot_list[0]
            for plot in self.trace_plot_list[1:]:
                plot.setYLink(first_plot)
                plot.setXLink(first_plot)
        self.hist_plot.addItem(self.line)
        if self.scatter is not None:
            self.hist_plot.addItem(self.scatter)

    def plot_element_reset(self):
        self.selected = 0
        if self.hist_plot is not None:
            [self.hist_plot.removeItem(item) for item in self.hist_plot.items[2:]]
        if self.trace_plot is not None:
           self.trace_panel.clear()
           self.trace_plot = None
        self.trace_plot_list = []

    def invalidate_output(self):
        self.hist_plot.clear()
        self.trace_panel.clear()
        self.hist_plot = None
        self.line = None
        self.scatter = None
        self.trace_plot = None
        self.trace_plot_list = []


class MatrixTracePlot(pg.GraphicsLayoutWidget):
    def __init__(self):
        pg.GraphicsLayoutWidget.__init__(self)
        self.setRenderHints(self.renderHints() | pg.QtGui.QPainter.Antialiasing)
