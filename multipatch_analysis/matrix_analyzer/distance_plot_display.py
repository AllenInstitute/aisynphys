"""
Distance vs connection probability plots for Matrix Analyzer.

"""

from __future__ import print_function, division
import pyqtgraph as pg
import numpy as np
from pyqtgraph.widgets.ColorMapWidget import ColorMapParameter
from pyqtgraph import parametertree as ptree
from pyqtgraph.parametertree import Parameter
from pyqtgraph import parametertree as ptree
from neuroanalysis.ui.plot_grid import PlotGrid
from multipatch_analysis.ui.graphics import distance_plot 


class DistancePlotTab(pg.QtGui.QWidget):
    def __init__(self):
        pg.QtGui.QWidget.__init__(self)
        self.layout = pg.QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.v_splitter = pg.QtGui.QSplitter()
        self.v_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.distance_plot = DistancePlot()
        self.ptree = ptree.ParameterTree(showHeader=False)
        self.ptree.setParameters(self.distance_plot.params)
        self.layout.addWidget(self.v_splitter)
        self.v_splitter.addWidget(self.ptree)
        self.v_splitter.addWidget(self.distance_plot.grid)
        self.v_splitter.setSizes([1, 5000])   

    
class DistancePlot(object):
    def __init__(self):
        self.grid = PlotGrid()
        self.grid.set_shape(2, 1)
        self.grid.grid.ci.layout.setRowStretchFactor(0, 5)
        self.grid.grid.ci.layout.setRowStretchFactor(1, 10)
        self.plots = (self.grid[1,0], self.grid[0,0])
        self.plots[0].grid = self.grid
        self.plots[0].addLegend()
        self.grid.show()
        self.plots[0].setLabels(bottom=('distance', 'm'), left='connection probability')
        self.plots[0].setXRange(0, 200e-6)
        self.params = Parameter.create(name='Distance binning window', type='float', value=40.e-6, step=10.e-6, suffix='m', siPrefix=True)
        self.element_plot = None
        self.elements = []
        self.element_colors = []
        self.results = None
        self.color = None
        self.name = None

        self.params.sigTreeStateChanged.connect(self.update_plot)

    def plot_distance(self, results, color, name, size=10, suppress_scatter=False):
        """Results needs to be a DataFrame or Series object with 'connected' and 'distance' as columns

        """
        connected = results['connected']
        distance = results['distance'] 
        dist_win = self.params.value()
        if self.results is None:
            self.name = name
            self.color = color
            self.results = results
        if suppress_scatter is True:
        #suppress scatter plot for all results (takes forever to plot)
            plots = list(self.plots)
            plots[1] = None
            self.dist_plot = distance_plot(connected, distance, plots=plots, color=color, name=name, size=size, window=dist_win, spacing=dist_win)
        else:
            self.dist_plot = distance_plot(connected, distance, plots=self.plots, color=color, name=name, size=size, window=dist_win, spacing=dist_win)
        
        return self.dist_plot

    def invalidate_output(self):
        self.grid.clear()

    def element_distance(self, element, color, add_to_list=True):
        if add_to_list is True:
            self.element_colors.append(color)
            self.elements.append(element)
        pre = element['pre_class'][0].name
        post = element['post_class'][0].name
        name = ('%s->%s' % (pre, post))
        self.element_plot = self.plot_distance(element, color=color, name=name, size=15)

    def element_distance_reset(self, results, color, name, suppress_scatter=False):
        self.elements = []
        self.element_colors = []
        self.grid.clear()
        self.dist_plot = self.plot_distance(results, color=color, name=name, size=10, suppress_scatter=suppress_scatter)

    def update_plot(self):
        self.invalidate_output()
        self.plot_distance(self.results, self.color, self.name, suppress_scatter=True)
        if self.element_plot is not None:
            for element, color in zip(self.elements, self.element_colors):
                self.element_distance(element, color, add_to_list=False)



