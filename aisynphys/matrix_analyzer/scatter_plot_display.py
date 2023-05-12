"""
ScatterPlotWidgets for Matrix Analyzer. 
Allows plotting of Matrix data on a per element basis and per pair basis.

"""

from __future__ import print_function, division
import pyqtgraph as pg
import pandas as pd
import numpy as np
from ..ui.actions import PairActions
from aisynphys.ui.ScatterPlotWidget import ScatterPlotWidget


class ScatterPlotTab(pg.QtWidgets.QWidget):
    def __init__(self):
        pg.QtWidgets.QWidget.__init__(self)
        self.layout = pg.QtWidgets.QGridLayout()
        self.setLayout(self.layout)
        self.v_splitter = pg.QtWidgets.QSplitter()
        self.v_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.layout.addWidget(self.v_splitter)
        self.element_scatter = ElementScatterPlot()
        self.pair_scatter = PairScatterPlot()#self.element_scatter.plot)
        self.v_splitter.addWidget(self.element_scatter)
        self.v_splitter.addWidget(self.pair_scatter)

class ScatterPlots(ScatterPlotWidget):
    def __init__(self):
        ScatterPlotWidget.__init__(self)

    def set_fields(self, fields):
        self.fields = [('pair_class', {'mode': 'enum'})]
        self.fields.extend([f for f in fields if f != ('None', {})])
        self.setFields(self.fields)

    def color_selected_element(self, color, pre_class, post_class):
        pair_name = '-'.join([pre_class.name, post_class.name])
        try:
            pair_map = self.colorMap.child('pair_class')
        except KeyError:
            pair_map = self.colorMap.addNew('pair_class')
        pair_map['Values', pair_name] = pg.mkColor(color)

    def reset_element_color(self):
        try:
            pair_class_map = self.colorMap.child('pair_class')
            self.colorMap.removeChild(pair_class_map)
        except:
            return

    def invalidate_output(self):
        self.data = None

class ElementScatterPlot(ScatterPlots):
    def __init__(self):
        ScatterPlots.__init__(self)
    
        header = pg.QtWidgets.QLabel()
        header.setText('<span style="font-weight: bold">Element-wise Scatter Plot</span>')
        self.ctrlPanel.insertWidget(0, header)
        self.selected_points = []

    def set_data(self, data):
        field_data = data.xs('metric_summary', axis='columns', level=1, drop_level=True)
        field_data.reset_index(inplace=True)
        field_data['pair_class'] = field_data.apply(lambda row: '-'.join([row.pre_class.name, row.post_class.name]), axis=1)
        field_data.drop(columns=['pre_class', 'post_class'])
        rec_data = field_data.to_records()
        names = tuple([str(name) for name in rec_data.dtype.names]) # for some reason the unicode throws off
        rec_data.dtype.names = names

        self.fields['pair_class']['values'] = list(field_data.pair_class)
        self.setData(rec_data)

    def plotClicked(self, plot, points):
        if len(self.selected_points) > 0:
            for pt, style in self.selected_points:
                brush, pen, size = style
                try:
                    pt.setBrush(brush)
                    pt.setPen(pen)
                    pt.setSize(size)
                except AttributeError:
                    pass
        self.selected_points = []
        for pt in points:
            style = (pt.brush(), pt.pen(), pt.size())
            self.selected_points.append([pt, style])
            data = pt.data()
            element = ('%s -> %s ' % (data.pre_class.name, data.post_class.name))
            print('Clicked:' '%s' % element)
            fields = self.fieldList.selectedItems()
            for field in fields:
                field_name = field.text()
                value = data[field_name]
                print('%s: %s' % (field_name, pg.siFormat(value)))
            pt.setBrush(pg.mkBrush('y'))
            pt.setSize(15)
        self.sigScatterPlotClicked.emit(self, points)

class PairScatterPlot(ScatterPlots):
    def __init__(self):
        ScatterPlots.__init__(self)
        
        header = pg.QtWidgets.QLabel()
        header.setText('<span style="font-weight: bold">Pair-wise Scatter Plot</span>')
        self.ctrlPanel.insertWidget(0, header)
        self.top_plot = None
        self.menu = pg.QtWidgets.QMenu()
        self.pair_actions = PairActions()
        for act in self.pair_actions.actions.values():
            self.menu.addAction(act)

    def contextMenuEvent(self, event):
        self.menu.popup(event.globalPos())


    def set_data(self, data):
        data['pair_class'] = data.apply(lambda row: '-'.join([row.pre_class.name, row.post_class.name]), axis=1)
        rec_data = data.to_records()
        names = tuple([str(name) for name in rec_data.dtype.names]) # for some reason the unicode throws off
        rec_data.dtype.names = names

        self.fields['pair_class']['values'] = list(data.pair_class)
        self.setData(rec_data)

    def filter_selected_element(self, pre_class, post_class):
        pair_name = '-'.join([pre_class.name, post_class.name])
        try:
            pair_filter = self.filter.child('pair_class')
        except KeyError:
            pair_filter = self.filter.addNew('pair_class')
            for child in pair_filter.children():
                child.setValue(False)
        pair_filter.child(pair_name).setValue(True)

    def reset_element_filter(self):
        try:
            pair_class_filter = self.filter.child('pair_class')
            self.filter.removeChild(pair_class_filter)
        except:
            return

    def plotClicked(self, plot, points):
        if self.top_plot is not None:
            self.plot.removeItem(self.top_plot)
        #     for pt, style in self.selected_points:
        #         brush, pen, size = style
        #         try:
        #             pt.setBrush(brush)
        #             pt.setPen(pen)
        #             pt.setSize(size)
        #         except AttributeError:
        #             pass
        selected_points = [[], []]
        for pt in points:
            # style = (pt.brush(), pt.pen(), pt.size())
            x, y = pt.pos()
            selected_points[0].append(x)
            selected_points[1].append(y)
            data = pt.data()
            pair = data.index
            print('Clicked:' '%s' % pair)
            print('pre-class: %s' % data['pre_class'])
            print('post-class: %s' % data['post_class'])
            fields = self.fieldList.selectedItems()
            for field in fields:
                field_name = field.text()
                value = data[field_name]
                print('%s: %s' % (field_name, pg.siFormat(value)))
            # pt.setBrush(pg.mkBrush('y'))
            # pt.setSize(15)
        self.top_plot = self.plot.plot(selected_points[0], selected_points[1], pen=None, symbol='o', symbolBrush='y', symbolSize=15)
        self.sigScatterPlotClicked.emit(self, points)
        self.pair_actions.pair = pair

