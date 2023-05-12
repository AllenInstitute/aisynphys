"""
ScatterPlotWidgets for Matrix Analyzer. 
Allows plotting of Matrix data on a per cell basis.

"""

from __future__ import print_function, division
import pyqtgraph as pg
import colorsys
import pandas as pd
import numpy as np
from aisynphys.ui.ScatterPlotWidget import ScatterPlotWidget

class CellScatterTab(pg.QtWidgets.QWidget):
    def __init__(self):
        pg.QtWidgets.QWidget.__init__(self)
        self.layout = pg.QtWidgets.QGridLayout()
        self.setLayout(self.layout)
        self.v_splitter = pg.QtWidgets.QSplitter()
        # self.v_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        # self.layout.addWidget(self.v_splitter)
        self.cell_scatter = CellScatterPlot()
        self.layout.addWidget(self.cell_scatter)
        # self.morpho_scatter = MorphoScatter()
        # self.v_splitter.add(self.morpho_scatter)

        # self.colorMap.sigColorMapChanged.connect(self.set_colors)

class CellScatterPlot(ScatterPlotWidget):
    def __init__(self):
        ScatterPlotWidget.__init__(self)
        self.selected_points = []

    def set_fields(self, fields):
        self.fields = [('CellClass', {'mode': 'enum'})]
        self.fields.extend([f for f in fields if f != ('None', {})])
        self.setFields(self.fields)

    def set_data(self, data):
        rec_data = data.to_records()
        names = tuple([str(name) for name in rec_data.dtype.names]) # for some reason the unicode throws off
        rec_data.dtype.names = names

        for field, options in self.fields.items():
            data_type = options.get('mode')
            values = options.get('values')
            defaults = options.get('defaults')
            if data_type == 'enum':
                if values is None:
                    unique_values = list(set(data.get(field)))
                    if field =='CellClass':
                        self.fields[field]['values'] =  sorted(unique_values, key=lambda x: (x.name is None, x.name))
                    else: 
                        unique_values = [v for v in unique_values if isinstance(v, str)]
                        self.fields[field]['values'] =  sorted(unique_values, key=lambda x: (x is None, x)) 
                if defaults is None:
                    n_colors = len(set(data.get(field))) if values is None else len(values)
                    self.fields[field]['defaults'] = {'colormap': [pg.intColor(n, n_colors) for n in np.arange(n_colors)]}
              
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
            print('Clicked:' '%s' % data.index)
            # fields = self.fieldList.selectedItems()
            for field, options in self.fields.items():
                value = data[field]
                if value is not None:
                    if options.get('mode') == 'range':
                        print('%s: %s' % (field, pg.siFormat(value)))
                    elif options.get('mode') == 'enum':
                        print('%s: %s' % (field, value))
            pt.setBrush(pg.mkBrush('y'))
            pt.setSize(15)
        self.sigScatterPlotClicked.emit(self, plot, points)

    def color_selected_element(self, color, pre_class, post_class):
        try:
            cell_class_map = self.colorMap.child('CellClass')
            cell_class_style = self.style.child('CellClass')
        except KeyError:
            cell_class_map = self.colorMap.addNew('CellClass')
            cell_class_style = self.style.addNew('CellClass')
        for cell_class in self.fields['CellClass']['values']:
            cell_class_map['Values', cell_class] = pg.mkColor((128, 128, 128))
        cell_class_map['Values', pre_class.name] = pg.mkColor(color)
        cell_class_map['Values', post_class.name] = pg.mkColor(color)
        pre_style = [c for c in cell_class_style.children() if c.name() == pre_class.name][0]
        pre_style.setValue(True)
        pre_style['Symbol'] = 't'
        post_style = [c for c in cell_class_style.children() if c.name() == post_class.name][0]
        post_style.setValue(True)
        post_style['Symbol'] = 't2'

    def reset_element_color(self):
        try:
            cell_class_map = self.colorMap.child('CellClass')
            self.colorMap.removeChild(cell_class_map)
            cell_class_style = self.style.child('CellClass')
            self.style.removeChild(cell_class_style)
        except:
            return

    def invalidate_output(self):
        self.data = None

# class PatchSeqScatter(CellScatterPlot):
#     def __init__(self):
#         CellScatterPlot.__init__(self)

# class MorphoScatter(CellScatterPlot):
#     def __init__(self):
#         CellScatterPlot.__init__(self)