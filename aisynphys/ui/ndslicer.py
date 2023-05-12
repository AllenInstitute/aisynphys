from collections import OrderedDict
import numpy as np
from pyqtgraph.Qt import QtWidgets, QtCore
import pyqtgraph as pg
import pyqtgraph.dockarea


class NDSlicer(QtWidgets.QWidget):
    """Tool for visualizing 1D and 2D slices from an ND array.
    
    Parameters
    ----------
    axes : ordered dict
        Description of array axes to expect. Format is::
        
            {'axis_name': {'values': array}}
    """
    
    selection_changing = QtCore.Signal(object)
    selection_changed = QtCore.Signal(object)
    
    def __init__(self, axes):
        QtWidgets.QWidget.__init__(self)
        self.layout = QtWidgets.QGridLayout()
        self.setLayout(self.layout)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        self.dockarea = pg.dockarea.DockArea()
        self.layout.addWidget(self.dockarea)

        self.viewers = []
        self.data = None
        self.axes = OrderedDict([(ax, AxisData(name=ax, **ax_info)) for ax, ax_info in axes.items()])
        
        self.params = pg.parametertree.Parameter(name='params', type='group', children=[
            {'name': 'index', 'type': 'group', 'children': [{'name': ax, 'type': 'int', 'value': self.axes[ax].selection} for ax in self.axes]},
            {'name': 'max project', 'type': 'group', 'children': [{'name': ax, 'type': 'bool', 'value': False} for ax in self.axes]},
            ColorAxisParam(name='color axis', slicer=self),
            {'name': 'optimize', 'type': 'group', 'children': [{'name': ax, 'type': 'bool', 'value': False} for ax in self.axes]},
            {'name': '1D views', 'type': 'group'},
            MultiAxisParam(ndim=2, slicer=self),
        ])

        pos = {'position': 'right', 'relativeTo': None}
        for ax in axes:
            ch = pg.parametertree.Parameter.create(name=ax, type='bool', value=True)
            self.params.child('1D views').addChild(ch)
            ch.sigValueChanged.connect(self.one_d_show_changed)
            ch.viewer, dock = self.add_view(axes=[ax], position=pos)
            pos = {'position': 'bottom', 'relativeTo': dock}
            
        self.params.child('index').sigTreeStateChanged.connect(self.index_param_changed)
        self.params.child('max project').sigTreeStateChanged.connect(self.max_project_changed)
        self.params.child('color axis').axis_color_changed.connect(self.axis_color_changed)
        
        self.ctrl_split = pg.QtWidgets.QSplitter(pg.QtCore.Qt.Vertical)
        
        self.ptree = pg.parametertree.ParameterTree(showHeader=False)
        self.ptree.setParameters(self.params, showTop=False)
        self.ctrl_split.addWidget(self.ptree)
        
        self.histlut = pg.HistogramLUTWidget()
        self.histlut.sigLevelsChanged.connect(self.histlut_changed)
        self.histlut.sigLookupTableChanged.connect(self.histlut_changed)
        self.ctrl_split.addWidget(self.histlut)
        
        self.ctrl_dock = pg.dockarea.Dock("view selection")
        self.ctrl_dock.addWidget(self.ctrl_split)
        self.dockarea.addDock(self.ctrl_dock, 'left')

    def set_data(self, data, axes=None):
        """Set the data to be displayed.
        
        Parameters
        ----------
        data : array
            Data array of any dimensionality to be displayed
        axes : dict
            Optional description of axes in *data*.
        """
        self.data = data
        axes = axes or {}
        for ax,info in axes.items():
            for k,v in info.items():
                setattr(self.axes[ax], k, v)
        for viewer in self.viewers:
            viewer.set_data(self.data, self.axes)
        data_lim = (self.data.min(), self.data.max())
        self.histlut.setLevels(*data_lim)
        self.histlut.setHistogramRange(*data_lim)
        
    def add_view(self, axes, position=None):
        """Add a new 1D or 2D view.

        Parameters
        ----------
        axes : list
            List of 1 or 2 axis names to show in the view.
        position : dict | None
            Extra arguments used to set the position of the new dock (see
            pyqtgraph.dockarea.DockArea.addDock). 
        """
        dock = pg.dockarea.Dock("viewer", area=self.dockarea)
        if len(axes) == 1:
            viewer = OneDViewer(axes)
        elif len(axes) == 2:
            viewer = TwoDViewer(axes)
        dock.addWidget(viewer)
        position = position or {'position': 'right'}
        self.dockarea.addDock(dock, **position)
        viewer.dock = dock
        viewer.selection_changed.connect(self.viewer_selection_changed)
        viewer.selection_changing.connect(self.viewer_selection_changing)
        self.viewers.append(viewer)
        viewer.set_data(self.data, self.axes)
        self.histlut_changed()
        return viewer, dock

    def set_selection(self, axes, emit=True):
        """Set the current selection along any set of axes.

        Parameters:
        -----------
        axes : dict
            Dictionary of {'axis_name': value} pairs specifying the new selection.
        """
        # update index parameters
        for ax,val in axes.items():
            self._set_axis_value(ax, val)

        # auto-optimize other parameters if requested
        optimize_axes = []
        for k in self.axes:
            # only optimize axes if they have been marked for optimizatin _and_ they are not being explicitly set here
            if k not in axes and self.params['optimize', k]:
                optimize_axes.append(k)
        if len(optimize_axes) > 0:
            # find best index for optimized axes
            index = self.index()
            opt_data = self.data
            for i,k in list(enumerate(self.axes.keys()))[::-1]:
                if k not in optimize_axes:
                    # take single index for axes that are not being optimized
                    opt_data = np.take(opt_data, index[k], axis=i)
            max_ind = np.unravel_index(opt_data.argmax(), opt_data.shape)
            for i,k in enumerate(optimize_axes):
                ax = self.axes[k]
                val = ax.values[max_ind[i]]
                self._set_axis_value(k, val)

        # process updates
        for viewer in self.viewers:
            viewer.update_selection()
        if emit:
            self.selection_changed.emit(self)

    def _set_axis_value(self, ax, val):
        with pg.SignalBlock(self.params.child('index').sigTreeStateChanged, self.index_param_changed):
            self.params['index', ax] = self.axes[ax].index_at(val)
        self.axes[ax].selection = val

    def set_index(self, index):
        """Set currently selected indices on any axes.

        Parameters:
        -----------
        index : list
            List of indices to select on each axis.
        """
        select_values = {}
        for i,x in enumerate(index):
            ax = list(self.axes.values())[i]
            select_values[ax.name] = ax.values[x]
        self.set_selection(select_values)

    def selection(self):
        """Return an ordered dictionary of the currently selected values::

            {'axis_name': value, ...}
        """
        vals = OrderedDict([(ax,val.selection) for ax,val in self.axes.items()])
        return vals

    def index(self):
        """Return an ordered dictionary of the currently selected indices::

            {'axis_name': index, ...}
        """
        return OrderedDict([(ax,val.index) for ax,val in self.axes.items()])

    def one_d_show_changed(self, param):
        param.viewer.dock.setVisible(param.value())
        
    def viewer_selection_changed(self, viewer, axes):
        self.set_selection(axes, emit=True)

    def viewer_selection_changing(self, viewer, axes):
        self.set_selection(axes, emit=False)
        self.selection_changing.emit(self)

    def index_param_changed(self, param, changes):
        sel = {}
        for param, change, value in changes:
            if change != 'value':
                continue            
            sel[param.name()] = self.axes[param.name()].value_at(value)
        self.set_selection(sel)

    def histlut_changed(self):
        for viewer in self.viewers:
            if not isinstance(viewer, TwoDViewer):
                continue
            viewer.set_image_params(self.histlut.getLevels(), self.histlut.getLookupTable(n=1024))

    def max_project_changed(self):
        for ax in self.axes.values():
            ax.max_project = self.params['max project', ax.name]
        for viewer in self.viewers:
            viewer.update_selection()

    def axis_color_changed(self):
        for ax in self.axes.values():
            if ax.name != self.params['color axis', 'axis']:
                ax.colors = None
            else:
                ax.colors = self.params.child('color axis').colors
        for viewer in self.viewers:
            viewer.update_selection()
            

class AxisData(object):
    def __init__(self, name, values):
        self.name = name
        self.values = values
        self.selection = values[0]
        self.max_project = False
        self.colors = None

    def index_at(self, x):
        return np.argmin(np.abs(self.values - x))
        
    def value_at(self, i):
        return self.values[i]

    @property
    def index(self):
        return self.index_at(self.selection)


class MultiAxisParam(pg.parametertree.types.GroupParameter):
    def __init__(self, ndim, slicer):
        self.ndim = ndim
        self.slicer = slicer
        pg.parametertree.types.GroupParameter.__init__(self, name="%dD views"%ndim, addText="Add new..")
        
    def addNew(self):
        axis_names = list(self.slicer.axes.keys())
        param = pg.parametertree.Parameter(name="%dD view" % self.ndim, autoIncrementName=True, type='group', children=[
            {'name': 'axis %d' % i, 'type': 'list', 'values': axis_names, 'value': axis_names[i]} for i in range(self.ndim)
        ], removable=True)
        self.addChild(param)
        viewer, dock = self.slicer.add_view(axes=list(self.slicer.axes.keys())[:self.ndim])
        param.viewer = viewer
        param.sigTreeStateChanged.connect(self.axes_changed)
        return param

    def axes_changed(self, param, changes):
        axes = [ch.value() for ch in param.children()]
        if len(axes) == 2 and axes[0] == axes[1]:
            return
        param.viewer.set_selected_axes(axes)


class ColorAxisParam(pg.parametertree.types.GroupParameter):
    
    axis_color_changed = QtCore.Signal(object)  # self
    
    def __init__(self, name, slicer):
        self.slicer = slicer
        pg.parametertree.types.GroupParameter.__init__(self, name=name, children=[
            {'name': 'axis', 'type': 'list', 'values': ['none']},
            {'name': 'mode', 'type': 'list', 'values': ['discrete', 'range'], 'visible': False},
            {'name': 'colors', 'type': 'group', 'visible': False},
        ])
        
        self.update_axes()
        self.child('axis').sigValueChanged.connect(self.axis_changed)
        self.child('colors').sigTreeStateChanged.connect(self.colors_changed)
        
    def update_axes(self):
        axis_names = list(self.slicer.axes.keys())
        self.child('axis').setLimits(['none'] + axis_names)

    def axis_changed(self):
        with pg.SignalBlock(self.child('colors').sigTreeStateChanged, self.colors_changed):
            if self['axis'] == 'none':
                self.child('mode').setOpts(visible=False)
                self.child('colors').setOpts(visible=False)
            else:
                self.child('mode').setOpts(visible=True)
                self.child('colors').setOpts(visible=True)
                values = self.slicer.axes[self['axis']].values
                for ch in self.child('colors').children():
                    self.child('colors').removeChild(ch)
                for i,v in enumerate(values):
                    color = pg.mkColor((i, int(len(values)*1.2)))
                    ch = pg.parametertree.types.SimpleParameter(name=str(v), type='color', value=color)
                    self.child('colors').addChild(ch)
            
        self.colors_changed()

    def colors_changed(self):
        self.colors = [ch.value() for ch in self.child('colors').children()]
        self.axis_color_changed.emit(self)


class Viewer(object):
    def __init__(self, ax):
        self.data = None
        self.data_axes = None
        self.set_selected_axes(ax)

    def set_data(self, data, axes):
        self.data = data
        self.data_axes = axes
        self.update_selection()
        
    def set_selected_axes(self, axes):
        self.selected_axes = axes
        if self.data is not None:        
            self.update_selection()

    def update_selection(self):
        self.update_display()
        
    def update_display(self):
        raise NotImplementedError()
        
    def get_data(self):
        data = self.data
        
        # slice or flatten non-visible axes
        axis_names = list(self.data_axes.keys())
        colormap_axis = None
        removed = 0
        for i in range(len(axis_names)):
            j = i - removed
            ax_name = axis_names[j]
            ax = self.data_axes[ax_name]
            if ax_name in self.selected_axes:
                continue
            if ax.colors is not None:
                colormap_axis = ax
                continue
            if ax.max_project:
                # max projection across this axis
                data = data.max(axis=j)
            else:
                # slice this axis
                data = data.take(ax.index, j)
            axis_names.pop(j)
            removed += 1
            
        # re-order visible axes
        order = [axis_names.index(ax) for ax in self.selected_axes]
        if colormap_axis is not None:
            order.append(axis_names.index(colormap_axis.name))
        # data = data.transpose(np.argsort(order))
        data = data.transpose(order)
        
        colormap = None if colormap_axis is None else colormap_axis.colors
        
        return data, colormap


class OneDViewer(Viewer, pg.PlotWidget):
    selection_changing = QtCore.Signal(object, object)  # self, {axis: value, ...}
    selection_changed = QtCore.Signal(object, object)  # self, {axis: value, ...}
    
    def __init__(self, axes):
        pg.PlotWidget.__init__(self)
        self.line = self.addLine(x=0, movable=True)
        self.curves = []

        Viewer.__init__(self, axes)

        self.line.sigDragged.connect(self.line_moved)
        self.line.sigPositionChangeFinished.connect(self.line_move_finished)

    def update_selection(self):
        axis = self.selected_axes[0]
        with pg.SignalBlock(self.line.sigPositionChangeFinished, self.line_move_finished):
            self.line.setValue(self.data_axes[axis].index)
        Viewer.update_selection(self)

    def line_moved(self):
        ax = self.selected_axes[0]
        val = self.data_axes[ax].value_at(int(np.round(self.line.value())))
        self.selection_changing.emit(self, {ax: val})
        
    def line_move_finished(self):
        ax = self.selected_axes[0]
        val = self.data_axes[ax].value_at(int(np.round(self.line.value())))
        self.selection_changed.emit(self, {ax: val})
        
    def update_display(self):
        if self.data is None:
            self.clear_curves()
            return
            
        axis = self.selected_axes[0]
        self.setLabels(bottom=axis)
        data, colors = self.get_data()
            
        if colors is None:
            data = data[..., np.newaxis]
            colors = ['w']
            
        new_curves = []
        for i in range(data.shape[-1]):
            c = pg.PlotCurveItem(data[...,i], pen=colors[i], antialias=True)
            new_curves.append(c)
            self.addItem(c)

        self.clear_curves()
        self.curves = new_curves
            
        axvals = self.data_axes[axis].values
        self.getAxis('bottom').setTicks([[(i, "%0.2g"%axvals[i]) for i in range(len(axvals))]])
        
    def clear_curves(self):
        for c in self.curves:
            self.removeItem(c)
            self.curves = []
        

class TwoDViewer(Viewer, pg.GraphicsLayoutWidget):
    selection_changing = QtCore.Signal(object, object)  # self, {axis: value, ...}
    selection_changed = QtCore.Signal(object, object)  # self, {axis: value, ...}

    def __init__(self, axes):
        self.data_bounds = (0, 1)
        self.levels = None
        self.lut = None
        
        pg.GraphicsLayoutWidget.__init__(self)
        self.plot = self.addPlot()
        self.plot.invertY(False)
        self.plot.setAspectLocked(False)
        
        self.image = pg.ImageItem()
        self.plot.addItem(self.image)
        self.image.setPos(-0.5, -0.5)
        
        self.lines = [self.plot.addLine(x=0, movable=True), self.plot.addLine(y=0, movable=True)]
        self.lines[0]._viewer_axis = 0
        self.lines[1]._viewer_axis = 1
        
        Viewer.__init__(self, axes)
        for line in self.lines:
            line.sigDragged.connect(self.line_moved)
        for line in self.lines:
            line.sigPositionChangeFinished.connect(self.line_move_finished)

    def set_data(self, data, axes):
        if data is not None:
            self.data_bounds = (data.min(), data.max())
        Viewer.set_data(self, data, axes)
        
    def set_image_params(self, levels, lut):
        self.image.setLevels(levels)
        self.image.setLookupTable(lut)
        self.levels = levels
        self.lut = lut
        
    def update_selection(self):
        for i,line in enumerate(self.lines):
            ax = self.selected_axes[i]
            with pg.SignalBlock(line.sigPositionChangeFinished, self.line_move_finished):
                line.setValue(self.data_axes[ax].index)
        Viewer.update_selection(self)

    def line_moved(self):
        line = self.sender()
        ax = self.selected_axes[line._viewer_axis]
        axes = {ax: self.data_axes[ax].value_at(int(np.round(line.value())))}
        self.selection_changing.emit(self, axes)

    def line_move_finished(self):
        line = self.sender()
        ax = self.selected_axes[line._viewer_axis]
        axes = {ax: self.data_axes[ax].value_at(int(np.round(line.value())))}
        self.selection_changed.emit(self, axes)

    def update_display(self):
        if self.data is None:
            self.image.setImage(np.zeros((1, 1)))
            return
        axes = self.selected_axes
        self.plot.setLabels(left=axes[1], bottom=axes[0])
        data, colors = self.get_data()
        
        levels = self.levels or self.data_bounds
        
        if colors is None:
            self.image.setImage(data, levels=levels)
        else:
            comp = np.zeros(data.shape[:-1] + (3,), dtype=data.dtype)
            for i,color in enumerate(colors):
                color = np.array(pg.colorTuple(color)[:3], dtype=float).reshape(1, 1, 3) / 255.
                comp += data[..., i:i+1] * color
            self.image.setImage(comp, levels=levels)

        xvals = self.data_axes[axes[0]].values
        yvals = self.data_axes[axes[1]].values
        self.plot.getAxis('bottom').setTicks([[(i, "%0.2g"%xvals[i]) for i in range(len(xvals))]])
        self.plot.getAxis('left').setTicks([[(i, "%0.2g"%yvals[i]) for i in range(len(yvals))]])


if __name__ == '__main__':
    pg.mkQApp()
    pg.dbg()
    
    data = np.random.normal(size=(10, 20, 40)) * 0.1
    data += np.sin(np.linspace(0, 2*np.pi, data.shape[0]))[:, None, None]
    data += np.sin(np.linspace(0, 4*np.pi, data.shape[1]))[None, :, None]
    data += np.sin(np.linspace(0, 8*np.pi, data.shape[2]))[None, None, :]
    
    axes = OrderedDict([
        ('X', {'values': np.arange(data.shape[0])}),
        ('Y', {'values': np.arange(data.shape[1]) * 200}),
        ('Z', {'values': np.arange(data.shape[2]) * 0.01 + 3.4}),
    ])
        
    nds = NDSlicer(axes)
    nds.show()
    
    nds.set_data(data)
    