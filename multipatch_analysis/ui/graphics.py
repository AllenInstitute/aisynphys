# *-* coding: utf-8 *-*
from __future__ import print_function, division
import numpy as np
import pyqtgraph as pg
from neuroanalysis.stats import binomial_ci
from neuroanalysis.ui.plot_grid import PlotGrid
from statsmodels.stats.proportion import proportion_confint



class MatrixItem(pg.QtGui.QGraphicsItemGroup):
    """GraphicsItem displaying a table with column / row labels and text in
    each cell.
    
    Parameters
    ----------
    text : 2d array or nested lists
        Strings to display inside each cell
    fgcolor : 2d array or nested lists
        Text colors for each cell
    bgcolor : 2d array or nested lists
        Background colors for each cell
    rows : 1d array or list
        Strings or tuples to display as row header. If tuple, then there are multiple headers,
        with rows grouped together. For example, ``[('l4', 'pv'), ('l4', 'sst'), ('l5', 'pv'), ...]``
        specifies rows grouped first by layer (l4, l5) then by cre type.
    cols : 1d array or list
        Strings to display as col header. See *rows*.
    size : float
        Width of each cell
    border_color : 2d array or nested lists
        Border colors for each cell
    header_color : str | tuple
        Color of header text
    """
    class SignalHandler(pg.QtCore.QObject):
        """Because we can't subclass from both QObject and QGraphicsRectItem at the same time
        """
        sigClicked = pg.QtCore.Signal(object, object, object, object) # self, event, row, col

    def __init__(self, text, fgcolor, bgcolor, rows, cols, size=50, border_color='k', header_color='w'):
        pg.QtGui.QGraphicsItemGroup.__init__(self)
        self._signalHandler = MatrixItem.SignalHandler()
        self.sigClicked = self._signalHandler.sigClicked
        self.cell_size = size
        self.header_color = header_color

        self.header_labels = {}
        self.group_items = {}
        self._make_header(rows, 'left')
        self._make_header(cols, 'top')

        self.cell_labels = []
        self.cells = []
        for i,row in enumerate(rows):
            self.cells.append([])
            self.cell_labels.append([])
            for j,col in enumerate(cols):
                x = j * size
                y = i * size
                rect = MatrixElementItem(size, parent=self)
                rect.setPos(x, y)
                rect.setBrush(pg.mkBrush(bgcolor[i][j]))
                rect.setPen(pg.mkPen(border_color[i][j]))
                rect.setZValue(-10)
                rect.row = i
                rect.col = j
                rect.sigClicked.connect(self.element_clicked)
                self.cells[-1].append(rect)

                txt = pg.QtGui.QGraphicsTextItem(text[i][j], parent=self)
                br = txt.boundingRect()
                txt.setTextWidth(br.width())
                txt.setPos(x + size/2 - br.center().x(), y + size/2 - br.center().y())
                txt.setDefaultTextColor(pg.mkColor(fgcolor[i][j]))
                self.cell_labels[-1].append(txt)

        br = pg.QtCore.QRectF()
        for item in self.childItems():
            br = br.united(self.mapRectFromItem(item, item.boundingRect()))
        self._bounding_rect = br

        font_size = size / 2
        for i in ('Presynaptic', 'Postsynaptic'):
            html = '<span style="font-size: %dpx; font-weight: bold">%s</span>' % (font_size, i)
            item = pg.QtGui.QGraphicsTextItem("", parent=self)
            item.setHtml(html)
            if i == 'Presynaptic':
                item.rotate(-90)
                x = self.boundingRect().left() - item.boundingRect().height()
                y = item.boundingRect().width() + ((self.boundingRect().bottom() - item.boundingRect().width()) / 2)
            elif i == 'Postsynaptic':
                x = (self.boundingRect().right() - item.boundingRect().width()) / 2
                y = self.boundingRect().top() - item.boundingRect().height()

            item.setPos(x, y)

    def element_clicked(self, rect, event):
        self.sigClicked.emit(self, event, rect.row, rect.col)  


    def _make_header(self, labels, side):
        padding = 5
        if isinstance(labels[0], tuple):
            # draw groups first
            grp_labels, labels = zip(*labels)
        else:
            grp_labels = None

        self.header_labels[side] = [self._make_header_text(label, i, side, padding=padding) for i,label in enumerate(labels)]
        
        if grp_labels is not None:
            width = max(map(lambda l: l.boundingRect().width(), self.header_labels[side]))
            # measure range for each group
            grps = []
            for i,l in enumerate(grp_labels):
                if len(grps) == 0 or grps[-1][0] != l:
                    grps.append([l, i, i])
                else:
                    grps[-1][2] = i

            self._make_header_groups(grps, side, width)
            padding = 3

    def _make_header_text(self, txt, rowcol, side, padding=5, font_size=None):
        size = self.cell_size
        if font_size is None:
            font_size = size / 3.
        small_font_size = font_size * 0.7
        align = pg.QtCore.Qt.AlignLeft if side == 'top' else pg.QtCore.Qt.AlignRight
        lines = str(txt).split('\n')
        html = '<div style="line-height: 70%%">'
        for i,line in enumerate(lines):
            fs = font_size if i == 0 else small_font_size
            if i > 0:
                html += '<br>'
            html += '<span style="font-size: %dpx;">%s</span>' % (fs, line)
        html += '</div>'
        item = pg.QtGui.QGraphicsTextItem("", parent=self)
        item.setHtml(html)
        item.setTextWidth(item.boundingRect().width())

        doc = item.document()
        opts = doc.defaultTextOption()
        opts.setAlignment(align)
        doc.setDefaultTextOption(opts)

        if side == 'top':
            item.rotate(-90)
        br = item.mapRectToParent(item.boundingRect())
        if side == 'top':
            item.setPos(rowcol * size + size/2 - br.center().x(), -br.bottom() - padding)
        elif side == 'left':
            item.setPos(-br.right() - padding, rowcol * size + size/2. - br.center().y())
        else:
            raise ValueError("side must be top or left")
        item.setDefaultTextColor(pg.mkColor(self.header_color))
        item.setTextWidth(item.boundingRect().width())  # needed to enable text alignment
        return item

    def _make_header_groups(self, grps, side, width):
        size = self.cell_size
        self.group_items[side] = []
        for label, start, stop in grps:
            n_cells = 1 + stop - start
            w = n_cells * size
            h = size
            path = pg.QtGui.QPainterPath()
            path.moveTo(0, 0)
            path.lineTo(w, 0)
            path.lineTo(w, -width)
            path.lineTo(w/2, -width-h*0.5)
            path.lineTo(0, -width)
            path.closeSubpath()
            item = pg.QtGui.QGraphicsPathItem(self)
            item.setPath(path)
            item.setBrush(pg.mkBrush(0.8))
            item.setPen(pg.mkPen(color='k', width=1))
            item.setZValue(-5)
            if side == 'top':
                br = item.mapRectToParent(item.boundingRect())
                item.setPos(start * size - br.left(), -br.bottom())
            elif side == 'left':
                item.rotate(-90)
                br = item.mapRectToParent(item.boundingRect())
                item.setPos(-br.right(), start * size - br.top())
            else:
                raise ValueError("side must be top or left")

            txt = self._make_header_text(label, (start+stop)/2., side, padding=(10 + width + size*0.5), font_size=size/2.)
            self.group_items[side].append((item, txt))

    def boundingRect(self):
        return self._bounding_rect

class MatrixElementItem(pg.QtGui.QGraphicsRectItem):
    class SignalHandler(pg.QtCore.QObject):
        """Because we can't subclass from both QObject and QGraphicsRectItem at the same time
        """
        sigClicked = pg.QtCore.Signal(object, object) # self, event

    def __init__(self, size, parent):
        self._signalHandler = MatrixElementItem.SignalHandler()
        self.sigClicked = self._signalHandler.sigClicked
        pg.QtGui.QGraphicsRectItem.__init__(self, 0, 0, size, size, parent=parent)

    def mouseClickEvent(self, event):
        self.sigClicked.emit(self, event)    
    
def distance_plot(connected, distance, plots=None, color=(100, 100, 255), size=10, window=40e-6, spacing=None, name=None, fill_alpha=30):
    """Draw connectivity vs distance profiles with confidence intervals.
    
    Parameters
    ----------
    connected : boolean array
        Whether a synaptic connection was found for each probe
    distance : array
        Distance between cells for each probe
    plots : list of PlotWidget | PlotItem
        (optional) Two plots used to display distance profile and scatter plot.
    color : tuple
        (R, G, B) color values for line and confidence interval. The confidence interval
        will be drawn with alpha=100
    size: int
        size of scatter plot symbol
    window : float
        Width of distance window over which proportions are calculated for each point on
        the profile line.
    spacing : float
        Distance spacing between points on the profile line


    Note: using a spacing value that is smaller than the window size may cause an
    otherwise smooth decrease over distance to instead look more like a series of downward steps.
    """
    color = pg.colorTuple(pg.mkColor(color))[:3]
    connected = np.array(connected).astype(float)
    distance = np.array(distance)

    # scatter plot connections probed
    if plots is None:
        grid = PlotGrid()
        grid.set_shape(2, 1)
        grid.grid.ci.layout.setRowStretchFactor(0, 5)
        grid.grid.ci.layout.setRowStretchFactor(1, 10)
        plots = (grid[1,0], grid[0,0])
        plots[0].grid = grid
        plots[0].addLegend()
        grid.show()
        plots[0].setLabels(bottom=('distance', 'm'), left='connection probability')

    if plots[1] is not None:
         # scatter points a bit
        pts = np.vstack([distance, connected]).T
        conn = pts[:,1] == 1
        unconn = pts[:,1] == 0
        if np.any(conn):
            cscat = pg.pseudoScatter(pts[:,0][conn], spacing=10e-6, bidir=False)
            mx = abs(cscat).max()
            if mx != 0:
                cscat = cscat * 0.2# / mx
            pts[:,1][conn] = -5e-5 - cscat
        if np.any(unconn):
            uscat = pg.pseudoScatter(pts[:,0][unconn], spacing=10e-6, bidir=False)
            mx = abs(uscat).max()
            if mx != 0:
                uscat = uscat * 0.2# / mx
            pts[:,1][unconn] = uscat
        
        plots[1].setXLink(plots[0])
        plots[1].hideAxis('bottom')
        plots[1].hideAxis('left')

        color2 = color + (100,)
        scatter = plots[1].plot(pts[:,0], pts[:,1], pen=None, symbol='o', labels={'bottom': ('distance', 'm')}, size=size, symbolBrush=color2, symbolPen=None, name=name)
        scatter.scatter.opts['compositionMode'] = pg.QtGui.QPainter.CompositionMode_Plus

    # use a sliding window to plot the proportion of connections found along with a 95% confidence interval
    # for connection probability

    if spacing is None:
        spacing = window / 4.0
        
    xvals = np.arange(window / 2.0, 500e-6, spacing)
    upper = []
    lower = []
    prop = []
    ci_xvals = []
    for x in xvals:
        minx = x - window / 2.0
        maxx = x + window / 2.0
        # select points inside this window
        mask = (distance >= minx) & (distance <= maxx)
        pts_in_window = connected[mask]
        # compute stats for window
        n_probed = pts_in_window.shape[0]
        n_conn = pts_in_window.sum()
        if n_probed == 0:
            prop.append(np.nan)
        else:
            prop.append(n_conn / n_probed)
            ci = proportion_confint(n_conn, n_probed, method='beta')
            lower.append(ci[0])
            upper.append(ci[1])
            ci_xvals.append(x)

    # plot connection probability and confidence intervals
    color2 = [c / 3.0 for c in color]
    mid_curve = plots[0].plot(xvals, prop, pen={'color': color, 'width': 3}, antialias=True, name=name)
    upper_curve = plots[0].plot(ci_xvals, upper, pen=(0, 0, 0, 0), antialias=True)
    lower_curve = plots[0].plot(ci_xvals, lower, pen=(0, 0, 0, 0), antialias=True)
    upper_curve.setVisible(False)
    lower_curve.setVisible(False)
    color2 = color + (fill_alpha,)
    fill = pg.FillBetweenItem(upper_curve, lower_curve, brush=color2)
    fill.setZValue(-10)
    plots[0].addItem(fill, ignoreBounds=True)
    
    return plots, ci_xvals, prop, upper, lower


        
