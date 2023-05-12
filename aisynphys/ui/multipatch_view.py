import numpy as np
from scipy.ndimage import gaussian_filter
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets, QtCore
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.miesnwb import MiesNwb
from neuroanalysis.filter import remove_artifacts


class MultipatchMatrixView(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.layout = QtWidgets.QGridLayout()
        self.setLayout(self.layout)
        self.layout.setContentsMargins(0, 0, 0, 0)

        self.plots = PlotGrid(parent=self)
        self.layout.addWidget(self.plots, 0, 0)
        self.plots.scene().sigMouseClicked.connect(self._plot_clicked)

        #self.pair_view = PairAnalyzer()

        self.params = pg.parametertree.Parameter(name='params', type='group', children=[
            {'name': 'show', 'type': 'list', 'values': ['sweep avg', 'sweep avg + sweeps', 'sweeps', 'pulse avg']},
            {'name': 'lowpass', 'type': 'bool', 'value': True, 'children': [
                {'name': 'sigma', 'type': 'float', 'value': 200e-6, 'step': 1e-5, 'limits': [0, None], 'suffix': 's', 'siPrefix': True},
            ]},
            {'name': 'first pulse', 'type': 'int', 'value': 0, 'limits': [0, None]},
            {'name': 'last pulse', 'type': 'int', 'value': 7, 'limits': [0, None]},
            {'name': 'window', 'type': 'float', 'value': 30e-3, 'step': 1e-3, 'limits': [0, None], 'suffix': 's', 'siPrefix': True},
            {'name': 'remove artifacts', 'type': 'bool', 'value': True, 'children': [
                {'name': 'window', 'type': 'float', 'suffix': 's', 'siPrefix': True, 'value': 1e-3, 'step': 1e-4, 'bounds': [0, None]},
            ]},
            {'name': 'remove baseline', 'type': 'bool', 'value': True},
            {'name': 'show ticks', 'type': 'bool', 'value': True},
        ])
        self.params.sigTreeStateChanged.connect(self._params_changed)

    def show_group(self, grp):
        self.show_sweeps(grp.sweeps)

    def data_selected(self, sweeps, channels):
        self.sweeps = sweeps
        self.channels = channels
        self._update_plots(auto_range=True)

    def _params_changed(self, *args):
        self._update_plots()

    def _plot_clicked(self, ev):
        item = self.plots.scene().itemAt(ev.scenePos())
        r,c = self.plots.item_index(item)
        for i in range(self.plots.rows):
            for j in range(self.plots.cols):
                color = None if (i, j) != (r, c) else pg.mkColor(30, 30, 50)
                self.plots[i,j].vb.setBackgroundColor(color)

    def _update_plots(self, auto_range=False):
        sweeps = self.sweeps
        chans = self.channels
        self.plots.clear()
        if len(sweeps) == 0 or len(chans) == 0:
            return
        
        # collect data
        data = MiesNwb.pack_sweep_data(sweeps)
        data, stim = data[...,0], data[...,1]  # unpack stim and recordings
        dt = sweeps[0].recordings[0]['primary'].dt

        # mask for selected channels
        mask = np.array([ch in chans for ch in sweeps[0].devices])
        data = data[:, mask]
        stim = stim[:, mask]
        chans = np.array(sweeps[0].devices)[mask]

        modes = [sweeps[0][ch].clamp_mode for ch in chans]
        
        # get pulse times for each channel
        stim = stim[0]
        diff = stim[:,1:] - stim[:,:-1]
        # note: the [1:] here skips the test pulse
        on_times = [np.argwhere(diff[i] > 0)[1:,0] for i in range(diff.shape[0])]
        off_times = [np.argwhere(diff[i] < 0)[1:,0] for i in range(diff.shape[0])]

        # remove capacitive artifacts from adjacent electrodes
        if self.params['remove artifacts']:
            npts = int(self.params['remove artifacts', 'window'] / dt)
            for i in range(stim.shape[0]):
                for j in range(stim.shape[0]):
                    if i == j:
                        continue
                    
                    # are these headstages adjacent?
                    hs1, hs2 = chans[i], chans[j]
                    if abs(hs2-hs1) > 3:
                        continue
                    
                    # remove artifacts
                    for k in range(len(on_times[i])):
                        on = on_times[i][k]
                        off = off_times[i][k]
                        data[:, j, on:on+npts] = data[:, j, max(0,on-npts):on].mean(axis=1)[:,None]
                        data[:, j, off:off+npts] = data[:, j, max(0,off-npts):off].mean(axis=1)[:,None]

        # lowpass filter
        if self.params['lowpass']:
            data = gaussian_filter(data, (0, 0, self.params['lowpass', 'sigma'] / dt))

        # prepare to plot
        window = int(self.params['window'] / dt)
        n_sweeps = data.shape[0]
        n_channels = data.shape[1]
        self.plots.set_shape(n_channels, n_channels)
        self.plots.setClipToView(True)
        self.plots.setDownsampling(True, True, 'peak')
        self.plots.enableAutoRange(False, False)

        show_sweeps = 'sweeps' in self.params['show']
        show_sweep_avg = 'sweep avg' in self.params['show']
        show_pulse_avg = self.params['show'] == 'pulse avg'

        for i in range(n_channels):
            for j in range(n_channels):
                plt = self.plots[i, j]
                start = on_times[j][self.params['first pulse']] - window
                if start < 0:
                    frontpad = -start
                    start = 0
                else:
                    frontpad = 0
                stop = on_times[j][self.params['last pulse']] + window

                # select the data segment to be displayed in this matrix cell
                # add padding if necessary
                if frontpad == 0:
                    seg = data[:, i, start:stop].copy()
                else:
                    seg = np.empty((data.shape[0], stop + frontpad), data.dtype)
                    seg[:, frontpad:] = data[:, i, start:stop]
                    seg[:, :frontpad] = seg[:, frontpad:frontpad+1]

                # subtract off baseline for each sweep
                if self.params['remove baseline']:
                    seg -= seg[:, :window].mean(axis=1)[:,None]

                if show_sweeps:
                    alpha = 100 if show_sweep_avg else 200
                    color = (255, 255, 255, alpha)
                    t = np.arange(seg.shape[1]) * dt
                    for k in range(n_sweeps):
                        plt.plot(t, seg[k], pen={'color': color, 'width': 1}, antialias=True)

                if show_sweep_avg or show_pulse_avg:
                    # average selected segments over all sweeps
                    segm = seg.mean(axis=0)

                    if show_pulse_avg:
                        # average over all selected pulses
                        pulses = []
                        for k in range(self.params['first pulse'], self.params['last pulse'] + 1):
                            pstart = on_times[j][k] - on_times[j][self.params['first pulse']]
                            pstop = pstart + (window * 2)
                            pulses.append(segm[pstart:pstop])
                        # for p in pulses:
                        #     t = np.arange(p.shape[0]) * dt
                        #     plt.plot(t, p)
                        segm = np.vstack(pulses).mean(axis=0)

                    t = np.arange(segm.shape[0]) * dt

                    if i == j:
                        color = (80, 80, 80)
                    else:
                        dif = segm - segm[:window].mean()
                        qe = 30 * np.clip(dif, 0, 1e20).mean() / segm[:window].std()
                        qi = 30 * np.clip(-dif, 0, 1e20).mean() / segm[:window].std()
                        if modes[i] == 'ic':
                            qi, qe = qe, qi  # invert color metric for current clamp 
                        g = 100
                        r = np.clip(g + max(qi, 0), 0, 255)
                        b = np.clip(g + max(qe, 0), 0, 255)
                        color = (r, g, b)

                    plt.plot(t, segm, pen={'color': color, 'width': 1}, antialias=True)

                if self.params['show ticks']:
                    vt = pg.VTickGroup((on_times[j]-start) * dt, [0, 0.15], pen=0.4)
                    plt.addItem(vt)

                # Link all plots along x axis
                plt.setXLink(self.plots[0, 0])

                if i == j:
                    # link y axes of all diagonal plots
                    plt.setYLink(self.plots[0, 0])
                else:
                    # link y axes of all plots within a row
                    plt.setYLink(self.plots[i, (i+1) % 2])  # (i+1)%2 just avoids linking to 0,0

                if i < n_channels - 1:
                    plt.getAxis('bottom').setVisible(False)
                if j > 0:
                    plt.getAxis('left').setVisible(False)

                if i == n_channels - 1:
                    plt.setLabels(bottom=('CH%d'%chans[j], 's'))
                if j == 0:
                    plt.setLabels(left=('CH%d'%chans[i], 'A' if modes[i] == 'vc' else 'V'))

        if auto_range:
            r = 14e-12 if modes[i] == 'vc' else 5e-3
            self.plots[0, 1].setYRange(-r, r)
            r = 2e-9 if modes[i] == 'vc' else 100e-3
            self.plots[0, 0].setYRange(-r, r)

            self.plots[0, 0].setXRange(t[0], t[-1])

        