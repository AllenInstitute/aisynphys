# coding: utf8
from __future__ import print_function, division
import numpy as np
import pyqtgraph as pg
import pyqtgraph.multiprocess
from pyqtgraph.Qt import QtGui, QtCore
from aisynphys.stochastic_release_model import StochasticReleaseModel
from aisynphys.ui.ndslicer import NDSlicer


class ModelDisplayWidget(QtGui.QWidget):
    """UI containing an NDSlicer for visualizing the complete model parameter space, and
    a ModelSingleResultWidget for showing more detailed results from individual points in the parameter space.
    """
    def __init__(self, model_runner):
        QtGui.QWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)
        self.splitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.layout.addWidget(self.splitter)
        
        self.slicer = NDSlicer(model_runner.param_space.axes())
        self.slicer.selection_changed.connect(self.selection_changed)
        self.splitter.addWidget(self.slicer)
        
        self.result_widget = ModelSingleResultWidget()
        self.splitter.addWidget(self.result_widget)

        # set up a few default 2D slicer views
        v1 = self.slicer.params.child('2D views').addNew()
        v1['axis 0'] = 'n_release_sites'
        v1['axis 1'] = 'base_release_probability'
        v2 = self.slicer.params.child('2D views').addNew()
        v2['axis 0'] = 'vesicle_recovery_tau'
        v2['axis 1'] = 'facilitation_recovery_tau'
        self.slicer.dockarea.moveDock(v2.viewer.dock, 'bottom', v1.viewer.dock)
        v3 = self.slicer.params.child('2D views').addNew()
        v3['axis 0'] = 'vesicle_recovery_tau'
        v3['axis 1'] = 'base_release_probability'
        v4 = self.slicer.params.child('2D views').addNew()
        v4['axis 0'] = 'facilitation_amount'
        v4['axis 1'] = 'facilitation_recovery_tau'
        self.slicer.dockarea.moveDock(v4.viewer.dock, 'bottom', v3.viewer.dock)
        
        # turn on max projection for all parameters by default
        for ch in self.slicer.params.child('max project'):
            if ch.name() == 'synapse':
                continue
            ch.setValue(True)
        
        self.model_runner = model_runner
        self.param_space = model_runner.param_space
        
        result_img = np.zeros(self.param_space.result.shape)
        for ind in np.ndindex(result_img.shape):
            result_img[ind] = self.param_space.result[ind]['likelihood']
        self.slicer.set_data(result_img)
        self.results = result_img
        
        # select best result
        best = np.unravel_index(np.argmax(result_img), result_img.shape)
        self.select_result(best)

        max_like = self.results.max()
        
        # if results are combined across synapses, set up colors
        if 'synapse' in self.param_space.params:
            self.slicer.params['color axis', 'axis'] = 'synapse'
            syn_axis = list(self.param_space.params.keys()).index('synapse')
            max_img = np.array([result_img.take(i, axis=syn_axis).max() for i in range(result_img.shape[syn_axis])])
            max_like = max_img.min()
            max_img = max_img.min() / max_img
            syns = self.param_space.params['synapse']
            for i in syns:
                c = pg.colorTuple(pg.intColor(i, len(syns)*1.2))
                c = pg.mkColor(c[0]*max_img[i], c[1]*max_img[i], c[2]*max_img[i])
                self.slicer.params['color axis', 'colors', str(i)] = c
            
        # set histogram range
        self.slicer.histlut.setLevels(max_like * 0.85, max_like)

    def selection_changed(self, slicer):
        index = self.selected_index()
        self.select_result(index, update_slicer=False)

    def selected_index(self):
        return tuple(self.slicer.index().values())

    def get_result(self, index=None):
        if index is None:
            index = self.selected_index()
        return self.param_space.result[index]

    def select_result(self, index, update_slicer=True):
        result = self.get_result(index)
        result['params'].update(self.param_space[index])
        
        # re-run the model to get the complete results
        full_result = self.model_runner.run_model(result['params'], full_result=True, show=True)
        self.result_widget.set_result(full_result)
        
        print("----- Selected result: -----")
        print("  model parameters:")
        for k,v in full_result['params'].items():
            print("    {:30s}: {}".format(k, v))
        if 'optimized_params' in full_result:
            print("  optimized parameters:")
            for k,v in full_result['optimized_params'].items():
                print("    {:30s}: {}".format(k, v))
        if 'optimization_init' in full_result:
            print("  initial optimization parameters:")
            for k,v in full_result['optimization_init'].items():
                print("    {:30s}: {}".format(k, v))
        if 'optimization_result' in full_result:
            opt = full_result['optimization_result']
            print("  optimization results:")
            print("    nfev:", opt.nfev)
            print("    message:", opt.message)
            print("    success:", opt.success)
            print("    status:", opt.status)
        if 'optimization_info' in full_result:
            print("  optimization info:")
            for k,v in full_result['optimization_info'].items():
                print("    {:30s}: {}".format(k, v))
        print("  likelihood: {}".format(full_result['likelihood']))
        
        if update_slicer:
            self.slicer.set_index(index)


class ModelSingleResultWidget(QtGui.QWidget):
    """Plots event amplitudes and distributions for a single stochastic model run.
    """
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)

        self.splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.splitter)
        
        self.btn_layout = QtGui.QVBoxLayout()
        self.btn_widget = QtGui.QWidget()
        self.btn_widget.setLayout(self.btn_layout)
        self.splitter.addWidget(self.btn_widget)
        
        self.opt_path_plot = pg.PlotWidget()

        self.panels = {
            'events': ModelEventPlot(self),
            'correlation': ModelEventCorrelationPlot(self),
            'induction': ModelInductionPlot(self),
            'optimization': ModelOptimizationPlot(self),
        }
        for name, panel in self.panels.items():
            btn = QtGui.QPushButton(name)
            self.btn_layout.addWidget(btn)
            btn.setCheckable(True)
            btn.toggled.connect(panel.set_visible)
            btn.setMaximumWidth(150)
            panel.show_btn = btn
        
        self.panels['events'].show_btn.setChecked(True)

    def set_result(self, result):
        self.result = result
        for p in self.panels.values():
            p.result_changed()


class ModelResultView(object):
    """Displays one aspect of a model result.
    """
    def __init__(self, parent, visible=False):
        """Responsible for attaching widgets to parent.
        """
        self._parent = parent
        self._visible = False
        self.widget = QtGui.QWidget()
        self.layout = QtGui.QGridLayout()
        self.widget.setLayout(self.layout)
        parent.splitter.addWidget(self.widget)
        self.layout.setSpacing(2)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self._need_update = False
        
        self.set_visible(visible)
        
    def set_visible(self, visible):
        self._visible = visible
        self.widget.setVisible(visible)
        if visible and self._need_update:
            self.update_display()

    @property
    def visible(self):
        return self._visible

    def result_changed(self):
        """Called when the model result has been set on parent
        """
        self._need_update = True
        if self.visible:
            self.update_display()
        
    def update_display(self):
        """Extend in subclass to update displayed results (taken from self._parent.result)
        """
        self._need_update = False    
        

class ModelEventPlot(ModelResultView):
    def __init__(self, parent):
        ModelResultView.__init__(self, parent)

        self.event_view = pg.GraphicsLayoutWidget()
        self.layout.addWidget(self.event_view)

        self.plots = {
            'likelihood': self.event_view.addPlot(0, 0, title="model likelihood vs compressed time"),
            'amplitude': self.event_view.addPlot(1, 0, title="event amplitude vs compressed time"),
        }        
        self.plots['amplitude'].setXLink(self.plots['likelihood'])
        self.state_keys = ['release_probability', 'facilitation', 'sensitization']
        for i,state_key in enumerate(self.state_keys):
            self.plots[state_key] = self.event_view.addPlot(2+i, 0, title=state_key + " vs compressed time")
            self.plots[state_key].setXLink(self.plots['likelihood'])
        
        self.amp_dist_plot = self.event_view.addPlot(0, 1, title="amplitude distributions", rowspan=3)
        self.amp_dist_plot.setMaximumWidth(500)
        self.amp_dist_plot.selected_items = []

        self.amp_sample_values = np.linspace(-0.005, 0.005, 800)

        self.ctrl = QtGui.QWidget()
        self.hl = QtGui.QHBoxLayout()
        self.hl.setSpacing(2)
        self.hl.setContentsMargins(0, 0, 0, 0)
        self.ctrl.setLayout(self.hl)
        self.layout.addWidget(self.ctrl)

        self.plot_checks = {
            'likelihood': QtGui.QCheckBox('likelihood'),
            'amplitude': QtGui.QCheckBox('amplitude'),
            'release_probability': QtGui.QCheckBox('release probability'),
            'facilitation': QtGui.QCheckBox('facilitation'),
            'sensitization': QtGui.QCheckBox('sensitization'),
        }
        self.plot_checks['amplitude'].setChecked(True)
        for name,c in self.plot_checks.items():
            self.hl.addWidget(c)
            c.toggled.connect(self.update_display)

    def update_display(self):
        ModelResultView.update_display(self)

        for k in self.plots:
            if self.plot_checks[k].isChecked():
                self.plots[k].setVisible(True)
                self.plots[k].setMaximumHeight(10000)
            else:
                self.plots[k].setVisible(False)
                self.plots[k].setMaximumHeight(0)

        full_result = self._parent.result
        model = full_result['model']
        result = full_result['result']
        pre_state = full_result['pre_spike_state']
        post_state = full_result['post_spike_state']
        
        # color events by likelihood
        cmap = pg.ColorMap([0, 1.0], [(0, 0, 0), (255, 0, 0)])
        threshold = 10
        err_colors = cmap.map((threshold - result['likelihood']) / threshold)
        brushes = [pg.mkBrush(c) if np.isfinite(result['likelihood'][i]) else pg.mkBrush(None) for i,c in enumerate(err_colors)]

        # log spike intervals to make visualization a little easier
        compressed_spike_times = np.empty(len(result['spike_time']))
        compressed_spike_times[0] = 0.0
        np.cumsum(np.diff(result['spike_time'])**0.25, out=compressed_spike_times[1:])

        if self.plot_checks['likelihood'].isChecked():
            self.plots['likelihood'].clear()
            self.plots['likelihood'].plot(compressed_spike_times, result['likelihood'], pen=None, symbol='o', symbolBrush=brushes)
        
        if self.plot_checks['amplitude'].isChecked():
            self.plots['amplitude'].clear()
            self.plots['amplitude'].plot(compressed_spike_times, result['expected_amplitude'], pen=None, symbol='x', symbolPen=0.5, symbolBrush=brushes)
            amp_sp = self.plots['amplitude'].plot(compressed_spike_times, result['amplitude'], pen=None, symbol='o', symbolBrush=brushes)
            amp_sp.scatter.sigClicked.connect(self.amp_sp_clicked)

        for k in self.state_keys:
            if not self.plot_checks[k].isChecked():
                continue
            self.plots[k].clear()
            self.plots[k].plot(compressed_spike_times, pre_state[k], pen=None, symbol='t', symbolBrush=brushes)
            self.plots[k].plot(compressed_spike_times, post_state[k], pen=None, symbol='o', symbolBrush=brushes)

        # plot full distribution of event amplitudes
        self.amp_dist_plot.clear()
        bins = np.linspace(np.nanmin(result['amplitude']), np.nanmax(result['amplitude']), 40)
        d_amp = bins[1] - bins[0]
        amp_hist = np.histogram(result['amplitude'], bins=bins)
        self.amp_dist_plot.plot(amp_hist[1], amp_hist[0] / (amp_hist[0].sum() * d_amp), stepMode=True, fillLevel=0, brush=0.3)

        # plot average model event distribution
        amps = self.amp_sample_values
        d_amp = amps[1] - amps[0]
        total_dist = np.zeros(len(amps))
        for i in range(result.shape[0]):
            state = pre_state[i]
            if not np.all(np.isfinite(tuple(state))):
                continue
            total_dist += model.likelihood(amps, state)
        total_dist /= total_dist.sum() * d_amp
        self.amp_dist_plot.plot(amps, total_dist, fillLevel=0, brush=(255, 0, 0, 50))
    
    def amp_sp_clicked(self, sp, pts):
        result = self._parent.result['result']
        i = pts[0].index()
        state = self.pre_state[i]
        expected_amp = result[i]['expected_amplitude']
        measured_amp = result[i]['amplitude']
        amps = self.amp_sample_values
        
        for item in self.amp_dist_plot.selected_items:
            self.amp_dist_plot.removeItem(item)
        l = self.model.likelihood(amps, state)
        p = self.amp_dist_plot.plot(amps, l / l.sum(), pen=(255, 255, 0, 100))
        l1 = self.amp_dist_plot.addLine(x=measured_amp)
        l2 = self.amp_dist_plot.addLine(x=expected_amp, pen='r')
        self.amp_dist_plot.selected_items = [p, l1, l2]


class ModelOptimizationPlot(ModelResultView):
    """Plots the mini amplitude optimization path
    """
    def __init__(self, parent):
        ModelResultView.__init__(self, parent)
        self.plot = pg.PlotWidget()
        self.layout.addWidget(self.plot)
        
    def update_display(self):
        ModelResultView.update_display(self)

        result = self._parent.result
        plt = self.plot
        x = result['optimization_path']['mini_amplitude']
        y = result['optimization_path']['likelihood']
        brushes = [pg.mkBrush((i, int(len(x)*1.2))) for i in range(len(x))]
        plt.clear()
        plt.plot(x, y, pen=None, symbol='o', symbolBrush=brushes)
        plt.addLine(x=result['optimized_params']['mini_amplitude'])
        plt.addLine(y=result['likelihood'])


class ModelEventCorrelationPlot(ModelResultView):
    """Show correlation in amplitude between adjacent events. 

    The motivation here is that if synaptic release causes vesicle depletion, then the amplitude of
    two adjacent events in a train should be anti-correlated (a large release on one event causes more vesicle
    depletion and therefore more depression; the following event should be smaller). 
    """
    def __init__(self, parent):
        ModelResultView.__init__(self, parent)
        self.plot = pg.PlotWidget(labels={'left': 'amp 2', 'bottom': 'amp 1'})
        self.plot.showGrid(True, True)
        self.plot.setAspectLocked()
        self.layout.addWidget(self.plot)

        self.ctrl = QtGui.QWidget()
        self.hl = QtGui.QHBoxLayout()
        self.hl.setSpacing(2)
        self.hl.setContentsMargins(0, 0, 0, 0)
        self.ctrl.setLayout(self.hl)
        self.layout.addWidget(self.ctrl)

        self.mode_radios = {
            'all_events': QtGui.QRadioButton('all events'),
            'first_in_train': QtGui.QRadioButton('first in train'),
        }
        self.mode_radios['all_events'].setChecked(True)
        for name,r in self.mode_radios.items():
            self.hl.addWidget(r)
            r.toggled.connect(self.update_display)
        
    def update_display(self):
        ModelResultView.update_display(self)
        self.plot.clear()

        result = self._parent.result
        spikes = result['result']['spike_time']
        amps = result['result']['amplitude']

        if self.mode_radios['all_events'].isChecked():
            cmap = pg.ColorMap(np.linspace(0, 1, 4), np.array([[255, 255, 255], [255, 255, 0], [255, 0, 0], [0, 0, 0]], dtype='ubyte'))
            cvals = [((np.log(dt) / np.log(10)) + 2) / 4. for dt in np.diff(spikes)]
            brushes = [pg.mkBrush(cmap.map(c)) for c in cvals]
            x = amps[:-1]
            y = amps[1:]
        else:
            brushes = pg.mkBrush('w')

            # find all first+second pulse amps that come from the same recording
            x_mask = result['event_meta']['pulse_number'] == 1
            y_mask = result['event_meta']['pulse_number'] == 2
            x_ids = result['event_meta']['sync_rec_ext_id'][x_mask]
            y_ids = result['event_meta']['sync_rec_ext_id'][y_mask]
            common_ids = list(set(x_ids) & set(y_ids))
            id_mask = np.isin(result['event_meta']['sync_rec_ext_id'], common_ids)

            x = amps[x_mask & id_mask]
            y = amps[y_mask & id_mask]

        self.plot.plot(x, y, pen=None, symbol='o', symbolBrush=brushes)


class ModelInductionPlot(ModelResultView):
    def __init__(self, parent):
        ModelResultView.__init__(self, parent)
        self.lw = pg.GraphicsLayoutWidget()
        self.ind_plots = [self.lw.addPlot(0, 0), self.lw.addPlot(1, 0), self.lw.addPlot(2, 0)]
        self.rec_plot = self.lw.addPlot(0, 1)
        self.corr_plot = self.lw.addPlot(1, 1, rowspan=2)
        # self.corr_plot.setAspectLocked()
        self.corr_plot.showGrid(True, True)
        self.layout.addWidget(self.lw)
        
    def update_display(self):
        ModelResultView.update_display(self)
        result = self._parent.result
        spikes = result['result']['spike_time']
        amps = result['result']['amplitude']
        meta = result['event_meta']
        
        self.corr_plot.clear()
        
        # generate a list of all trains sorted by stimulus
        trains = {}  # {ind_f: {rec_d: [[a1, a2, ..a12], [b1, b2, ..b12], ...], ...}, ...}
        current_sweep = None
        skip_sweep = False
        current_train = []
        for i in range(len(amps)):
            sweep_id = meta['sync_rec_ext_id'][i]
            ind_f = meta['induction_frequency'][i]
            rec_d = meta['recovery_delay'][i]
            if sweep_id != current_sweep:
                skip_sweep = False
                current_sweep = sweep_id
                current_train = []
                ind_trains = trains.setdefault(ind_f, {})
                rec_trains = ind_trains.setdefault(rec_d, [])
                rec_trains.append(current_train)
            if skip_sweep:
                continue
            if not np.isfinite(amps[i]) or not np.isfinite(spikes[i]):
                skip_sweep = True
                continue
            current_train.append(amps[i])

        # scatter plots of event amplitudes sorted by pulse number
        for ind_i, ind_f in enumerate([20, 50, 100]):
            ind_trains = trains.get(ind_f, {})
            
            # collect all induction events by pulse number
            ind_pulses = [[] for i in range(12)]
            for rec_d, rec_trains in ind_trains.items():
                for train in rec_trains:
                    for i,amp in enumerate(train):
                        ind_pulses[i].append(amp)
                        
            x = []
            y = []
            for i in range(12):
                if len(ind_pulses[i]) == 0:
                    continue
                y.extend(ind_pulses[i])
                xs = pg.pseudoScatter(np.array(ind_pulses[i]), bidir=True, shuffle=True)
                xs /= np.abs(xs).max() * 4
                x.extend(xs + i)

            self.ind_plots[ind_i].clear()
            self.ind_plots[ind_i].plot(x, y, pen=None, symbol='o')
            
            # re-model based on mean amplitudes
            mean_times = np.arange(12) / ind_f
            mean_times[8:] += 0.25
            model = result['model']
            params = result['params'].copy()
            params.update(result['optimized_params'])
            mean_result = model.measure_likelihood(mean_times, amplitudes=None, params=params)
            
            expected_amps = mean_result['result']['expected_amplitude']
            self.ind_plots[ind_i].plot(expected_amps, pen='w', symbol='d', symbolBrush='y')
        
            # normalize events by model prediction
            x = []
            y = []
            for rec_d, rec_trains in ind_trains.items():
                for train in rec_trains:
                    train = [t - expected_amps[i] for i,t in enumerate(train)]
                    for i in range(1, len(train)):
                        x.append(train[i-1])
                        y.append(train[i])
            x = np.array(x)
            y = np.array(y)
            
            y1 = y[x<0]
            y2 = y[x>0]
            x1 = pg.pseudoScatter(y1, bidir=True)
            x2 = pg.pseudoScatter(y2, bidir=True)
            x1 = 0.25 * x1 / x1.max()
            x2 = 0.25 * x2 / x2.max()
            self.corr_plot.plot(x1, y1, pen=None, symbol='o')
            self.corr_plot.plot(x2 + 1, y2, pen=None, symbol='o')
            # self.corr_plot.plot(x, y, pen=None, symbol='o', symbolBrush=(ind_i, 4))
            
        # scatter plot of event pairs normalized by model expectation
        
        