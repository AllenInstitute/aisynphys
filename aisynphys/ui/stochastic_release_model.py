# coding: utf8
from __future__ import print_function, division
import numpy as np
import scipy.stats
import pyqtgraph as pg
import pyqtgraph.multiprocess
from pyqtgraph.Qt import QtWidgets, QtCore
from aisynphys.stochastic_release_model import StochasticReleaseModel
from aisynphys.ui.ndslicer import NDSlicer


data_color = (0, 128, 255)
model_color = (255, 128, 0)


class ModelDisplayWidget(QtWidgets.QWidget):
    """UI containing an NDSlicer for visualizing the complete model parameter space, and
    a ModelSingleResultWidget for showing more detailed results from individual points in the parameter space.
    """
    def __init__(self, model_runner):
        QtWidgets.QWidget.__init__(self)
        self.layout = QtWidgets.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)
        
        self.slicer = NDSlicer(model_runner.param_space.axes())
        self.slicer.selection_changed.connect(self.selection_changed)
        self.layout.addWidget(self.slicer)
        
        # set up a few default 2D slicer views
        v1 = self.slicer.params.child('2D views').addNew()
        v1['axis 0'] = 'n_release_sites'
        v1['axis 1'] = 'base_release_probability'
        v2 = self.slicer.params.child('2D views').addNew()
        v2['axis 0'] = 'depression_amount'
        v2['axis 1'] = 'base_release_probability'
        self.slicer.dockarea.moveDock(v2.viewer.dock, 'bottom', v1.viewer.dock)
        v3 = self.slicer.params.child('2D views').addNew()
        v3['axis 0'] = 'depression_tau'
        v3['axis 1'] = 'base_release_probability'
        v4 = self.slicer.params.child('2D views').addNew()
        v4['axis 0'] = 'facilitation_amount'
        v4['axis 1'] = 'base_release_probability'
        self.slicer.dockarea.moveDock(v4.viewer.dock, 'bottom', v3.viewer.dock)

        self.result_widget = ModelSingleResultWidget()
        self.result_dock = pg.dockarea.Dock('model results')
        self.result_dock.addWidget(self.result_widget)
        self.slicer.dockarea.addDock(self.result_dock, 'bottom')

        # turn on max projection for all parameters by default
        for ch in self.slicer.params.child('max project'):
            if ch.name() == 'synapse':
                continue
            ch.setValue(True)
        
        self.model_runner = model_runner
        self.param_space = model_runner.param_space
        
        # result_img = np.zeros(self.param_space.result.shape)
        # for ind in np.ndindex(result_img.shape):
        #     result_img[ind] = self.param_space.result[ind].likelihood
        result_img = self.param_space.result['likelihood']
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
        params = self.param_space.params_at_index(index)
        
        # re-run the model to get the complete results
        full_result = self.model_runner.run_model(params, show=True)
        self.result_widget.set_result(full_result)
        
        print("----- Selected result: -----")
        print("  model parameters:")
        for k,v in full_result.params.items():
            print("    {:30s}: {}".format(k, v))
        if hasattr(full_result, 'optimized_params'):
            print("  optimized parameters:")
            for k,v in full_result.optimized_params.items():
                print("    {:30s}: {}".format(k, v))
        if hasattr(full_result, 'optimization_init'):
            print("  initial optimization parameters:")
            for k,v in full_result.optimization_init.items():
                print("    {:30s}: {}".format(k, v))
        if hasattr(full_result, 'optimization_result'):
            opt = full_result.optimization_result
            print("  optimization results:")
            print("    nfev:", opt.nfev)
            print("    message:", opt.message)
            print("    success:", opt.success)
            print("    status:", opt.status)
        if hasattr(full_result, 'optimization_info'):
            print("  optimization info:")
            for k,v in full_result.optimization_info.items():
                print("    {:30s}: {}".format(k, v))
        print("  likelihood: {}".format(full_result.likelihood))
        
        if update_slicer:
            self.slicer.set_index(index)


class ModelSingleResultWidget(QtWidgets.QWidget):
    """Plots event amplitudes and distributions for a single stochastic model run.
    """
    def __init__(self):
        self.result = None
        self._random_model_result = None

        QtWidgets.QWidget.__init__(self)
        self.layout = QtWidgets.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)

        self.splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.splitter)
        
        self.btn_layout = QtWidgets.QVBoxLayout()
        self.btn_widget = QtWidgets.QWidget()
        self.btn_widget.setLayout(self.btn_layout)
        self.splitter.addWidget(self.btn_widget)
        
        self.opt_path_plot = pg.PlotWidget()

        self.panels = {
            'events': ModelEventPlot(self),
            'induction': ModelInductionPlot(self),
            'correlation': ModelEventCorrelationPlot(self),
            'optimization': ModelOptimizationPlot(self),
        }
        for name, panel in self.panels.items():
            btn = QtWidgets.QPushButton(name)
            self.btn_layout.addWidget(btn)
            btn.setCheckable(True)
            btn.toggled.connect(panel.set_visible)
            btn.setMaximumWidth(150)
            panel.show_btn = btn
        
        self.panels['events'].show_btn.setChecked(True)

    def set_result(self, result):
        self.result = result
        self._random_model_result = None
        for p in self.panels.values():
            p.result_changed()

    def random_model_result(self):
        """Re-run the model with the same parameters, generating a random sample of event amplitudes.

        The random sample is cached until set_result is called again.
        """
        if self._random_model_result is None:
            self._random_model_result = self.result.model.run_model(
                self.result.result['spike_time'], 
                amplitudes='random', 
                params=self.result.all_params)
        return self._random_model_result


class ModelResultView(object):
    """Displays one aspect of a model result.
    """
    def __init__(self, parent, visible=False):
        """Responsible for attaching widgets to parent.
        """
        self._parent = parent
        self._visible = False
        self.widget = QtWidgets.QWidget()
        self.layout = QtWidgets.QGridLayout()
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
    def parent(self):
        return self._parent

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
        """Extend in subclass to update displayed results (taken from self.parent.result)
        """
        self._need_update = False    
        

class ModelEventPlot(ModelResultView):
    """Plot event amplitudes, likelihoods, and model state vs time.

    Clicking on individual points shows a plot of the model predicted amplitude distribution for that event.
    """
    def __init__(self, parent):
        ModelResultView.__init__(self, parent)

        self.event_view = pg.GraphicsLayoutWidget()
        self.layout.addWidget(self.event_view)

        self.plots = {
            'amplitude': self.event_view.addPlot(0, 0, title="event amplitude vs compressed time"),
            'likelihood': self.event_view.addPlot(1, 0, title="model likelihood vs compressed time"),
        }
        self.state_keys = ['vesicle_pool', 'release_probability', 'depression', 'facilitation']
        for i,state_key in enumerate(self.state_keys):
            self.plots[state_key] = self.event_view.addPlot(2+i, 0, title=state_key + " vs compressed time")
        
        self.amp_dist_plot = self.event_view.addPlot(0, 1, title="amplitude distributions", rowspan=10)
        self.amp_dist_plot.setMaximumWidth(500)
        self.amp_dist_plot.selected_items = []

        self.amp_sample_values = np.linspace(-0.005, 0.005, 800)

        self.stim_regions = []

        self.ctrl = QtWidgets.QWidget()
        self.hl = QtWidgets.QHBoxLayout()
        self.hl.setSpacing(2)
        self.hl.setContentsMargins(0, 0, 0, 0)
        self.ctrl.setLayout(self.hl)
        self.layout.addWidget(self.ctrl)

        self.plot_checks = {
            'amplitude': QtWidgets.QCheckBox('amplitude'),
            'likelihood': QtWidgets.QCheckBox('likelihood'),
            'vesicle_pool': QtWidgets.QCheckBox('vesicle pool'),
            'release_probability': QtWidgets.QCheckBox('release probability'),
            'depression': QtWidgets.QCheckBox('depression'),
            'facilitation': QtWidgets.QCheckBox('facilitation'),
        }
        self.plot_checks['amplitude'].setChecked(True)
        for name,c in self.plot_checks.items():
            self.hl.addWidget(c)
            c.toggled.connect(self.update_display)

        for name,plot in self.plots.items():
            plot.setXLink(self.plots['amplitude'])

    def update_display(self):
        """Called when we need to update displayed results (from self.parent.result)
        """
        ModelResultView.update_display(self)

        for k in self.plots:
            if self.plot_checks[k].isChecked():
                self.plots[k].setVisible(True)
                self.plots[k].setMaximumHeight(10000)
            else:
                self.plots[k].setVisible(False)
                self.plots[k].setMaximumHeight(0)

        full_result = self.parent.result
        model = full_result.model
        result = full_result.result
        pre_state = full_result.pre_spike_state
        post_state = full_result.post_spike_state
        
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

            # show regions for each type of stimulus
            last_stim_name = None
            rgns = []
            for i, stim_name in enumerate(full_result.event_meta['stim_name']):
                ev_time = compressed_spike_times[i]
                if stim_name != last_stim_name:
                    rgns.append([stim_name, ev_time, ev_time])
                    last_stim_name = stim_name
                else:
                    rgns[-1][1] = ev_time
            for stim_name, start_time, stop_time in rgns:
                rgn = pg.LinearRegionItem(values=[start_time, stop_time], orientation='vertical', brush=(len(self.stim_regions), 10), pen=None, movable=False)
                rgn.lines[0].label = pg.InfLineLabel(rgn.lines[1], stim_name, movable=False, rotateAxis=(1, 0), position=0, anchors=[(0, 0), (0, 0)])
                rgn.setZValue(-10)
                rgn.setOpacity(0.3)
                self.plots['amplitude'].addItem(rgn)
                self.stim_regions.append(rgn)

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
        # user clicked an event in the scatter plot; show: 
        # - the amplitude of the selected event (red line)
        # - expectation value from model distribution (yellow line)
        # - the full model distribution (yellow histogram)

        result = self.parent.result.result
        pre_state = self.parent.result.pre_spike_state
        model = self.parent.result.model
        i = pts[0].index()
        state = pre_state[i]
        expected_amp = result[i]['expected_amplitude']
        measured_amp = result[i]['amplitude']
        amps = self.amp_sample_values
        
        for item in self.amp_dist_plot.selected_items:
            self.amp_dist_plot.removeItem(item)
        l = model.likelihood(amps, state)
        p = self.amp_dist_plot.plot(amps, l, pen=(255, 255, 0, 100))
        l1 = self.amp_dist_plot.addLine(x=measured_amp, pen='r')
        l2 = self.amp_dist_plot.addLine(x=expected_amp, pen='y')
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
        x = result.optimization_path['mini_amplitude']
        y = result.optimization_path['likelihood']
        brushes = [pg.mkBrush((i, int(len(x)*1.2))) for i in range(len(x))]
        plt.clear()
        plt.plot(x, y, pen=None, symbol='o', symbolBrush=brushes)
        plt.addLine(x=result.optimized_params['mini_amplitude'])
        plt.addLine(y=result.likelihood)


class ModelEventCorrelationPlot(ModelResultView):
    """Show correlation in amplitude between adjacent events. 

    The motivation here is that if synaptic release causes vesicle depletion, then the amplitude of
    two adjacent events in a train should be anti-correlated (a large release on one event causes more vesicle
    depletion and therefore more depression; the following event should be smaller). 
    """
    def __init__(self, parent):
        ModelResultView.__init__(self, parent)
        self.view = pg.GraphicsLayoutWidget()
        self.layout.addWidget(self.view)

        self.plots = [[self.view.addPlot(row=i, col=0, labels={'left': '%d:%d'%(2**i, 2**(i+1))}), self.view.addPlot(row=i, col=1)] for i in range(3)]
        for row in self.plots:
            for plot in row:
                plot.showGrid(True, True)
                plot.label = pg.TextItem()
                plot.label.setParentItem(plot.vb)
            row[1].setXLink(row[0])
            row[1].setYLink(row[0])
        self.plots[0][0].setTitle('data')
        self.plots[0][1].setTitle('model')

        # self.ctrl = QtWidgets.QWidget()
        # self.hl = QtWidgets.QHBoxLayout()
        # self.hl.setSpacing(2)
        # self.hl.setContentsMargins(0, 0, 0, 0)
        # self.ctrl.setLayout(self.hl)
        # self.layout.addWidget(self.ctrl)

        # self.mode_radios = {
        #     'all_events': QtWidgets.QRadioButton('all events'),
        #     'first_in_train': QtWidgets.QRadioButton('first in train'),
        # }
        # self.mode_radios['all_events'].setChecked(True)
        # for name,r in self.mode_radios.items():
        #     self.hl.addWidget(r)
        #     r.toggled.connect(self.update_display)
        
    def update_display(self):
        ModelResultView.update_display(self)
        for row in self.plots:
            for plot in row:
                plot.clear()

        result = self.parent.result
        spikes = result.result['spike_time']
        amps = result.result['amplitude']

        # re-simulate one random trial 
        model_result = self.parent.random_model_result()
        model_amps = model_result.result['amplitude']

        # this analysis relies on repeated structures in the stimulus, so we need to reconstruct
        # these from the event metadata
        recs = result.events_by_recording(require_contiguous=True)

        for i in range(3):
            n = (i + 1) * 2
            # find all recordings containing all of the pulses we want to analyze
            event_inds = []
            for rec in recs:
                if len(rec) < n:
                    continue
                event_inds.append(rec[:n])

            for j,this_amps in enumerate([amps, model_amps]):
                selected_amps = this_amps[np.array(event_inds)]

                x = selected_amps[:, :n//2].mean(axis=1)
                y = selected_amps[:, n//2:].mean(axis=1)

                self.plots[i][j].plot(x, y, pen=None, symbol='o', symbolPen=None, symbolBrush=[data_color, model_color][j])

                # plot regression
                slope, intercept, r, p, stderr = scipy.stats.linregress(x, y)
                assert np.isfinite(p)
                x = np.array([x.min(), x.max()])
                y = slope * x + intercept
                self.plots[i][j].plot(x, y, pen='y')
                self.plots[i][j].label.setText(f'r={r:0.2f} p={p:0.3f}')


class ModelInductionPlot(ModelResultView):
    def __init__(self, parent):
        ModelResultView.__init__(self, parent)
        self.lw = pg.GraphicsLayoutWidget()
        self.induction_freqs = [10, 20, 50, 100]
        self.ind_plots = [self.lw.addPlot(i, 0, labels={'left': f'{freq:d}Hz'}) for i,freq in enumerate(self.induction_freqs)]
        self.layout.addWidget(self.lw)
        
    def update_display(self):
        ModelResultView.update_display(self)
        result = self._parent.result
        spikes = result.result['spike_time']
        amps = result.result['amplitude']
        meta = result.event_meta

        # re-simulate one random trial 
        model_result = self.parent.random_model_result()
        model_amps = model_result.result['amplitude']
        
        # generate a list of all trains sorted by stimulus
        trains = result.events_by_stimulus()

        # scatter plots of event amplitudes sorted by pulse number
        for ind_i, ind_f in enumerate(self.induction_freqs):
            ind_trains = trains.get(ind_f, {})
            
            # collect all induction events by pulse number
            ind_pulses = [[] for i in range(12)]
            for rec_d, rec_trains in ind_trains.items():
                for train in rec_trains:
                    for i,ev_ind in enumerate(train):
                        ind_pulses[i].append(ev_ind)
            
            real_x = []
            real_y = []
            real_avg_y = []
            model_x = []
            model_y = []
            model_avg_y = []
            for i in range(12):
                if len(ind_pulses[i]) == 0:
                    real_avg_y = [np.nan] * 12
                    continue

                inds = np.array(ind_pulses[i])
                for this_amp, x, avg_y, y, sign in ((amps, real_x, real_avg_y, real_y, -1), (model_amps, model_x, model_avg_y, model_y, 1)):
                    amp = this_amp[inds]
                    y.extend(amp)
                    if len(amp) == 0:
                        avg_y.append(np.nan)
                    else:
                        avg_y.append(amp.mean())                    
                    xs = pg.pseudoScatter(np.array(amp), bidir=False, shuffle=True)
                    xs /= np.abs(xs).max() * 4
                    x.extend(i + xs * sign)

            self.ind_plots[ind_i].clear()
            self.ind_plots[ind_i].plot(real_x, real_y, pen=None, symbol='o', symbolPen=None, symbolBrush=data_color+(200,), symbolSize=5)
            self.ind_plots[ind_i].plot(np.arange(12), real_avg_y, pen=data_color+(100,), symbol='d', symbolPen=None, symbolBrush=data_color+(100,), symbolSize=5)
            self.ind_plots[ind_i].plot(np.array(model_x)+0.1, model_y, pen=None, symbol='o', symbolPen=None, symbolBrush=model_color+(200,), symbolSize=5, zValue=-1)
            
            # re-model based on mean amplitudes
            mean_times = np.arange(12) / ind_f
            mean_times[8:] += 0.25
            model = result.model
            params = result.all_params
            mean_result = model.run_model(mean_times, amplitudes='expected', params=params)
            
            # plot model distribution expectation values
            expected_amps = mean_result.result['expected_amplitude']
            self.ind_plots[ind_i].plot(expected_amps, pen=model_color+(100,), symbol='d', symbolBrush=model_color+(100,), symbolPen=None, zValue=-10)
