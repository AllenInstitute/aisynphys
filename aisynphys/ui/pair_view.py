import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets, QtCore
from collections import OrderedDict
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.ui.filter import SignalFilter, ArtifactRemover
from neuroanalysis.ui.baseline import BaselineRemover
from neuroanalysis.ui.fitting import FitExplorer
from neuroanalysis.data import TSeries
from neuroanalysis.spike_detection import detect_evoked_spikes
from neuroanalysis import fitting
from neuroanalysis.baseline import float_mode
from neuroanalysis.stats import ragged_mean


class PairView(QtWidgets.QWidget):
    """For analyzing pre/post-synaptic pairs.
    """
    def __init__(self, parent=None):
        self.sweeps = []
        self.channels = []

        self.current_event_set = None
        self.event_sets = []

        QtWidgets.QWidget.__init__(self, parent)

        self.layout = QtWidgets.QGridLayout()
        self.setLayout(self.layout)
        self.layout.setContentsMargins(0, 0, 0, 0)

        self.vsplit = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.layout.addWidget(self.vsplit, 0, 0)
        
        self.pre_plot = pg.PlotWidget()
        self.post_plot = pg.PlotWidget()
        for plt in (self.pre_plot, self.post_plot):
            plt.setClipToView(True)
            plt.setDownsampling(True, True, 'peak')
        
        self.post_plot.setXLink(self.pre_plot)
        self.vsplit.addWidget(self.pre_plot)
        self.vsplit.addWidget(self.post_plot)
        
        self.response_plots = PlotGrid()
        self.vsplit.addWidget(self.response_plots)

        self.event_splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.vsplit.addWidget(self.event_splitter)
        
        self.event_table = pg.TableWidget()
        self.event_splitter.addWidget(self.event_table)
        
        self.event_widget = QtWidgets.QWidget()
        self.event_widget_layout = QtWidgets.QGridLayout()
        self.event_widget.setLayout(self.event_widget_layout)
        self.event_splitter.addWidget(self.event_widget)
        
        self.event_splitter.setSizes([600, 100])
        
        self.event_set_list = QtWidgets.QListWidget()
        self.event_widget_layout.addWidget(self.event_set_list, 0, 0, 1, 2)
        self.event_set_list.addItem("current")
        self.event_set_list.itemSelectionChanged.connect(self.event_set_selected)
        
        self.add_set_btn = QtWidgets.QPushButton('add')
        self.event_widget_layout.addWidget(self.add_set_btn, 1, 0, 1, 1)
        self.add_set_btn.clicked.connect(self.add_set_clicked)
        self.remove_set_btn = QtWidgets.QPushButton('del')
        self.event_widget_layout.addWidget(self.remove_set_btn, 1, 1, 1, 1)
        self.remove_set_btn.clicked.connect(self.remove_set_clicked)
        self.fit_btn = QtWidgets.QPushButton('fit all')
        self.event_widget_layout.addWidget(self.fit_btn, 2, 0, 1, 2)
        self.fit_btn.clicked.connect(self.fit_clicked)

        self.fit_plot = PlotGrid()

        self.fit_explorer = None

        self.artifact_remover = ArtifactRemover(user_width=True)
        self.baseline_remover = BaselineRemover()
        self.filter = SignalFilter()
        
        self.params = pg.parametertree.Parameter(name='params', type='group', children=[
            {'name': 'pre', 'type': 'list', 'values': []},
            {'name': 'post', 'type': 'list', 'values': []},
            self.artifact_remover.params,
            self.baseline_remover.params,
            self.filter.params,
            {'name': 'time constant', 'type': 'float', 'suffix': 's', 'siPrefix': True, 'value': 10e-3, 'dec': True, 'minStep': 100e-6}
            
        ])
        self.params.sigTreeStateChanged.connect(self._update_plots)

    def data_selected(self, sweeps, channels):
        self.sweeps = sweeps
        self.channels = channels
        
        self.params.child('pre').setLimits(channels)
        self.params.child('post').setLimits(channels)
        
        self._update_plots()

    def _update_plots(self):
        sweeps = self.sweeps
        self.current_event_set = None
        self.event_table.clear()
        
        # clear all plots
        self.pre_plot.clear()
        self.post_plot.clear()

        pre = self.params['pre']
        post = self.params['post']
        
        # If there are no selected sweeps or channels have not been set, return
        if len(sweeps) == 0 or pre == post or pre not in self.channels or post not in self.channels:
            return

        pre_mode = sweeps[0][pre].clamp_mode
        post_mode = sweeps[0][post].clamp_mode
        for ch, mode, plot in [(pre, pre_mode, self.pre_plot), (post, post_mode, self.post_plot)]:
            units = 'A' if mode == 'vc' else 'V'
            plot.setLabels(left=("Channel %d" % ch, units), bottom=("Time", 's'))
        
        # Iterate over selected channels of all sweeps, plotting traces one at a time
        # Collect information about pulses and spikes
        pulses = []
        spikes = []
        post_traces = []
        for i,sweep in enumerate(sweeps):
            pre_trace = sweep[pre]['primary']
            post_trace = sweep[post]['primary']
            
            # Detect pulse times
            stim = sweep[pre]['command'].data
            sdiff = np.diff(stim)
            on_times = np.argwhere(sdiff > 0)[1:, 0]  # 1: skips test pulse
            off_times = np.argwhere(sdiff < 0)[1:, 0]
            pulses.append(on_times)

            # filter data
            post_filt = self.artifact_remover.process(post_trace, list(on_times) + list(off_times))
            post_filt = self.baseline_remover.process(post_filt)
            post_filt = self.filter.process(post_filt)
            post_traces.append(post_filt)
            
            # plot raw data
            color = pg.intColor(i, hues=len(sweeps)*1.3, sat=128)
            color.setAlpha(128)
            for trace, plot in [(pre_trace, self.pre_plot), (post_filt, self.post_plot)]:
                plot.plot(trace.time_values, trace.data, pen=color, antialias=False)

            # detect spike times
            spike_times = []
            spike_info = []
            for on, off in zip(on_times, off_times):
                spikes = detect_evoked_spikes(sweep[pre], [on, off])
                spike_info.append(spikes)
                for spike in spikes:
                    if spike['max_slope_time'] is None:
                        continue
                    spike_times.append(spike['max_slope_time'])
            spikes.append(spike_info)
                    
            dt = pre_trace.dt
            vticks = pg.VTickGroup(spike_times, yrange=[0.0, 0.2], pen=color)
            self.pre_plot.addItem(vticks)

        # Iterate over spikes, plotting average response
        all_responses = []
        avg_responses = []
        fits = []
        fit = None
        
        npulses = max(map(len, pulses))
        self.response_plots.clear()
        self.response_plots.set_shape(1, npulses+1) # 1 extra for global average
        self.response_plots.setYLink(self.response_plots[0,0])
        for i in range(1, npulses+1):
            self.response_plots[0,i].hideAxis('left')
        units = 'A' if post_mode == 'vc' else 'V'
        self.response_plots[0, 0].setLabels(left=("Averaged events (Channel %d)" % post, units))
        
        fit_pen = {'color':(30, 30, 255), 'width':2, 'dash': [1, 1]}
        for i in range(npulses):
            # get the chunk of each sweep between spikes
            responses = []
            all_responses.append(responses)
            for j, sweep in enumerate(sweeps):
                # get the current spike
                if i >= len(spikes[j]):
                    continue
                spike = spikes[j][i]
                if spike is None:
                    continue
                
                # find next spike
                next_spike = None
                for sp in spikes[j][i+1:]:
                    if sp is not None:
                        next_spike = sp
                        break
                    
                # determine time range for response
                max_len = int(40e-3 / dt)  # don't take more than 50ms for any response
                start = spike['rise_index']
                if next_spike is not None:
                    stop = min(start + max_len, next_spike['rise_index'])
                else:
                    stop = start + max_len
                    
                # collect data from this trace
                trace = post_traces[j]
                d = trace.data[start:stop].copy()
                responses.append(d)

            if len(responses) == 0:
                continue
                
            # extend all responses to the same length and take nanmean
            avg = ragged_mean(responses, method='clip')
            avg -= float_mode(avg[:int(1e-3/dt)])
            avg_responses.append(avg)
            
            # plot average response for this pulse
            start = np.median([sp[i]['rise_index'] for sp in spikes if sp[i] is not None]) * dt
            t = np.arange(len(avg)) * dt
            self.response_plots[0,i].plot(t, avg, pen='w', antialias=True)

            # fit!
            fit = self.fit_psp(avg, t, dt, post_mode)
            fits.append(fit)
            
            # let the user mess with this fit
            curve = self.response_plots[0,i].plot(t, fit.eval(), pen=fit_pen, antialias=True).curve
            curve.setClickable(True)
            curve.fit = fit
            curve.sigClicked.connect(self.fit_curve_clicked)
            
        # display global average
        global_avg = ragged_mean(avg_responses, method='clip')
        t = np.arange(len(global_avg)) * dt
        self.response_plots[0,-1].plot(t, global_avg, pen='w', antialias=True)
        global_fit = self.fit_psp(global_avg, t, dt, post_mode)
        self.response_plots[0,-1].plot(t, global_fit.eval(), pen=fit_pen, antialias=True)
            
        # display fit parameters in table
        events = []
        for i,f in enumerate(fits + [global_fit]):
            if f is None:
                continue
            if i >= len(fits):
                vals = OrderedDict([('id', 'avg'), ('spike_time', np.nan), ('spike_stdev', np.nan)])
            else:
                spt = [s[i]['peak_index'] * dt for s in spikes if s[i] is not None]
                vals = OrderedDict([('id', i), ('spike_time', np.mean(spt)), ('spike_stdev', np.std(spt))])
            vals.update(OrderedDict([(k,f.best_values[k]) for k in f.params.keys()]))
            events.append(vals)
            
        self.current_event_set = (pre, post, events, sweeps)
        self.event_set_list.setCurrentRow(0)
        self.event_set_selected()

    def fit_psp(self, data, t, dt, clamp_mode):
        mode = float_mode(data[:int(1e-3/dt)])
        sign = -1 if data.mean() - mode < 0 else 1
        params = OrderedDict([
            ('xoffset', (2e-3, 3e-4, 5e-3)),
            ('yoffset', data[0]),
            ('amp', sign * 10e-12),
            #('k', (2e-3, 50e-6, 10e-3)),
            ('rise_time', (2e-3, 50e-6, 10e-3)),
            ('decay_tau', (4e-3, 500e-6, 50e-3)),
            ('rise_power', (2.0, 'fixed')),
        ])
        if clamp_mode == 'ic':
            params['amp'] = sign * 10e-3
            #params['k'] = (5e-3, 50e-6, 20e-3)
            params['rise_time'] = (5e-3, 50e-6, 20e-3)
            params['decay_tau'] = (15e-3, 500e-6, 150e-3)
        
        fit_kws = {'xtol': 1e-3, 'maxfev': 100}
        
        psp = fitting.Psp()
        return psp.fit(data, x=t, fit_kws=fit_kws, params=params)

    def add_set_clicked(self):
        if self.current_event_set is None:
            return
        ces = self.current_event_set
        self.event_sets.append(ces)
        item = QtWidgets.QListWidgetItem("%d -> %d" % (ces[0], ces[1]))
        self.event_set_list.addItem(item)
        item.event_set = ces
        
    def remove_set_clicked(self):
        if self.event_set_list.currentRow() == 0:
            return
        sel = self.event_set_list.takeItem(self.event_set_list.currentRow())
        self.event_sets.remove(sel.event_set)
        
    def event_set_selected(self):
        sel = self.event_set_list.selectedItems()[0]
        if sel.text() == "current":
            self.event_table.setData(self.current_event_set[2])
        else:
            self.event_table.setData(sel.event_set[2])
    
    def fit_clicked(self):
        from synaptic_release import ReleaseModel
        
        self.fit_plot.clear()
        self.fit_plot.show()
        
        n_sets = len(self.event_sets)
        self.fit_plot.set_shape(n_sets, 1)
        l = self.fit_plot[0, 0].legend
        if l is not None:
            l.scene.removeItem(l)
            self.fit_plot[0, 0].legend = None
        self.fit_plot[0, 0].addLegend()
        
        for i in range(n_sets):
            self.fit_plot[i,0].setXLink(self.fit_plot[0, 0])
        
        spike_sets = []
        for i,evset in enumerate(self.event_sets):
            evset = evset[2][:-1]  # select events, skip last row (average)
            x = np.array([ev['spike_time'] for ev in evset])
            y = np.array([ev['amp'] for ev in evset])
            x -= x[0]
            x *= 1000
            y /= y[0]
            spike_sets.append((x, y))
            
            self.fit_plot[i,0].plot(x/1000., y, pen=None, symbol='o')
            
        model = ReleaseModel()
        dynamics_types = ['Dep', 'Fac', 'UR', 'SMR', 'DSR']
        for k in dynamics_types:
            model.Dynamics[k] = 0
        
        fit_params = []
        with pg.ProgressDialog("Fitting release model..", 0, len(dynamics_types)) as dlg:
            for k in dynamics_types:
                model.Dynamics[k] = 1
                fit_params.append(model.run_fit(spike_sets))
                dlg += 1
                if dlg.wasCanceled():
                    return
        
        max_color = len(fit_params)*1.5

        for i,params in enumerate(fit_params):
            for j,spikes in enumerate(spike_sets):
                x, y = spikes
                t = np.linspace(0, x.max(), 1000)
                output = model.eval(x, params.values())
                y = output[:,1]
                x = output[:,0]/1000.
                self.fit_plot[j,0].plot(x, y, pen=(i,max_color), name=dynamics_types[i])

    def fit_curve_clicked(self, curve):
        if self.fit_explorer is None:
            self.fit_explorer = FitExplorer(curve.fit)
        else:
            self.fit_explorer.set_fit(curve.fit)
        self.fit_explorer.show()