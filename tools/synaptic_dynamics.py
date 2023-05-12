import pyqtgraph as pg
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import TSeries
from neuroanalysis.fitting import StackedPsp
from aisynphys.database import default_db as db
from aisynphys.ui.experiment_browser import ExperimentBrowser
from aisynphys.dynamics import pulse_response_query, sorted_pulse_responses


class DynamicsWindow(pg.QtWidgets.QSplitter):
    def __init__(self):
        self.loaded_pair = None
        
        pg.QtWidgets.QSplitter.__init__(self, pg.QtCore.Qt.Horizontal)
        self.ctrl_split = pg.QtWidgets.QSplitter(pg.QtCore.Qt.Vertical)
        self.addWidget(self.ctrl_split)
        
        self.browser = ExperimentBrowser()
        self.ctrl_split.addWidget(self.browser)
        
        self.ptree = pg.parametertree.ParameterTree()
        self.ctrl_split.addWidget(self.ptree)
        
        self.params = pg.parametertree.Parameter.create(name='params', type='group', children=[
            {'name': 'show spikes', 'type': 'bool', 'value': True},
            {'name': 'subtract baseline', 'type': 'bool', 'value': True},
            {'name': 'stimulus filter', 'type': 'group'},
        ])
        self.ptree.setParameters(self.params)
        
        self.scroll_area = pg.QtWidgets.QScrollArea()
        self.addWidget(self.scroll_area)
        
        self.view = pg.GraphicsLayoutWidget()
        self.scroll_area.setWidget(self.view)
        
        self.resize(1600, 1000)
        
        self.plots = []
        
        self.browser.itemSelectionChanged.connect(self.browser_item_selected)
        self.params.sigTreeStateChanged.connect(self.plot_all)

    def clear(self):
        for plt in self.plots:
            self.view.removeItem(plt)
        self.plots = []
    
    def browser_item_selected(self):
        with pg.BusyCursor():
            selected = self.browser.selectedItems()
            if len(selected) != 1:
                return
            item = selected[0]
            if not hasattr(item, 'pair'):
                return
            pair = item.pair

            self.load_pair(pair)
        
    def load_pair(self, pair):
        if pair is not self.loaded_pair:
            print("Loading:", pair)
            q = pulse_response_query(pair, data=True, spike_data=True)
            self.sorted_recs = sorted_pulse_responses(q.all())
            self.stim_keys = sorted(list(self.sorted_recs.keys()))
            self.update_params()
            self.loaded_pair = pair
        
        self.plot_all()
        
    def update_params(self):
        with pg.SignalBlock(self.params.sigTreeStateChanged, self.plot_all):
            stim_param = self.params.child('stimulus filter')
            for ch in stim_param.children():
                stim_param.removeChild(ch)
            
            for k in self.stim_keys:
                param = pg.parametertree.Parameter.create(name=str(k), type="bool", value="True")
                stim_param.addChild(param)
        
    def plot_all(self):
        with pg.BusyCursor():
            self.clear()
            show_spikes = self.params['show spikes']

            for i,stim_key in enumerate(self.stim_keys):
                if self.params['stimulus filter', str(stim_key)] is False:
                    continue
                
                plt = DynamicsPlot()
                self.plots.append(plt)
                self.view.addItem(plt)
                self.view.nextRow()
                plt.set_title("%s  %0.0f Hz  %0.2f s" % stim_key)
                prs = self.sorted_recs[stim_key]
                plt.set_data(prs, show_spikes=show_spikes, subtract_baseline=self.params['subtract baseline'])

            plt_height = max(400 if show_spikes else 250, self.scroll_area.height() / len(self.stim_keys))
            self.view.setFixedHeight(plt_height * len(self.plots))
            self.view.setFixedWidth(self.scroll_area.width() - self.scroll_area.verticalScrollBar().width())
        

class DynamicsPlot(pg.GraphicsLayout):
    def __init__(self):
        pg.GraphicsLayout.__init__(self)
        self.show_spikes = False
        
        self.label = pg.TextItem()
        self.label.setParentItem(self)
        self.label.setPos(200, 0)
        self.spike_plot = self.addPlot()
        self.spike_plot.setVisible(False)
        self.nextRow()
        self.data_plot = self.addPlot()
        self.spike_plot.setXLink(self.data_plot)
        
    def set_title(self, title):
        self.label.setText(title)
        
    def set_data(self, data, show_spikes=False, subtract_baseline=True):
        self.spike_plot.setVisible(show_spikes)
        self.spike_plot.enableAutoRange(False, False)
        self.data_plot.enableAutoRange(False, False)
        psp = StackedPsp()
        
        for recording in data:
            pulses = sorted(list(data[recording].keys()))
            for pulse_n in pulses:
                rec = data[recording][pulse_n]
                # spike-align pulse + offset for pulse number
                spike_t = rec.StimPulse.first_spike_time
                if spike_t is None:
                    spike_t = rec.StimPulse.onset_time + 1e-3
                    
                qc_pass = rec.PulseResponse.in_qc_pass if rec.Synapse.synapse_type == 'in' else rec.PulseResponse.ex_qc_pass
                pen = (255, 255, 255, 100) if qc_pass else (200, 50, 0, 100)
                
                t0 = rec.PulseResponse.data_start_time - spike_t
                ts = TSeries(data=rec.data, t0=t0, sample_rate=db.default_sample_rate)
                c = self.data_plot.plot(ts.time_values, ts.data, pen=pen)
                
                # arrange plots nicely
                y0 = 0 if not subtract_baseline else ts.time_slice(None, 0).median()
                shift = (pulse_n * 35e-3 + (30e-3 if pulse_n > 8 else 0), -y0)
                zval = 0 if qc_pass else -10
                c.setPos(*shift)
                c.setZValue(zval)

                if show_spikes:
                    t0 = rec.spike_data_start_time - spike_t
                    spike_ts = TSeries(data=rec.spike_data, t0=t0, sample_rate=db.default_sample_rate)
                    c = self.spike_plot.plot(spike_ts.time_values, spike_ts.data, pen=pen)
                    c.setPos(*shift)
                    c.setZValue(zval)
                    
                # evaluate recorded fit for this response
                fit_par = rec.PulseResponseFit
                if fit_par.fit_amp is None:
                    continue
                fit = psp.eval(
                    x=ts.time_values, 
                    exp_amp=fit_par.fit_exp_amp,
                    exp_tau=fit_par.fit_decay_tau,
                    amp=fit_par.fit_amp,
                    rise_time=fit_par.fit_rise_time,
                    decay_tau=fit_par.fit_decay_tau,
                    xoffset=fit_par.fit_latency,
                    yoffset=fit_par.fit_yoffset,
                    rise_power=2,
                )
                pen = (0, 255, 0, 100) if qc_pass else (50, 150, 0, 100)
                c = self.data_plot.plot(ts.time_values, fit, pen=pen)
                c.setZValue(10)
                c.setPos(*shift)

                if not qc_pass:
                    print("qc fail: ", rec.PulseResponse.meta.get('qc_failures', 'no qc failures recorded'))
                    
        self.spike_plot.enableAutoRange(True, True)
        self.data_plot.enableAutoRange(True, True)


if __name__ == '__main__':
    import sys, argparse
    from aisynphys import config
    
    parser = argparse.ArgumentParser(parents=[config.parser])
    parser.add_argument('experiment_id', type=str, nargs='?')
    parser.add_argument('pre_cell_id', type=str, nargs='?')
    parser.add_argument('post_cell_id', type=str, nargs='?')
    args = parser.parse_args()
    
    app = pg.mkQApp()
    if sys.flags.interactive == 1:
        pg.dbg()
    
    win = DynamicsWindow()
    win.show()
    
    if args.post_cell_id is not None:
        expt = db.experiment_from_ext_id(args.experiment_id)
        win.browser.populate([expt], synapses=True)
        pair = expt.pairs[args.pre_cell_id, args.post_cell_id]
        win.browser.select_pair(pair.id)
    else:
        win.browser.populate(synapses=True)
    
    if sys.flags.interactive == 0:
        app.exec_()
