import numpy as np
import pyqtgraph as pg
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.fitting import StackedPsp
from aisynphys.data import PulseResponseList


class AvgResponseFitUi(object):
    def __init__(self):
        self.widget = pg.QtWidgets.QSplitter(pg.QtCore.Qt.Horizontal)
        self.plots = {}
        for clamp_mode in ('ic', 'vc'):
            units = {'ic': 'V', 'vc': 'A'}[clamp_mode]
            for holding in (-70, -55):
                prp = PulseResponsePlot("%s %d" % (clamp_mode, holding), units)
                self.plots[clamp_mode, holding] = prp
                self.widget.addWidget(prp.grid)

    def show_pulse_responses(self, prs):
        for (clamp_mode, holding), plot in self.plots.items():
            plot.show_pulse_responses(prs[clamp_mode, holding]['qc_pass'], prs[clamp_mode, holding]['qc_fail'])

    def show_data_notes(self, notes_rec):
        if notes_rec is None:
            return
        notes = notes_rec.notes
        for (clamp_mode, holding), plot in self.plots.items():
            params = notes['fit_parameters']['fit'][clamp_mode][str(holding)].copy()
            if len(params) == 0:
                continue
            params.pop('nrmse')
            params.setdefault('exp_tau', params['decay_tau'])
            params.setdefault('exp_amp', 0)
            qc_pass = notes['fit_pass'][clamp_mode][str(holding)]
            plot.show_expected_fit(params, qc_pass)

    def show_fit_results(self, clamp_mode, holding, results, average, qc_pass):
        plot = self.plots[clamp_mode, holding]
        plot.show_fit_results(results, average, qc_pass)


class PulseResponsePlot(object):
    def __init__(self, title, units):
        self.grid = PlotGrid()
        self.grid.set_shape(2, 1)
        self.grid.grid.ci.layout.setRowStretchFactor(0, 3)
        self.grid.grid.ci.layout.setRowStretchFactor(1, 8)
        self.spike_plot = self.grid[0, 0]
        self.response_plot = self.grid[1, 0]

        self.spike_plot.hideAxis('bottom')
        self.spike_plot.setTitle(title)
        self.spike_plot.setLabel('left', text="presynaptic spike")
        self.spike_plot.addLine(x=0)

        self.response_plot.setLabel('bottom', text='Time from spike', units='s')
        self.response_plot.setXLink(self.spike_plot)
        # self.response_plot.setLabel('left', text="%d holding" % int(holding), units=units)

        self.response_plot.enableAutoRange(False, False)
        
        self.items = []

    def show_pulse_responses(self, passed_prs, failed_prs):
        for prs, pen in [(failed_prs, (255, 0, 0, 40)), (passed_prs, (255, 255, 255, 40))]:
            if len(prs) == 0:
                continue

            prl = PulseResponseList(prs)

            post_ts = prl.post_tseries(align='spike', bsub=True)
            for ts in post_ts:
                item = self.response_plot.plot(ts.time_values, ts.data, pen=pen)
                self.items.append(item)

            pre_ts = prl.pre_tseries(align='spike', bsub=True)
            for ts in pre_ts:
                item = self.spike_plot.plot(ts.time_values, ts.data, pen=pen)
                self.items.append(item)

        self.response_plot.autoRange()

    def show_expected_fit(self, fit_params, qc_pass):
        psp = StackedPsp()
        t = np.linspace(-10e-3, 20e-3, 1000)
        v = psp.eval(x=t, **fit_params)
        pen = {'color': (0, 150, 0) if qc_pass else (255, 100, 0), 'dash': [5, 5]} 
        self.response_plot.plot(t, v, pen=pen, zValue=10)

    def show_fit_results(self, results, average, qc_pass):
        self.response_plot.plot(average.time_values, average.data, pen='b')

        psp = StackedPsp()
        t = average.time_values
        v = psp.eval(x=t, **results.best_values)
        pen = {'color':  (0, 150, 0) if qc_pass else (255, 100, 0), 'width': 2}
        self.response_plot.plot(t, v, pen=pen)

    def plot_fit(self, trace, fit, holding, fit_pass=False):
        if holding == '-55':
            if self.fit_item_55 is not None:
                self.trace_plots[0].removeItem(self.fit_item_55)
            self.fit_item_55 = pg.PlotDataItem(trace.time_values, fit, name='-55 holding', pen={'color': self.fit_color[fit_pass], 'width': 3})
            self.trace_plots[0].addItem(self.fit_item_55)
        
        elif holding == '-70':
            if self.fit_item_70 is not None:
                self.trace_plots[1].removeItem(self.fit_item_70)
            self.fit_item_70 = pg.PlotDataItem(trace.time_values, fit, name='-70 holding', pen={'color': self.fit_color[fit_pass], 'width': 3})
            self.trace_plots[1].addItem(self.fit_item_70)

    def color_fit(self, name, value):
        if '-55' in name:
            if self.fit_item_55 is not None:
                self.fit_item_55.setPen({'color': self.fit_color[value], 'width': 3})
        if '-70' in name:
            if self.fit_item_70 is not None:
                self.fit_item_70.setPen({'color': self.fit_color[value], 'width': 3})

    def clear_plots(self):
        for item in self.items + [self.fit_item_55, self.fit_item_70]:
            if item is None:
                continue
            item.scene().removeItem(item)
        self.items = []
        
        self.plots[-1].autoRange()
        self.plots[-1].setXRange(-5e-3, 10e-3)
        self.fit_item_70 = None
        self.fit_item_55 = None
        
