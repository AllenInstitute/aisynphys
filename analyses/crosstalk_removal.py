# coding: utf8
"""
Question: can we digitally remove pipette electrical / capacitive crosstalk somehow?

"""
from __future__ import division, print_function
import numpy as np
import pyqtgraph as pg

from neuroanalysis.data import TraceList
from neuroanalysis import filter
from multipatch_analysis.ui.experiment_browser import ExperimentBrowser
from multipatch_analysis.connection_strength import get_amps, get_baseline_amps
from multipatch_analysis import database as db
from multipatch_analysis import data


class StimResponseList(object):
    """A collection of stimulus-response objects that provides some common analysis tools
    """
    def __init__(self, srs):
        self.srs = srs

    def get_tseries(self, series, bsub=True, align='stim', bsub_window=(-3e-3, 0)):
        """Return a TraceList of timeseries, optionally baseline-subtracted and time-aligned.
        
        Parameters
        ----------
        series : str
            "stim", "pre", or "post"
        """
        assert series in ('stim', 'pre', 'post'), "series must be one of 'stim', 'pre', or 'post'"
        tseries = []
        for i,sr in enumerate(self.srs):
            ts = getattr(sr, series + '_tseries')
            if bsub:
                bstart = sr.stim_pulse.onset_time + bsub_window[0]
                bstop = sr.stim_pulse.onset_time + bsub_window[1]
                baseline = np.median(ts.time_slice(bstart, bstop).data)
                ts = ts - baseline
            if align is not None:
                if align == 'stim':
                    t_align = sr.stim_pulse.onset_time 
                elif align == 'pre':
                    t_align = sr.stim_pulse.spikes[0].max_dvdt_time
                elif align == 'post':
                    raise NotImplementedError()
                else:
                    raise ValueError("invalid time alignment mode %r" % align)
                t_align = t_align or 0
                ts = ts.copy(t0=ts.t0-t_align)
            tseries.append(ts)
        return TraceList(tseries)

    def __iter__(self):
        for sr in self.srs:
            yield sr

    def __getitem__(self, index):
        item = self.srs.__getitem__(index)
        if isinstance(item, list):
            return StimResponseList(item)
        else:
            return item

    def __len__(self):
        return len(self.srs)


class CrosstalkAnalyzer(object):
    """User interface for exploring crosstalk removal methods
    """
    def __init__(self, expt_browser=None):
        if expt_browser is None:
            expt_browser = ExperimentBrowser()
        self.expt_browser = expt_browser
        expt_browser.itemSelectionChanged.connect(self.expt_selection_changed)
        self.session = db.Session()
        self.pair = None
        self.response_list = None
        self.corrected_response_list = None
        self._plot_items = []

        self.params = pg.parametertree.Parameter(name="crosstalk removal", type="group", children=[
            dict(name='data', type='group', children=[
                dict(name='pre mode', type='list', values=['ic', 'vc']),
                dict(name='post mode', type='list', values=['ic', 'vc']),
                dict(name='limit', type='int', value=100),
            ]),
            dict(name='correction', type='group', children=[
                dict(name='stimulus', type='bool', value=True, children=[
                    dict(name='auto', type='action'),
                    dict(name='shift', type='int', value=0),
                    dict(name='scale', type='float', value=30e3, dec=True, step=0.5),
                    dict(name='lowpass', type='float', value=7e3, dec=True, step=0.5, suffix='Hz', siPrefix=True),
                ]),
                dict(name='charging', type='bool', value=True, children=[
                    dict(name='auto', type='action'),
                    dict(name='capacitance', type='float', value=10e-12, suffix='F', dec=True, step=0.5, siPrefix=True),
                    dict(name='resistance', type='float', value=10e6, suffix=u'Ω', dec=True, step=0.5, siPrefix=True),
                    dict(name='scale', type='float', value=30e3, dec=True, step=0.5),
                    dict(name='lowpass', type='float', value=7e3, dec=True, step=0.5, suffix='Hz', siPrefix=True),
                ]),
                dict(name='spike', type='bool', value=False, children=[
                    dict(name='scale', type='float', value=0.001, dec=True, step=0.5),
                    dict(name='lowpass', type='float', value=10e3, dec=True, step=0.5, suffix='Hz', siPrefix=True),
                ]),
                dict(name='spike dv/dt', type='bool', value=True, children=[
                    dict(name='plot only', type='bool', value=False),
                    dict(name='scale', type='float', value=0.03, dec=True, step=0.5),
                    dict(name='lowpass', type='float', value=550., dec=True, step=0.5, suffix='Hz', siPrefix=True),
                ]),
            ]),
            dict(name='display', type='group', children=[
                dict(name='limit', type='int', value=20),
                dict(name='plot spike aligned', type='bool', value=False),
                dict(name='plot dv/dt', type='bool', value=False),
                dict(name='plot lowpass', type='float', value=10e3, dec=True, step=0.5, suffix='Hz', siPrefix=True),
            ])
        ])
        self.ptree = pg.parametertree.ParameterTree()
        self.ptree.setParameters(self.params)
        self.params.child('data').sigTreeStateChanged.connect(self.data_filter_changed)
        self.params.child('correction').sigTreeStateChanged.connect(self.update_analysis)
        self.params.child('display').sigTreeStateChanged.connect(self.update_display)
        self.params.child('correction', 'stimulus', 'auto').sigActivated.connect(self.auto_stim_clicked)

        self.pw = pg.GraphicsLayoutWidget()
        
        self.plt1 = self.pw.addPlot(0, 0)
        self.plt2 = self.pw.addPlot(1, 0)
        self.plt3 = self.pw.addPlot(2, 0)
        self.plt2.setXLink(self.plt1)
        self.plt3.setXLink(self.plt1)

    def show(self):
        self.win = pg.QtGui.QSplitter(pg.QtCore.Qt.Horizontal)
        self.vsplit = pg.QtGui.QSplitter(pg.QtCore.Qt.Vertical)
        self.win.addWidget(self.vsplit)
        self.vsplit.addWidget(self.expt_browser)
        self.vsplit.addWidget(self.ptree)
        self.win.addWidget(self.pw)
        self.win.show()
        return self.win
        
    def expt_selection_changed(self):
        sel = self.expt_browser.selectedItems()
        if len(sel) == 0:
            return
        sel = sel[0]
        if not hasattr(sel, 'pair'):
            return
        self.load_pair(sel.pair)

    def select_pair(self, pair_id):
        self.expt_browser.select_pair(pair_id)

    def load_pair(self, pair=None):
        pair = pair or self.pair
        if pair is None:
            return
            
        limit = self.params['data', 'limit']
        with pg.BusyCursor():
            
            pre_prec = db.aliased(db.PatchClampRecording)
            post_prec = db.aliased(db.PatchClampRecording)
            q = db.query(db.PulseResponse).join(db.StimPulse).join(db.Pair)
            q = q.join(pre_prec, db.PulseResponse.recording_id==pre_prec.recording_id)
            q = q.join(post_prec, db.StimPulse.recording_id==post_prec.recording_id)
            q = q.filter(db.Pair.id==pair.id)
            q = q.filter(pre_prec.clamp_mode==self.params['data', 'pre mode'])
            q = q.filter(post_prec.clamp_mode==self.params['data', 'post mode'])
            prs = q.limit(limit).all()
            
            # preload pulse response data
            ids = [pr.id for pr in prs]
            db.query(db.PulseResponse.data).filter(db.PulseResponse.id.in_(ids)).all()
            
            self.response_list = StimResponseList(prs)
            self.update_analysis()

    def data_filter_changed(self):
        self.load_pair()

    def update_analysis(self):
        with pg.BusyCursor():
            selected_responses = self.response_list
            self.corrected_response_list = StimResponseList([self.process_response(sr) for sr in selected_responses])
            self.update_display()

    def clear_plots(self):
        for item, plt in self._plot_items:
            plt.removeItem(item)
        self._plot_items = []

    def update_display(self):
        with pg.BusyCursor():
            align = 'pre' if self.params['display', 'plot spike aligned'] else 'stim'
            
            self.clear_plots()

            rl = self.corrected_response_list
            if len(rl) == 0:
                return
            
            stim_ts = rl.get_tseries('stim', align=align)
            self._plot_ts(stim_ts, self.plt1)
        
            pre_ts = rl.get_tseries('pre', align=align)
            self._plot_ts(pre_ts, self.plt2)

            post_ts = rl.get_tseries('post', align=align)
            dvdt = self.params['display', 'plot dv/dt']
            lowpass = self.params['display', 'plot lowpass']
            self._plot_ts(post_ts, self.plt3, dvdt=dvdt, lowpass=lowpass)

    def display_filter(self, ts, dvdt, lowpass):
        if dvdt is True:
            ts = ts.copy(data=np.diff(ts.data))
        if lowpass is not None:
            try:
                ts = filter.bessel_filter(ts, lowpass, bidir=True)
            except ValueError:
                pass
        return ts
        
    def _plot_ts(self, ts_list, plt, dvdt=False, lowpass=None):
        limit = self.params['display', 'limit']
        for ts in ts_list[:limit]:
            ts = self.display_filter(ts, dvdt, lowpass)
            item = plt.plot(ts.time_values, ts.data, pen=(255, 255, 255, 80), antialias=True)
            self._plot_items.append((item, plt))
        avg = ts_list.mean()
        avg = self.display_filter(avg, dvdt, lowpass)
        item = plt.plot(avg.time_values, avg.data, pen={'color':'g', 'width':2}, shadowPen={'color':'k', 'width':3}, antialias=True)
        self._plot_items.append((item, plt))

    def process_response(self, sr, plt=None):
        stim = sr.stim_tseries
        pre = sr.pre_tseries
        post = sr.post_tseries.copy()

        if self.params['correction', 'stimulus']:
            pulse_start = sr.stim_pulse.onset_time
            pulse_stop = pulse_start + sr.stim_pulse.duration
            stim_correction = post.copy(data=np.zeros(len(post)))
            stim_correction.t0 = stim_correction.t0 + stim_correction.dt * self.params['correction', 'stimulus', 'shift']
            stim_correction.time_slice(pulse_start, pulse_stop).data[:] = self.params['correction', 'stimulus', 'scale'] * sr.stim_pulse.amplitude
            stim_correction = filter.bessel_filter(stim_correction, self.params['correction', 'stimulus', 'lowpass'], bidir=False)
            post = post - stim_correction.data

        if self.params['correction', 'charging']:
            tau = self.params['correction', 'charging', 'capacitance'] * self.params['correction', 'charging', 'resistance']
            scale = self.params['correction', 'charging', 'scale']
            pulse_start = sr.stim_pulse.onset_time
            pulse_stop = pulse_start + sr.stim_pulse.duration
            charging_correction = post.copy(data=np.zeros(len(post)))

            during_pulse = charging_correction.time_slice(pulse_start, pulse_stop)
            t = during_pulse.time_values - during_pulse.t0
            during_pulse.data[:] = scale * sr.stim_pulse.amplitude * (1.0 - np.exp(-t / tau))

            after_pulse = charging_correction.time_slice(pulse_stop, None)
            t = after_pulse.time_values - after_pulse.t0
            after_pulse.data[:] = during_pulse.data[-1] * np.exp(-t / tau)

            charging_correction = filter.bessel_filter(charging_correction, self.params['correction', 'charging', 'lowpass'], bidir=False)
            post = post - charging_correction.data

        if self.params['correction', 'spike']:
            spike = pre.copy(data=pre.data * self.params['correction', 'spike', 'scale'])
            spike_correction = filter.bessel_filter(spike, self.params['correction', 'spike', 'lowpass'], bidir=False)
            start = max(spike.t0, post.t0)
            stop = min(spike.t_end, post.t_end)

            post = post.time_slice(start, stop) - spike_correction.time_slice(start, stop).data

        if self.params['correction', 'spike dv/dt']:
            spike_diff = np.diff(pre.data) * self.params['correction', 'spike dv/dt', 'scale']
            spike = pre.copy(data=spike_diff)
            spike_correction = filter.bessel_filter(spike, self.params['correction', 'spike dv/dt', 'lowpass'], bidir=False)
            start = max(spike.t0, post.t0)
            stop = min(spike.t_end, post.t_end)

            post = post.time_slice(start, stop) - spike_correction.time_slice(start, stop).data

        return data.PulseResponse(stim_pulse=sr.stim_pulse, post_tseries=post)

    def auto_stim_clicked(self):
        print("hi")

        

if __name__ == '__main__':
    import sys
    app = pg.mkQApp()
    pg.dbg()
    
    cta = CrosstalkAnalyzer()

    if len(sys.argv) > 1:
        pair_id = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
        cta.select_pair(pair_id)

    win = cta.show()
    win.resize(1400, 1000)
    
    if sys.flags.interactive == 0:
        app.exec_()
