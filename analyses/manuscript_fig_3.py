# coding: utf8
"""
2018 E-E manuscript fig 3:
Analysis of detection limits vs synaptic strength, kinetics, and background noise
"""
from __future__ import print_function, division
from datetime import datetime
import os, pickle
import multiprocessing
import pyqtgraph as pg
import numpy as np
from scipy import stats, ndimage

from neuroanalysis.data import TSeries, TSeriesList
import strength_analysis
from aisynphys.database import database as db




def write_csv(fh, data, description, units='membrane voltage (V)'):
    """Used to generate csv file accompanying figure.
    """
    if isinstance(data, TSeries):
        write_csv(fh, data.time_values, description + " time (s)")
        write_csv(fh, data.data, description + " %s" % units)
    else:
        cols = ['"' + description + '"'] + list(data)
        line = ','.join(map(str, cols))
        fh.write(line)
        fh.write('\n')


if __name__ == '__main__':
    csv_file = open("manuscript_fig_3.csv", 'wb')

    # silence warnings about fp issues
    np.seterr(all='ignore')

    # Three connections selected for analysis
    show_conns = [
        # (expt_uid, pre_cell, post_cell)
        ("high signal, low noise", (1499725138.07, 7, 4)),
        ("low signal, high noise", (1494969844.93, 6, 1)),
        ("high noise, no detected connection", (1489009391.46, 3, 5)),
    ]
    
    pg.mkQApp()
    pg.dbg()
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')

    win = pg.GraphicsLayoutWidget()
    win.show()
    win.resize(900, 900)

    # set up scatter plot
    scatter_plot = win.addPlot(0, 0, rowspan=len(show_conns))
    scatter_plot.setLogMode(True, True)
    scatter_plot.setAspectLocked()
    scatter_plot.setFixedWidth(350)
    scatter_plot.showGrid(True, True, alpha=0.5)

    # read all pair records from DB
    classifier = strength_analysis.get_pair_classifier(seed=0, use_vc_features=False)
    conns = strength_analysis.query_all_pairs(classifier)

    # filter
    mask = np.isfinite(conns['ic_deconv_amp_mean'])
    filtered = conns[mask]

    # remove recordings with gain errors
    mask = filtered['ic_deconv_amp_mean'] < 0.02

    # remove recordings with high crosstalk
    mask &= abs(filtered['ic_crosstalk_mean']) < 60e-6

    # remove recordings with low sample count
    mask &= filtered['ic_n_samples'] > 50

    typs = filtered['pre_cre_type']
    mask &= typs == filtered['post_cre_type']

    typ_mask = ((typs == 'sim1') | (typs == 'tlx3') | (typs == 'unknown') | (typs == 'rorb') | (typs == 'ntsr1'))
    mask &= typ_mask

    # This is the final set of experiments we will include in the analysis here
    filtered = filtered[mask]


    # plot signal vs background for all pairs
    brushes = [pg.mkBrush('y') if c['synapse'] else pg.mkBrush(0.5) for c in filtered]

    c_mask = filtered['synapse'] == True
    u_mask = ~c_mask

    signal = filtered['ic_fit_amp']
    background = filtered['minimum_amplitude']

    # plot connected pairs
    x, y = background[c_mask], signal[c_mask]
    c_plot = scatter_plot.plot(x, y, pen=None, symbol='d', symbolPen='k', symbolBrush=(0, 255, 255), symbolSize=10, data=filtered[c_mask])
    write_csv(csv_file, x, "Figure 3A connected pairs; minimum detectable amplitude (V)")
    write_csv(csv_file, y, "Figure 3A connected pairs; PSP amplitude (V)")

    # plot unconnected pairs
    x, y = background[u_mask], signal[u_mask]
    u_plot = scatter_plot.plot(x, y, pen=None, symbol='o', symbolPen=None, symbolBrush=(50, 50, 50, 80), symbolSize=4, data=filtered[u_mask])
    write_csv(csv_file, x, "Figure 3A unconnected pairs; minimum detectable amplitude (V)")
    write_csv(csv_file, y, "Figure 3A unconnected pairs; PSP amplitude (V)")


    # plot for showing distribution of response amplitudes and detection limits
    profile_plot = win.addPlot(3, 0, colspan=5)
    profile_plot.setMinimumHeight(300)
    profile_plot.setXRange(2e-6, 2e-3, padding=0)
    profile_plot.setLogMode(True, False)

    
    # Select pairs to use for amplitude profile analysis
    prof_conns = filtered[c_mask]
    prof_amps = prof_conns['ic_fit_amp']

    # Plot distribution of connection amplitudes
    n_bins = 30
    bins = 10e-6 * ((2e-3/10e-6)**(1./n_bins)) ** np.arange(n_bins+1)
    amp_hist = list(np.histogram(prof_amps, bins=bins))
    amp_hist[0] = ndimage.gaussian_filter(amp_hist[0].astype(float), 1)
    x, y = amp_hist[1], amp_hist[0]
    profile_plot.plot(x, y, stepMode=True, fillLevel=0, brush=(200, 200, 200), pen='k')
    write_csv(csv_file, x, "Figure 3E measured PSP amplitude distribution bin edges (V)")
    write_csv(csv_file, y, "Figure 3E measured PSP amplitude distribution connections per bin")

    # Plot probability of detection vs PSP amplitude 
    min_amps = prof_conns['minimum_amplitude']
    min_amps = min_amps[np.isfinite(min_amps)]
    limit_hist = list(np.histogram(min_amps, bins=bins))
    limit_hist[0] = limit_hist[0].astype(float) / limit_hist[0].sum()
    bin_centers = (bins[1:] * bins[:-1]) ** 0.5
    limit_prof = np.cumsum(limit_hist[0])
    x, y = bin_centers, limit_prof * amp_hist[0].max()
    profile_plot.plot(x, y, pen='r')
    write_csv(csv_file, x, "Figure 3E detection probability x values (V)")
    write_csv(csv_file, y, "Figure 3E detection probability y values")

    # Plot corrected amplitude distribution
    corrected_prof = ndimage.gaussian_filter(amp_hist[0], 0) / limit_prof
    x, y = bins, corrected_prof
    profile_plot.plot(x, y, stepMode=True, fillLevel=0, brush=(120, 120, 120), pen='k').setZValue(-10)
    write_csv(csv_file, x, "Figure 3E estimated PSP amplitude distribution bin edges (V)")
    write_csv(csv_file, y, "Figure 3E estimated PSP amplitude distribution connections per bin")
    
    print("Global connectivity correction factor:", corrected_prof.sum() / amp_hist[0].sum())

    trace_plots = []
    deconv_plots = []
    hist_plots = []

    session = db.session()

    def add_connection_plots(i, name, timestamp, pre_id, post_id):
        global session, win, filtered
        p = pg.debug.Profiler(disabled=True, delayed=False)
        trace_plot = win.addPlot(i, 1)
        trace_plots.append(trace_plot)
        trace_plot.setYRange(-1.4e-3, 2.1e-3)
        # deconv_plot = win.addPlot(i, 2)
        # deconv_plots.append(deconv_plot)
        # deconv_plot.hide()
        
        hist_plot = win.addPlot(i, 2)
        hist_plots.append(hist_plot)
        limit_plot = win.addPlot(i, 3)
        limit_plot.addLegend()
        limit_plot.setLogMode(True, False)
        limit_plot.addLine(y=classifier.prob_threshold)

        # Find this connection in the pair list
        idx = np.argwhere((abs(filtered['acq_timestamp'] - timestamp) < 1) & (filtered['pre_cell_id'] == pre_id) & (filtered['post_cell_id'] == post_id))
        if idx.size == 0:
            print("not in filtered connections")
            return
        idx = idx[0,0]
        p()

        # Mark the point in scatter plot
        scatter_plot.plot([background[idx]], [signal[idx]], pen='k', symbol='o', size=10, symbolBrush='r', symbolPen=None)
        
        # Plot example traces and histograms
        for plts in [trace_plots]:#, deconv_plots]:
            plt = plts[-1]
            plt.setXLink(plts[0])
            plt.setYLink(plts[0])
            plt.setXRange(-10e-3, 17e-3, padding=0)
            plt.hideAxis('left')
            plt.hideAxis('bottom')
            plt.addLine(x=0)
            plt.setDownsampling(auto=True, mode='peak')
            plt.setClipToView(True)
            hbar = pg.QtWidgets.QGraphicsLineItem(0, 0, 2e-3, 0)
            hbar.setPen(pg.mkPen(color='k', width=5))
            plt.addItem(hbar)
            vbar = pg.QtWidgets.QGraphicsLineItem(0, 0, 0, 100e-6)
            vbar.setPen(pg.mkPen(color='k', width=5))
            plt.addItem(vbar)


        hist_plot.setXLink(hist_plots[0])
        
        pair = session.query(db.Pair).filter(db.Pair.id==filtered[idx]['pair_id']).all()[0]
        p()
        amps = strength_analysis.get_amps(session, pair)
        p()
        base_amps = strength_analysis.get_baseline_amps(session, pair, amps=amps, clamp_mode='ic')
        p()
        
        q = strength_analysis.response_query(session)
        p()
        q = q.join(strength_analysis.PulseResponseStrength)
        q = q.filter(strength_analysis.PulseResponseStrength.id.in_(amps['id']))
        q = q.join(db.MultiPatchProbe)
        q = q.filter(db.MultiPatchProbe.induction_frequency < 100)
        # pre_cell = db.aliased(db.Cell)
        # post_cell = db.aliased(db.Cell)
        # q = q.join(db.Pair).join(db.Experiment).join(pre_cell, db.Pair.pre_cell_id==pre_cell.id).join(post_cell, db.Pair.post_cell_id==post_cell.id)
        # q = q.filter(db.Experiment.id==filtered[idx]['experiment_id'])
        # q = q.filter(pre_cell.ext_id==pre_id)
        # q = q.filter(post_cell.ext_id==post_id)

        fg_recs = q.all()
        p()

        traces = []
        deconvs = []
        for i,rec in enumerate(fg_recs[:100]):
            result = strength_analysis.analyze_response_strength(rec, source='pulse_response', lpf=True, lowpass=2000,
                                                remove_artifacts=False, bsub=True)
            trace = result['raw_trace']
            trace.t0 = -result['spike_time']
            trace = trace - np.median(trace.time_slice(-0.5e-3, 0.5e-3).data)
            traces.append(trace)
            trace_plot.plot(trace.time_values, trace.data, pen=(0, 0, 0, 20))
            write_csv(csv_file, trace, "Figure 3B; {name}; trace {trace_n}".format(name=name, trace_n=i))

            # trace = result['dec_trace']
            # trace.t0 = -result['spike_time']
            # trace = trace - np.median(trace.time_slice(-0.5e-3, 0.5e-3).data)
            # deconvs.append(trace)            
            # # deconv_plot.plot(trace.time_values, trace.data, pen=(0, 0, 0, 20))

        # plot average trace
        mean = TSeriesList(traces).mean()
        trace_plot.plot(mean.time_values, mean.data, pen={'color':'g', 'width': 2}, shadowPen={'color':'k', 'width': 3}, antialias=True)
        write_csv(csv_file, mean, "Figure 3B; {name}; average".format(name=name))
        # mean = TSeriesList(deconvs).mean()
        # # deconv_plot.plot(mean.time_values, mean.data, pen={'color':'g', 'width': 2}, shadowPen={'color':'k', 'width': 3}, antialias=True)

        # add label
        label = pg.LabelItem(name)
        label.setParentItem(trace_plot)


        p("analyze_response_strength")

        # bins = np.arange(-0.0005, 0.002, 0.0001) 
        # field = 'pos_amp'
        bins = np.arange(-0.001, 0.015, 0.0005) 
        field = 'pos_dec_amp'
        n = min(len(amps), len(base_amps))
        hist_y, hist_bins = np.histogram(base_amps[:n][field], bins=bins)
        hist_plot.plot(hist_bins, hist_y, stepMode=True, pen=None, brush=(200, 0, 0, 150), fillLevel=0)
        write_csv(csv_file, hist_bins, "Figure 3C; {name}; background noise amplitude distribution bin edges (V)".format(name=name))
        write_csv(csv_file, hist_y, "Figure 3C; {name}; background noise amplitude distribution counts per bin".format(name=name))
        
        hist_y, hist_bins = np.histogram(amps[:n][field], bins=bins)
        hist_plot.plot(hist_bins, hist_y, stepMode=True, pen='k', brush=(0, 150, 150, 100), fillLevel=0)
        write_csv(csv_file, hist_bins, "Figure 3C; {name}; PSP amplitude distribution bin edges (V)".format(name=name))
        write_csv(csv_file, hist_y, "Figure 3C; {name}; PSP amplitude distribution counts per bin".format(name=name))
        p()

        pg.QtWidgets.QApplication.processEvents()


        # Plot detectability analysis
        q = strength_analysis.baseline_query(session)
        q = q.join(strength_analysis.BaselineResponseStrength)
        q = q.filter(strength_analysis.BaselineResponseStrength.id.in_(base_amps['id']))
        # q = q.limit(100)
        bg_recs = q.all()

        def clicked(sp, pts):
            data = pts[0].data()
            print("-----------------------\nclicked:", data['rise_time'], data['amp'], data['prediction'], data['confidence'])
            for r in data['results']:
                print({k:r[k] for k in classifier.features})
            traces = data['traces']
            plt = pg.plot()
            bsub = [t.copy(data=t.data - np.median(t.time_slice(0, 1e-3).data)) for t in traces]
            for t in bsub:
                plt.plot(t.time_values, t.data, pen=(0, 0, 0, 50))
            mean = TSeriesList(bsub).mean()
            plt.plot(mean.time_values, mean.data, pen='g')


        # def analyze_response_strength(recs, source, dtype):
        #     results = []
        #     for i,rec in enumerate(recs):
        #         result = strength_analysis.analyze_response_strength(rec, source)
        #         results.append(result)
        #     return str_analysis_result_table(results)



        # measure background connection strength
        bg_results = [strength_analysis.analyze_response_strength(rec, 'baseline') for rec in bg_recs]
        bg_results = strength_analysis.str_analysis_result_table(bg_results, bg_recs)

        # for this example, we use background data to simulate foreground
        # (but this will be biased due to lack of crosstalk in background data)
        fg_recs = bg_recs

        # now measure foreground simulated under different conditions
        amps = 2e-6 * 2**np.arange(9)
        amps[0] = 0
        rtimes = [1e-3, 2e-3, 4e-3, 6e-3]
        dt = 1 / db.default_sample_rate
        results = np.empty((len(amps), len(rtimes)), dtype=[('results', object), ('predictions', object), ('confidence', object), ('traces', object), ('rise_time', float), ('amp', float)])
        print("  Simulating synaptic events..")

        cachefile = 'fig_3_cache.pkl'
        if os.path.exists(cachefile):
            cache = pickle.load(open(cachefile, 'rb'))
        else:
            cache = {}
        pair_key = (timestamp, pre_id, post_id)
        pair_cache = cache.setdefault(pair_key, {})

        for j,rtime in enumerate(rtimes):
            new_results = False
            for i,amp in enumerate(amps):
                print("---------------------------------------    %d/%d  %d/%d      \r" % (i,len(amps),j,len(rtimes)),)
                result = pair_cache.get((rtime, amp))
                if result is None:
                    result = strength_analysis.simulate_connection(fg_recs, bg_results, classifier, amp, rtime)
                    pair_cache[rtime, amp] = result
                    new_results = True

                for k,v in result.items():
                    results[i,j][k] = v

            x, y = amps, [np.mean(x) for x in results[:,j]['confidence']]
            c = limit_plot.plot(x, y, pen=pg.intColor(j, len(rtimes)*1.3, maxValue=150), symbol='o', antialias=True, name="%dus"%(rtime*1e6), data=results[:,j], symbolSize=4)
            write_csv(csv_file, x, "Figure 3D; {name}; {rise_time:0.3g} ms rise time; simulated PSP amplitude (V)".format(name=name, rise_time=rtime*1000))
            write_csv(csv_file, y, "Figure 3D; {name}; {rise_time:0.3g} ms rise time; classifier decision probability".format(name=name, rise_time=rtime*1000))
            c.scatter.sigClicked.connect(clicked)
            pg.QtWidgets.QApplication.processEvents()

            if new_results:
                pickle.dump(cache, open(cachefile, 'wb'))

        pg.QtWidgets.QApplication.processEvents()

    # Handle selected individual connections
    next_row = 0
    for name, key in show_conns:
        ts, pre_id, post_id = key
        add_connection_plots(next_row, name, ts, pre_id, post_id)
        next_row += 1

    def clicked(sp, pts):
        global next_row
        print("-----------")
        for pt in pts:
            d = pt.data()
            print(d['acq_timestamp'], d['pre_cell_id'], d['post_cell_id'])
            add_connection_plots(next_row, "", d['acq_timestamp'], d['pre_cell_id'], d['post_cell_id'])
            next_row += 1

    c_plot.scatter.sigClicked.connect(clicked)
    u_plot.scatter.sigClicked.connect(clicked)


csv_file.close()