import sys, copy, argparse, datetime
from collections import OrderedDict
from random import shuffle, seed
import numpy as np
import pyqtgraph as pg
from pyqtgraph import parametertree as ptree
from pyqtgraph.parametertree import Parameter
from pyqtgraph.widgets.DataFilterWidget import DataFilterParameter

from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import TSeriesList
from neuroanalysis.fitting import Psp, StackedPsp
from neuroanalysis.ui.fitting import FitExplorer

from aisynphys.ui.experiment_browser import ExperimentBrowser
from aisynphys.avg_response_fit import get_pair_avg_fits, response_query, sort_responses
from aisynphys.database import default_db as db
import aisynphys.data.data_notes_db as notes_db
from aisynphys.data import PulseResponseList
from aisynphys.avg_response_fit import fit_avg_pulse_response


default_latency = 1e-3
comment_hashtag = [
    '#doublespike',
    '#doublepsp',
    '#badspikes',
    '#fixable',
    '#needsecondopinion',
    '#lostcause',
    '#MVP',
    '#crosstalk',
    '#badqc',
    '#risetime',
    '#baseline',
    '#nodatatofit',
    '#polysynaptic',
    '#doublechecked',
    '#syncretypmismatch',
    '#badfit',
]

comment_hashtag.sort(key=lambda x:x[1])
comment_hashtag = [''] + comment_hashtag

modes = ['vc', 'ic']
holdings = [-55, -70]


class PairAnalysisWindow(pg.QtWidgets.QWidget):
    def __init__(self, default_session, notes_session):
        pg.QtWidgets.QWidget.__init__(self)
        self.default_session = default_session
        self.notes_session = notes_session
        self.layout = pg.QtWidgets.QGridLayout()
        self.setLayout(self.layout)
        self.h_splitter = pg.QtWidgets.QSplitter()
        self.h_splitter.setOrientation(pg.QtCore.Qt.Horizontal)
        self.pair_analyzer = PairAnalysis()
        self.fit_explorer = None
        self.ctrl_panel = self.pair_analyzer.ctrl_panel
        self.user_params = self.ctrl_panel.user_params
        self.output_params = self.ctrl_panel.output_params
        self.ptree = ptree.ParameterTree(showHeader=False)
        self.pair_param = Parameter.create(name='Current Pair', type='str', readonly=True)
        self.ptree.addParameters(self.pair_param)
        self.ptree.addParameters(self.user_params, showTop=False)
        self.fit_ptree = ptree.ParameterTree(showHeader=False)
        self.fit_ptree.addParameters(self.output_params, showTop=False)
        self.save_btn = pg.FeedbackButton('Save Analysis')
        self.expt_btn = pg.QtWidgets.QPushButton('Set Experiments')
        self.fit_btn = pg.QtWidgets.QPushButton('Fit Responses')
        self.ic_plot = self.pair_analyzer.ic_plot
        self.vc_plot = self.pair_analyzer.vc_plot
        self.select_ptree = ptree.ParameterTree(showHeader=False)
        self.hash_select = Parameter.create(name='Hashtags', type='group', children=
            [{'name': 'With multiple selected:', 'type': 'list', 'values': ['Include if any appear', 'Include if all appear'], 'value': 'Include if any appear'}]+
            [{'name': '#', 'type': 'bool'}] +
            [{'name': ht, 'type': 'bool'} for ht in comment_hashtag[1:]])
        self.rigs = db.query(db.Experiment.rig_name).distinct().all()
        self.operators = db.query(db.Experiment.operator_name).distinct().all()
        self.rig_select = Parameter.create(name='Rig', type='group', children=[{'name': rig[0], 'type': 'bool'} for rig in self.rigs])
        self.operator_select = Parameter.create(name='Operator', type='group', children=[{'name': operator[0], 'type': 'bool'} for operator in self.operators])
        self.data_type = Parameter.create(name='Reduce data to:', type='group', children=[
            {'name': 'Pairs with data', 'type': 'bool', 'value': True},
            {'name': 'Synapse is None', 'type': 'bool'}])
        [self.select_ptree.addParameters(param) for param in [self.data_type, self.rig_select, self.operator_select, self.hash_select]]
        self.experiment_browser = self.pair_analyzer.experiment_browser
        self.v_splitter = pg.QtWidgets.QSplitter()
        self.v_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.expt_splitter = pg.QtWidgets.QSplitter()
        self.expt_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.h_splitter.addWidget(self.expt_splitter)
        self.expt_splitter.addWidget(self.select_ptree)
        self.expt_splitter.addWidget(self.expt_btn)
        self.h_splitter.addWidget(self.v_splitter)
        self.v_splitter.addWidget(self.experiment_browser)
        self.v_splitter.addWidget(self.ptree)
        self.v_splitter.addWidget(self.fit_btn)
        self.v_splitter.addWidget(self.fit_ptree)
        self.v_splitter.addWidget(self.save_btn)
        self.v_splitter.setSizes([200, 20, 20,400, 20])
        # self.next_pair_button = pg.QtWidgets.QPushButton("Load Next Pair")
        # self.v_splitter.addWidget(self.next_pair_button)
        self.h_splitter.addWidget(self.vc_plot.grid)
        self.h_splitter.addWidget(self.ic_plot.grid)
        self.fit_compare = self.pair_analyzer.fit_compare
        self.meta_compare = self.pair_analyzer.meta_compare
        self.v2_splitter = pg.QtWidgets.QSplitter()
        self.v2_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.v2_splitter.addWidget(self.fit_compare)
        self.v2_splitter.addWidget(self.meta_compare)
        self.h_splitter.addWidget(self.v2_splitter)
        self.h_splitter.setSizes([100, 200, 200, 200, 500])
        self.layout.addWidget(self.h_splitter)
        self.fit_compare.hide()
        self.meta_compare.hide()
        self.setGeometry(280, 130, 1500, 900)
        self.show()

        # self.next_pair_button.clicked.connect(self.load_next_pair)
        self.experiment_browser.itemSelectionChanged.connect(self.selected_pair)
        self.save_btn.clicked.connect(self.save_to_db)
        self.expt_btn.clicked.connect(self.get_expts)
        self.fit_btn.clicked.connect(self.pair_analyzer.fit_response_update)

    def save_to_db(self):
        try:
            self.pair_analyzer.save_to_db()
            self.save_btn.success()
        except:
            self.save_btn.failure()
            raise

    def get_expts(self):
        expt_query = db.query(db.Experiment)
        synapse_none = self.data_type['Synapse is None']
        if synapse_none:
            subquery = db.query(db.Pair.experiment_id).filter(db.Pair.has_synapse==None).subquery()
            expt_query = expt_query.filter(db.Experiment.id.in_(subquery))
        selected_rigs = [rig.name() for rig in self.rig_select.children() if rig.value() is True]
        if len(selected_rigs) != 0:
            expt_query = expt_query.filter(db.Experiment.rig_name.in_(selected_rigs))
        selected_operators = [operator.name() for operator in self.operator_select.children() if operator.value()is True]
        if len(selected_operators) != 0:
            expt_query = expt_query.filter(db.Experiment.operator_name.in_(selected_operators))
        selected_hashtags = [ht.name() for ht in self.hash_select.children()[1:] if ht.value() is True]
        if len(selected_hashtags) != 0:
            timestamps = self.get_expts_hashtag(selected_hashtags)
            expt_query = expt_query.filter(db.Experiment.ext_id.in_(timestamps))
        expts = expt_query.all()
        self.set_expts(expts)

    def get_expts_hashtag(self, selected_hashtags):
        q = self.notes_session.query(notes_db.PairNotes)
        pairs_to_include = []
        note_pairs = q.all()
        note_pairs.sort(key=lambda p: p.expt_id)
        for p in note_pairs:
            comments = p.notes.get('comments')
            if comments is None:
                continue    
            if len(selected_hashtags) == 1:
                hashtag = selected_hashtags[0]
                if hashtag == '#':
                    if hashtag in comments and all([ht not in comments for ht in comment_hashtag[1:]]):
                        print(p.expt_id, p.pre_cell_id, p.post_cell_id, comments)
                        pairs_to_include.append(p)
                else:
                    if hashtag in comments:
                        print(p.expt_id, p.pre_cell_id, p.post_cell_id, comments)
                        pairs_to_include.append(p)
                
            if len(selected_hashtags) > 1:
                hashtag_present = [ht in comments for ht in selected_hashtags]
                or_expts = self.hash_select['With multiple selected:'] == 'Include if any appear'
                and_expts = self.hash_select['With multiple selected:'] == 'Include if all appear'
                if or_expts and any(hashtag_present):
                    print(p.expt_id, p.pre_cell_id, p.post_cell_id, comments)
                    pairs_to_include.append(p)
                if and_expts and all(hashtag_present):
                    print(p.expt_id, p.pre_cell_id, p.post_cell_id, comments)
                    pairs_to_include.append(p)

        return set([pair.expt_id for pair in pairs_to_include])

    def set_expts(self, expts):
        with pg.BusyCursor():
            self.experiment_browser.clear()
            has_data = self.data_type['Pairs with data']
            if not has_data:
                self.experiment_browser.populate(experiments=expts, all_pairs=True)
            else:
                self.experiment_browser.populate(experiments=expts)

    def selected_pair(self):
        with pg.BusyCursor():
            self.fit_compare.hide()
            self.meta_compare.hide()
            selected = self.experiment_browser.selectedItems()
            if len(selected) != 1:
                return
            item = selected[0]
            if hasattr(item, 'pair') is False:
                return
            pair = item.pair
            ## check to see if the pair has already been analyzed
            expt_id = pair.experiment.ext_id
            pre_cell_id = pair.pre_cell.ext_id
            post_cell_id = pair.post_cell.ext_id
            record = notes_db.get_pair_notes_record(expt_id, pre_cell_id, post_cell_id, session=self.notes_session)
            
            self.pair_param.setValue(pair)
            if record is None:
                self.pair_analyzer.load_pair(pair, self.default_session)
                self.pair_analyzer.analyze_responses()
                # self.pair_analyzer.fit_responses()
            else:
                self.pair_analyzer.load_pair(pair, self.default_session, record=record)
                self.pair_analyzer.analyze_responses()
                self.pair_analyzer.load_saved_fit(record)

    def explore_fit(self, mode, holding):
        fit = self.pair_analyzer.last_fit[mode, holding]
        self.fit_explorer = FitExplorer(fit)
        self.fit_explorer.show()


class ControlPanel(object):
    def __init__(self):
        self.user_latency = Parameter.create(name='User Latency', type='float', suffix='s', siPrefix=True, dec=True, value=default_latency)
        self.synapse = Parameter.create(name='Synapse call', type='list', values={'Excitatory': 'ex', 'Inhibitory': 'in', 'None': None})
        self.polysynaptic = Parameter.create(name='Polysynaptic call', type='list', values={'Excitatory': 'ex', 'Inhibitory': 'in', 'Mixed': 'mix', 'None': None})
        self.gap = Parameter.create(name='Gap junction call', type='bool')
        fit_cat = ['-55 VC', '-70 VC', '-55 IC', '-70 IC']
        fit_param = [
            {'name': c, 'type': 'group', 'children': [
                {'name': 'Parameter:', 'type': 'str', 'readonly': True, 'value': 'Amplitude, Latency, Rise time, Decay tau, NRMSE'},
                {'name': 'Value:', 'type': 'str', 'readonly': True},
                {'name': 'Fit Pass', 'type': 'bool'},
            ]} for c in fit_cat
        ]
        self.fit_params = Parameter.create(name='Fit parameters', type='group', children=fit_param)
        self.warn_param = Parameter.create(name='Warnings', type='text', readonly=True)
        self.comments = Parameter.create(name='Comments', type='group', children=[
            {'name': 'Hashtag', 'type': 'list', 'values': comment_hashtag, 'value': ''},
            {'name': '', 'type': 'text'}
        ])
        self.user_params = Parameter.create(name='user_params', type='group', children=[
            self.user_latency,
            self.synapse,
            self.polysynaptic,
            self.gap,
        ])
        self.output_params = Parameter.create(name='output_params', type='group', children=[
            self.fit_params,
            self.warn_param,
            self.comments,
        ])

        self.comments.child('Hashtag').sigValueChanged.connect(self.add_text_to_comments)

    def add_text_to_comments(self):
        text = self.comments['Hashtag']
        comments = self.comments['']
        update_comments = comments + text + '\n'
        self.comments.child('').setValue(update_comments)
        
    def update_user_params(self, **kargs):
        for k, v in kargs.items():
            self.user_params.child(k).setValue(v)

    def update_fit_params(self, fit_params, fit_pass=False):
        param_names = ['amp', 'xoffset', 'rise_time', 'decay_tau', 'nrmse']
        for mode in modes:
            if mode == 'vc':
                suffix = ['A', 's', 's', 's', '']
            elif mode == 'ic':
                suffix = ['V', 's', 's', 's', '']
            for holding in holdings:
                group = self.output_params.child('Fit parameters', str(holding) + ' ' + mode.upper())
                values = fit_params[mode][str(holding)]
                format_values = self.format_fit_output(values, param_names, suffix)
                group.child('Value:').setValue(format_values)
                if fit_pass:
                    group.child('Fit Pass').setValue(False)


    def format_fit_output(self, values, name, suffix):
        format_list = []
        for p in zip(name, suffix):
            if values.get(p[0]) is None:
                format_list.append('nan')
            else:
                value = values[p[0]]
                if p[0] == 'nrmse':
                    p_format = ('%0.2f' % value)
                else:
                    p_format = pg.siFormat(value, suffix=p[1])
                format_list.append(p_format)
        output = ", ".join(format_list)
        return output

    def set_latency(self, superline):
        value = superline.pos()
        self.user_params.child('User Latency').setValue(value)



class TSeriesPlot(pg.GraphicsLayoutWidget):
    def __init__(self, title, units):
        pg.GraphicsLayoutWidget.__init__(self)
        self.grid = PlotGrid()
        self.grid.set_shape(4, 1)
        self.grid.grid.ci.layout.setRowStretchFactor(0, 3)
        self.grid.grid.ci.layout.setRowStretchFactor(1, 8)
        self.grid.grid.ci.layout.setRowStretchFactor(2, 5)
        self.grid.grid.ci.layout.setRowStretchFactor(3, 10)
        self.grid.show()
        self.trace_plots = (self.grid[1, 0], self.grid[3, 0])
        self.spike_plots = (self.grid[0, 0], self.grid[2, 0])
        self.plots = self.spike_plots + self.trace_plots
        for plot in self.plots[:-1]:
            plot.hideAxis('bottom')
        self.plots[-1].setLabel('bottom', text='Time from spike', units='s')
        self.fit_item_55 = None
        self.fit_item_70 = None
        self.fit_color = {True: 'g', False: 'r'}
        self.qc_color = {'qc_pass': (255, 255, 255, 100), 'qc_fail': (255, 0, 0, 100)}

        self.plots[0].setTitle(title)
        for (plot, holding) in zip(self.trace_plots, holdings):
            plot.setXLink(self.plots[-1])
            plot.setLabel('left', text="%d holding" % int(holding), units=units)
        for plot in self.spike_plots:
            plot.setXLink(self.plots[-1])
            plot.setLabel('left', text="presynaptic spike")
            plot.addLine(x=0)
        self.plots[-1].setXRange(-5e-3, 10e-3)
        
        self.items = []

    def plot_responses(self, pulse_responses):
        self.plot_traces(pulse_responses)
        self.plot_spikes(pulse_responses)

    def plot_traces(self, pulse_responses):  
        for i, holding in enumerate(pulse_responses.keys()):
            for qc, prs in pulse_responses[holding].items():
                if len(prs) == 0:
                    continue
                prl = PulseResponseList(prs)
                
                for mode in ['spike', 'peak', 'stim']:
                    try:
                        post_ts = prl.post_tseries(align=mode, bsub=True, alignment_failure_mode='average')
                        break
                    except Exception:
                        post_ts = None
                        continue
                if post_ts is None:
                    continue        
                for trace in post_ts:
                    item = self.trace_plots[i].plot(trace.time_values, trace.data, pen=self.qc_color[qc])
                    if qc == 'qc_fail':
                        item.setZValue(-10)
                    self.items.append(item)
                if qc == 'qc_pass':
                    grand_trace = post_ts.mean()
                    item = self.trace_plots[i].plot(grand_trace.time_values, grand_trace.data, pen={'color': 'b', 'width': 2})
                    self.items.append(item)
            self.trace_plots[i].autoRange()
            self.trace_plots[i].setXRange(-5e-3, 10e-3)
            # y_range = [grand_trace.data.min(), grand_trace.data.max()]
            # self.plots[i].setYRange(y_range[0], y_range[1], padding=1)

    def plot_spikes(self, pulse_responses):
        for i, holding in enumerate(pulse_responses.keys()):
            for prs in pulse_responses[holding].values():
                if len(prs) == 0:
                    continue
                prl = PulseResponseList(prs)
                for mode in ['spike', 'peak', 'stim']:
                    try:
                        pre_ts = prl.pre_tseries(align=mode, bsub=True, alignment_failure_mode='average')
                        break
                    except Exception:
                        pre_ts = None
                        continue
                if pre_ts is None:
                    continue
                for pr, spike in zip(prl, pre_ts):
                    # pr.stim_pulse.n_spikes can == 1 but the spike time (ie max slope) is None, failing
                    # the postsynaptic responses. Consider using pr.stim_pulse.first_spike_time != None
                    # and qc failing these spikes as the traces are as well.
                    qc = 'qc_pass' if pr.stim_pulse.first_spike_time is not None else 'qc_fail'
                    
                    item = self.spike_plots[i].plot(spike.time_values, spike.data, pen=self.qc_color[qc])
                    if qc == 'qc_fail':
                        item.setZValue(-10)
                    self.items.append(item)

    def plot_fit(self, trace, holding, fit_pass=False):
        if holding == -55:
            if self.fit_item_55 is not None:
                self.trace_plots[0].removeItem(self.fit_item_55)
            self.fit_item_55 = pg.PlotDataItem(trace.time_values, trace.data, name='-55 holding', pen={'color': self.fit_color[fit_pass], 'width': 3})
            self.trace_plots[0].addItem(self.fit_item_55)
        
        elif holding == -70:
            if self.fit_item_70 is not None:
                self.trace_plots[1].removeItem(self.fit_item_70)
            self.fit_item_70 = pg.PlotDataItem(trace.time_values, trace.data, name='-70 holding', pen={'color': self.fit_color[fit_pass], 'width': 3})
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
        

class SuperLine(pg.QtCore.QObject):
    sigPositionChanged = pg.QtCore.Signal(object)
    sigPositionChangeFinished = pg.QtCore.Signal(object)

    def __init__(self):
        pg.QtCore.QObject.__init__(self)
        self.lines = []

    def new_line(self, x_pos):
        self.line = pg.InfiniteLine(x_pos, pen={'color': 'y', 'width': 3}, movable=True)
        self.line.setZValue(100)
        self.lines.append(self.line)
        self.line.sigPositionChanged.connect(self.line_sync)
        self.line.sigPositionChangeFinished.connect(self.move_finished)
        return self.line

    def move_finished(self):
        self.sigPositionChangeFinished.emit(self)

    def line_sync(self, moved_line):
        for line in self.lines:
            with pg.SignalBlock(line.sigPositionChanged, self.line_sync):
                line.setValue(moved_line.value())
        self.sigPositionChanged.emit(self)

    def pos(self):
        line_positions = []
        for line in self.lines:
            line_positions.append(line.value())
        n_positions = set(line_positions)
        if len(n_positions) > 1:
            raise Exception("Lines are out of sync and reporting different positions")
        position = list(n_positions)[0]

        return position
        
    def set_value(self, value, block_fit=False):
        for line in self.lines:
            with pg.SignalBlock(line.sigPositionChanged, self.line_sync):
                line.setValue(value or 1e-3)
        self.sigPositionChanged.emit(self)
        if block_fit is False:
            self.sigPositionChangeFinished.emit(self)

    def set_value_from_ctrl_panel(self, ctrl_panel):
        value = ctrl_panel.value()
        self.set_value(value, block_fit=True)


class PairAnalysis(object):
    def __init__(self):
        self.ctrl_panel = ControlPanel()
        
        self.latency_superline = SuperLine()
        self.latency_superline.sigPositionChanged.connect(self.ctrl_panel.set_latency)
       
        self.ic_plot = TSeriesPlot('Current Clamp', 'V')
        for plot in self.ic_plot.trace_plots:
            plot.addItem(self.latency_superline.new_line(default_latency))

        self.vc_plot = TSeriesPlot('Voltage Clamp', 'A')
        for plot in self.vc_plot.trace_plots:
            plot.addItem(self.latency_superline.new_line(default_latency))

        self.user_latency = self.ctrl_panel.user_latency
        self.user_latency.sigValueChanged.connect(self.latency_superline.set_value_from_ctrl_panel)
            
        self.ctrl_panel.output_params.child('Fit parameters').sigTreeStateChanged.connect(self.colorize_fit)
        self.experiment_browser = ExperimentBrowser()
        self.fit_compare = pg.DiffTreeWidget()
        self.meta_compare = pg.DiffTreeWidget()
        self.nrmse_thresh = 4
        self.sorted_responses = None
        self.signs = {
            'vc': {
                -55: {'ex': -1, 'in': 1}, 
                -70: {'ex': -1, 'in': 0},
            },
            'ic':{
                -55: {'ex': 1, 'in': -1},
                -70: {'ex': 1, 'in': 0},
            },
        }

        self.fit_precision = {
            'amp': {'vc': 14, 'ic': 8},
            'exp_amp': {'vc': 14, 'ic': 8},
            'exp_tau': {'vc': 8, 'ic': 8},
            'decay_tau': {'vc': 8, 'ic': 8},
            'nrmse': {'vc': 2, 'ic': 2},
            'rise_time': {'vc': 7, 'ic': 6},
            'rise_power': {'vc': 0, 'ic': 0},
            'xoffset': {'vc': 7, 'ic': 7},
            'yoffset': {'vc': 5, 'ic': 14},
        }

    def colorize_fit(self, param, changes):
        for c in changes:
            p, change, info = c
            if p.name() != 'Fit Pass':
                continue
            if 'VC' in p.parent().name():
                self.vc_plot.color_fit(p.parent().name(), info)
            if 'IC' in p.parent().name():
                self.ic_plot.color_fit(p.parent().name(), info)

    def reset_display(self):
        self.vc_plot.clear_plots()
        self.ic_plot.clear_plots()
        self.latency_superline.set_value(default_latency, block_fit=True)
        self.ctrl_panel.output_params.child('Comments', 'Hashtag').setValue('')
        self.ctrl_panel.output_params.child('Comments', '').setValue('')
        self.ctrl_panel.output_params.child('Warnings').setValue('')
        
    def load_pair(self, pair, default_session, record=None):
        self.record = record
        self.initial_fit_parameters = OrderedDict([
            ('vc', {'-55': {}, '-70': {}}), 
            ('ic', {'-55': {}, '-70': {}}),
        ])
        self.output_fit_parameters = OrderedDict([
            ('vc', {'-55': {}, '-70': {}}), 
            ('ic', {'-55': {}, '-70': {}}),
        ])
        self.fit_params = {'initial': self.initial_fit_parameters, 'fit': self.output_fit_parameters}

        with pg.BusyCursor():
            self.reset_display()
            self.pair = pair
            print ('loading responses for %s...' % pair)
            q = response_query(default_session, pair)
            self.pulse_responses = [q.PulseResponse for q in q.all()]
            print('got %d pulse responses' % len(self.pulse_responses))
                
            if pair.has_synapse is True and pair.synapse is not None:
                synapse_type = pair.synapse.synapse_type
            else:
                synapse_type = None
            pair_params = {'Synapse call': synapse_type, 'Gap junction call': pair.has_electrical}
            self.ctrl_panel.update_user_params(**pair_params)
            self.ctrl_panel.update_fit_params(self.fit_params['fit'], fit_pass=True)

    def analyze_responses(self):
        self.sorted_responses = sort_responses(self.pulse_responses)
        got_fitable_responses = False
        for mode in modes:
            for holding in holdings:
                got_fitable_responses = got_fitable_responses or len(self.sorted_responses[mode, holding]['qc_pass']) > 0
        if not got_fitable_responses:
            print('No fitable responses, bailing out')
        
        self.vc_plot.plot_responses({holding: self.sorted_responses['vc', holding] for holding in holdings})
        self.ic_plot.plot_responses({holding: self.sorted_responses['ic', holding] for holding in holdings})    
      
    def fit_response_update(self):
        latency = self.ctrl_panel.user_params['User Latency']
        self.fit_responses(latency=latency)

    def fit_responses(self, latency=None):
        if latency is None:
            latency_window = [0.5e-3, 10e-3]
        else:
            latency_window = [latency-100e-6, latency+100e-6]

        with pg.ProgressDialog("curve fitting..", maximum=len(modes)*len(holdings)) as dlg:
            self.last_fit = {}
            for mode in modes:
                for holding in holdings:
                    self.fit_pass = False
                    sign = self.signs[mode][holding].get(self.ctrl_panel.user_params['Synapse call'], 
                            self.signs[mode][holding].get(self.ctrl_panel.user_params['Polysynaptic call'], 0))
                    
                    # ofp, x_offset, best_fit = fit_avg_response(self.traces, mode, holding, latency, sign)
                    prs = self.sorted_responses[mode, holding]['qc_pass']
                    if len(prs) == 0:
                        dlg += 1
                        continue
                    
                    fit, avg = fit_avg_pulse_response(prs, latency_window=latency_window, sign=sign)
                    fit_ts = avg.copy(data=fit.best_fit)
                    self.last_fit[mode, holding] = fit
                    
                    self.initial_fit_parameters[mode][str(holding)]['xoffset'] = latency
                    self.output_fit_parameters[mode][str(holding)]['nrmse'] = fit.nrmse()
                    self.output_fit_parameters[mode][str(holding)].update(fit.best_values)
                    self.fit_pass = fit.nrmse() < self.nrmse_thresh
                    self.ctrl_panel.output_params.child('Fit parameters', str(holding) + ' ' + mode.upper(), 'Fit Pass').setValue(self.fit_pass)
                    if mode == 'vc':
                        self.vc_plot.plot_fit(fit_ts, holding, self.fit_pass)
                    elif mode == 'ic':
                        self.ic_plot.plot_fit(fit_ts, holding, self.fit_pass)
                    dlg += 1
                    if dlg.wasCanceled():
                        raise Exception("User canceled fit")
        self.fit_params['initial'] = self.initial_fit_parameters
        self.fit_params['fit'] = self.output_fit_parameters

        self.ctrl_panel.update_fit_params(self.fit_params['fit'])
        self.generate_warnings() 

    def generate_warnings(self):
        self.warnings = []
        latency_mode = []
        for mode in modes:
            latency_holding = []
            for holding in holdings:
                initial_latency = self.fit_params['initial'][mode][str(holding)].get('xoffset')
                fit_latency = self.fit_params['fit'][mode][str(holding)].get('xoffset')
                if fit_latency is None or initial_latency is None:
                    continue
                if abs(np.diff([fit_latency, initial_latency])) > 0.1e-3:
                    warning = 'Initial latency and fit latency differ by %s' % pg.siFormat(abs(np.diff([fit_latency, initial_latency])), suffix='s')
                latency_holding.append(fit_latency)
            if len(latency_holding) == 2 and abs(np.diff(latency_holding)) > 0.01e-3:
                warning = 'Fit latencies for %s mode do not match' % mode
                self.warnings.append(warning)
            latency_mode.append(np.mean(latency_holding))
        latency_diff = np.diff(latency_mode)[0]
        if  abs(latency_diff) > 0.2e-3:
            warning = 'Latency across modes differs by %s' % pg.siFormat(latency_diff, suffix='s')
            self.warnings.append(warning)

        min_latency = np.nanmin(latency_mode) < 0.4e-3
        max_latency = np.nanmax(latency_mode) > 3.1e-3
        
        if min_latency and self.ctrl_panel.user_params['Gap junction call'] is False:
            self.warnings.append("Short latency; is this a gap junction?") 


        guess_sign = []
        guess = None
        for mode, fits1 in self.output_fit_parameters.items():
            for holding, fit in fits1.items():
                if 'amp' not in fit:
                    continue
                if mode == 'ic':
                    guess = 1 if fit['amp'] > 0 else -1
                elif mode == 'vc':
                    guess = -1 if fit['amp'] > 0 else 1
                guess_sign.append(guess)
        if np.all(np.array(guess) == 1):
            guess = "ex"
        elif np.all(np.array(guess) == -1):
            guess = "in"
        if guess is None:
            self.warnings.append("Mixed amplitude signs; pick ex/in carefully.")
        elif guess != self.ctrl_panel.user_params['Synapse call'] and not min_latency and not max_latency:
            self.warnings.append("Looks like an %s synapse??" % guess)
        elif self.ctrl_panel.user_params['Polysynaptic call'] != guess and max_latency:
            self.warnings.append("Looks like a %s polysynaptic synpase??" % guess)

        print_warning = '\n'.join(self.warnings)
        self.ctrl_panel.output_params.child('Warnings').setValue(print_warning)

    def save_to_db(self):
        fit_pass = {}
        for mode in modes:
            fit_pass[mode] = {}
            for holding in holdings:
                fit_pass[mode][str(holding)] = self.ctrl_panel.output_params['Fit parameters', str(holding) + ' ' + mode.upper(), 'Fit Pass']

        expt_id = self.pair.experiment.ext_id
        pre_cell_id = self.pair.pre_cell.ext_id
        post_cell_id = self.pair.post_cell.ext_id
        meta = {
            'expt_id': expt_id,
            'pre_cell_id': pre_cell_id,
            'post_cell_id': post_cell_id,
            'synapse_type': self.ctrl_panel.user_params['Synapse call'],
            'polysynaptic_type': self.ctrl_panel.user_params['Polysynaptic call'],
            'gap_junction': self.ctrl_panel.user_params['Gap junction call'],
            'comments': self.ctrl_panel.output_params['Comments', ''],
        }

        if self.ctrl_panel.user_params['Synapse call'] is not None or self.ctrl_panel.user_params['Polysynaptic call'] is not None:
            meta.update({
            'fit_parameters': self.fit_params,
            'fit_pass': fit_pass,
            'fit_warnings': self.warnings,
            })
        
        session = notes_db.db.session(readonly=False)
        record = notes_db.get_pair_notes_record(expt_id, pre_cell_id, post_cell_id, session=session)

        if record is None:
            entry = notes_db.PairNotes(
                expt_id=expt_id,
                pre_cell_id=pre_cell_id,
                post_cell_id=post_cell_id, 
                notes=meta,
                modification_time=datetime.datetime.now(),
            )
            session.add(entry)
            session.commit()
        else:
            self.print_pair_notes(meta, record)
            msg = pg.QtWidgets.QMessageBox.question(None, "Pair Analysis", 
                "The record you are about to save conflicts with what is in the Pair Notes database.\nYou can see the differences highlighted in red.\nWould you like to overwrite?",
                pg.QtWidgets.QMessageBox.Yes | pg.QtWidgets.QMessageBox.No)
            if msg == pg.QtWidgets.QMessageBox.Yes:
                record.notes = meta
                record.modification_time = datetime.datetime.now()
                session.commit() 
            else:
                raise Exception('Save Cancelled')
        session.close()

    def print_pair_notes(self, meta, saved_rec):
        meta_copy = copy.deepcopy(meta)
        current_fit = {k:v for k, v in meta_copy.get('fit_parameters', {}).get('fit', {}).items()}
        saved_fit = {k:v for k, v in saved_rec.notes.get('fit_parameters',{}).get('fit', {}).items()}
        
        for mode in modes:
            for holding in holdings:
                if not current_fit:
                    continue
                current_fit[mode][str(holding)] = {k:round(v, self.fit_precision[k][mode]) for k, v in current_fit[mode][str(holding)].items()}
                if not saved_fit:
                    continue
                saved_fit[mode][str(holding)] = {k:round(v, self.fit_precision[k][mode]) for k, v in saved_fit[mode][str(holding)].items()}

        self.fit_compare.setData(current_fit, saved_fit)
        self.fit_compare.trees[0].setHeaderLabels(['Current Fit Parameters', 'type', 'value'])
        self.fit_compare.trees[1].setHeaderLabels(['Saved Fit Parameters', 'type', 'value'])

        current_meta = {k:v for k, v in meta_copy.items() if k != 'fit_parameters'} 
        saved_meta = {k:v for k, v in saved_rec.notes.items() if k != 'fit_parameters'} 
        saved_meta.update({k:str(saved_meta[k]) for k in ['comments', 'expt_id', 'pre_cell_id', 'post_cell_id'] if k in saved_meta.keys()})

        self.meta_compare.setData(current_meta, saved_meta)
        self.meta_compare.trees[0].setHeaderLabels(['Current Metadata', 'type', 'value'])
        self.meta_compare.trees[1].setHeaderLabels(['Saved Metadata', 'type', 'value'])

        self.fit_compare.show()
        self.meta_compare.show()

    def load_saved_fit(self, record):
        data = record.notes        
        pair_params = {'Synapse call': data.get('synapse_type', None), 'Polysynaptic call': data.get('polysynaptic_type', None), 'Gap junction call': data.get('gap_junction', False)}
        self.ctrl_panel.update_user_params(**pair_params)
        self.warnings = data.get('fit_warnings', [])
        self.ctrl_panel.output_params.child('Warnings').setValue('\n'.join(self.warnings))
        self.ctrl_panel.output_params.child('Comments', '').setValue(data.get('comments', ''))

        # some records may be stored with no fit if a synapse is not present.
        if 'fit_parameters' not in data:
            return
            
        self.fit_params = data['fit_parameters']
        self.ctrl_panel.update_fit_params(data['fit_parameters']['fit'])
        self.output_fit_parameters = data['fit_parameters']['fit']
        self.initial_fit_parameters = data['fit_parameters']['initial']

        initial_vc_latency = (
            data['fit_parameters']['initial']['vc']['-55'].get('xoffset') or
            data['fit_parameters']['initial']['vc']['-70'].get('xoffset') or
            1e-3
        )
        initial_ic_latency = (
            data['fit_parameters']['initial']['ic']['-55'].get('xoffset') or
            data['fit_parameters']['initial']['ic']['-70'].get('xoffset') or
            1e-3
        )

        latency_diff = np.diff([initial_vc_latency, initial_ic_latency])[0]
        if abs(latency_diff) < 100e-6:
            self.latency_superline.set_value(initial_vc_latency, block_fit=True)
        else:
            fit_pass_vc = [data['fit_pass']['vc'][str(h)] for h in holdings]
            fit_pass_ic = [data['fit_pass']['ic'][str(h)] for h in holdings]
            if any(fit_pass_vc):
                self.latency_superline.set_value(initial_vc_latency, block_fit=True)
            elif any(fit_pass_ic):
                self.latency_superline.set_value(initial_ic_latency, block_fit=True)
            else:
                self.latency_superline.set_value(initial_vc_latency, block_fit=True)

        for mode in modes:
            for holding in holdings:
                fit_pass = data['fit_pass'][mode][str(holding)]
                self.ctrl_panel.output_params.child('Fit parameters', str(holding) + ' ' + mode.upper(), 'Fit Pass').setValue(fit_pass)
                fit_params = copy.deepcopy(data['fit_parameters']['fit'][mode][str(holding)])
                
                if fit_params:
                    fit_params.pop('nrmse', None)
                    if mode == 'ic':
                        p = StackedPsp() 
                    if mode == 'vc':
                        p = Psp()
                        
                    # make a list of spike-aligned postsynaptic tseries
                    tsl = PulseResponseList(self.sorted_responses[mode, holding]['qc_pass']).post_tseries(align='spike', bsub=True)
                    if len(tsl) == 0:
                        continue
                    
                    # average all together
                    avg = tsl.mean()
                    
                    fit_params.setdefault('exp_tau', fit_params['decay_tau'])
                    fit_psp = p.eval(x=avg.time_values, **fit_params)
                    fit_tseries = avg.copy(data=fit_psp)
                    
                    if mode == 'vc':
                        self.vc_plot.plot_fit(fit_tseries, holding, fit_pass=fit_pass)
                    if mode == 'ic':
                        self.ic_plot.plot_fit(fit_tseries, holding, fit_pass=fit_pass)            



if __name__ == '__main__':
    app = pg.mkQApp()
    parser = argparse.ArgumentParser()
    parser.add_argument('--timestamps', type=str, nargs='*')
    parser.add_argument('--dbg', default=False, action='store_true')
    parser.add_argument('expt_id', type=str, nargs='?', default=None)
    parser.add_argument('pre_cell_id', type=str, nargs='?', default=None)
    parser.add_argument('post_cell_id', type=str, nargs='?', default=None)

    args = parser.parse_args(sys.argv[1:])

    if args.dbg:
        pg.dbg()

    default_session = db.session()
    notes_session = notes_db.db.session()
    # timestamps = [r.acq_timestamp for r in db.query(db.Experiment.acq_timestamp).all()]
    timestamps = []
    
    mw = MainWindow(default_session, notes_session)
    if args.timestamps is not None:
        timestamps = args.timestamps
    elif args.expt_id is not None:
        timestamps = [args.expt_id]
    
    q = default_session.query(db.Experiment).filter(db.Experiment.acq_timestamp.in_(timestamps))
    expts = q.all()
    mw.set_expts(expts)

    if None not in (args.expt_id, args.pre_cell_id, args.post_cell_id):
        expt = db.experiment_from_ext_id(args.expt_id)
        pair = expt.pairs[args.pre_cell_id, args.post_cell_id]
        mw.experiment_browser.select_pair(pair.id)

    if sys.flags.interactive == 0:
        app.exec_()
