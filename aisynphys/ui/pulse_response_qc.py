from __future__ import division, print_function

import time, pdb
import pyqtgraph as pg
import pyqtgraph.console
from pyqtgraph import parametertree
from pyqtgraph.parametertree import Parameter
from neuroanalysis.ui.user_test import UserTestUi

class PulseResponseQCUI(object):
    def __init__(self, title=None):
        self.pw = pg.GraphicsLayoutWidget()
        self.pw.addLabel(title)
        self.in_plt = self.pw.addPlot(row=1, col=0, title='Inhibitory Pulse QC')
        self.ex_plt = self.pw.addPlot(row=2, col=0, title='Excitatory Pulse QC')
        self.qc_text = pg.DataTreeWidget()
        self.qc_text.setHeaderLabels(['QC failures', 'type', 'value'])
        self.console = pg.console.ConsoleWidget()

        self.widget = pg.QtWidgets.QSplitter(pg.QtCore.Qt.Vertical)
        self.widget.addWidget(self.pw)        
        self.widget.addWidget(self.qc_text)
        self.widget.addWidget(self.console)
        self.widget.resize(1000, 900)
        self.widget.setSizes([800, 200])
        self.widget.show()

        self.trace_color = {True: 'g', False: 'r', None: 'w'}

    def clear(self):
        self.in_plt.clear()
        self.ex_plt.clear()

    def show_result(self, pulse_response, clamp_mode, ex_pass, in_pass, failures):
        pr = pulse_response['primary']
        pre_pulse = pr.time_slice(pr.t0, pr.t0 + 5e-3)
        base = pre_pulse.median()
        base_std = pre_pulse.std()
        pen = pg.mkPen('y', style=pg.QtCore.Qt.DotLine)
        for plt, qc_pass in zip([self.in_plt, self.ex_plt], [in_pass, ex_pass]):
            color = self.trace_color[qc_pass]
            plt.plot(pr.time_values, pr.data, pen=color)
            plt.addLine(x=pr.t0, pen=pen)
            plt.addLine(x=pr.t0 + 5e-3, pen=pen)
            # pdb.set_trace()
            if clamp_mode == 'vc':
                plt.addLine(y=base + base_std, pen=pen)
                plt.addLine(y=base - base_std, pen=pen)
                plt.addLine(y=base + 15e-12, pen='y')
                plt.addLine(y=base - 15e-12, pen='y')
                plt.setLabel('left', text='Voltage Clamp', units='A')
            if clamp_mode == 'ic':
                plt.addLine(y=base + base_std, pen=pen)
                plt.addLine(y=base - base_std, pen=pen)
                plt.addLine(y=base + 1.5e-3, pen='y')
                plt.addLine(y=base - 1.5e-3, pen='y')
                plt.addLine(y=-40e-3, pen='m')
                plt.addLine(y=-45e-3, pen='c')
                plt.setLabel('left', text='Current Clamp', units='V')
                plt.setYRange(-80e-3, -60e-3)
        if clamp_mode == 'ic':
            self.ex_plt.addLine(y=-85e-3, pen='c')
            self.in_plt.addLine(y=-60e-3, pen='c')

        # qc_fail_text = '\n'.join(failures)
        # self.qc_text.setValue(qc_fail_text)
        self.qc_text.setData(failures)


class PulseResponseQCTestUi(UserTestUi):
    """UI for manually pass/failing pulse response QC unit tests.
    """
    def __init__(self):
        expected_display = PulseResponseQCUI('expected result')
        current_display = PulseResponseQCUI('current result')
        UserTestUi.__init__(self, expected_display, current_display)

    def show_inputs(self, **inputs):
        self.inputs = inputs

    def show_results(self, expected, current):
        self.diff_widget.setData(expected, current)
        post_rec = self.inputs['post_rec']
        window = self.inputs['window']
        pulse_response = post_rec.time_slice(window[0], window[1])
        clamp_mode = post_rec.clamp_mode
        if expected is None:
            expected = (None, None, None)
        if current is None:
            current = (None, None, None)
        self.display2.show_result(pulse_response, clamp_mode, *current)
        self.display1.show_result(pulse_response, clamp_mode, *expected)
       
