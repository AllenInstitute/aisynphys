import pickle, sys
import numpy as np
import pyqtgraph as pg
from aisynphys.data import MultiPatchSyncRecAnalyzer, MultiPatchProbe
from aisynphys.ui.pulse_response_qc import PulseResponseQCUI
from aisynphys.qc import PulseResponseQCTestCase, pulse_response_qc_pass
from aisynphys.database import default_db as db

pg.dbg()

expt_id = sys.argv[1]
pre_cell_id = sys.argv[2]
post_cell_id = sys.argv[3]

ui = PulseResponseQCUI()
skip_btn = pg.QtWidgets.QPushButton('skip')
ui.widget.addWidget(skip_btn)
save_btn = pg.QtWidgets.QPushButton('save')
ui.widget.addWidget(save_btn)

def iter_resps():
    expt = db.experiment_from_timestamp(expt_id)
    pre_cell = expt.cells[pre_cell_id]
    pre_dev = pre_cell.electrode.device_id
    post_cell = expt.cells[post_cell_id]
    post_dev = post_cell.electrode.device_id
    sweeps = expt.data.contents

    for srec in sweeps:
        try:
            post_rec = srec[post_dev]
        except KeyError:
            continue
        if not isinstance(post_rec, MultiPatchProbe):
            continue
        pre_rec = srec[pre_dev]
        mpa = MultiPatchSyncRecAnalyzer(srec)
        responses = mpa.get_spike_responses(pre_rec, post_rec, align_to='pulse', require_spike=False)
        for resp in responses:
            window = [resp['rec_start'], resp['rec_stop']]
            spike = None if len(resp['spikes']) == 0 else resp['spikes'][0]
            n_spikes = 0 if spike is None else 1        
            yield (srec, post_rec, window, n_spikes) #this will feed into qc.pulse_response_qc_pass

all_resps = iter_resps()
last_result = None

def load_next():
    global all_resps, last_result, ui
    try:
        (srec, post_rec, window, n_spikes) = next(all_resps)
    except StopIteration:
        ui.widget.hide()
        return

    ex_qc_pass, in_qc_pass, qc_failures = pulse_response_qc_pass(post_rec=post_rec, window=window, n_spikes=n_spikes, adjacent_pulses=[])
    pulse_response = post_rec.time_slice(window[0], window[1])
    clamp = post_rec.clamp_mode
    ui.clear()
    ui.show_result(pulse_response, clamp, ex_qc_pass, in_qc_pass, qc_failures)

    tc = PulseResponseQCTestCase()
    tc._meta = {
        'expt_id': expt_id,
        'post_cell_id': post_cell_id,
        'sweep_id': srec.key,
    }
    tc._input_args = {
        'post_rec': post_rec,
        'window': window,
        'n_spikes': n_spikes,
        'adjacent_pulses': [],
    }

    last_result = tc

def save_and_load_next():
    global last_result

    # write results out to test file
    test_file = 'test_data/pulse_response_qc/%s.pkl' % (last_result.name)
    last_result.save_file(test_file)

    load_next()


skip_btn.clicked.connect(load_next)
save_btn.clicked.connect(save_and_load_next)
load_next()