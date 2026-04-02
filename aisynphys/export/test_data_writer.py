import os
import h5py
import numpy as np
import pytest
from pathlib import Path
from unittest.mock import MagicMock
from aisynphys.export.data_writer import ExperimentWriter

class MockRecording:
    def __init__(self, data, units):
        self.data = data
        self.units = units

class MockSweep:
    def __init__(self, key, devices, recordings, meta):
        self.key = key
        self.devices = devices
        self.recordings = recordings
        self.recording_dict = {d: MagicMock(meta=meta) for d in devices}
        # Mocking the recording object to have 'primary' attribute
        for d in devices:
            self.recordings[devices.index(d)] = {
                'primary': MockRecording(recordings[devices.index(d)], 'V' if 'V' in str(recordings[devices.index(d)]) else 'A')
            }

@pytest.fixture
def mock_expt():
    expt = MagicMock()
    expt.id = 12345
    expt.ext_id = "expt_ext_12345"
    expt.date = "2023-01-01"
    expt.target_region = "V1"

    # Mock electrodes
    elec1 = MagicMock()
    elec1.device_id = 0
    elec1.cell.id = 101
    elec1.cell.ext_id = "cell_ext_101"
    elec1.cell.electrode = elec1

    expt.electrodes = {0: elec1}
    expt.cells = {"1": MagicMock()} # Electrode 0 corresponds to cell_key "1"

    # Mock slice
    expt.slice.species = "mouse"
    expt.slice.age = "P30"
    expt.slice.sex = "M"
    expt.slice.hemisphere = "left"

    # Mock data
    expt.data.__enter__.return_value = MagicMock()
    data_file = expt.data.__enter__.return_value

    # Mock sweep
    recording = MagicMock()
    recording.__getitem__.side_effect = lambda k: MockRecording(np.zeros(100), 'A') if k == 'primary' else MagicMock()
    recording.clamp_mode = 'vc'
    recording.holding_potential = -0.070
    recording.stimulus.description = 'test'
    recording.stimulus.items = []

    sweep = MagicMock()
    sweep.key = 0
    sweep.devices = [0]
    sweep.recordings = [recording]
    sweep.recording_dict = {0: MagicMock()}
    sweep.recording_dict[0].meta = {'notebook': {'Sampling interval': 0.1, 'Async AD 1: Bath Temperature': 34.0}}

    data_file.contents = [sweep]

    expt.pair_list = []

    return expt

def test_experiment_writer_init(tmp_path):
    writer = ExperimentWriter(output_dir=tmp_path)
    assert writer.output_dir == Path(tmp_path)
    assert writer.output_dir.exists()

def test_scale_recording():
    rec_v = {'primary': MockRecording(np.array([1.0], dtype='float32'), 'V')}
    scaled_v = ExperimentWriter.scale_recording(rec_v)
    assert scaled_v[0] == 1000.0

    rec_a = {'primary': MockRecording(np.array([1.0], dtype='float32'), 'A')}
    scaled_a = ExperimentWriter.scale_recording(rec_a)
    assert scaled_a[0] == 1e12

def test_index_sweeps():
    recording = MagicMock()
    recording.clamp_mode = 'vc'
    recording.holding_potential = -0.070
    recording.stimulus.description = 'test'

    sweep = MagicMock()
    sweep.key = 10
    sweep.devices = [0]
    sweep.recordings = [recording]

    data_file = MagicMock()
    data_file.contents = [sweep]

    index = ExperimentWriter.index_sweeps(data_file)
    assert index[0][-70.0] == [10]

def test_get_temperature():
    sweep = MagicMock()
    sweep.devices = [0]
    sweep.recording_dict = {0: MagicMock()}
    sweep.recording_dict[0].meta = {'notebook': {'Async AD 1: Bath Temperature': 35.5}}

    data_file = MagicMock()
    data_file.contents = [sweep]

    temp = ExperimentWriter.get_temperature(data_file)
    assert temp == 35.5
