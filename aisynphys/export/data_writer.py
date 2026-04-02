import os
import gc
import h5py
import numpy as np
from pathlib import Path
from .constants import INTRINSIC_FIELDS

class ExperimentWriter:
    CONDITIONS = {0: -55, 1: -70}
    BLANK_ADD_TIME = 0.02
    CHUNK_SIZE = 50_000
    COMPRESSION = "gzip"
    COMPRESSION_OPTS = 4
    SIGNAL_KEYS = ["signal_-55mV", "signal_-70mV"]

    def __init__(self, output_dir="output"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self._h5f = None
        self._datasets = {}
        self._metadata_store = {}

    def _get_or_create_dataset(self, out_file, cell_name, signal_key, dtype, electrode, expt):
        if self._h5f is None:
            self._h5f = h5py.File(out_file, "w")
            # Disable HDF5 chunk cache to prevent memory accumulation
            self._h5f.id.set_chunk_cache(0, 0, 0)

        key = (cell_name, signal_key)
        if key in self._datasets:
            return self._datasets[key]

        grp = self._h5f.require_group(cell_name)
        if "electrode" not in grp.attrs:
            grp.attrs["cell_id"] = cell_name
            grp.attrs["expt_ext_id"] = self.safe_attr(getattr(expt, "ext_id", None))
            grp.attrs["species"] = self.safe_attr(getattr(expt.slice, "species", None))
            grp.attrs["age"] = self.safe_attr(getattr(expt.slice, "age", None))
            grp.attrs["electrode"] = electrode
            grp.attrs["cell_ext_id"] = self.safe_attr(expt.electrodes[electrode].cell.ext_id)
            grp.attrs["holding_potentials"] = list(self.CONDITIONS.values())
            grp.attrs["sampling_rate_unit"] = "seconds"
            grp.attrs["signal_unit"] = "pA"

        subgrp = grp.require_group(signal_key)
        data_ds = subgrp.create_dataset(
            "data",
            shape=(0,),
            maxshape=(None,),
            dtype=dtype,
            chunks=(self.CHUNK_SIZE,),
            compression=self.COMPRESSION,
            compression_opts=self.COMPRESSION_OPTS,
        )

        mask_ds = subgrp.create_dataset(
            "mask",
            shape=(0,),
            maxshape=(None,),
            dtype=bool,
            chunks=(self.CHUNK_SIZE,),
            compression=self.COMPRESSION,
            compression_opts=self.COMPRESSION_OPTS,
        )

        self._datasets[key] = [data_ds, mask_ds, 0]
        self._metadata_store[key] = {
            "n_sweeps": 0,
            "sweep_lengths": [],
            "sampling_intervals": [],
            "n_samples": 0,
        }
        return self._datasets[key]

    @staticmethod
    def _append_data(ds_tuple, new_data, new_mask):
        data_ds, mask_ds, size = ds_tuple
        n = len(new_data)
        data_ds.resize(size + n, axis=0)
        mask_ds.resize(size + n, axis=0)
        data_ds[size:size + n] = new_data
        mask_ds[size:size + n] = new_mask
        ds_tuple[2] += n  # update size

    def write(self, expt):
        if expt.data is None:
            return None

        exp_name = f"exp_{expt.id:05d}"
        out_file = self.output_dir / f"{exp_name}.h5"

        self._h5f = None
        self._datasets = {}
        self._metadata_store = {}
        metadata_rows = []

        # Build cell_id to electrode mapping
        cell_to_electrode = {}
        for elec_id, elec in expt.electrodes.items():
            if elec.cell is not None:
                cell_to_electrode[elec.cell.id] = elec_id

        # Build postsynaptic cell to presynaptic electrodes mapping
        presynaptic_electrodes = {}
        for pair in expt.pair_list:
            if pair.synapse:
                post_cell_id = pair.post_cell_id
                pre_cell_id = pair.pre_cell_id
                pre_electrode = cell_to_electrode.get(pre_cell_id)
                if pre_electrode is not None:
                    if post_cell_id not in presynaptic_electrodes:
                        presynaptic_electrodes[post_cell_id] = []
                    presynaptic_electrodes[post_cell_id].append(pre_electrode)

        with expt.data as data_file:
            syn_index = self.build_synapse_index(expt)
            sweep_index = self.index_sweeps(data_file)

            # Pre-fetch temperature once safely
            temp = self.get_temperature(data_file)

            for electrode in expt.data.contents[0].devices:
                if expt.electrodes[electrode].cell is None:
                    continue

                cell_id = expt.electrodes[electrode].cell.id
                cell_name = f"cell_{cell_id:06d}"

                passive_cache = {"cap": [], "rin": [], "rseries": []}

                for cond_idx, v_hold in self.CONDITIONS.items():
                    cell_key = str(electrode + 1)

                    if cell_key not in expt.cells:
                        continue

                    sweep_indices = sweep_index.get(electrode, {}).get(v_hold, [])
                    if len(sweep_indices) == 0:
                        continue

                    signal_key = self.SIGNAL_KEYS[cond_idx]

                    for sweep_ix in sweep_indices:
                        sweep = data_file.contents[sweep_ix]

                        devices = sweep.devices
                        if electrode not in devices:
                            continue

                        device_map = {d: i for i, d in enumerate(devices)}
                        device_ix = device_map[electrode]
                        recording = sweep.recordings[device_ix]

                        # Ensure data is cast to float32 and scaled to standard units
                        sweep_data = self.scale_recording(recording)

                        meta = sweep.recording_dict[electrode].meta['notebook']
                        if 'Sampling interval' in meta:
                            sample_interval = np.round(meta['Sampling interval'] * 1e-3, 6)
                        else:
                            sample_interval = np.round(meta['Minimum Sampling interval'] * 1e-3, 6)

                        mask = np.zeros_like(sweep_data, dtype=bool)

                        # --- stimulus blanking (current electrode) ---
                        self._blank_recording(recording, mask, sample_interval, self.BLANK_ADD_TIME)

                        # --- stimulus blanking (neighboring electrodes) ---
                        for neighbor in (electrode - 1, electrode + 1):
                            ix = device_map.get(neighbor)
                            if ix is None:
                                continue
                            try:
                                self._blank_recording(sweep.recordings[ix], mask, sample_interval, 0)
                            except Exception:
                                continue

                        # --- stimulus blanking (presynaptic electrodes) ---
                        pre_elecs = presynaptic_electrodes.get(cell_id, [])
                        for pre_elec in pre_elecs:
                            ix = device_map.get(pre_elec)
                            if ix is None:
                                continue
                            try:
                                self._blank_recording(sweep.recordings[ix], mask, sample_interval, 0)
                            except Exception:
                                continue

                        # create dataset ONLY when first real data appears
                        ds_tuple = self._get_or_create_dataset(
                            out_file,
                            cell_name,
                            signal_key,
                            sweep_data.dtype,
                            electrode,
                            expt
                        )

                        # append immediately (no buffering)
                        self._append_data(ds_tuple, sweep_data, mask)

                        meta_dict = self._metadata_store[(cell_name, signal_key)]
                        meta_dict["n_sweeps"] += 1
                        meta_dict["sweep_lengths"].append(len(sweep_data))
                        meta_dict["sampling_intervals"].append(sample_interval)
                        meta_dict["n_samples"] += len(sweep_data)

                        tp = getattr(recording, "nearest_test_pulse", None)
                        if tp:
                            if tp.capacitance is not None:
                                passive_cache["cap"].append(tp.capacitance)
                            if tp.input_resistance is not None:
                                passive_cache["rin"].append(tp.input_resistance)
                            if tp.access_resistance is not None:
                                passive_cache["rseries"].append(tp.access_resistance)

                row = self.extract_electrode_metadata(expt, electrode, syn_index, temp)

                # overwrite passive with cached version
                row.update({
                        'capacitance': np.nanmedian(passive_cache["cap"]) * 1e12 if passive_cache["cap"] else np.nan,
                        'resistance': np.nanmedian(passive_cache["rin"]) * 1e-6 if passive_cache["rin"] else np.nan,
                        'access_resistance': np.nanmedian(passive_cache["rseries"]) * 1e-6 if passive_cache["rseries"] else np.nan,
                    })

                metadata_rows.append(row)

        # --- finalize ---
        if self._h5f is None:
            return {
                "expt_id": expt.id,
                "ext_id": getattr(expt, "ext_id", None),
                "n_cells": 0,
                "output_file": None,
                "status": "no_data",
            }

        self._finalize_metadata(metadata_rows)

        self._h5f.flush()
        self._h5f.close()
        self._h5f = None

        # Explicit garbage collection to prevent memory accumulation
        gc.collect()

        cells_with_data = {
          cell_name for (cell_name, _) in self._metadata_store.keys()
        }
        n_cells = len(cells_with_data)

        result = {
            "expt_id": expt.id,
            "ext_id": getattr(expt, "ext_id", None),
            "n_cells": n_cells,
            "output_file": str(out_file),
            "status": "OK"
        }

        return result

    def _finalize_metadata(self, metadata_rows):
        for (cell_name, signal_key), meta in self._metadata_store.items():
            grp = self._h5f[cell_name][signal_key]

            sweep_lengths = np.array(meta["sweep_lengths"], dtype=np.int32)
            sampling_intervals = np.array(meta["sampling_intervals"], dtype=np.float32)

            duration_s = float(np.sum(sweep_lengths * sampling_intervals))
            cond_idx = self.SIGNAL_KEYS.index(signal_key)

            grp.attrs["n_sweeps"] = meta["n_sweeps"]
            grp.attrs["sweep_lengths"] = sweep_lengths
            grp.attrs["sweep_starts"] = np.cumsum(sweep_lengths) - sweep_lengths
            grp.attrs["sampling_intervals"] = sampling_intervals
            grp.attrs["n_samples"] = meta["n_samples"]
            grp.attrs["duration_s"] = duration_s
            grp.attrs["voltage_mV"] = self.CONDITIONS.get(cond_idx, np.nan)

        if metadata_rows:
            keys = sorted(metadata_rows[0].keys())
            grp = self._h5f.require_group("metadata_table")

            for key in keys:
                col = [row.get(key) for row in metadata_rows]

                # --- detect if column is numeric ---
                is_numeric = True
                for v in col:
                    if v is None:
                        continue
                    if isinstance(v, (int, float, np.number)):
                        continue
                    is_numeric = False
                    break

                if is_numeric:
                    data = np.array(
                        [np.nan if v is None else v for v in col],
                        dtype=np.float32
                    )
                else:
                    data = np.array(
                        [str(v) if v is not None else "" for v in col],
                        dtype="S"
                    )
                grp.create_dataset(key, data=data)

    @staticmethod
    def get_temperature(data_file):
        """Safely extract temperature from the first sweep in the data file."""
        try:
            # Try to get temperature from the first sweep
            # Using the first available electrode
            first_sweep = data_file.contents[0]
            first_elec = first_sweep.devices[0]
            temp = first_sweep.recording_dict[first_elec].meta['notebook']['Async AD 1: Bath Temperature']
            return temp
        except Exception:
            return None

    @staticmethod
    def index_sweeps(data_file):
        """Pre-index sweeps in the data file for faster lookup.

        Returns a nested dictionary: {electrode_id: {holding_potential: [sweep_indices]}}
        """
        index = {}
        for sweep in data_file.contents:
            for electrode_id in sweep.devices:
                device_ix = sweep.devices.index(electrode_id)
                rec = sweep.recordings[device_ix]

                if rec.clamp_mode != 'vc':
                    continue

                stim = rec.stimulus.description
                if stim in ['MIES_Blowout_DA_0', 'Chirp_DA_0', 'Mixedf_DA_0']:
                    holding = 'NaN'
                else:
                    holding = np.round(rec.holding_potential * 1e3)

                if electrode_id not in index:
                    index[electrode_id] = {}
                if holding not in index[electrode_id]:
                    index[electrode_id][holding] = []
                index[electrode_id][holding].append(sweep.key)

        return index

    @staticmethod
    def get_blank_segments_from_cmd_trace(stim_data, extra_points):
        baseline = stim_data[0]
        sdiff = np.diff(stim_data)
        changes = np.argwhere(sdiff != 0)[:, 0] + 1
        pulses = []
        for i, start in enumerate(changes):
            if (stim_data[start] - baseline) > 0:
                stop = changes[i+1] + int(extra_points) if (i+1 < len(changes)) else len(stim_data)
                pulses.append((start, stop))

        return pulses

    @staticmethod
    def get_blank_segments(stim_item, sample_interval, extra_points):
        stim_start = stim_item.start_time
        pulses = []
        for pulse in range(stim_item.n_pulses):
            start = stim_start + pulse * stim_item.interval
            stop = stim_start + pulse * stim_item.interval + stim_item.pulse_duration

            pulses.append((int(start/sample_interval), int(stop/sample_interval) + int(extra_points)))

        return pulses

    @staticmethod
    def scale_recording(recording):
        """Scale recording data to standard units (pA or mV) and cast to float32."""
        scaling = 1e3 if recording['primary'].units == 'V' else 1e12
        return (recording['primary'].data * scaling).astype('float32')

    @staticmethod
    def safe_attr(value):
        return "" if value is None else value

    @staticmethod
    def build_synapse_index(expt):
        syn_index = {}

        for pair in expt.pair_list:
            if pair.synapse:
                syn_index[pair.post_cell_id] = pair

        return syn_index

    @staticmethod
    def get_synaptic_properties(syn_index, cell_id):
        results = {
            'psc_amplitude': None,
            'psc_rise_time': None,
            'psc_decay_tau': None,
            'psc_fit_amplitude': None,
            'psp_amplitude': None,
            'psp_rise_time': None,
            'psp_decay_tau': None,
            'psp_fit_amplitude': None,
        }

        pair = syn_index.get(cell_id)
        if pair is None:
            return results

        syn = pair.synapse

        vc_fit = None
        ic_fit = None

        for fit in getattr(syn, "avg_response_fits", []):
            mode = getattr(fit, "clamp_mode", None)
            if mode == 'vc':
                vc_fit = getattr(fit, "fit_amp", None)
            elif mode == 'ic':
                ic_fit = getattr(fit, "fit_amp", None)

        results.update({
            'psc_amplitude': getattr(syn, 'psc_amplitude', None),
            'psc_rise_time': getattr(syn, 'psc_rise_time', None),
            'psc_decay_tau': getattr(syn, 'psc_decay_tau', None),
            'psc_fit_amplitude': vc_fit,
            'psp_amplitude': getattr(syn, 'psp_amplitude', None),
            'psp_rise_time': getattr(syn, 'psp_rise_time', None),
            'psp_decay_tau': getattr(syn, 'psp_decay_tau', None),
            'psp_fit_amplitude': ic_fit,
        })

        return results

    def _blank_recording(self, recording, mask, sample_interval, blank_add_time):
        stimulus_info = False
        for item in recording.stimulus.items:
            if item.description == 'test pulse':
                start = int(item.start_time / sample_interval) - 2
                end = int((item.start_time + item.duration * 2) / sample_interval)
                mask[max(start, 0):min(end, len(mask))] = True

            elif item.type == 'SquarePulseTrain':
                stimulus_info = True
                segments = self.get_blank_segments(
                    item,
                    sample_interval,
                    blank_add_time / sample_interval
                )
                for start, stop in segments:
                    mask[max(start - 2, 0):min(stop + 5, len(mask))] = True

        if not stimulus_info:
            try:
                cmd = recording['command'].data
                segments = self.get_blank_segments_from_cmd_trace(
                    cmd,
                    self.BLANK_ADD_TIME / sample_interval
                )
                for start, stop in segments:
                    mask[max(start - 2, 0):min(stop + 5, len(mask))] = True
            except (KeyError, AttributeError, ValueError):
                pass

    def extract_electrode_metadata(self, expt, electrode, syn_index, temperature=None):
        cell = expt.electrodes[electrode].cell
        morph = getattr(cell, "morphology", None)
        loc = getattr(cell, "cortical_location", None)
        intrinsic = getattr(cell, "intrinsic", None)
        elec = cell.electrode

        row = {
            # --- Cell ---
            "cell_id": cell.id,
            "cell_ext_id": cell.ext_id,
            "cortical_layer": getattr(loc, "cortical_layer", None),
            "fractional_depth": getattr(loc, "fractional_depth", None),
            "fractional_layer_depth": getattr(loc, "fractional_layer_depth", None),

            "dendrite_type": getattr(morph, "dendrite_type", None),
            "pyramidal": getattr(morph, "pyramidal", None),
            "cell_class": getattr(morph, "cell_class", None),
            "cell_class_nonsynaptic": getattr(morph, "cell_class_nonsynaptic", None),
            "cre_type": getattr(morph, "cre_type", None),

            "apical_trunc_distance": getattr(morph, "apical_trunc_distance", None),
            "qual_morpho_type": getattr(morph, "qual_morpho_type", None),

            # --- Experiment ---
            "expt_id": expt.id,
            "expt_ext_id": expt.ext_id,
            "target_region": getattr(expt, "target_region", None),
            "date": expt.date,

            # --- Slice ---
            "species": getattr(expt.slice, "species", None),
            "age": getattr(expt.slice, "age", None),
            "sex": getattr(expt.slice, "sex", None),
            "hemisphere": getattr(expt.slice, "hemisphere", None),
            "temperature_c": temperature,

            # --- Electrode ---
            "device_id": elec.device_id,
        }

        # --- Intrinsic ---
        if intrinsic is not None:
            for field in INTRINSIC_FIELDS:
                row[f"intrinsic_{field}"] = getattr(intrinsic, field, None)
        else:
            for field in INTRINSIC_FIELDS:
                row[f"intrinsic_{field}"] = None

        # --- Synaptic ---
        row.update(self.get_synaptic_properties(syn_index, cell.id))

        return row