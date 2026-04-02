import os
import shutil
import gc
import numpy as np
import h5py
from pathlib import Path

class ExperimentWriter:
    INTRINSIC_FIELDS = [
        "rheobase", "fi_slope", "input_resistance", "input_resistance_ss",
        "sag", "tau", "sag_peak_t", "sag_depol", "sag_peak_t_depol",
        "ap_upstroke_downstroke_ratio", "ap_width", "ap_upstroke",
        "ap_downstroke", "ap_threshold_v", "ap_peak_deltav",
        "ap_fast_trough_deltav", "firing_rate_rheo", "latency_rheo",
        "firing_rate_40pa", "latency_40pa", "adaptation_index", "isi_cv",
        "chirp_peak_freq", "chirp_3db_freq", "chirp_peak_ratio",
        "chirp_peak_impedance", "chirp_sync_freq", "chirp_inductive_phase",
        "isi_adapt_ratio", "upstroke_adapt_ratio", "downstroke_adapt_ratio",
        "width_adapt_ratio", "threshold_v_adapt_ratio",
    ]

    CONDITIONS = {0: -55, 1: -70}
    SIGNAL_KEYS = ["signal_-55mV", "signal_-70mV"]
    BLANK_ADD_TIME = 0.02
    COMPRESSION = "gzip"
    COMPRESSION_OPTS = 4

    def __init__(self, out_file, chunk_size=50_000):
        self.out_file = Path(out_file)
        self.out_file.parent.mkdir(exist_ok=True, parents=True)
        self.chunk_size = chunk_size

        self.h5f = h5py.File(
            str(self.out_file),
            "w",
            rdcc_nbytes=2 * 1024**2,
            rdcc_nslots=500,
        )

        self.datasets = {} # (cell_name, signal_key) -> [data_ds, mask_ds, current_size]
        self.metadata_store = {} # (cell_name, signal_key) -> metadata dict
        self.metadata_rows = []

    def _get_dataset(self, cell_name, signal_key, dtype, electrode, expt):
        key = (cell_name, signal_key)
        if key in self.datasets:
            return self.datasets[key]

        grp = self.h5f.require_group(f"cells/{cell_name}")
        if "electrode" not in grp.attrs:
            grp.attrs["cell_id"] = cell_name
            grp.attrs["expt_ext_id"] = self._safe_attr(getattr(expt, "ext_id", None))
            grp.attrs["species"] = self._safe_attr(getattr(expt.slice, "species", None))
            grp.attrs["age"] = self._safe_attr(getattr(expt.slice, "age", None))
            grp.attrs["electrode"] = electrode
            grp.attrs["cell_ext_id"] = self._safe_attr(expt.electrodes[electrode].cell.ext_id)
            grp.attrs["holding_potentials"] = list(self.CONDITIONS.values())
            grp.attrs["sampling_rate_unit"] = "seconds"
            grp.attrs["signal_unit"] = "pA"

        sig_grp = grp.require_group(f"signals/{signal_key}")

        data_ds = sig_grp.create_dataset(
            "data",
            shape=(0,),
            maxshape=(None,),
            dtype=dtype,
            chunks=(self.chunk_size,),
            compression=self.COMPRESSION,
            compression_opts=self.COMPRESSION_OPTS,
        )
        # Use low-level interface to disable chunk cache for these datasets if needed,
        # but the suggestion had it.
        data_ds.id.set_chunk_cache(0, 0, 0)

        mask_ds = sig_grp.create_dataset(
            "mask",
            shape=(0,),
            maxshape=(None,),
            dtype=bool,
            chunks=(self.chunk_size,),
            compression=self.COMPRESSION,
            compression_opts=self.COMPRESSION_OPTS,
        )
        mask_ds.id.set_chunk_cache(0, 0, 0)

        self.datasets[key] = [data_ds, mask_ds, 0]
        self.metadata_store[key] = {
            "n_sweeps": 0,
            "sweep_lengths": [],
            "sampling_intervals": [],
            "n_samples": 0,
        }
        return self.datasets[key]

    def _append(self, ds_tuple, data, mask):
        data_ds, mask_ds, size = ds_tuple
        n = len(data)

        data_ds.resize(size + n, axis=0)
        mask_ds.resize(size + n, axis=0)

        data_ds[size:size + n] = data
        mask_ds[size:size + n] = mask

        ds_tuple[2] += n

    def _safe_attr(self, value):
        return "" if value is None else value

    def _get_vc_sweep_indices(self, data_file, electrode_id, v_hold=None):
        sweep_indices = []
        for sweep in data_file.contents:
            devices = sweep.devices
            if electrode_id not in devices:
                continue

            device_ix = devices.index(electrode_id)
            rec = sweep.recordings[device_ix]
            if rec.clamp_mode != 'vc':
                continue

            stim = rec.stimulus.description
            if stim in ['MIES_Blowout_DA_0', 'Chirp_DA_0', 'Mixedf_DA_0']:
                holding = 'NaN'
            else:
                holding = np.round(rec.holding_potential * 1e3)

            if holding == v_hold:
                sweep_indices.append(sweep.key)
        return sweep_indices

    def _get_blank_segments_from_cmd_trace(self, stim_data, extra_points):
        baseline = stim_data[0]
        sdiff = np.diff(stim_data)
        changes = np.argwhere(sdiff != 0)[:, 0] + 1
        pulses = []
        for i, start in enumerate(changes):
            if (stim_data[start] - baseline) > 0:
                stop = changes[i+1] + int(extra_points) if (i+1 < len(changes)) else len(stim_data)
                pulses.append((start, stop))
        return pulses

    def _get_blank_segments(self, stim_item, sample_interval, extra_points):
        stim_start = stim_item.start_time
        pulses = []
        for pulse in range(stim_item.n_pulses):
            start = stim_start + pulse * stim_item.interval
            stop = stim_start + pulse * stim_item.interval + stim_item.pulse_duration
            pulses.append((int(start/sample_interval), int(stop/sample_interval) + int(extra_points)))
        return pulses

    def _build_synapse_index(self, expt):
        syn_index = {}
        for pair in expt.pair_list:
            if pair.synapse:
                syn_index[pair.post_cell_id] = pair
        return syn_index

    def _get_synaptic_properties(self, syn_index, cell_id):
        results = {
            'psc_amplitude': None, 'psc_rise_time': None, 'psc_decay_tau': None,
            'psc_fit_amplitude': None, 'psp_amplitude': None, 'psp_rise_time': None,
            'psp_decay_tau': None, 'psp_fit_amplitude': None,
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

    def _extract_electrode_metadata(self, expt, electrode, syn_index):
        cell = expt.electrodes[electrode].cell
        morph = cell.morphology
        loc = cell.cortical_location
        intrinsic = getattr(cell, "intrinsic", None)
        elec = cell.electrode

        try:
            temp = expt.data.contents[0].recording_dict[electrode].meta['notebook']['Async AD 1: Bath Temperature']
        except Exception:
            temp = None

        row = {
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
            "expt_id": expt.id,
            "expt_ext_id": expt.ext_id,
            "target_region": getattr(expt, "target_region", None),
            "date": expt.date,
            "species": getattr(expt.slice, "species", None),
            "age": getattr(expt.slice, "age", None),
            "sex": getattr(expt.slice, "sex", None),
            "hemisphere": getattr(expt.slice, "hemisphere", None),
            "temperatature_c": temp,
            "device_id": elec.device_id,
        }

        for field in self.INTRINSIC_FIELDS:
            row[f"intrinsic_{field}"] = getattr(intrinsic, field, None) if intrinsic is not None else None

        row.update(self._get_synaptic_properties(syn_index, cell.id))
        return row

    def process_experiment(self, expt):
        if expt.data is None:
            return {"expt_id": expt.id, "status": "no_data"}

        syn_index = self._build_synapse_index(expt)
        with expt.data as data_file:
            for electrode in data_file.contents[0].devices:
                if expt.electrodes[electrode].cell is None:
                    continue

                cell_id = expt.electrodes[electrode].cell.id
                cell_name = f"cell_{cell_id:06d}"
                print(f"\tProcessing {cell_name}")

                passive_cache = {"cap": [], "rin": [], "rseries": []}

                for cond_idx, v_hold in self.CONDITIONS.items():
                    cell_key = str(electrode + 1)
                    if cell_key not in expt.cells:
                        continue

                    sweep_indices = self._get_vc_sweep_indices(data_file, electrode, v_hold)
                    if not sweep_indices:
                        continue

                    signal_key = self.SIGNAL_KEYS[cond_idx]

                    for sweep_ix in sweep_indices:
                        sweep = data_file.contents[sweep_ix]
                        if electrode not in sweep.devices:
                            continue

                        device_ix = sweep.devices.index(electrode)
                        recording = sweep.recordings[device_ix]

                        # --- data ---
                        sweep_data = recording['primary'].data.astype(np.float32, copy=False)
                        scaling = 1e3 if recording['primary'].units == 'V' else 1e12
                        sweep_data *= scaling

                        # --- sampling ---
                        meta = sweep.recording_dict[electrode].meta['notebook']
                        sample_interval = np.round(
                            meta.get('Sampling interval',
                                     meta.get('Minimum Sampling interval', 0)) * 1e-3,
                            6
                        )

                        # --- mask ---
                        mask = np.zeros_like(sweep_data, dtype=bool)
                        stimulus_info = False

                        for item in recording.stimulus.items:
                            if item.description == 'test pulse':
                                start = int(item.start_time / sample_interval) - 2
                                end = int((item.start_time + item.duration * 2) / sample_interval)
                                mask[max(start, 0):min(end, len(mask))] = True
                            elif item.type == 'SquarePulseTrain':
                                stimulus_info = True
                                segments = self._get_blank_segments(
                                    item, sample_interval, self.BLANK_ADD_TIME / sample_interval
                                )
                                for start, stop in segments:
                                    mask[max(start - 2, 0):min(stop + 5, len(mask))] = True

                        if not stimulus_info:
                            cmd = recording['command'].data
                            segments = self._get_blank_segments_from_cmd_trace(
                                cmd, self.BLANK_ADD_TIME / sample_interval
                            )
                            for start, stop in segments:
                                mask[max(start - 2, 0):min(stop + 5, len(mask))] = True

                        sweep_data[mask] = np.nan

                        # --- append ---
                        ds_tuple = self._get_dataset(cell_name, signal_key, sweep_data.dtype, electrode, expt)
                        self._append(ds_tuple, sweep_data, mask)

                        # --- metadata store ---
                        mstore = self.metadata_store[(cell_name, signal_key)]
                        mstore["n_sweeps"] += 1
                        mstore["sweep_lengths"].append(len(sweep_data))
                        mstore["sampling_intervals"].append(sample_interval)
                        mstore["n_samples"] += len(sweep_data)

                        # --- passive ---
                        tp = getattr(recording, "nearest_test_pulse", None)
                        if tp:
                            if tp.capacitance is not None: passive_cache["cap"].append(tp.capacitance)
                            if tp.input_resistance is not None: passive_cache["rin"].append(tp.input_resistance)
                            if tp.access_resistance is not None: passive_cache["rseries"].append(tp.access_resistance)

                        del sweep, recording, sweep_data, mask

                # --- metadata row ---
                row = self._extract_electrode_metadata(expt, electrode, syn_index)
                row.update({
                    'capacitance': np.nanmedian(passive_cache["cap"]) * 1e12 if passive_cache["cap"] else np.nan,
                    'resistance': np.nanmedian(passive_cache["rin"]) * 1e-6 if passive_cache["rin"] else np.nan,
                    'access_resistance': np.nanmedian(passive_cache["rseries"]) * 1e-6 if passive_cache["rseries"] else np.nan,
                })
                self.metadata_rows.append(row)
                gc.collect()

        self._finalize_attributes()
        self._write_cell_table()

        raw_file_path = data_file.filename
        self.h5f.flush()
        self.h5f.close()
        gc.collect()

        if os.path.exists(raw_file_path):
            os.remove(raw_file_path)
            shutil.rmtree(os.path.dirname(raw_file_path), ignore_errors=True)

        cells_with_data = {key[0] for key in self.metadata_store.keys()}
        return {
            "expt_id": expt.id,
            "ext_id": getattr(expt, "ext_id", None),
            "n_cells": len(cells_with_data),
            "output_file": str(self.out_file),
            "status": "OK"
        }

    def _finalize_attributes(self):
        for (cell_name, signal_key), meta in self.metadata_store.items():
            grp = self.h5f[f"cells/{cell_name}/signals/{signal_key}"]
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

    def _write_cell_table(self):
        if not self.metadata_rows:
            return

        keys = sorted(self.metadata_rows[0].keys())
        grp = self.h5f.require_group("cell_table")

        for key in keys:
            col = [row.get(key) for row in self.metadata_rows]
            is_numeric = all(v is None or isinstance(v, (int, float, np.number)) for v in col)

            if is_numeric:
                data = np.array([np.nan if v is None else v for v in col], dtype=np.float32)
            else:
                data = np.array([str(v) if v is not None else "" for v in col], dtype="S")
            grp.create_dataset(key, data=data)
