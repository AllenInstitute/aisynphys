import os
import gc
import time
import h5py
import psutil
import shutil
import numpy as np
import multiprocessing as mp
from pathlib import Path



INTRINSIC_FIELDS = [
    "rheobase",
    "fi_slope",
    "input_resistance",
    "input_resistance_ss",
    "sag",
    "tau",
    "sag_peak_t",
    "sag_depol",
    "sag_peak_t_depol",
    "ap_upstroke_downstroke_ratio",
    "ap_width",
    "ap_upstroke",
    "ap_downstroke",
    "ap_threshold_v",
    "ap_peak_deltav",
    "ap_fast_trough_deltav",
    "firing_rate_rheo",
    "latency_rheo",
    "firing_rate_40pa",
    "latency_40pa",
    "adaptation_index",
    "isi_cv",
    "chirp_peak_freq",
    "chirp_3db_freq",
    "chirp_peak_ratio",
    "chirp_peak_impedance",
    "chirp_sync_freq",
    "chirp_inductive_phase",
    "isi_adapt_ratio",
    "upstroke_adapt_ratio",
    "downstroke_adapt_ratio",
    "width_adapt_ratio",
    "threshold_v_adapt_ratio",
]


def get_vc_sweep_indices(data_file, electrode_id, v_hold=None):
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


def get_blank_segments(stim_item, sample_interval, extra_points):
  stim_start = stim_item.start_time
  pulses = []
  for pulse in range(stim_item.n_pulses):
    start = stim_start + pulse * stim_item.interval
    stop = stim_start + pulse * stim_item.interval + stim_item.pulse_duration

    pulses.append((int(start/sample_interval), int(stop/sample_interval) + int(extra_points)))

  return pulses


def safe_attr(value):
  return "" if value is None else value


def build_synapse_index(expt):
    syn_index = {}

    for pair in expt.pair_list:
        if pair.synapse:
            syn_index[pair.post_cell_id] = pair

    return syn_index


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


def extract_electrode_metadata(expt, electrode, syn_index):
    cell = expt.electrodes[electrode].cell
    morph = cell.morphology
    loc = cell.cortical_location
    intrinsic = getattr(cell, "intrinsic", None)
    elec = cell.electrode

    # Pre-fetch temperature once safely
    try:
        temp = expt.data.contents[0].recording_dict[electrode].meta['notebook']['Async AD 1: Bath Temperature']
    except Exception:
        temp = None

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
        "temperatature_c": temp,

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
    row.update(get_synaptic_properties(syn_index, cell.id))

    return row



    OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)

CONDITIONS = {0: -55, 1: -70}
BLANK_ADD_TIME = 0.02
CHUNK_SIZE = 50_000
COMPRESSION = "gzip"
COMPRESSION_OPTS = 4

SIGNAL_KEYS = ["signal_-55mV", "signal_-70mV"]


def save_data_from_expt(expt):
    if expt.data is None:
        return None

    exp_name = f"exp_{expt.id:05d}"
    out_file = OUTPUT_DIR / f"{exp_name}.h5"

    h5f = None  # lazy open
    datasets = {}  # (cell, signal) -> (data_ds, mask_ds, current_size)
    metadata_store = {}  # (cell_name, signal_key) -> metadata dict
    metadata_rows = []

    def get_or_create_dataset(cell_name, signal_key, dtype, electrode):
        nonlocal h5f

        if h5f is None:
            h5f = h5py.File(out_file, "w")

        key = (cell_name, signal_key)

        if key in datasets:
            return datasets[key]

        grp = h5f.require_group(cell_name)
        if "electrode" not in grp.attrs:
            grp.attrs["cell_id"] = cell_name
            grp.attrs["expt_ext_id"] = safe_attr(getattr(expt, "ext_id", None))
            grp.attrs["species"] = safe_attr(getattr(expt.slice, "species", None))
            grp.attrs["age"] = safe_attr(getattr(expt.slice, "age", None))
            grp.attrs["electrode"] = electrode
            grp.attrs["cell_ext_id"] = safe_attr(expt.electrodes[electrode].cell.ext_id)
            grp.attrs["holding_potentials"] = list(CONDITIONS.values())
            grp.attrs["sampling_rate_unit"] = "seconds"
            grp.attrs["signal_unit"] = "pA"

        subgrp = grp.require_group(signal_key)

        data_ds = subgrp.create_dataset(
            "data",
            shape=(0,),
            maxshape=(None,),
            dtype=dtype,
            chunks=(CHUNK_SIZE,),
            compression=COMPRESSION,
            compression_opts=COMPRESSION_OPTS,
        )

        mask_ds = subgrp.create_dataset(
            "mask",
            shape=(0,),
            maxshape=(None,),
            dtype=bool,
            chunks=(CHUNK_SIZE,),
            compression=COMPRESSION,
            compression_opts=COMPRESSION_OPTS,
        )

        datasets[key] = [data_ds, mask_ds, 0]
        metadata_store[key] = {
            "n_sweeps": 0,
            "sweep_lengths": [],
            "sampling_intervals": [],
            "n_samples": 0,
        }
        return datasets[key]

    def append_data(ds_tuple, new_data, new_mask):
        data_ds, mask_ds, size = ds_tuple
        n = len(new_data)

        data_ds.resize(size + n, axis=0)
        mask_ds.resize(size + n, axis=0)

        data_ds[size:size + n] = new_data
        mask_ds[size:size + n] = new_mask

        ds_tuple[2] += n  # update size

    with expt.data as data_file:
        syn_index = build_synapse_index(expt)
        for electrode in expt.data.contents[0].devices:
            if expt.electrodes[electrode].cell is None:
                continue

            cell_id = expt.electrodes[electrode].cell.id
            cell_name = f"cell_{cell_id:06d}"
            print(f"\tProcessing {cell_name}")

            passive_cache = {"cap": [], "rin": [], "rseries": []}

            for cond_idx, v_hold in CONDITIONS.items():
                cell_key = str(electrode + 1)

                if cell_key not in expt.cells:
                    continue

                sweep_indices = get_vc_sweep_indices(data_file, electrode, v_hold)
                if len(sweep_indices) == 0:
                    continue

                signal_key = SIGNAL_KEYS[cond_idx]

                for sweep_ix in sweep_indices:
                    sweep = data_file.contents[sweep_ix]

                    devices = sweep.devices
                    if electrode not in devices:
                        continue

                    device_map = {d: i for i, d in enumerate(devices)}
                    device_ix = device_map[electrode]
                    recording = sweep.recordings[device_ix]

                    scaling = 1e3 if recording['primary'].units == 'V' else 1e12
                    sweep_data = recording['primary'].data * scaling

                    meta = sweep.recording_dict[electrode].meta['notebook']
                    if 'Sampling interval' in meta:
                        sample_interval = np.round(meta['Sampling interval'] * 1e-3, 6)
                    else:
                        sample_interval = np.round(meta['Minimum Sampling interval'] * 1e-3, 6)

                    mask = np.zeros_like(sweep_data, dtype=bool)
                    stimulus_info = False

                    # --- stimulus blanking ---
                    for item in recording.stimulus.items:
                        if item.description == 'test pulse':
                            start = int(item.start_time / sample_interval) - 2
                            end = int((item.start_time + item.duration * 2) / sample_interval)
                            mask[max(start, 0):min(end, len(mask))] = True

                        elif item.type == 'SquarePulseTrain':
                            stimulus_info = True
                            segments = get_blank_segments(
                                item,
                                sample_interval,
                                BLANK_ADD_TIME / sample_interval
                            )
                            for start, stop in segments:
                                mask[max(start - 2, 0):min(stop + 5, len(mask))] = True

                    if not stimulus_info:
                        cmd = sweep.recordings[device_ix]['command'].data
                        segments = get_blank_segments_from_cmd_trace(
                            cmd,
                            BLANK_ADD_TIME / sample_interval
                        )
                        for start, stop in segments:
                            mask[max(start - 2, 0):min(stop + 5, len(mask))] = True

                    # create dataset ONLY when first real data appears
                    ds_tuple = get_or_create_dataset(
                        cell_name,
                        signal_key,
                        sweep_data.dtype,
                        electrode
                    )

                    # append immediately (no buffering)
                    append_data(ds_tuple, sweep_data, mask)

                    meta = metadata_store[(cell_name, signal_key)]
                    meta["n_sweeps"] += 1
                    meta["sweep_lengths"].append(len(sweep_data))
                    meta["sampling_intervals"].append(sample_interval)
                    meta["n_samples"] += len(sweep_data)

                    tp = getattr(recording, "nearest_test_pulse", None)
                    if tp:
                        if tp.capacitance is not None:
                            passive_cache["cap"].append(tp.capacitance)
                        if tp.input_resistance is not None:
                            passive_cache["rin"].append(tp.input_resistance)
                        if tp.access_resistance is not None:
                            passive_cache["rseries"].append(tp.access_resistance)

            row = extract_electrode_metadata(expt, electrode, syn_index)

            # overwrite passive with cached version
            row.update({
                    'capacitance': np.nanmedian(passive_cache["cap"]) * 1e12 if passive_cache["cap"] else np.nan,
                    'resistance': np.nanmedian(passive_cache["rin"]) * 1e-6 if passive_cache["rin"] else np.nan,
                    'access_resistance': np.nanmedian(passive_cache["rseries"]) * 1e-6 if passive_cache["rseries"] else np.nan,
                })

            metadata_rows.append(row)

    # --- finalize ---
    if h5f is None:
        print(f"No valid data for {exp_name}, skipping file.")
        return {
            "expt_id": expt.id,
            "ext_id": getattr(expt, "ext_id", None),
            "n_cells": 0,
            "output_file": None,
            "status": "no_data",
        }

    for (cell_name, signal_key), meta in metadata_store.items():
        grp = h5f[cell_name][signal_key]

        sweep_lengths = np.array(meta["sweep_lengths"], dtype=np.int32)
        sampling_intervals = np.array(meta["sampling_intervals"], dtype=np.float32)

        duration_s = float(np.sum(sweep_lengths * sampling_intervals))

        cond_idx = SIGNAL_KEYS.index(signal_key)

        grp.attrs["n_sweeps"] = meta["n_sweeps"]
        grp.attrs["sweep_lengths"] = sweep_lengths
        grp.attrs["sweep_starts"] = np.cumsum(sweep_lengths) - sweep_lengths
        grp.attrs["sampling_intervals"] = sampling_intervals
        grp.attrs["n_samples"] = meta["n_samples"]
        grp.attrs["duration_s"] = duration_s
        grp.attrs["voltage_mV"] = CONDITIONS.get(cond_idx, np.nan)


    if metadata_rows:
      keys = sorted(metadata_rows[0].keys())

      grp = h5f.require_group("metadata_table")

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

    h5f.flush()
    h5f.close()
    raw_file_path = data_file.filename
    gc.collect()
    data_file.close()

    cells_with_data = {
      cell_name for (cell_name, _) in metadata_store.keys()
    }
    n_cells = len(cells_with_data)

    result = {
        "expt_id": expt.id,
        "ext_id": getattr(expt, "ext_id", None),
        "n_cells": n_cells,
        "output_file": str(out_file),
        "status": "OK"
    }

    try:
      del data_file, h5f, sweep_data, mask, sweep, recording, datasets, metadata_store, sweep_indices, devices, meta, cmd
    except:
      pass

    if os.path.exists(raw_file_path):
      os.remove(raw_file_path)
      # remove folder
      folder = os.path.dirname(raw_file_path)
      shutil.rmtree(folder)

    gc.collect()
    print(f"Saved: {out_file}")

    return result



if __name__ == "__main__":
    from aisynphys.database import SynphysDatabase
    from datetime import timedelta


    def print_mem():
        process = psutil.Process(os.getpid())
        print(f"RAM: {process.memory_info().rss / 1e9:.2f} GB")


    def worker(expt, queue):
        try:
            result = save_data_from_expt(expt)
        except Exception as e:
            result = {
                "expt_id": expt.id,
                "ext_id": getattr(expt, "ext_id", None),
                "n_cells": 0,
                "output_file": None,
                "status": f"error: {e}",
            }

        queue.put(result)

    db = SynphysDatabase.load_current('small')
    experiments = db.list_experiments()
    expts_with_data = [e for e in experiments if e.ephys_file]

    cells = sum(
        1
        for e in expts_with_data
        for el in db.experiment_from_ext_id(e.ext_id).electrodes
        if el.stop_time
        and (el.stop_time - el.start_time) > timedelta(minutes=2)
    )

    print(f"Number of experiments: {len(experiments)}")
    print(f"Number of experiments with data: {len(expts_with_data)}")
    print(f"Number of cells with data: {cells}")


    log_file = OUTPUT_DIR / "processing_log.txt"
    results = []


    with open(log_file, "a") as f:

        np.random.seed(5632)
        random_indices = np.random.choice(len(expts_with_data), size=10, replace=False)

        for idx in random_indices:

            expt = db.experiment_from_ext_id(expts_with_data[idx].ext_id)

            print(f"Experiment: {expt.id} ({expt.ext_id})")

            queue = mp.Queue()

            start = time.time()

            p = mp.Process(target=worker, args=(expt, queue))
            p.start()
            p.join()

            end = time.time()

            result = queue.get()
            result["time_s"] = end - start

            results.append(result)

            log_line = (
                f"expt_id={result['expt_id']}, "
                f"ext_id={result['ext_id']}, "
                f"cells={result['n_cells']}, "
                f"time={result['time_s']:.2f}s, "
                f"status={result['status']}"
            )

            print(log_line)
            f.write(log_line + "\n")

            print_mem()
            print("")

    total_expts = len(results)
    successful = [r for r in results if r["status"] == "OK"]

    total_cells = sum(r["n_cells"] for r in successful)

    summary = f"""
    SUMMARY
    -------
    Experiments processed: {total_expts}
    Successful: {len(successful)}
    Failed: {total_expts - len(successful)}
    Total cells with data: {total_cells}
    Avg cells per experiment: {total_cells / max(len(successful),1):.2f}
    """

    print(summary)