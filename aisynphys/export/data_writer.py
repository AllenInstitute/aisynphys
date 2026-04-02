import numpy as np
import h5py
import gc


class ExperimentWriter:

    def __init__(self, out_file, chunk_size=50_000):
        self.out_file = out_file
        self.chunk_size = chunk_size

        self.h5f = h5py.File(
            out_file,
            "w",
            rdcc_nbytes=2 * 1024**2,
            rdcc_nslots=500,
        )

        self.datasets = {}
        self.metadata_rows = []

    # -----------------------------
    # Dataset handling
    # -----------------------------
    def _get_dataset(self, cell_name, signal_key, dtype):

        key = (cell_name, signal_key)

        if key in self.datasets:
            return self.datasets[key]

        grp = self.h5f.require_group(f"cells/{cell_name}")
        sig_grp = grp.require_group(f"signals/{signal_key}")

        data_ds = sig_grp.create_dataset(
            "data",
            shape=(0,),
            maxshape=(None,),
            dtype=dtype,
            chunks=(self.chunk_size,),
            compression="gzip",
            compression_opts=4,
        )
        data_ds.id.set_chunk_cache(0, 0, 0)

        mask_ds = sig_grp.create_dataset(
            "mask",
            shape=(0,),
            maxshape=(None,),
            dtype=bool,
            chunks=(self.chunk_size,),
            compression="gzip",
            compression_opts=4,
        )
        mask_ds.id.set_chunk_cache(0, 0, 0)

        self.datasets[key] = [data_ds, mask_ds, 0]

        return self.datasets[key]

    def _append(self, ds_tuple, data, mask):
        data_ds, mask_ds, size = ds_tuple
        n = len(data)

        data_ds.resize(size + n, axis=0)
        mask_ds.resize(size + n, axis=0)

        data_ds[size:size+n] = data
        mask_ds[size:size+n] = mask

        ds_tuple[2] += n

    # -----------------------------
    # Main processing
    # -----------------------------
    def process_experiment(self, expt):

        syn_index = build_synapse_index(expt)

        with expt.data as data_file:

            devices = expt.data.contents[0].devices

            for electrode in devices:

                cell = expt.electrodes[electrode].cell
                if cell is None:
                    continue

                cell_name = f"cell_{cell.id:06d}"

                # --- passive cache ---
                passive_cache = {"cap": [], "rin": [], "rseries": []}

                for cond_idx, v_hold in CONDITIONS.items():

                    signal_key = SIGNAL_KEYS[cond_idx]

                    sweep_indices = get_vc_sweep_indices(data_file, electrode, v_hold)
                    if not sweep_indices:
                        continue

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
                        dt = np.round(
                            meta.get('Sampling interval',
                                     meta.get('Minimum Sampling interval')) * 1e-3,
                            6
                        )

                        # --- mask ---
                        mask = np.zeros_like(sweep_data, dtype=bool)

                        # (your masking logic here...)

                        sweep_data[mask] = np.nan

                        # --- append ---
                        ds = self._get_dataset(cell_name, signal_key, sweep_data.dtype)
                        self._append(ds, sweep_data, mask)

                        # --- passive ---
                        tp = getattr(recording, "nearest_test_pulse", None)
                        if tp:
                            if tp.capacitance is not None:
                                passive_cache["cap"].append(tp.capacitance)
                            if tp.input_resistance is not None:
                                passive_cache["rin"].append(tp.input_resistance)
                            if tp.access_resistance is not None:
                                passive_cache["rseries"].append(tp.access_resistance)

                        del sweep, recording, sweep_data, mask

                # --- metadata row ---
                row = extract_electrode_metadata(expt, electrode, syn_index)

                # overwrite passive with cached version
                row.update({
                    'capacitance': np.nanmedian(passive_cache["cap"]) * 1e12 if passive_cache["cap"] else np.nan,
                    'resistance': np.nanmedian(passive_cache["rin"]) * 1e-6 if passive_cache["rin"] else np.nan,
                    'access_resistance': np.nanmedian(passive_cache["rseries"]) * 1e-6 if passive_cache["rseries"] else np.nan,
                })

                self.metadata_rows.append(row)

                gc.collect()

        # --- write metadata table ---
        self._write_cell_table()

        self.h5f.flush()
        self.h5f.close()

        del self.h5f
        gc.collect()

    # -----------------------------
    # Metadata table
    # -----------------------------
    def _write_cell_table(self):

        if not self.metadata_rows:
            return

        keys = sorted(self.metadata_rows[0].keys())

        grp = self.h5f.require_group("cell_table")

        for key in keys:
            col = [row.get(key) for row in self.metadata_rows]

            # convert safely
            if isinstance(col[0], (int, float, np.number)) or col[0] is None:
                data = np.array([np.nan if v is None else v for v in col], dtype=np.float32)
            else:
                data = np.array([str(v) if v is not None else "" for v in col], dtype="S")

            grp.create_dataset(key, data=data)