import h5py
import numpy as np
from pathlib import Path



def _decode(arr):
    if arr.dtype.kind == "S":
        return np.array([x.decode() for x in arr])
    return arr


class ExperimentFile:

    def __init__(self, path):
        self.path = path
        self._meta_cache = None

    def load_metadata(self):
        if self._meta_cache is not None:
            return self._meta_cache

        with h5py.File(self.path, "r") as f:
            table = f["metadata_table"]

            meta = {
                k: _decode(table[k][:])
                for k in table.keys()
            }

        self._meta_cache = meta
        return meta

    def iter_cells(self):
        meta = self.load_metadata()
        n = len(next(iter(meta.values())))

        for i in range(n):
            yield CellView(self, i)


class CellView:

    def __init__(self, expt_file, index):
        self.expt_file = expt_file
        self.index = index

    def get(self, key):
        return self.expt_file._meta_cache[key][self.index]

    @property
    def cell_id(self):
        return self.get("cell_id")

    def get_signal(self, signal_key):
        with h5py.File(self.expt_file.path, "r") as f:
            ds = f[f"cell_{int(self.cell_id):06d}/{signal_key}/data"]
            return ds[:]  # lazy load here

    def get_mask(self, signal_key):
        with h5py.File(self.expt_file.path, "r") as f:
            return f[f"cell_{int(self.cell_id):06d}/{signal_key}/mask"][:]

    def get_time(self, signal_key):
        with h5py.File(self.expt_file.path, "r") as f:
            grp = f[f"cell_{int(self.cell_id):06d}/{signal_key}"]

            lengths = grp.attrs["sweep_lengths"]
            starts = grp.attrs["sweep_starts"]
            dt = grp.attrs["sampling_intervals"]

        return reconstruct_time_axis(starts, lengths, dt)


def reconstruct_time_axis(starts, lengths, sampling_intervals):
    time = np.empty(np.sum(lengths), dtype=np.float32)

    for i, (start, length, dt) in enumerate(zip(starts, lengths, sampling_intervals)):
        t = np.arange(length, dtype=np.float32) * dt
        time[start:start+length] = t

    return time


class Dataset:

    def __init__(self, folder):
        self.files = [ExperimentFile(p) for p in Path(folder).glob("*.h5")]

        self._global_index = None

    def _build_index(self):

        rows = []

        for f in self.files:
            meta = f.load_metadata()

            n = len(next(iter(meta.values())))

            for i in range(n):
                row = {k: meta[k][i] for k in meta}
                row["_file"] = f
                row["_index"] = i
                rows.append(row)

        self._global_index = rows

    def query(self, **filters):

        if self._global_index is None:
            self._build_index()

        results = []

        for row in self._global_index:

            match = True
            for k, v in filters.items():
                if row.get(k) != v:
                    match = False
                    break

            if match:
                results.append(
                    CellView(row["_file"], row["_index"])
                )

        return results