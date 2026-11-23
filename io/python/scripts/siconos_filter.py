"""
Description: Filter the contents of a Siconos mechanics-IO HDF5 simulation file.
"""

# Lighter imports before command line parsing
import argparse
import os
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "fns_in", metavar="input", type=str, nargs="+", help="input file(s) (HDF5)"
)
parser.add_argument(
    "fn_out", metavar="output", type=str, nargs=1, help="output file (HDF5)"
)
parser.add_argument(
    "--start",
    metavar="time",
    type=float,
    help="time in seconds to cut the start of the simulation",
)
parser.add_argument(
    "--end",
    metavar="time",
    type=float,
    help="time in seconds to cut the end of the simulation",
)
parser.add_argument(
    "--interval", metavar="time", type=float, help="time between frames to preserve"
)
parser.add_argument("--gzip", action="store_true", help="enable compression in copy")
parser.add_argument(
    "--single", action="store_true", help="use single-precision floats in copy"
)
parser.add_argument(
    "--exclude", type=str, help="specify objects to exclude from copy (comma-separated)"
)
parser.add_argument(
    "--attr",
    type=str,
    action="append",
    help="specify attributes to change during copy (obj.name=value)",
)
parser.add_argument("-V", "--version", action="version", version="@SICONOS_VERSION@")

if __name__ == "__main__":
    args = parser.parse_args()

# Heavier imports after command line parsing
import h5py
import numpy as np


class CopyVisitor:
    """The CopyVisitor is called for each group and dataset in the HDF5
    file, and is responsible for copying the structure to the new
    HDF5 file."""

    def __init__(self, time_filter=None, object_filter=None, attr_filter=None):
        self.time_filter = None
        if time_filter is not None:
            self.time_filter = np.vectorize(time_filter, otypes=[bool])
        self.object_filter = object_filter
        self.time_idx = None
        self.times = None
        self.excluded_objects = None
        # rows of the input kept in each table that has a step index
        self.kept = {}

        def copy_attrs(obj_to, obj_from):
            for a in obj_from.attrs:
                value = None
                if attr_filter:
                    value = attr_filter(obj_to, obj_from, a)
                if value is None:
                    value = obj_from.attrs[a]
                obj_to.attrs[a] = value

        self.copy_attrs = copy_attrs

    def kept_steps(self, index):
        """Steps (t, first_row, n_rows) of a step index whose time is kept."""
        steps = index[:]
        return steps[np.isin(steps[:, 0].astype(self.times.dtype), self.times)]

    def visitor(self, path, obj):
        gr = io_out
        names = path.split("/")

        # Fill list of excluded objects if it has not been done yet
        if self.excluded_objects is None and self.object_filter is not None:
            if self.object_filter is not None:
                inp = obj.file["data/input"]
                self.excluded_objects = [
                    x.attrs["id"]
                    for name, x in inp.items()
                    if not self.object_filter(name, x)
                ]
                print(f"Excluding object IDs {self.excluded_objects}")

        # If we are copying an excluded object, return early
        if self.excluded_objects is not None and "data/input/" in path:
            if "id" not in obj.attrs:
                id = obj.parent.attrs["id"]
            else:
                id = obj.attrs["id"]
            if id in self.excluded_objects:
                return

        # Get indexes of filtered times in data/dynamic
        if self.time_filter is not None and (
            self.time_idx is None or self.times is None
        ):
            dyn = obj.file["data/dynamic"]
            self.time_idx = self.time_filter(dyn[:, 0]).nonzero()[0]
            self.times = dyn[self.time_idx, 0]

        # Create parent groups
        if len(names) > 1:
            for i, name in enumerate(names[:-1]):
                if name in gr:
                    gr = gr[name]
                else:
                    gr = gr.create_group(name)
                    gr_in = io_in["/".join(names[: i + 1])]
                    self.copy_attrs(gr, gr_in)

        # Determine chunks argument
        if obj.__class__ == h5py.Dataset:
            chunks = None
            if (
                hasattr(obj, "chunks")
                and obj.chunks is not None
                and np.all(np.array(obj.chunks) <= np.array(obj.shape))
            ):
                chunks = obj.chunks

            shape = obj.shape
            time_idx = None

            # Tables with a step index (data/step_index/<table>)
            name = path[len("data/") :]
            index = None
            if path.startswith("data/") and "/" not in name:
                index = io_in.get("data/step_index/" + name)

            # Filter current indexes
            tables = ["data/cf", "data/dynamic", "data/velocities", "data/static"]
            if index is not None or path in tables:
                if self.time_idx is not None:
                    # Get indexes of corresponding times in current dataset
                    if index is not None:
                        # rows of the kept steps, also for the tables without time column
                        steps = self.kept_steps(index)[:, 1:].astype(np.int64)
                        time_idx = np.concatenate(
                            [np.empty(0, np.int64)]
                            + [np.arange(first, first + n) for first, n in steps]
                        )
                    elif path == "data/dynamic":
                        time_idx = self.time_idx
                    # Time-filter all but static objects
                    elif path != "data/static":
                        time_idx = np.isin(obj[:, 0], self.times).nonzero()[0]

                # Additionally remove any lines referencing excluded objects
                if self.excluded_objects is not None and path in tables:
                    # ds_id in column 1, in cf: ds1, ds2 the last 2 columns
                    ids = obj[:, -2:] if path == "data/cf" else obj[:, 1:2]
                    excluded = np.isin(ids, self.excluded_objects).any(axis=1)
                    ex_idx = excluded.nonzero()[0]
                    if time_idx is None:
                        time_idx = np.arange(obj.shape[0])
                    time_idx = np.setdiff1d(time_idx, ex_idx)

                if index is not None and time_idx is not None:
                    self.kept[name] = time_idx

            # Shape of filtered dataset
            if time_idx is not None:
                shape = (len(time_idx),) + tuple(shape[1:])

            # Create the dataset, supply compression and dtype
            # arguments, possibly overridden by command line arguments
            comp = (obj.compression is True or args.gzip) and chunks is not None
            if comp:
                chunks = (4000,) + tuple(obj.shape[1:])
            # --single: floats only, the step index (row numbers) stays f8
            dtype = obj.dtype
            if args.single and dtype.kind == "f" and "data/step_index/" not in path:
                dtype = np.dtype("f4")
            ds = gr.create_dataset(
                obj.name,
                dtype=dtype,
                shape=shape,
                maxshape=obj.maxshape,
                chunks=chunks,
                compression=[obj.compression, True][comp],
                compression_opts=[obj.compression_opts, 9][comp],
                shuffle=comp,
                fletcher32=obj.fletcher32,
            )
            self.copy_attrs(ds, obj)

            # Copy the filtered or unfiltered dataset
            if obj.shape == ():  # scalar datasets: /format, data/script/source
                ds[()] = obj[()]
            elif time_idx is not None:
                if len(time_idx) == 0:
                    pass
                elif len(time_idx) == 1:
                    ds[0, :] = obj[time_idx[0], :]
                else:
                    ds[: len(time_idx), :] = obj[time_idx, :]
            else:
                ds[:] = obj

        # Some groups might be empty but we copy them anyway for their attributes
        elif obj.__class__ == h5py.Group:
            if path not in io_out:
                gr = io_out.create_group(path)
                for a in obj.attrs:
                    gr.attrs[a] = obj.attrs[a]

        else:
            print(f'Unknown type "{path}": {obj.__class__!s}')

    def rebuild_step_index(self):
        """After the copy: data/step_index/<table> of the output for the rows kept in <table>."""
        for name, rows in self.kept.items():
            steps = self.kept_steps(io_in["data/step_index/" + name])
            first = steps[:, 1].astype(np.int64)
            end = first + steps[:, 2].astype(np.int64)
            n_rows = np.searchsorted(rows, end) - np.searchsorted(rows, first)
            index = io_out["data/step_index/" + name]
            index.resize(len(steps), 0)
            index[:] = np.column_stack([steps[:, 0], np.cumsum(n_rows) - n_rows, n_rows])


if __name__ == "__main__":
    if os.path.exists(args.fn_out[0]):
        print(f'Output file "{args.fn_out[0]}" already exists!')
        sys.exit(1)

    re_exclude = lambda _: False
    if args.exclude is not None:
        # re_exclude = re.compile(args.exclude).match
        ex = args.exclude.split(",")
        re_exclude = lambda x: x in ex

    class TimeFilter:
        def __init__(self):
            self.marker = None
            self.last = None

        def __call__(self, t):
            res = True
            if args.end is not None:
                res = res and t <= args.end
            if args.start is not None:
                res = res and t >= args.start

            # Interval filter is true only for every arg.interval time index
            if res is True and self.marker is None:
                self.marker = t
            if args.interval is not None and self.marker is not None:
                if self.last is None:
                    if t >= self.marker:
                        self.last = t
                        self.marker += args.interval
                    else:
                        res = False
                else:
                    if t > self.last:
                        res = False
                        self.last = None
            return res

    class AttrFilter:
        """Return a value to assign to a given attribute,
        or None to copy the original value"""

        def __init__(self, attr_values):
            def pair(av):
                a, v = av.split("=")
                try:
                    v = int(v)
                except ValueError:
                    try:
                        v = float(v)
                    except ValueError:
                        pass
                return a.split("."), v

            self.attrs = [pair(v) for v in attr_values]

        def __call__(self, obj_to, obj_from, attr):
            if "nslaw" not in obj_from.name:
                return None
            for a, v in self.attrs:
                if obj_from.name == "/data/" + "/".join(a[:-1]) and attr == a[-1]:
                    return v

    with h5py.File(args.fns_in[0], mode="r") as io_in:
        # a file without version attributes is a v1.0 file
        version = (io_in.attrs.get("version_major", 1), io_in.attrs.get("version_minor", 0))
        if version[0] != 1:
            print(f'"{args.fns_in[0]}": format version {version[0]}.{version[1]}, only v1 is filtered')
            sys.exit(1)
        with h5py.File(args.fn_out[0], mode="w") as io_out:
            # attributes of the file: version_major, version_minor, dimension
            for a in io_in.attrs:
                io_out.attrs[a] = io_in.attrs[a]
            visitor = CopyVisitor(
                time_filter=TimeFilter(),
                object_filter=lambda name, obj: not re_exclude(name),
                attr_filter=args.attr and AttrFilter(args.attr),
            )
            io_in.visititems(visitor.visitor)
            visitor.rebuild_step_index()
