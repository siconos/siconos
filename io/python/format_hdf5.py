#!/usr/bin/env @Python_EXECUTABLE@
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2025 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import h5py
import numpy as np
import json
import os
import re

vlen_dtype = h5py.vlen_dtype(str) # variable-length string, used to create "vlen" datasets

def check_datatype(data, format, symbols={}):
    # Examples
    # "dimension": {"enum": [2, 3], "default": 3}
    # "stops": {"dtype": "float", "shape": [null, 3], "opt": true, "note": "[axis, position, direction]"},
    # symbols: value of each symbol of the shapes ("dim": 3, ...); None (null) in a shape = free size
    data = np.asarray(data)
    if "enum" in format:
        enum = format["enum"]
        if data not in enum:
            raise ValueError(f"{data} not in {enum}")
    if "dtype" in format:
        dtypes = format["dtype"] if isinstance(format["dtype"], list) else [format["dtype"]]
        for dtype in dtypes: # first dtype of the list that accepts data
            if dtype == "vlen" or dtype == "str":
                if data.dtype.kind in "USO": # str, bytes (h5py 3 read) or object
                    break
            elif np.can_cast(data.dtype, dtype):
                break
        else:
            raise ValueError(f"{data} can not be casted to the expected type : {data.dtype} -- expected : {dtypes}")
    if "shape" in format:
        shapes = format["shape"] if isinstance(format["shape"][0], list) else [format["shape"]]
        for shape in shapes: # first shape of the list that matches data
            if data.ndim == len(shape) and all(s_format is None or s_data == symbols.get(s_format, s_format)
                                               for s_data, s_format in zip(data.shape, shape)):
                break
        else:
            raise ValueError(f"{data} does not have the expected shapes : {data.shape} -- expected : {shapes}")

def info(layout):
    """attrs["info"] of a table, from the layout of its columns ('x[k]' = k columns, 'x[*]' = the rest)."""
    lines, col = [], 0
    for column in layout:
        name, _, size = column.partition("[")
        size = size.rstrip("]")
        if size == "*":
            lines.append(f"[{col}:] : {name}")
            break
        n = int(size) if size else 1
        lines.append(f"[{col}] : {name}" if n == 1 else f"[{col}:{col + n}] : {name}")
        col += n
    return ",\n ".join(lines)

class FormatAdaptor:
    def __init__(self, filename, mode="r", version=(1,0), compression=False):
        if version[0] == 1:
            with open(os.path.join(os.path.dirname(__file__), "format_v1.json"), "r", encoding="utf-8") as f:
                self.format = json.load(f)
        self.attributes = self.format["attributes"]
        self.symbols_format = self.format["symbols"]
        self.data_format = self.format["data"]
        self.compression = compression
        self.hdf5 = h5py.File(filename, mode)
        self.writable = self.hdf5.mode != "r"
        self.data = self.hdf5["data"] if "data" in self.hdf5 else self.hdf5.create_group("data") if self.writable else None
        self.dimension = int(self.hdf5.attrs.get("dimension", self.attributes["dimension"]["default"]))
        check_datatype(self.dimension, self.attributes["dimension"])

    @property
    def symbols(self):
        """Value of each symbol for the current dimension of the file (it may change while writing)."""
        dimension = str(int(self.hdf5.attrs.get("dimension", self.attributes["dimension"]["default"])))
        return {name: values[dimension] for name, values in self.symbols_format.items()}

    def open(self, exclude=()):
        """Every entry of the format in /data: created if missing (writable file), else None if
        the format says it is optional ("opt"), else reading fails. Entries named at write
        time ('<...>') and those of exclude are left to the caller."""
        if self.data is None:
            raise ValueError(f"no /data group in {self.hdf5.filename}")
        for name, spec in self.data_format.items():
            if "<" in name or name in exclude or name in self.data:
                continue
            if self.writable:
                self.get_group(name) if "dtype" not in spec else self.get_dataset(name)
            elif spec.get("opt"):
                print(f"[io.mechanics] Warning - '{name}' not in {self.hdf5.filename}")
            else:
                raise ValueError(f"'{name}' not in {self.hdf5.filename}")

    def __getitem__(self, name):
        """Entry of /data, None if absent."""
        return self.data.get(name) if self.data is not None else None

    def flush(self):
        self.hdf5.flush()

    def close(self):
        self.hdf5.close()

    def _spec(self, path):
        """Format entry of a path relative to /data ('<...>' keys match any name)."""
        spec = self.data_format
        for part in path.strip("/").split("/"):
            if part in spec:
                spec = spec[part]
                continue
            for key, sub in spec.items():
                if "<" in key and re.fullmatch(re.sub(r"<[^>]*>", ".+", key), part):
                    spec = sub
                    break
            else:
                raise KeyError(f"'{path}' not described by the format")
        return spec

    def add_into_group(self, group, dataset, /, data=None, compression=False, **dataset_kwargs):
        #add interaction, group : permanent interaction, dataset : <interaction>, ...
        current_group = self.get_group(group)
        if dataset not in current_group:
            dataset_kwargs = {k: v for k, v in dataset_kwargs.items() if v is not None} # None = not given
            spec = self._spec(f"{group}/{dataset}")
            attrs = spec.get("attrs", {})
            unknown = dataset_kwargs.keys() - attrs.keys()
            if unknown:
                raise KeyError(f"{group}/{dataset}: attributes {sorted(unknown)} not in the format")
            for name, fmt in attrs.items():
                if name in dataset_kwargs: check_datatype(dataset_kwargs[name], fmt, self.symbols)
                elif not fmt.get("opt"):
                    raise ValueError(f"{group}/{dataset}: required attribute '{name}' missing")
            if "dtype" not in spec: # no dtype in the format: group (input/<object>, external_functions/<function>)
                data_hdf5 = current_group.create_group(dataset)
            else:
                if data is not None: # data already carries its types: written as is
                    check_datatype(data, spec, self.symbols)
                    data_hdf5 = current_group.create_dataset(
                        dataset, data=data, compression="gzip" if compression else None)
                else: # empty structure built from the format
                    data_hdf5 = current_group.create_dataset(
                        dataset, shape=spec["shape"], dtype=vlen_dtype if spec["dtype"] == "vlen" else spec["dtype"],
                        compression="gzip" if compression else None)
            data_hdf5.attrs.update(dataset_kwargs)
        return current_group[dataset]

    def get_dataset(self, dataset):
        if dataset in self.data:
            return self.data[dataset]
        if not self.writable:
            return None
        spec = self._spec(dataset)
        ncols = spec["shape"][-1]
        comp = self.compression and ncols > 0
        table = self.data.create_dataset(dataset, (0, ncols), maxshape=(None, ncols), dtype=spec["dtype"],
                                         chunks=(4000, ncols) if comp else None,
                                         compression="gzip" if comp else None, compression_opts=9 if comp else None)
        if "layout" in spec:
            table.attrs["info"] = info(spec["layout"])
        return table

    def get_group(self, group):
        if group in self.data:
            return self.data[group]
        if not self.writable:
            return None
        return self.data.require_group(group)
