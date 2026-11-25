# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# SPDX-FileCopyrightText: INRIA
# SPDX-License-Identifier: Apache-2.0
#

"""Reference trajectories for the python scripts of the vkernel tests, as
ref_check.hpp for the C++ ones.

The trajectory of a run is the rows of data/dynamic (t, id, positions) of
its hdf5 file, followed by the velocities of data/velocities (the same rows:
one per dynamic object and output step). It is compared to <name>.ref, next
to this file, by the relative Frobenius error |d - r| / max(|r|, 1).

The name of a reference carries the parameters of the run it is for: a run
with other parameters has none, and is not checked, unless
SICONOS_CHECK_REF=1 (the ctest tests), where a missing reference fails.
With SICONOS_WRITE_REF=1 in the environment, the reference is written
instead: run the script once that way, then check the .ref in.
"""

import os

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def trajectory(hdf5, every=1):
    """the rows of every every-th output time, and of the last one"""
    with h5py.File(hdf5, "r") as f:
        dynamic = np.asarray(f["data/dynamic"][:], dtype=np.float64)
        velocities = np.asarray(f["data/velocities"][:], dtype=np.float64)
    if dynamic.shape[0] != velocities.shape[0]:
        raise ValueError(
            f"{hdf5}: {dynamic.shape[0]} positions, {velocities.shape[0]} velocities"
        )
    times, k = np.unique(dynamic[:, 0], return_inverse=True)
    keep = (k % every == 0) | (k == len(times) - 1)
    return np.hstack([dynamic, velocities[:, 2:]])[keep]


def check(hdf5, name, tol, every=1):
    """True if the trajectory of hdf5 (every every-th output time) is the one
    of <name>.ref within tol (prints the error); writes <name>.ref with
    SICONOS_WRITE_REF=1."""
    data = trajectory(hdf5, every)
    ref_file = os.path.join(HERE, name + ".ref")
    if os.environ.get("SICONOS_WRITE_REF") == "1":
        np.savetxt(ref_file, data, fmt="%.17g")
        print(f"{ref_file}: written, {data.shape[0]} rows")
        return True
    if not os.path.exists(ref_file):
        if os.environ.get("SICONOS_CHECK_REF") == "1":
            print(f"FAILED: no reference {ref_file}")
            return False
        print(f"no reference {ref_file} for these parameters: not checked")
        return True
    ref = np.loadtxt(ref_file, ndmin=2)
    if ref.shape != data.shape:
        print(f"FAILED: {ref_file}: shape {data.shape}, reference {ref.shape}")
        return False
    error = np.linalg.norm(data - ref) / max(np.linalg.norm(ref), 1.0)
    row, col = np.unravel_index(np.argmax(np.abs(data - ref)), data.shape)
    print(
        f"{ref_file}: relative error {error:.3e} (tolerance {tol:.1e}), worst at row "
        f"{row + 1} column {col + 1}: {data[row, col]!r} vs {ref[row, col]!r}"
    )
    if not error <= tol:  # NaN fails too
        print(f"FAILED: {ref_file}")
        return False
    return True
