# SPDX-FileCopyrightText: INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
import numpy as np
import siconos.numerics as sn


def test_projection():
    r = np.array([1.0, 1.0, 1.0])
    status = sn.soc_projection_3d(r, 0.5)
    rr = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
    status_r = sn.rolling_cone_projection(rr, 0.5, 0.5)
    return status, status_r
