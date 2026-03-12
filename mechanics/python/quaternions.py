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
#

"""Utilities to handle quaternions in the siconos.mechanics package

usage example:

    >>> import siconos.mechanics.quaternions as quat
    >>> phi, theta, psi = quat.euler_from_quaternion(q0, q1, q2, q3)

"""
import numpy as np
from math import cos, sin


def euler_from_quaternion(q0, q1, q2, q3):
    """
    Convert a quaternion into Euler angles (phi, theta, psi).

    Args:
        q0 (float or ndarray): Composante scalaire du quaternion.
        q1 (float or ndarray): Première composante vectorielle.
        q2 (float or ndarray): Deuxième composante vectorielle.
        q3 (float or ndarray): Troisième composante vectorielle.

    Returns:
        tuple of ndarray: (phi, theta, psi), angles d'Euler en radians.

    Example:
        >>> euler_from_quaternion(1, 0, 0, 0)
        (0.0, 0.0, 0.0)
    """
    phi = np.arctan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1**2 + q2**2))
    theta = np.arcsin(2 * (q0 * q2 - q3 * q1))
    psi = np.arctan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2**2 + q3**2))

    return phi, theta, psi


def quaternion_get(orientation):
    """
    Converts the orientation (either axis-angle or quaternion) into a
    quaternion representation.

    If the input is already a quaternion, it is returned as-is.

    Parameters
    ----------
    orientation : np array or tuple
        The orientation in either of the following formats:
        - Axis-angle: ([x, y, z], angle) where [x, y, z] is a 3D vector and angle is a float.
        - Quaternion: [w, x, y, z] where w, x, y, z are the quaternion components.
    """

    if isinstance(orientation, np.ndarray) or (
        isinstance(orientation, list) and len(orientation) == 4
    ):
        orientation = np.asarray(orientation, dtype=np.float64)

        if orientation.shape[0] != 4:
            raise AssertionError("quaternion_get. The quaternion must be of size 4")
        return orientation

    if isinstance(orientation, tuple) or isinstance(orientation, list):
        if len(orientation) != 2:  # orientation is a tuple
            raise AssertionError("quaternion_get. Wrong input format")
        # axis + angle
        axis, angle = orientation
        if len(axis) != 3:
            raise ValueError("Axis must be a 3D vector.")
        if not isinstance(angle, float):
            raise ValueError("Angle must be a float.")
        axis_norm = np.linalg.norm(axis)
        if axis_norm == 0:
            raise ValueError("Axis vector cannot be the zero vector.")
        n = sin(angle / 2.0) / axis_norm
        return np.array([cos(angle / 2.0), axis[0] * n, axis[1] * n, axis[2] * n])

    else:
        raise ValueError(
            "Orientation must be either an axis-angle pair (3D vector + float)"
            + " or a quaternion."
        )


def quaternion_multiply(q1, q0):
    w0, x0, y0, z0 = q0
    w1, x1, y1, z1 = q1
    return np.array(
        [
            -x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
            x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
            -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
            x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0,
        ],
        dtype=np.float64,
    )
