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
import vtk

vtkMath = vtk.vtkMath()


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
    Get quaternion from orientation
    """
    if len(orientation) == 2:
        # axis + angle
        axis = orientation[0]
        if not (len(axis) == 3):
            raise AssertionError("quaternion_get. axis in not a 3D vector")
        angle = orientation[1]
        if not isinstance(angle, float):
            raise AssertionError("quaternion_get. angle must be a float")
        n = sin(angle / 2.0) / np.linalg.norm(axis)
        ori = [cos(angle / 2.0), axis[0] * n, axis[1] * n, axis[2] * n]
    else:
        if not (len(orientation) == 4):
            raise AssertionError("quaternion_get. The quaternion must be of size 4")
        # a given quaternion
        ori = orientation
    return ori


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


class Quaternion:

    def __init__(self, *args):
        self._data = vtk.vtkQuaternion[float](*args)

    def __mul__(self, q):
        r = Quaternion()
        vtkMath.MultiplyQuaternion(self._data, q._data, r._data)
        return r

    def __getitem__(self, i):
        return self._data[i]

    def conjugate(self):
        r = Quaternion((self[0], self[1], self[2], self[3]))
        r._data.Conjugate()
        return r

    def rotate(self, v):
        pv = Quaternion((0, v[0], v[1], v[2]))
        rv = self * pv * self.conjugate()
        # assert(rv[0] == 0)
        return [rv[1], rv[2], rv[3]]

    def axisAngle(self):
        r = [0, 0, 0]
        a = self._data.GetRotationAngleAndAxis(r)
        return r, a
