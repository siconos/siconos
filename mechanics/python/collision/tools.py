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
"""
tools.py - Geometry and contact tools for Siconos mechanics module

This module provides classes to describe shapes, volumes, and contactors
used in mechanical simulations, with optional transformation, material, and
collision metadata.

Classes
-------
- MovedShape : Base class for geometry with translation and orientation.
- Shape : A named shape with optional instance identifier.
- Volume : A Shape with mass and density parameters.
- Contactor : A Shape belonging to a collision group, with contact metadata.
"""
from dataclasses import dataclass, field
from typing import Optional, Union, Tuple, List
from math import cos, sin
from numpy.linalg import norm
import numpy as np


@dataclass
class Material:
    """
    Some material properties that may be associated to shapes.

    Parameters
    ----------
    density : float, optional
        The material density (default is None).
    """

    density: Optional[float] = None


@dataclass
class MovedShape:
    """
    A geometrical shape moved by a translation and an orientation relative
    to its parent frame.

    The orientation can be provided either as a quaternion [w, x, y, z]
    or as an (axis, angle) pair.

    Parameters
    ----------
    shape_name : str
        Reference name or identifier of the shape.

    data : object, optional
        Actual shape geometry object, typically from an external CAD or OCC representation.

    relative_translation : tuple of 3 floats, optional
        Translation vector in the local body frame. Defaults to (0.0, 0.0, 0.0).

    relative_orientation : tuple, optional
        Either a quaternion [w, x, y, z], or (axis, angle) where
        axis is a list of 3 floats and angle is a float.
        Defaults to (1.0, 0.0, 0.0, 0.0), i.e., identity rotation.
    """

    shape_name: str
    data: Optional[object] = None
    relative_translation: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    relative_orientation: Union[
        Tuple[float, float, float, float], Tuple[Tuple[float, float, float], float]
    ] = (1.0, 0.0, 0.0, 0.0)
    translation: np.ndarray = field(init=False)
    orientation: np.ndarray = field(init=False)

    def __post_init__(self):
        # Convert the relative translation into a NumPy array (modifiable)
        # and easier to use in mechanics_run
        self.translation = np.array(self.relative_translation, dtype=np.float64)

        ori = self.relative_orientation

        if isinstance(ori, tuple) and len(ori) == 2:
            axis, angle = ori
            assert len(axis) == 3
            assert isinstance(angle, float)
            n = sin(angle / 2) / norm(axis)
            # Maybe it would be better to normalize axis ?
            ori = np.array(
                [
                    cos(angle / 2),
                    axis[0] * n,
                    axis[1] * n,
                    axis[2] * n,
                ],
                dtype=np.float64,
            )
        else:
            assert len(ori) == 4
            ori = np.array(ori, dtype=np.float64)

        self.orientation = ori


@dataclass
class Shape(MovedShape):
    """
    A MovedShape with an optional instance name.

    Parameters
    ----------
    shape_name : str
        Reference name or identifier of the shape.

    data : object, optional
        Actual shape geometry object, typically from an external CAD or OCC representation.

    instance_name : str, optional
        Optional identifier for this shape instance.

    relative_translation : tuple of 3 floats, optional
        Translation vector in the local body frame. Defaults to (0.0, 0.0, 0.0).

    relative_orientation : tuple, optional
        Either a quaternion [w, x, y, z], or (axis, angle) where
        axis is a list of 3 floats and angle is a float.
        Defaults to (1.0, 0.0, 0.0, 0.0), i.e., identity rotation.
    """

    instance_name: Optional[str] = None


@dataclass
class Volume(Shape):
    """
    A Shape with mass properties and material parameters.

    Inherits from Shape and adds physical characteristics like mass and density

    Parameters
    ----------
    shape_name : str
        Reference name of the shape.

    shape_data : object, optional
        Actual shape geometry object.

    instance_name : str, optional
        Optional unique name for this instance.

    mass : float, optional
        Total mass. If None, can be computed from density and volume.

    parameters : Material, optional
        Material properties (e.g., density).

    relative_translation : tuple of 3 floats, optional
        Translation in local coordinates. Defaults to (0, 0, 0).

    relative_orientation : tuple, optional
        Orientation as quaternion [w, x, y, z] or (axis, angle).
    """

    mass: Optional[float] = None
    parameters: Material = field(default_factory=lambda: Material(density=1.0))


@dataclass
class Contactor(Shape):
    """
    A Shape that belongs to a collision group and may carry contact-specific data.

    Depending on the geometrical engine used, some information may be added
    to the contactor, such as the kind of contact and the concerned part of the shape.
    Contact laws must then be defined between collision groups.


    Parameters
    ----------
    shape_name : str
        Reference name of the shape.

    shape_data : object, optional
        Actual shape geometry.

    instance_name : str, optional
        Identifier for this instance of the shape.

    collision_group : int, optional
        Group ID used in collision detection. Default is 0.

    parameters : object, optional
        Extra parameters attached to the contactor (e.g., friction, tags).

    contact_type : object, optional
        Type of contact (e.g., face, edge), engine-specific.

    contact_index : int, optional
        Index of the concerned sub-entity.

    relative_translation : tuple of 3 floats, optional
        Translation in local coordinates. Defaults to (0, 0, 0).

    relative_orientation : tuple, optional
        Orientation as quaternion [w, x, y, z] or (axis, angle).
    """

    collision_group: int = 0
    parameters: Optional[object] = None
    contact_type: Optional[object] = None
    contact_index: Optional[int] = None

    def __post_init__(self):
        super().__post_init__()
        # We keep 'group' for compat with legacy version
        # but it's better to use collision_group everywhere
        self.group = self.collision_group
