# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2026 INRIA.
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
#

"""
Tools to generate random convex blocks (mostly used in rockfall simulations)
"""

import numpy as np
import trimesh
from dataclasses import dataclass

from siconos.mechanics.collision.convexhull import ConvexHull
from siconos.mechanics.collision.tools import Contactor

# Note FP: all this stuff should be in mechanics but we need io input ...


@dataclass
class RockShapeConfig:
    """
    Parameters to define a 'random' shape block

    Attributes
    ----------
    nb_pts : int
        number of vertices of the generated convex hull
        Must be >= 4. Default=50.
    y_aspect_ratio : float
        ratio (relative to x axis) along y-axis. Default=1
        0.8 means the block is 80% as wide in y as in x.
    z_aspect_ratio : float
        same as y_ratio in z dir. Default = 1.
    volume_min : float
        minimum volume for the block
    volume_max : float
        max volume for the block
    """

    nb_pts: int = 50
    y_aspect_ratio: float = 1.0
    z_aspect_ratio: float = 1.0
    volume_min: float = 1.0
    volume_max: float = 2.0


@dataclass
class RocksDropConfig:
    """
    Parameters to define

    Attributes
    ----------
    drop_zone : string
        input file (stl) defining the drop zone
    number_of_rocks : int
        number of rocks to be dropped. Default=10
    height_fall_min : float
        mininmal drop height
    height_fall_max : float
        maximal drop height

        dest_vol : float
        Required volume of the generated block. Default = 1.
    """

    drop_zone: str
    number_of_rocks: int = 10
    height_fall_min: float = 2.0
    height_fall_max: float = 2.0


def normalize_shape(vertices, y_aspect_ratio, z_aspect_ratio, dest_vol):
    """
    Transform a raw set of 3D points into a 'block', convex and normalized
    to match aspect ratios and a given volume.

    Parameters
    ----------
    vertices : array like object of shape (nb_points, 3)
        The 3D points set.

    y_aspect_ratio : float
        normalisation ratio in y direction relative to x.
        For example, 0.8 means the block is 80% as wide in y as in z.
    z_aspect_ratio : float
        same as y_aspect_ration for z dir.
    dest_vol : float
        Targeted volume of the normalized shape.

    Returns
    -------
    vertices : ndarray of shape (new_nb_points, 3)
        Convex hull vertices, centered at the origin, normalized
        to fit with aspect ratios and volume.
        Interior (concave) points of initial set have been removed.

    """
    vertices = np.array(vertices)
    # 1- Cleanup vertices (= remove concave points):
    # NOTE: phch.vertices are not really vertices, but rather a list of index of points
    # (into self.__base_vertices__)
    # that forms the faces of the convex shell :
    # [[i_face1_pt1, i_face1_pt2, i_face3_pt3], [i_face2_pt1, i_face2_pt2, i_face2_pt3], ...]
    # So using numpy "flatten" then "unique" is a way to select only the points
    # into self.__base_vertices__
    #  that are involved in the convex shell (concave points are filtered out).

    siconos_convex_hull = ConvexHull(vertices)
    phch = siconos_convex_hull.hull  # "PyHull Convex Hull"
    vertices = phch.points[np.sort(np.unique(np.array(phch.vertices).flatten()))]

    # 2- transform the shape to fit y_aspect_ratio and z_aspect_ratio:

    #        2.1- transform x coordinates to exactly fit [0, 1] interval.
    vertices[:, 0] -= vertices[:, 0].min()
    vertices[:, 0] /= vertices[:, 0].max()

    #        2.2- transform y coordinates to exactly fit [0, y_aspect_ratio] interval.
    vertices[:, 1] -= vertices[:, 1].min()
    vertices[:, 1] /= vertices[:, 1].max()
    vertices[:, 1] *= y_aspect_ratio

    #        2.3- transform z coordinates to exactly fit [0, z_aspect_ratio] interval.
    vertices[:, 2] -= vertices[:, 2].min()
    vertices[:, 2] /= vertices[:, 2].max()
    vertices[:, 2] *= z_aspect_ratio
    # new convex hull, cleaned and with correct aspect ratio.
    siconos_convex_hull = ConvexHull(vertices)

    # 3- Move the shape center of mass to origin [0, 0, 0]:
    vertices -= siconos_convex_hull.centroid()
    siconos_convex_hull = ConvexHull(vertices)  # center the shape

    # 4- Set volume
    # assume density:=1 for now, so inertia will be recomputed later as
    # it is proportionnal to density
    base_inertia, start_volume = siconos_convex_hull.inertia([0.0, 0.0, 0.0])
    vol_factor = dest_vol / start_volume
    coord_factor = vol_factor ** (1 / 3)
    vertices *= coord_factor  # homothetic transform -> now vol=dest_vol

    return vertices


def generate_shape(nbPts, y_aspect_ratio, z_aspect_ratio, dest_vol):
    """
    Generate a random convex 3D rock block shape.

    Parameters
    ----------
    nbPts : int
        Required number of vertices of the output shape. Must be >= 4.
    y_aspect_ratio : float
        Expected ratio along Y (see normalize_shape).
    z_aspect_ratio : float
        Expected ratio along Z (see normalize_shape).
    dest_vol : float
        Expected volume of the generated shape.

    Returns
    -------
    norm_vertices : ndarray of shape (nbPts, 3)
        Vertices of the generated convex block, centered at the origin,
        normalized to fit with aspect ratios and volume.
    """
    # np.random.seed(0)

    # NOTE: generate way too much points as pyhull (qhull) will have to ignore a lot
    # of them to get a CONVEX hull. So at first we don't have a good control on nbPts...
    base_vertices = np.random.rand(nbPts * 100, 3)
    norm_vertices = normalize_shape(
        base_vertices, y_aspect_ratio, z_aspect_ratio, dest_vol
    )

    # NOTE: after normalize_shape, the shape is now convex,
    # but we don't respect nbPts vertices,
    # we hope to have more so we just have to remove some of them.
    nb_points_to_remove = len(norm_vertices) - nbPts

    if nb_points_to_remove < 0:
        raise ValueError(
            "In random shape generation, nb_points_to_remove was negative:",
            nb_points_to_remove,
        )

    del_index_vector = np.random.randint(
        norm_vertices.shape[0], size=nb_points_to_remove
    )
    norm_vertices = np.delete(norm_vertices, del_index_vector, axis=0)

    # renormalize shape as we likely changed its volume and size above.
    norm_vertices = normalize_shape(
        norm_vertices, y_aspect_ratio, z_aspect_ratio, dest_vol
    )

    return norm_vertices


def generate_random_blocks(io, drop_config, rock_config):

    # Load the area above which the blocks will be dropped
    mesh = trimesh.load(drop_config.drop_zone)
    areas = mesh.area_faces
    prob = areas / areas.sum()

    # creation the set of blocks
    for i in range(0, drop_config.number_of_rocks):
        # print("n_rock", i)
        nameDs = "block" + str(i)
        nameShape = "block" + str(i) + "-shape"
        initial_angles = np.random.uniform(0, 2 * np.pi, 3)
        (c1, c2, c3), (s1, s2, s3) = np.cos(initial_angles), np.sin(initial_angles)
        initial_orientation = [
            c1 * c2 * c3 - s1 * s2 * s3,
            s1 * s2 * c3 + c1 * c2 * s3,
            s1 * c2 * c3 + c1 * s2 * s3,
            c1 * s2 * c3 - s1 * c2 * s3,
        ]

        # Position initiale
        height_fall = drop_config.height_fall_min + (
            drop_config.height_fall_max - drop_config.height_fall_min
        ) ** np.random.uniform(0, 1, 1)
        triangle_index = np.random.choice(len(mesh.faces), p=prob)
        triangle = mesh.triangles[triangle_index]
        A, B, C = triangle
        r1 = np.sqrt(np.random.rand())
        r2 = np.random.rand()
        point = (1 - r1) * A + r1 * (1 - r2) * B + r1 * r2 * C
        Xpos = point[0]
        Ypos = point[1]
        Zpos = point[2] + rock_config.volume_max**0.33 + height_fall[0]

        # posInit=[-0.5*heightmap.shape[0]+indx_dep[0][0],-0.5*heightmap.shape[1]+indx_dep[0][1],315]
        posInit = [Xpos, Ypos, Zpos]
        dest_vol = rock_config.volume_min + np.random.rand() * (
            rock_config.volume_max - rock_config.volume_min
        )
        vertices = generate_shape(
            rock_config.nb_pts,
            rock_config.y_aspect_ratio,
            rock_config.z_aspect_ratio,
            dest_vol,
        )
        # create the convexhull

        ch = ConvexHull(vertices)
        cm = ch.centroid()
        # move the vertices to center the center of mass at 0.0
        vertices = np.array(vertices)[:] - cm[:]
        # ch = ConvexHull(vertices)
        # cm = ch.centroid()
        inertia, area = ch.inertia(cm)
        mass_block = 2650.0 * 2.0

        inertia = inertia * mass_block
        # random block shape
        io.add_convex_shape(nameShape, vertices, outsideMargin=0.0)
        io.add_object(
            nameDs,
            [Contactor(nameShape, collision_group=100)],
            translation=posInit,
            velocity=[0, 0, 0, 0, 0, 0],
            orientation=initial_orientation,
            mass=mass_block,
            inertia=inertia,
        )
    print(f"Generation of {drop_config.number_of_rocks} blocks done.")
