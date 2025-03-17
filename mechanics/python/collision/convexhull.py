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
import math
import numpy
import scipy.spatial  # For convexHull


class Simplex(object):
    def __init__(self, coordinates):
        if not len(coordinates) == 4:
            raise RuntimeError("You must provide only 4 coordinates!")
        self._coordinates = numpy.array(coordinates)

    def volume(self):
        """
        volume: Return volume of simplex.
        Formula from http://de.wikipedia.org/wiki/Tetraeder
        """
        vA = numpy.array(self._coordinates[1]) - numpy.array(self._coordinates[0])
        vB = numpy.array(self._coordinates[2]) - numpy.array(self._coordinates[0])
        vC = numpy.array(self._coordinates[3]) - numpy.array(self._coordinates[0])

        return numpy.abs(numpy.dot(numpy.cross(vA, vB), vC)) / 6.0

    def det(self):
        # import numpy
        # J = numpy.zeros((3,3))
        # J[:,0]= numpy.array(self._coordinates[1])-numpy.array(self._coordinates[0])
        # J[:,1]= numpy.array(self._coordinates[2])-numpy.array(self._coordinates[0])
        # J[:,2]= numpy.array(self._coordinates[3])-numpy.array(self._coordinates[0])
        J = self._coordinates[:-1] - self._coordinates[-1]
        # print('J',J)
        det = abs(numpy.linalg.det(J))
        # print('det', det, 'volume',det/6.0)
        return det

    def centroid(self):

        # compute centroid
        cm = numpy.zeros(3)
        for coords in self._coordinates:
            cm = cm + coords
        cm = cm / len(self._coordinates)
        # print('centroid',cm)
        return cm

    def inertia(self, G):
        """
        inertia : Return the inertia w.r.t the global axis and the point G
        """
        #
        det = self.det()

        # change of variable
        v1 = self._coordinates[0] - numpy.array(G)
        v2 = self._coordinates[1] - numpy.array(G)
        v3 = self._coordinates[2] - numpy.array(G)
        v4 = self._coordinates[3] - numpy.array(G)

        imat = numpy.zeros((3, 3))
        imat[0, 0] = (
            det
            * (
                v1[1] * v1[1]
                + v1[1] * v2[1]
                + v2[1] * v2[1]
                + v1[1] * v3[1]
                + v2[1] * v3[1]
                + v3[1] * v3[1]
                + v1[1] * v4[1]
                + v2[1] * v4[1]
                + v3[1] * v4[1]
                + v4[1] * v4[1]
                + v1[2] * v1[2]
                + v1[2] * v2[2]
                + v2[2] * v2[2]
                + v1[2] * v3[2]
                + v2[2] * v3[2]
                + v3[2] * v3[2]
                + v1[2] * v4[2]
                + v2[2] * v4[2]
                + v3[2] * v4[2]
                + v4[2] * v4[2]
            )
            / 60.0
        )
        imat[1, 1] = (
            det
            * (
                v1[0] * v1[0]
                + v1[0] * v2[0]
                + v2[0] * v2[0]
                + v1[0] * v3[0]
                + v2[0] * v3[0]
                + v3[0] * v3[0]
                + v1[0] * v4[0]
                + v2[0] * v4[0]
                + v3[0] * v4[0]
                + v4[0] * v4[0]
                + v1[2] * v1[2]
                + v1[2] * v2[2]
                + v2[2] * v2[2]
                + v1[2] * v3[2]
                + v2[2] * v3[2]
                + v3[2] * v3[2]
                + v1[2] * v4[2]
                + v2[2] * v4[2]
                + v3[2] * v4[2]
                + v4[2] * v4[2]
            )
            / 60.0
        )
        imat[2, 2] = (
            det
            * (
                v1[0] * v1[0]
                + v1[0] * v2[0]
                + v2[0] * v2[0]
                + v1[0] * v3[0]
                + v2[0] * v3[0]
                + v3[0] * v3[0]
                + v1[0] * v4[0]
                + v2[0] * v4[0]
                + v3[0] * v4[0]
                + v4[0] * v4[0]
                + v1[1] * v1[1]
                + v1[1] * v2[1]
                + v2[1] * v2[1]
                + v1[1] * v3[1]
                + v2[1] * v3[1]
                + v3[1] * v3[1]
                + v1[1] * v4[1]
                + v2[1] * v4[1]
                + v3[1] * v4[1]
                + v4[1] * v4[1]
            )
            / 60.0
        )
        # a'
        imat[1, 2] = (
            -det
            * (
                2.0 * v1[1] * v1[2]
                + v2[1] * v1[2]
                + v3[1] * v1[2]
                + v4[1] * v1[2]
                + v1[1] * v2[2]
                + 2.0 * v2[1] * v2[2]
                + v3[1] * v2[2]
                + v4[1] * v2[2]
                + v1[1] * v3[2]
                + v2[1] * v3[2]
                + 2.0 * v3[1] * v3[2]
                + v4[1] * v3[2]
                + v1[1] * v4[2]
                + v2[1] * v4[2]
                + v3[1] * v4[2]
                + 2.0 * v4[1] * v4[2]
            )
            / 120.0
        )
        imat[2, 1] = imat[1, 2]
        # b'
        imat[0, 2] = (
            -det
            * (
                2.0 * v1[0] * v1[2]
                + v2[0] * v1[2]
                + v3[0] * v1[2]
                + v4[0] * v1[2]
                + v1[0] * v2[2]
                + 2.0 * v2[0] * v2[2]
                + v3[0] * v2[2]
                + v4[0] * v2[2]
                + v1[0] * v3[2]
                + v2[0] * v3[2]
                + 2.0 * v3[0] * v3[2]
                + v4[0] * v3[2]
                + v1[0] * v4[2]
                + v2[0] * v4[2]
                + v3[0] * v4[2]
                + 2.0 * v4[0] * v4[2]
            )
            / 120.0
        )
        imat[2, 0] = imat[0, 2]
        # c'
        imat[0, 1] = (
            -det
            * (
                2.0 * v1[0] * v1[1]
                + v2[0] * v1[1]
                + v3[0] * v1[1]
                + v4[0] * v1[1]
                + v1[0] * v2[1]
                + 2.0 * v2[0] * v2[1]
                + v3[0] * v2[1]
                + v4[0] * v2[1]
                + v1[0] * v3[1]
                + v2[0] * v3[1]
                + 2.0 * v3[0] * v3[1]
                + v4[0] * v3[1]
                + v1[0] * v4[1]
                + v2[0] * v4[1]
                + v3[0] * v4[1]
                + 2.0 * v4[0] * v4[1]
            )
            / 120.0
        )
        imat[1, 0] = imat[0, 1]
        [p, v] = numpy.linalg.eig(imat)
        # print('Principal inertia:')
        # print(p)
        # print('Principal direction :')
        # print(v)
        return imat


class ConvexHull(object):

    def __init__(self, coordinates):
        """
        Constructor
        """

        if len(coordinates) < 4:
            raise RuntimeError("You must provide at least 4 coordinates!")
        self._coordinates = numpy.array(coordinates)
        self.hull = scipy.spatial.ConvexHull(self._coordinates)

    def centroid(self):
        cm = numpy.zeros(3)
        volume = self.volume_divergence_theorem()
        for vertices in self.hull.vertices:
            a = self._coordinates[vertices[0]]
            b = self._coordinates[vertices[1]]
            c = self._coordinates[vertices[2]]
            n = -numpy.cross(b - a, c - a)
            cm[:] += (
                1.0
                / (48.0 * volume)
                * n[:]
                * ((a[:] + b[:]) ** 2 + (c[:] + b[:]) ** 2 + (c[:] + a[:]) ** 2)
            )
        # print('cm ---------- ' ,cm)
        # print('barycenter ---------', self.barycenter())
        return cm

    __centroid = centroid

    def barycenter(self):
        # compute barycenter
        b = numpy.zeros(3)
        for coords in self._coordinates:
            b = b + coords
        b = b / len(self._coordinates)
        # print('barycenter',cm)
        return b

    def volume_divergence_theorem(self):
        volume = 0
        for vertices in self.hull.vertices:
            a = self._coordinates[vertices[0]]
            b = self._coordinates[vertices[1]]
            c = self._coordinates[vertices[2]]
            n = -numpy.cross(b - a, c - a)
            # bary = self.barycenter()
            # print "numpy.dot(bary-a,n)", numpy.dot(bary-a,n)
            # assert (numpy.dot(bary-a,n)<0)
            # we assume that the facet are oriented by qhull
            # n= n/numpy.linalg.norm(n)
            volume += 1 / 6.0 * numpy.dot(a, n)
        return volume

    def volume(self):

        volume = 0
        # compute centroid
        c = self.__centroid()
        for vertices in self.hull.vertices:
            coords = [list(c)]
            for i in vertices:
                coords.append(self._coordinates[i])
            simplex = Simplex(coords)

            # print(vertices)
            # print(coords)
            # print(simplex.volume())
            volume += simplex.volume()

        return volume

    def inertia(self, G):
        volume = 0
        imat = numpy.zeros((3, 3))
        # compute centroid
        c = self.__centroid()
        for vertices in self.hull.vertices:
            coords = [list(c)]
            for i in vertices:
                coords.append(self._coordinates[i])
            simplex = Simplex(coords)
            volume += simplex.volume()
            imat += simplex.inertia(G)
        # print('inertia of convexHull:')
        # print(imat)
        # [p,v]=numpy.linalg.eig(imat)
        # print('Principal inertia:')
        # print(p)
        # print('Principal direction :')
        # print(v)

        return imat, volume


class ConvexHull2d(ConvexHull):

    def __init__(self, coordinates):
        """
        Constructor
        """

        if len(coordinates) < 3:
            raise RuntimeError("You must provide at least 3 coordinates!")

        # construct a list of 3D coordinates for a extruded polydron
        nb_vertices = len(coordinates)
        coord_3d = list(coordinates)
        coord_3d.extend(coordinates)
        for i, v in enumerate(coord_3d):
            if i < nb_vertices:
                coord_3d[i] = numpy.append(v, 0.0)
            else:
                coord_3d[i] = numpy.append(v, 1.0)

        super(ConvexHull2d, self).__init__(coord_3d)

    def centroid(self):
        cm = super(ConvexHull2d, self).centroid()
        return cm[0:2]

    def barycenter(self):
        b = super(ConvexHull2d, self).barycenter()
        return b[0:2]

    def area(self):

        area = super(ConvexHull2d, self).volume()
        return area

    def inertia(self, G):

        G = list(G)
        G.append(0.5)
        imat, area = super(ConvexHull2d, self).inertia(G)
        return imat[2, 2], area


if __name__ == "__main__":
    print("####### first example #########")
    coords = []

    coords.append([8.33220, -11.86875, 0.93355])
    coords.append([0.75523, 5.0000, 16.37072])
    coords.append([52.61236, 5.0000, -5.38580])
    coords.append([2.000, 5.0000, 3.0000])
    s = Simplex(coords)
    print("volume of simplex", s.volume())

    imat = s.inertia(s.centroid())
    print("inertia:")
    print(imat)
    I_check = numpy.zeros((3, 3))
    I_check[0, 0] = 43520.33257
    I_check[1, 1] = 194711.28938
    I_check[2, 2] = 191168.76173
    I_check[1, 2] = -4417.66150
    I_check[2, 1] = I_check[1, 2]
    I_check[0, 2] = 46343.16662
    I_check[2, 0] = I_check[0, 2]
    I_check[0, 1] = -11996.20119
    I_check[1, 0] = I_check[0, 1]

    print("I_check", I_check)
    print("check", I_check - imat)
    print("error", numpy.linalg.norm(I_check - imat) / numpy.linalg.norm(I_check))

    print("####### second example #########")

    coords = []

    coords.append([0, 0, 0])
    coords.append([1, 0, 0])
    coords.append([0, 1, 0])
    coords.append([0, 0, 1])

    s = Simplex(coords)
    print("volume of simplex", s.volume())
    print(s.inertia(s.centroid()))
    # t = Polyeder(coords)
    # #print t.volume()

    # print ('volume', t.volume_by_cm())

    print("####### third example unit cube #########")

    coords.append([0, 1, 1])
    coords.append([1, 0, 1])
    coords.append([1, 1, 0])
    coords.append([1, 1, 1])

    p = ConvexHull(coords)
    # #print p.volume()

    print("volume", p.volume())
    print("volume_divergence_theorem", p.volume_divergence_theorem())
    imat, v = p.inertia(p.centroid())
    print("inertia:")
    print(imat)

    I_check = numpy.zeros((3, 3))
    I_check[0, 0] = (1 + 1) / 12.0
    I_check[1, 1] = (1 + 1) / 12.0
    I_check[2, 2] = (1 + 1) / 12.0
    # print('I_check',I_check)
    # print('check',I_check-imat)
    print("error", numpy.linalg.norm(I_check - imat) / numpy.linalg.norm(I_check))

    print("####### Fourth example octahedron of side #########")

    side = 1.0
    s = side / math.sqrt(2.0)
    coords = []
    coords.append([0, s, 0])
    coords.append([0, -s, 0])
    coords.append([s, 0, 0])
    coords.append([-s, 0, 0])
    coords.append([0, 0, s])
    coords.append([0, 0, -s])

    print(coords)
    p = ConvexHull(coords)
    # #print p.volume()

    print("volume", p.volume())
    print("volume_divergence_theorem", p.volume_divergence_theorem())

    centro = p.centroid()
    print("centroid", centro)

    imat, volume = p.inertia(centro)
    print("volume", volume)
    print("inertia:")
    print(imat)

    I_check = numpy.zeros((3, 3))
    I_check[0, 0] = volume * (side * side) / 10.0
    I_check[1, 1] = volume * (side * side) / 10.0
    I_check[2, 2] = volume * (side * side) / 10.0
    print("I_check", I_check)
    print("check", I_check - imat)
    print("error", numpy.linalg.norm(I_check - imat) / numpy.linalg.norm(I_check))
