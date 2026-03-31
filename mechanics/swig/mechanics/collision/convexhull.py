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
        volume: Return volume of simplex. Formula from http://de.wikipedia.org/wiki/Tetraeder
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
        return abs(numpy.linalg.det(J))

    def centroid(self):
        return numpy.mean(self.coordinates, axis=0)

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

        I = numpy.zeros((3, 3))
        I[0, 0] = (
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
        I[1, 1] = (
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
        I[2, 2] = (
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
        I[1, 2] = (
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
        I[2, 1] = I[1, 2]
        # b'
        I[0, 2] = (
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
        I[2, 0] = I[0, 2]
        # c'
        I[0, 1] = (
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
        I[1, 0] = I[0, 1]
        [p, v] = numpy.linalg.eig(I)
        # print('Principal inertia:')
        # print(p)
        # print('Principal direction :')
        # print(v)
        return I


class ConvexHull(object):

    def __init__(self, coordinates):
        """
        Constructor
        """

        if len(coordinates) < 4:
            raise RuntimeError("You must provide at least 4 coordinates!")
        self.hull = scipy.spatial.ConvexHull(coordinates)

    def centroid_3d(self):
        cm = numpy.zeros(3)
        vol = self.hull.volume  # self.volume_divergence_theorem()
        # volume is a method of Scipy ConvexHull
        inside = numpy.mean(self.hull.points, axis=0)
        for vertices in self.hull.simplices:
            a = self.hull.points[vertices[0]]
            b = self.hull.points[vertices[1]]
            c = self.hull.points[vertices[2]]
            n = numpy.cross(b - a, c - a)
            # ensure proper normal orientation
            if numpy.dot(n, a - inside) < 0:
                n = -n
            cm[:] += (
                1.0
                / (48.0 * vol)
                * n[:]
                * ((a[:] + b[:]) ** 2 + (c[:] + b[:]) ** 2 + (c[:] + a[:]) ** 2)
            )
        return cm

    def centroid(self):
        return self.centroid_3d()


    def barycenter(self):
        # compute barycenter
        b = numpy.zeros(3)
        for coords in self.hull.points:
            b = b + coords
        b = b / len(self.hull.points)
        # print('barycenter',cm)
        return b

    def volume_divergence_theorem(self):
        volume = 0
        inside = numpy.mean(self.hull.points, axis=0)
        for vertices in self.hull.simplices:
            a = self.hull.points[vertices[0]]
            b = self.hull.points[vertices[1]]
            c = self.hull.points[vertices[2]]
            n = numpy.cross(b - a, c - a)
            # ensure proper normal orientation
            if numpy.dot(n, a - inside) < 0:
                n = -n
            volume += 1 / 6.0 * numpy.dot(a, n)
        return volume

    def volume(self):

        volume = 0
        # compute centroid
        c = self.centroid()
        for vertices in self.hull.simplices:
            coords = [list(c)]
            for i in vertices:
                coords.append(self.hull.points[i])
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
        c = self.centroid_3d()
        for vertices in self.hull.simplices:
            coords = [list(c)]
            for i in vertices:
                coords.append(self.hull.points[i])
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
        I, area = super(ConvexHull2d, self).inertia(G)
        return I[2, 2], area


if __name__ == "__main__":
    print("####### first example #########")
    coords = []

    coords.append([8.33220, -11.86875, 0.93355])
    coords.append([0.75523, 5.0000, 16.37072])
    coords.append([52.61236, 5.0000, -5.38580])
    coords.append([2.000, 5.0000, 3.0000])
    s = Simplex(coords)
    print("volume of simplex", s.volume())

    I = s.inertia(s.centroid())
    print("inertia:")
    print(I)

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
    print("check", I_check - I)
    print("error", numpy.linalg.norm(I_check - I) / numpy.linalg.norm(I_check))

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
    I, v = p.inertia(p.centroid())
    print("inertia:")
    print(I)

    I_check = numpy.zeros((3, 3))
    I_check[0, 0] = (1 + 1) / 12.0
    I_check[1, 1] = (1 + 1) / 12.0
    I_check[2, 2] = (1 + 1) / 12.0
    # print('I_check',I_check)
    # print('check',I_check-I)
    print("error", numpy.linalg.norm(I_check - I) / numpy.linalg.norm(I_check))

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

    I, volume = p.inertia(centro)
    print("volume", volume)
    print("inertia:")
    print(I)

    I_check = numpy.zeros((3, 3))
    I_check[0, 0] = volume * (side * side) / 10.0
    I_check[1, 1] = volume * (side * side) / 10.0
    I_check[2, 2] = volume * (side * side) / 10.0
    print("I_check", I_check)
    print("check", I_check - I)
    print("error", numpy.linalg.norm(I_check - I) / numpy.linalg.norm(I_check))
