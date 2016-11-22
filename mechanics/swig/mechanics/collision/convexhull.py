import math
import numpy
class Simplex(object):
    def __init__(self,coordinates):
        if not len(coordinates) == 4:
            raise RuntimeError('You must provide only 4 coordinates!')
        self._coordinates = numpy.array(coordinates)

    def volume(self):
        '''
        volume: Return volume of simplex. Formula from http://de.wikipedia.org/wiki/Tetraeder
        '''
        import numpy

        vA = numpy.array(self._coordinates[1]) - numpy.array(self._coordinates[0])
        vB = numpy.array(self._coordinates[2]) - numpy.array(self._coordinates[0])
        vC = numpy.array(self._coordinates[3]) - numpy.array(self._coordinates[0])

        return numpy.abs(numpy.dot(numpy.cross(vA,vB),vC)) / 6.0

    def det(self):
        # import numpy
        # J = numpy.zeros((3,3))
        # J[:,0]= numpy.array(self._coordinates[1])-numpy.array(self._coordinates[0])
        # J[:,1]= numpy.array(self._coordinates[2])-numpy.array(self._coordinates[0])
        # J[:,2]= numpy.array(self._coordinates[3])-numpy.array(self._coordinates[0])
        J=self._coordinates[:-1]-self._coordinates[-1]
        #print('J',J)
        det=abs(numpy.linalg.det(J))
        #print('det', det, 'volume',det/6.0)
        return det

    def centroid(self):
        import numpy
        # compute centroid
        cm = numpy.zeros(3)
        for coords in self._coordinates:
            cm = cm + coords
        cm = cm /len(self._coordinates)
        #print('centroid',cm)
        return cm

    def inertia(self, G):
        '''
        inertia : Return the inertia w.r.t the global axis and the point G
        '''
        #
        det=self.det()

        # change of variable
        v1=  self._coordinates[0] - numpy.array(G)
        v2 = self._coordinates[1] - numpy.array(G)
        v3 = self._coordinates[2] - numpy.array(G)
        v4 = self._coordinates[3] - numpy.array(G)


        I = numpy.zeros((3,3))
        I[0,0] = det*(
            v1[1]*v1[1]  + v1[1]*v2[1]+ v2[1]*v2[1]+ v1[1]*v3[1]+ v2[1]*v3[1]
            + v3[1]*v3[1]+ v1[1]*v4[1]+ v2[1]*v4[1]+ v3[1]*v4[1]+ v4[1]*v4[1]
            + v1[2]*v1[2]+ v1[2]*v2[2]+ v2[2]*v2[2]+ v1[2]*v3[2]+ v2[2]*v3[2]+ v3[2]*v3[2]
            + v1[2]*v4[2]+ v2[2]*v4[2]+ v3[2]*v4[2]+ v4[2]*v4[2] )/ 60.0
        I[1,1] = det*(
            v1[0]*v1[0]  + v1[0]*v2[0]+ v2[0]*v2[0]+ v1[0]*v3[0]+ v2[0]*v3[0]
            + v3[0]*v3[0]+ v1[0]*v4[0]+ v2[0]*v4[0]+ v3[0]*v4[0]+ v4[0]*v4[0]
            + v1[2]*v1[2]+ v1[2]*v2[2]+ v2[2]*v2[2]+ v1[2]*v3[2]+ v2[2]*v3[2]+v3[2]*v3[2]
            + v1[2]*v4[2]+ v2[2]*v4[2]+ v3[2]*v4[2]+ v4[2]*v4[2] )/ 60.0
        I[2,2] = det*(
            v1[0]*v1[0]  + v1[0]*v2[0]+ v2[0]*v2[0]+ v1[0]*v3[0]+ v2[0]*v3[0]
            + v3[0]*v3[0]+ v1[0]*v4[0]+ v2[0]*v4[0]+ v3[0]*v4[0]+ v4[0]*v4[0]
            + v1[1]*v1[1]+ v1[1]*v2[1]+ v2[1]*v2[1]+ v1[1]*v3[1]+ v2[1]*v3[1]+v3[1]*v3[1]
            + v1[1]*v4[1]+ v2[1]*v4[1]+ v3[1]*v4[1]+ v4[1]*v4[1] )/ 60.0
        # a'
        I[1,2] = - det*(
              2.0*v1[1]*v1[2]+ v2[1]*v1[2] + v3[1]*v1[2] +  v4[1]*v1[2] + v1[1]*v2[2]
            + 2.0*v2[1]*v2[2]+ v3[1]*v2[2] + v4[1]*v2[2] +  v1[1]*v3[2] + v2[1]*v3[2]
            + 2.0*v3[1]*v3[2]+ v4[1]*v3[2] + v1[1]*v4[2] +  v2[1]*v4[2] + v3[1]*v4[2]
            + 2.0*v4[1]*v4[2]) /120.0
        I[2,1] = I[1,2]
        # b'
        I[0,2] = - det*(
              2.0*v1[0]*v1[2]+ v2[0]*v1[2] + v3[0]*v1[2] +  v4[0]*v1[2] + v1[0]*v2[2]
            + 2.0*v2[0]*v2[2]+ v3[0]*v2[2] + v4[0]*v2[2] +  v1[0]*v3[2] + v2[0]*v3[2]
            + 2.0*v3[0]*v3[2]+ v4[0]*v3[2] + v1[0]*v4[2] +  v2[0]*v4[2] + v3[0]*v4[2]
            + 2.0*v4[0]*v4[2]) /120.0
        I[2,0] = I[0,2]
        # c'
        I[0,1] = - det*(
              2.0*v1[0]*v1[1]+ v2[0]*v1[1] + v3[0]*v1[1] +  v4[0]*v1[1] + v1[0]*v2[1]
            + 2.0*v2[0]*v2[1]+ v3[0]*v2[1] + v4[0]*v2[1] +  v1[0]*v3[1] + v2[0]*v3[1]
            + 2.0*v3[0]*v3[1]+ v4[0]*v3[1] + v1[0]*v4[1] +  v2[0]*v4[1] + v3[0]*v4[1]
            + 2.0*v4[0]*v4[1]) /120.0
        I[1,0] = I[0,1]
        [p,v]=numpy.linalg.eig(I)
        # print('Principal inertia:')
        # print(p)
        # print('Principal direction :')
        # print(v)
        return I

class ConvexHull(object):

    def __init__(self,coordinates):
        '''
        Constructor
        '''

        if len(coordinates) < 4:
            raise RuntimeError('You must provide at least 4 coordinates!')
        self._coordinates = numpy.array(coordinates)
        from pyhull.convex_hull import ConvexHull
        self.hull = ConvexHull(self._coordinates)

    # def volume(self):
    #     from pyhull.delaunay import DelaunayTri

    #     delaunay = DelaunayTri(self._coordinates,joggle=True)
    #     volume = 0
    #     #print(delaunay.points)
    #     #print(delaunay.vertices)
    #     #print(delaunay.simplices)
    #     for vertices in delaunay.vertices:

    #         coords = [self._coordinates[i] for i in vertices]
    #         simplex = Simplex(coords)
    #         #print(vertices)
    #         #print(coords)
    #         #print(simplex.volume())
    #         volume += simplex.volume()

    #     return volume

    def centroid(self):
        cm = numpy.zeros(3)
        volume = self.volume_divergence_theorem()
        for vertices in self.hull.vertices:
            a = self._coordinates[vertices[0]]
            b = self._coordinates[vertices[1]]
            c = self._coordinates[vertices[2]]
            n = - numpy.cross(b-a,c-a)
            cm[:] += 1.0/(48.0*volume)*n[:]*((a[:]+b[:])**2 +(c[:]+b[:])**2+(c[:]+a[:])**2)
        #print('cm ---------- ' ,cm)
        #print('barycenter ---------', self.barycenter())
        return cm

    def barycenter(self):
        # compute barycenter
        b = numpy.zeros(3)
        for coords in self._coordinates:
            b = b + coords
        b = b /len(self._coordinates)
        #print('barycenter',cm)
        return b

    def volume_divergence_theorem(self):
        volume = 0
        for vertices in self.hull.vertices:
            a = self._coordinates[vertices[0]]
            b = self._coordinates[vertices[1]]
            c = self._coordinates[vertices[2]]
            n = - numpy.cross(b-a,c-a)
            # bary = self.barycenter()
            # print "numpy.dot(bary-a,n)", numpy.dot(bary-a,n)
            # assert (numpy.dot(bary-a,n)<0)
            # we assume that the facet are oriented by qhull
            #n= n/numpy.linalg.norm(n)
            volume += 1/6.0 *numpy.dot(a,n)
        return volume

    def volume(self):

        volume = 0
        # compute centroid
        c = self.centroid()
        for vertices in self.hull.vertices:
            coords=[list(c)]
            for i in vertices:
                coords.append( self._coordinates[i])
            simplex = Simplex(coords)

            #print(vertices)
            #print(coords)
            #print(simplex.volume())
            volume += simplex.volume()

        return volume

    def inertia(self,G):
        volume = 0
        I = numpy.zeros((3,3))
        # compute centroid
        c = self.centroid()
        for vertices in self.hull.vertices:
            coords=[list(c)]
            for i in vertices:
                coords.append( self._coordinates[i])
            simplex = Simplex(coords)
            volume += simplex.volume()
            I += simplex.inertia(G)
        #print('inertia of convexHull:')
        #print(I)
        #[p,v]=numpy.linalg.eig(I)
        #print('Principal inertia:')
        #print(p)
        #print('Principal direction :')
        #print(v)

        return I,volume

if __name__ == '__main__':
    print('####### first example #########')
    coords = []

    coords.append([8.33220,-11.86875,0.93355])
    coords.append([0.75523,5.0000,16.37072])
    coords.append([52.61236,5.0000,-5.38580])
    coords.append([2.000,5.0000,3.0000])
    s = Simplex(coords)
    print('volume of simplex',s.volume())

    I=s.inertia(s.centroid())
    print('inertia:')
    print(I)

    import numpy
    I_check= numpy.zeros((3,3))
    I_check[0,0] = 43520.33257
    I_check[1,1] = 194711.28938
    I_check[2,2] = 191168.76173
    I_check[1,2] = -4417.66150
    I_check[2,1] = I_check[1,2]
    I_check[0,2] = 46343.16662
    I_check[2,0] = I_check[0,2]
    I_check[0,1] = -11996.20119
    I_check[1,0] = I_check[0,1]

    print('I_check',I_check)
    print('check',I_check-I)
    print('error',numpy.linalg.norm(I_check-I)/numpy.linalg.norm(I_check))

    print('####### second example #########')


    coords = []

    coords.append([0,0,0])
    coords.append([1,0,0])
    coords.append([0,1,0])
    coords.append([0,0,1])

    s = Simplex(coords)
    print('volume of simplex',s.volume())
    print(s.inertia(s.centroid()))
    # t = Polyeder(coords)
    # #print t.volume()

    # print ('volume', t.volume_by_cm())

    print('####### third example unit cube #########')

    coords.append([0,1,1])
    coords.append([1,0,1])
    coords.append([1,1,0])
    coords.append([1,1,1])

    p = ConvexHull(coords)
    # #print p.volume()

    print ('volume', p.volume())
    print ('volume_divergence_theorem', p.volume_divergence_theorem())
    I,v=p.inertia(p.centroid())
    print('inertia:')
    print(I)

    I_check= numpy.zeros((3,3))
    I_check[0,0] = (1+1)/12.0
    I_check[1,1] = (1+1)/12.0
    I_check[2,2] = (1+1)/12.0
    #print('I_check',I_check)
    #print('check',I_check-I)
    print('error',numpy.linalg.norm(I_check-I)/numpy.linalg.norm(I_check))

    print('####### Fourth example octahedron of side #########')

    side = 1.0
    s = side/math.sqrt(2.0)
    coords=[]
    coords.append([0,s,0])
    coords.append([0,-s,0])
    coords.append([s,0,0])
    coords.append([-s,0,0])
    coords.append([0,0,s])
    coords.append([0,0,-s])

    print(coords)
    p = ConvexHull(coords)
    # #print p.volume()

    print ('volume', p.volume())
    print ('volume_divergence_theorem', p.volume_divergence_theorem())

    centro = p.centroid()
    print( 'centroid', centro)

    I,volume=p.inertia(centro)
    print ('volume', volume)
    print('inertia:')
    print(I)

    I_check= numpy.zeros((3,3))
    I_check[0,0] = volume*(side*side)/10.0
    I_check[1,1] = volume*(side*side)/10.0
    I_check[2,2] = volume*(side*side)/10.0
    print('I_check',I_check)
    print('check',I_check-I)
    print('error',numpy.linalg.norm(I_check-I)/numpy.linalg.norm(I_check))

