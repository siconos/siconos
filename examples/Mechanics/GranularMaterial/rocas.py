#!/usr/bin/env python

__all__ = ['create_rocas', 'una_roca']

import numpy, random, math
from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.convexhull import ConvexHull
import sys
def una_roca(io, name, cname, roca_size=0.05, density=1, trans=None, tob=None):
    # Definition of an irregular polyhedron as a convex shape

    rd = [math.pi/2 * random.gauss(0.5,0.2) for _ in range(16)]

    def vert(id1, id2, a, b, c):
        return (a*math.cos(rd[id1])*math.cos(rd[id2]),
                b*math.sin(rd[id1])*math.cos(rd[id2]),
                c*math.sin(rd[id2]))

    vertices = [ vert( 0,  1,   1,  1,  1),
                 vert( 2,  3,   1, -1,  1),
                 vert( 4,  5,  -1,  1,  1),
                 vert( 6,  7,  -1, -1,  1),
                 vert( 8,  9,   1,  1, -1),
                 vert(10, 11,   1, -1, -1),
                 vert(12, 13,  -1,  1, -1),
                 vert(14, 15,  -1, -1, -1) ]

    scale = roca_size / max(numpy.array(vertices).max(axis=0)
                            - numpy.array(vertices).min(axis=0))

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # correction of vertices such that 0 is the centroid
    vertices = (numpy.array(vertices)[:] - cm[:]) * scale

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # Definition of a polyhedron as a convex shape
    io.add_convex_shape(cname, vertices, insideMargin=0.1*roca_size)

    # computation of inertia and volume
    inertia,volume=ch.inertia(ch.centroid())

    # print('geometric inertia:', inertia)
    # print('volume:', volume)
    # print('mass:', volume*density)
    # print('inertia:', inertia*density)


    io.add_object(name,
                 [Contactor(cname)],
                 translation=trans,
                 #velocity=veloci,
                 mass=volume*density,
                 time_of_birth=tob,
                 inertia=inertia*density)


def un_cubo(io, name, cname, roca_size=0.05, density=1, trans=None, tob=None):
    # Definition of a cube as a convex shape

    vertices = numpy.array([[  1,  1,  1],
                            [  1, -1,  1],
                            [ -1,  1,  1],
                            [ -1, -1,  1],
                            [  1,  1, -1],
                            [  1, -1, -1],
                            [ -1,  1, -1],
                            [ -1, -1, -1]]) * roca_size

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # correction of vertices such that 0 is the centroid
    vertices = numpy.array(vertices)[:] - cm[:]

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # Definition of a polyhedron as a convex shape
    io.add_convex_shape(cname, vertices)

    # computation of inertia and volume
    inertia,volume=ch.inertia(ch.centroid())




    io.add_object(name,
                  [Contactor(cname)],
                  translation=trans,
                  #velocity=veloci,
                  mass=volume*density,
                  time_of_birth=tob,
                  inertia=inertia*density)


def create_rocas(io, n_layer=5, n_row=5, n_col=5, x_shift=3.0,
                 roca_size=0.05, top=0, rate=0.01, density=1,
                 distribution = ('uniform', 0.1)):

    N = n_layer*n_row*n_col

    dist, rng = distribution

    if dist == 'uniform':
        sizes = numpy.random.uniform(low = roca_size - rng/2,
                                     high = roca_size + rng/2,
                                     size = N)
    elif dist == 'double':
        sizes = numpy.hstack(
            (numpy.random.normal(scale = rng*0.2,
                                 loc   = roca_size - rng/2,
                                 size  = N/2),
             numpy.random.normal(scale = rng*0.2,
                                 loc   = roca_size + rng/2,
                                 size  = N/2)))
        numpy.random.shuffle(sizes)
        # Gives a rock size distribution with two sizes of rocks, with
        # the mean between both at roca_size
    elif dist == 'exp':
        # In a normal distribution, 10- and 90-percentiles are loc +/- rng*1.28.
        # Let's arrange the 10- and 90-percentiles of the distribution
        # to be in the desired range.
        sizes = numpy.random.exponential(1, N)
        bottom = numpy.percentile(sizes, 10)
        top = numpy.percentile(sizes, 90)
        scale = (rng*1.28) / (top - bottom)
        sizes = (sizes - bottom)*scale + roca_size - rng/2*1.28

    k=0
    print('Creation of the rocks')
    for n in range(n_layer):
        for i in range(n_row):
            for j in range(n_col):
                # initial translation
                if (k%100 == 0):
                    print('.', end='', flush=True)
                trans = [(i-n_row/2.0)*x_shift*roca_size,
                         (j-n_col/2.0)*x_shift*roca_size,
                         top]
                name = 'rock'+str(n)+'_'+str(i)+'_'+str(j)
                cname = 'RockCS'+str(n)+'_'+str(i)+'_'+str(j)
                una_roca(io, name, cname, sizes[k], density, trans,
                         tob = n*rate + random.random()*rate)
                k += 1
