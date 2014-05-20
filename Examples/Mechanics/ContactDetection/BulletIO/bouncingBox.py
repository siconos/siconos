#!/usr/bin/env python

from Siconos.Mechanics.ContactDetection import Contactor
from Siconos.Mechanics.ContactDetection.Bullet import IO

with IO.Hdf5(io_filename='bouncingBox.hdf5', mode='w') as io:

    io.insertConvexShape('Cube', [
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (-1.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
        (1.0, 1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, -1.0, 1.0)])

    #    io.insertPrimitiveShape('Cube', 'Box', (.5, .5, .5))

    io.insertPrimitiveShape('Ground', 'Box', (30, 30, .5))

    io.insertNewtonImpactFrictionNSL('contact', mu=0.7)

    io.insertObject('ground', [Contactor('Ground')],
                    position=[0, 0, 0])

    io.insertObject('cube', [Contactor('Cube')], position=[0, 0, 2],
                    velocity=[10, 0, 0, 1, 1, 1],
                    mass=1)


with IO.Hdf5(io_filename='bouncingBox.hdf5', mode='r+') as io:
    
    io.run()
