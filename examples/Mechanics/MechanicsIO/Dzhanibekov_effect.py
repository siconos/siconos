#!/usr/bin/env python

# The Dzhanibekov effect
# http://mathoverflow.net/questions/81960/the-dzhanibekov-effect-an-exercise-in-mechanics-or-fiction-explain-mathemat

import siconos
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5

options = siconos.mechanics.collision.bullet.SiconosBulletOptions()
options.worldScale = 0.01

with Hdf5() as io:

    # Definition of a tetrahedron as a convex shape
    io.addPrimitiveShape('Body1', 'Cylinder', (1, 6))
    io.addPrimitiveShape('Body2', 'Box', (2, .4, 11))
    
    io.addObject('roo', [Contactor('Body1'),
                         Contactor('Body2',
                                   relative_translation=[0, 3, 0])],
                 translation=[0, 0, 4],
                 # a small perturbation on z axis
                 velocity=[0, 0, 0, 0, 2, 0.0001],
                 mass=1,
                 inertia=[1, 10, 11])

with Hdf5(mode='r+') as io:

    io.run(with_timer=True,
           options=options,
           t0=0,
           T=90,
           h=0.005,
           Newton_max_iter=1,
           set_external_forces=lambda x: None)
