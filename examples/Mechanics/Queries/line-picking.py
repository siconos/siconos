#!/usr/bin/env python

#
# Example of how to "pick" objects, i.e., identifying objects based on
# a line intersection test.  Can be used to implement mouse
# interaction.  Here we just detect objects that pass through the
# vertical line at (x,y)=(0,0).
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.io.mechanics_io
from siconos.mechanics.collision.bullet import SiconosBulletOptions

options = SiconosBulletOptions()
options.worldScale = 1.0
options.contactBreakingThreshold = 0.04

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    # Definition of a sphere
    io.addPrimitiveShape('Sphere', 'Sphere', (0.5,))

    # Definition of the ground shape
    io.addPrimitiveShape('Ground', 'Box', (50, 10, 0.1))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.addNewtonImpactFrictionNSL('contact', mu=0.5, e=0.1)

    # Define a series of 10 spheres with an initial velocity to throw
    # them towards 0,0 so that we can detect them as they roll by.
    for i in range(10):
        io.addObject('sphere%d'%i, [Contactor('Sphere')],
                     translation=[i*2, 0, 2],
                     velocity=[-5, 0, 0, 0, 0, 0],
                     mass=1)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.addObject('ground', [Contactor('Ground')],
                 translation=[0, 0, -0.1])

class LinePicker():
    def initialize(self, io):
        self.io = io
        self.nsds = io._model.nonSmoothDynamicalSystem()
        self.manager = io._broadphase

    def step(self):
        query = self.manager.lineIntersectionQuery([0,0,10],[0,0,-10],False)
        names = [self.nsds.name(q.body) for q in query if q.body is not None]
        if len(names) > 0:
            print('time %f, picked %s'%(self.io.currentTime(), ', '.join(names)))

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(options=options,
           t0=0,
           T=10,
           h=0.005,
           theta=0.50001,
           controller = LinePicker(),
           Newton_max_iter=2,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=1000,
           tolerance=1e-5,
           verbose=False, verbose_progress=False)
