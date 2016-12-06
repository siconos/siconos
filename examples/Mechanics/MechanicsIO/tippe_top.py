#!/usr/bin/env python

# A tippe-top with Coulomb friction only & JeanMoreau time stepping.

import siconos
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import thetav, Hdf5
from math import pi
from matplotlib import pyplot as plt

siconos.io.mechanics_io.set_implementation('original')

# Cf A set-valued force law for spatial Coulomb-Contensou friction
# R.I. Leine, Ch. Glocker, European Journal of Mechanics, 2003

# Note: Here there is no Contensou friction (we are in the case where
# R --> 0) and with a Jean-Moreau timestepping method the time when the
# reverse equilibrium is reached depends on the time step.

# cf the conclusion of Leine & Glocker about the introduction of a
# higher order scheme: "The need for such a higher-order integration
# method became apparent during the numerical analysis of the
# Tippe-Top. The Tippe-Top experiences almost no damping and exhibits
# high-frequency oscillation in the nutation. Numerical simulation
# with the classical (low-order) time-stepping method yielded fastly
# diverging solutions, or a fast deadening of the nutational
# oscillation if integrated with a fully implicit version of the
# classical time-stepping method. Reducing the step-size led to
# different results as the divergence or deadening was weakened."

mu = 0.3
m = 6*1e-3     # kg
r1 = 1.5       # cm
r2 = 0.5       # cm
a1 = 0.3       # cm
a2 = 1.6       # cm
I1 = 8*1e-3    # kg . cm^2
I3 = 7*1e-3    # kg . cm^2

with Hdf5() as io:

    io.addPrimitiveShape('Body1', 'Sphere', (r1,))
    io.addPrimitiveShape('Body2', 'Cylinder', (r2, a2))
    io.addPrimitiveShape('Body3', 'Sphere', (r2,))
    io.addPrimitiveShape('Ground', 'Box', (100, 100, .5))

    io.addNewtonImpactFrictionNSL('contact', mu=mu)
    io.addObject('ground', [Contactor('Ground')], translation=[0, 0, 0])

    io.addObject('tippe-top', [Contactor('Body1',
                                         relative_translation=[0, 0, a1]),
                               Contactor('Body2',
                                         relative_orientation=([1, 0, 0],
                                                               pi/2),
                                         relative_translation=[0, 0, a2/2.]),
                               Contactor('Body3',
                                         relative_orientation=([1, 0, 0],
                                                               pi/2),
                                         relative_translation=[0, 0, a2])],
                 # we need to avoid contact at first step, so we let the top
                 # fall. This is not what is done in Leine & Glocker.
                 translation=[0, 0, r1-a1 + r1-a1],
                 orientation=([0, 1, 0], 0.1),
                 velocity=[0, 0, 0, .0, .0, 180],
                 mass=m)

with Hdf5(mode='r+') as io:

    io.run(with_timer=True,
           multipoints_iterations=False,
           gravity_scale=1./100.,
           t0=0,
           T=20,
           h=0.0001,
           Newton_max_iter=20)

    # plot of theta Euler angle.
    # to be compared with fig 12
    p = io.dynamic_data()[0:5000, :]

    plt.plot(p[:, 0], thetav(p[:, 5], p[:, 6], p[:, 7], p[:, 8]))
    plt.show()
