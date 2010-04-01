#!/usr/bin/env python

# Siconos-sample, Copyright INRIA 2005-2010.
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
# Siconos is a free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# Siconos is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Siconos; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Contact: Vincent ACARY, siconos-team@lists.gforge.fr
#

from numpy import *

from enthought.mayavi import mlab
from enthought.tvtk.tools import visual

from Siconos.Kernel import *

t0 = 0      # start time
T = 10      # end time
h = 0.005   # time step
r = 0.1     # ball radius
g = 9.81    # gravity
m = 1       # ball mass
e = 0.9     # restitution coeficient
theta = 0.5 # theta scheme



#
# dynamical system
#
x = array([1,0,0]) # initial position
v = array([0,0,0]) # initial velocity
mass = eye(3)      # mass matrix
mass[2,2]=3./5 * r * r

# the dynamical system
ball = LagrangianLinearTIDS(x,v,mass)

# set external forces 
weight = array([-m * g, 0, 0])
ball.setFExtPtr(weight)

#
# Interactions
#

# ball-floor
H = array([[1,0,0]])

nslaw = NewtonImpactNSL(e)
relation = LagrangianLinearTIR(H)
inter = Interaction(1, nslaw, relation)

#
# Model
#
bouncingBall = Model(t0,T)

# add the dynamical system to the non smooth dynamical system
bouncingBall.nonSmoothDynamicalSystem().insertDynamicalSystem(ball)

# link the interaction and the dynamical system
bouncingBall.nonSmoothDynamicalSystem().link(inter,ball);


#
# Simulation
#

# (1) OneStepIntegrators
OSI = Moreau(theta)
OSI.insertDynamicalSystem(ball)

# (2) Time discretisation --
t = TimeDiscretisation(t0,h)

# (3) one step non smooth problem
osnspb = LCP()

# (4) Simulation setup with (1) (2) (3)
s = TimeStepping(t)
s.insertIntegrator(OSI)
s.insertNonSmoothProblem(osnspb)

# end of model definition

#
# computation
#

# simulation initialization
bouncingBall.initialize(s)


#
# Animation
#

ground = visual.box(pos=(-r, 0., 0.), 
                    size=(0,5, 5), 
                    color=visual.color.white )

vball = visual.sphere(radius=r, color=visual.color.red)

def anime():
    s.computeOneStep()
    vball.pos = ball.q()
    s.nextStep()

iter = visual.iterate(10,anime)

# Choose a view angle, and display the figure
v = visual.get_viewer()
v.scene.y_plus_view()
v.scene.camera.position = 1,   1,   10
v.scene.render()

visual.show()

