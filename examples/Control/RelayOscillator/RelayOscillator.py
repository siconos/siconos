#!/usr/bin/env python

# /* Siconos-sample , Copyright INRIA 2005-2011.
#  * Siconos is a program dedicated to modeling, simulation and control
#  * of non smooth dynamical systems.
#  * Siconos is a free software; you can redistribute it and/or modify
#  * it under the terms of the GNU General Public License as published by
#  * the Free Software Foundation; either version 2 of the License, or
#  * (at your option) any later version.
#  * Siconos is distributed in the hope that it will be useful,
#  * but WITHOUT ANY WARRANTY; without even the implied warranty of
#  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  * GNU General Public License for more details.
#  *
#  * You should have received a copy of the GNU General Public License
#  * along with Siconos; if not, write to the Free Software
#  * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#  *
#  * Contact: Vincent ACARY vincent.acary@inrialpes.fr
# */
# //-----------------------------------------------------------------------
# //-----------------------------------------------------------------------

import math

t0 = 0.0
T = 100   # Total simulation time
h_step = 1.0e-2  # Time step
xinit = 3.0*math.sqrt(2.0)/(2.0*math.pi)    # initial voltage
Modeltitle = "RelayOscillator"

withPlot=True
if withPlot:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.pyplot import subplot, title, plot, grid, savefig

from siconos.kernel import FirstOrderLinearDS, FirstOrderLinearTIR, \
                           RelayNSL, Interaction,\
                           Model, EulerMoreauOSI, TimeDiscretisation, Relay,  \
                           TimeStepping,NonSmoothDynamicalSystem

#
# dynamical system
#
init_state = [0.0,xinit,0]

A = [[0,          1.0,     0.0],
     [0.0,        0.0,     1.0],
     [0.0,       -3.0,    -2.0]]

LSRelayOscillator=FirstOrderLinearDS(init_state, A)

#
# Interactions
#

C = [[1.0, 0.0,   0.0]]

D = [[0.0 ]]

B = [[0.0],
     [0.0],
     [1.0]]

LTIRRelayOscillator=FirstOrderLinearTIR(C,B)
LTIRRelayOscillator.setDPtr(D)

nslaw=RelayNSL(1)
InterRelayOscillator=Interaction(1, nslaw,LTIRRelayOscillator,1)


#
# Model
#
RelayOscillator=Model(t0,T,Modeltitle)

#   add the dynamical system in the non smooth dynamical system
myNSDS = NonSmoothDynamicalSystem()
myNSDS.insertDynamicalSystem(LSRelayOscillator)


#   link the interaction and the dynamical system
myNSDS.link(InterRelayOscillator,LSRelayOscillator)

RelayOscillator.setNonSmoothDynamicalSystemPtr(myNSDS)


#
# Simulation
#

# (1) OneStepIntegrators
theta = 0.5
aOSI = EulerMoreauOSI(theta)
# (2) Time discretisation
aTiDisc = TimeDiscretisation(t0,h_step)

# (3) Non smooth problem
aRelay = Relay()

# (4) Simulation setup with (1) (2) (3)
aTS = TimeStepping(aTiDisc,aOSI,aRelay)

# end of model definition

#
# computation
#

# simulation initialization
RelayOscillator.initialize(aTS)

k = 0
h = aTS.timeStep();
print("Timestep : ",h)
# Number of time steps
N = (T-t0)/h
print("Number of steps : ",N)

# Get the values to be plotted
# ->saved in a matrix dataPlot

from numpy import empty
dataPlot = empty([N+1,8])

x = LSRelayOscillator.x()
print("Initial state : ",x)
y = InterRelayOscillator.y(0)
print("First y : ",y)
lambda_ = InterRelayOscillator.lambda_(0)

while (k < N):
    aTS.computeOneStep()
    #aLCP.display()
    dataPlot[k, 0] = aTS.nextTime()
    #  inductor voltage
    dataPlot[k, 1] = x[0]
    dataPlot[k, 2] = x[1]
    dataPlot[k, 3] = x[2]
    dataPlot[k, 4] = y[0]
    dataPlot[k, 5] = lambda_[0]


    k += 1
    if k%1000==0:
        print("step =", k, " < ", N)
    aTS.nextStep()

if (withPlot) :
    #
    # plots
    #
    subplot(511)
    title('x1')
    plot(dataPlot[0:k-1,0], dataPlot[0:k-1,1])
    grid()
    subplot(512)
    title('x2')
    plot(dataPlot[0:k-1,0], dataPlot[0:k-1,2])
    grid()
    subplot(513)
    title('x3')
    plot(dataPlot[0:k-1,0], dataPlot[0:k-1,3])
    subplot(514)
    title('y')
    plot(dataPlot[0:k-1,0], dataPlot[0:k-1,4])
    grid()
    subplot(515)
    title('lambda')
    plot(dataPlot[0:k-1,0], dataPlot[0:k-1,5])
    grid()
    savefig("relay_oscillator.png")

