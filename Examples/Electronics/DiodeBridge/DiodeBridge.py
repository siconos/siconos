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
# //
# //  DiodeBridge  : sample of an electrical circuit involving :
# //	- a linear dynamical system consisting of an LC oscillator
# //	- a non smooth system (a 1000 Ohm resistor supplied through a 4 diodes bridge) in parallel
# //	  with the oscillator
# //
# //  Expected behavior : 
# //	The initial state (Vc = 10 V , IL = 0) of the oscillator provides an initial energy.
# //	The period is 2 Pi sqrt(LC) ~ 0,628 ms.
# //      The non smooth system is a full wave rectifier :
# //	each phase (positive and negative) of the oscillation allows current to flow
# //	through the resistor in a constant direction, resulting in an energy loss :
# //	the oscillation damps.
# //
# //  State variables : 
# //	- the voltage across the capacitor (or inductor)
# //	- the current through the inductor
# //
# //  Since there is only one dynamical system, the interaction is defined by :
# //	- complementarity laws between diodes current and voltage. Depending on
# //        the diode position in the bridge, y stands for the reverse voltage across the diode 
# //	  or for the diode current (see figure in the template file) 
# //	- a linear time invariant relation between the state variables and
# //	  y and lambda (derived from Kirchhoff laws)
# //
# //-----------------------------------------------------------------------


t0 = 0.0
T = 5.0e-3       # Total simulation time
h_step = 1.0e-6  # Time step
Lvalue = 1e-2  # inductance
Cvalue = 1e-6   # capacitance
Rvalue = 1e3    # resistance 
Vinit = 10.0    # initial voltage
Modeltitle = "DiodeBridge"

from matplotlib.pyplot import subplot, title, plot, grid, show

from Siconos.Kernel import FirstOrderLinearDS, FirstOrderLinearTIR, ComplementarityConditionNSL, Interaction, Model, Moreau, TimeDiscretisation, LCP, TimeStepping

from numpy import array, eye, empty
import numpy as np
#
# dynamical system
#
init_state = array([Vinit,0]) 
print init_state

A = array([[0,-1.0/Cvalue],
           [1.0/Lvalue,0]], order='FORTRAN')

LSDiodeBridge=FirstOrderLinearDS(init_state, A)

#
# Interactions
#

C = array([[0,0],[0,0],[-1.0,0],[1.0,0]],order='FORTRAN')
D = array([[1.0/Rvalue,1.0/Rvalue,-1,0],[1.0/Rvalue,1.0/Rvalue,0,-1],[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0]],order='FORTRAN')
B = array([[0.0,0.0,-1.0/Cvalue,1.0/Cvalue],[0.0,0.0,0.0,0.0]],order='FORTRAN')

# C = np.transpose(C)
# D = np.transpose(D)
# B = np.transpose(B)
print A
print B
print C
print D

LTIRDiodeBridge=FirstOrderLinearTIR(C,B)
LTIRDiodeBridge.setDPtr(D)

nslaw=ComplementarityConditionNSL(4)
InterDiodeBridge=Interaction(4, nslaw,LTIRDiodeBridge,1)
InterDiodeBridge.insert(LSDiodeBridge)

#
# Model
#

DiodeBridge=Model(t0,T,Modeltitle)
#   add the dynamical system in the non smooth dynamical system
DiodeBridge.nonSmoothDynamicalSystem().insertDynamicalSystem(LSDiodeBridge)
#   link the interaction and the dynamical system
DiodeBridge.nonSmoothDynamicalSystem().link(InterDiodeBridge,LSDiodeBridge)


#
# Simulation
#

# (1) OneStepIntegrators
theta = 0.5
aOSI = Moreau(LSDiodeBridge,theta)
 
# (2) Time discretisation
aTiDisc = TimeDiscretisation(t0,h_step)

# (3) Non smooth problem
aLCP = LCP()

# (4) Simulation setup with (1) (2) (3)
aTS = TimeStepping(aTiDisc,aOSI,aLCP)

# end of model definition

#
# computation
#

# simulation initialization
DiodeBridge.initialize(aTS)

k = 0
h = aTS.timeStep();
# Number of time steps
N = (T-t0)/h


# Get the values to be plotted 
# ->saved in a matrix dataPlot

dataPlot = empty((N+1,7))

x = LSDiodeBridge.x()
y = InterDiodeBridge.y(0)
lambda_ = InterDiodeBridge.lambda_(0)
k=0
# For the initial time step: 
#  time
dataPlot[0, 0] = t0
dataPlot[k, 0] = aTS.nextTime()
#  inductor voltage
dataPlot[k, 1] = LSDiodeBridge.x()[0]
# inductor current
dataPlot[k, 2] = LSDiodeBridge.x()[1]
# diode R1 current
    #print InterDiodeBridge.y(0)[0]
    #print y[0]
dataPlot[k, 3] = InterDiodeBridge.y(0)[0]
# diode R1 voltage
dataPlot[k, 4] = - InterDiodeBridge.lambda_(0)[0]
# diode F2 voltage 
dataPlot[k, 5] = - InterDiodeBridge.lambda_(0)[1]
# diode F1 current
dataPlot[k, 6] = - InterDiodeBridge.lambda_(0)[2]

print dataPlot

# time loop
k = 1
#while(aTS.nextTime() < T):
while (k < N+1):
    
    aTS.computeOneStep()
    #aLCP.display()
    dataPlot[k, 0] = aTS.nextTime()
    #  inductor voltage
    dataPlot[k, 1] = LSDiodeBridge.x()[0]
    # inductor current
    dataPlot[k, 2] = LSDiodeBridge.x()[1]
    # diode R1 current
    #print InterDiodeBridge.y(0)[0]
    #print y[0]
    dataPlot[k, 3] = InterDiodeBridge.y(0)[0]
    # diode R1 voltage
    dataPlot[k, 4] = - InterDiodeBridge.lambda_(0)[0]
    # diode F2 voltage 
    dataPlot[k, 5] = - InterDiodeBridge.lambda_(0)[1]
    # diode F1 current
    dataPlot[k, 6] = - InterDiodeBridge.lambda_(0)[2]

    k += 1
    aTS.nextStep()
    print aTS.nextTime()
    
#
# plots
#
subplot(411)
title('inductor voltage')
plot(dataPlot[:,0], dataPlot[:,1])
grid()
subplot(412)
title('inductor current')
plot(dataPlot[:,0], dataPlot[:,2])
grid()
subplot(413)
title('diode R1 current')
plot(dataPlot[:,0], dataPlot[:,3])
grid()
subplot(414)
title('diode R1 voltage')
plot(dataPlot[:,0], dataPlot[:,4])
grid()
show()

