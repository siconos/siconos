#!/usr/bin/env python

# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2016 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#

from siconos.kernel import NewtonEulerDS, NewtonImpactNSL,\
     NewtonEulerR, NewtonEulerFrom1DLocalFrameR, Interaction,\
     MoreauJeanOSI, TimeDiscretisation, LCP, TimeStepping,\
     changeFrameAbsToBody,changeFrameBodyToAbs,\
     rotationVectorFromQuaternion, quaternionFromRotationVector,\
     rotateAbsToBody,\
     SiconosVector, NonSmoothDynamicalSystem


import numpy as np

import math
        
t0 = 0.0     # start time
h = 0.001   # time step
N= 10000
T = h*N
theta = 0.5  # theta scheme

class UnstableRotation(NewtonEulerDS):

    def __init__(self,x, v):
        I = np.zeros((3, 3))
        I[0, 0] = 5.0
        I[1, 1] = 10.0
        I[2, 2] = 1.0
        m=1.0
        NewtonEulerDS.__init__(self,x, v, m, I)
        # Allocation of _MExt
        self.setMExtPtr(SiconosVector(3))
        # specify that MExt is expressed in the inertial frame.
        self.setIsMextExpressedInInertialFrame(True)

    def computeMExt(self, time, mExt=None):
        td = 2.0 - h
        if mExt is None :
            mExt = self._mExt

        if isinstance(mExt,SiconosVector):
            mExt.zero()
            if (0 <= time < td):
                mExt.setValue(0, 20.0)
            elif (td <= time <= td + h):
                mExt.setValue(1, 1.0 / (5.0 * h))
        else:
            mExt[:]=0
            if (0 <= time < td):
                mExt[0] =20.
            elif (td <= time <= td + h):
                mExt[1]= 1.0 / (5.0 * h)
        

#
# dynamical system
#
x = [0, 0, 0, 1.0, 0, 0, 0]  # initial configuration
v = [0, 0, 0, 0, 0, 0]  # initial velocity
unstableRotation = UnstableRotation(x, v)





class HeavyTop(NewtonEulerDS):

    def __init__(self,x, v):
        I = np.zeros((3, 3))
        I[0, 0] = 5.0
        I[1, 1] = 5.0
        I[2, 2] = 1.0
        m=1.0
        NewtonEulerDS.__init__(self,x, v, m, I)
        self._Mg=20
        self._l=1.0
        self.setComputeJacobianMIntqByFD(True)
        # Allocation of _mInt
        self._mInt = SiconosVector(3)


    def centermass(self,q):
        r= np.zeros(3)
        E3 = SiconosVector(3)
        E3.zero()
        E3.setValue(2,1.0)
        rotateAbsToBody(q,E3)
        r[0] = E3.getValue(0)
        r[1] = E3.getValue(1)
        r[2] = E3.getValue(2)
        return r
        
    def computeMInt(self, time, q, v, mInt=None):
        if mInt is None :
            mInt = self._mInt
        if isinstance(mInt,SiconosVector):
            r = self.centermass(q)
            m =  self._Mg*self._l*np.cross(r,[0,0,1.0])
            mInt.setValue(0,m[0])
            mInt.setValue(1,m[1])
            mInt.setValue(2,m[2])
            changeFrameAbsToBody(q,mInt)
            #print("mInt========")
            mInt.display()
        else:
            r = self.centermass(q)
            m =  self._Mg*self._l*np.cross(r,[0,0,1.0])
            m_sv = SiconosVector(m)
            changeFrameAbsToBody(q,m_sv)
            m_sv.display()
            mInt[0] = m_sv.getValue(0) 
            mInt[1] = m_sv.getValue(1) 
            mInt[2] = m_sv.getValue(2) 
            print("mInt", mInt)



rotationVector_init= SiconosVector(3)
rotationVector_init.zero()
rotationVector_init.setValue(0,0.3)
x=SiconosVector(7)
quaternionFromRotationVector(rotationVector_init,x)

#x = [0, 0, 0, 1.0, 0, 0, 0]  # initial configuration
v = [0, 0, 0, 0, 0, 50]  # initial velocity
heavytop = HeavyTop(x, v)



ds = unstableRotation

#ds = heavytop

ds.display()

# test swig director
# ds.computeMInt(1,x,v)
# ds._mInt.display()
# m=SiconosVector(3)
# ds.computeMInt(1,x,v,m)
# m.display()
# m=np.zeros(3)
# ds.computeMInt(1,x,v,m)
# print m
# raw_input()


# Non-Smooth Dynamical System
#
nsds = NonSmoothDynamicalSystem(t0, T)

# add the dynamical system to the non smooth dynamical system
nsds.insertDynamicalSystem(ds)

#
# Simulation
#

# (1) OneStepIntegrators
OSI = MoreauJeanOSI(theta)

# (2) Time discretisation --
t = TimeDiscretisation(t0, h)

# (3) one step non smooth problem
osnspb = LCP()

# (4) Simulation setup with (1) (2) (3)
s = TimeStepping(nsds, t, OSI, osnspb)
#s.setDisplayNewtonConvergence(True)
s.setNewtonTolerance(1e-10)
#s.setNewtonMaxIteration(1)

# end of model definition

#
# computation
#

# Get the values to be plotted
# ->saved in a matrix dataPlot
dataPlot = np.empty((N+1, 26))

#
# numpy pointers on dense Siconos vectors
#
q = ds.q()
v = ds.twist()
p = ds.p(1)

#
# initial data
#
k=0
dataPlot[k, 1] = q[0]
dataPlot[k, 2] = q[1]
dataPlot[k, 3] = q[2]
dataPlot[k, 4] = q[3]
dataPlot[k, 5] = q[4]
dataPlot[k, 6] = q[5]
dataPlot[k, 7] = q[6]

dataPlot[k, 8] = v[0]
dataPlot[k, 9] = v[1]
dataPlot[k, 10] = v[2]
dataPlot[k, 11] = v[3]
dataPlot[k, 12] = v[4]
dataPlot[k, 13] = v[5]

omega = v[3:6]
print("omega", omega)
angular_momentum = np.dot(ds.inertia(),omega)
am= SiconosVector(angular_momentum)
changeFrameBodyToAbs(q,am)

dataPlot[k, 14] = am.getValue(0)
dataPlot[k, 15] = am.getValue(1)
dataPlot[k, 16] = am.getValue(2)
dataPlot[k, 17] = am.norm2()

rotationVector = SiconosVector(3)
rotationVectorFromQuaternion(q[3],q[4],q[5],q[6], rotationVector)
dataPlot[k, 18] = rotationVector.getValue(0)
dataPlot[k, 19] = rotationVector.getValue(1)
dataPlot[k, 20] = rotationVector.getValue(2)


dataPlot[k, 22] = h* omega[0]
dataPlot[k, 23] = h* omega[1]
dataPlot[k, 24] = h* omega[2]
dataPlot[k, 25] = np.linalg.norm(h*omega)




k = 1

# time loop
while(s.hasNextEvent() and k < N):
    print(' ' )
    print (
        '------- k = ',
        k,
        '-----------------------------------------')
    print(' ' )
    s.computeOneStep()
    dataPlot[k, 0] = s.nextTime()    
    dataPlot[k, 1] = q[0]
    dataPlot[k, 2] = q[1]
    dataPlot[k, 3] = q[2]
    dataPlot[k, 4] = q[3]
    dataPlot[k, 5] = q[4]
    dataPlot[k, 6] = q[5]
    dataPlot[k, 7] = q[6]

    dataPlot[k, 8] = v[0]
    dataPlot[k, 9] = v[1]
    dataPlot[k, 10] = v[2]
    dataPlot[k, 11] = v[3]
    dataPlot[k, 12] = v[4]
    dataPlot[k, 13] = v[5]

    omega = v[3:6]
    angular_momentum = np.dot(ds.inertia(),omega)
    am= SiconosVector(angular_momentum)
    changeFrameBodyToAbs(q,am)
    a = np.zeros(1)
    a[0] = am.getValue(0)
    #a[1] = am.getValue(1)
    # print "omega", omega
    # print "angular_momentum", angular_momentum,
    # print "q=", q
    # print " norm(a[1:2])", np.linalg.norm(a) 
    #raw_input()
    dataPlot[k, 14] = am.getValue(0)
    dataPlot[k, 15] = am.getValue(1)
    dataPlot[k, 16] = am.getValue(2)
    dataPlot[k, 17] = am.norm2()

    
    rotationVector = SiconosVector(3)
    rotationVectorFromQuaternion(q[3],q[4],q[5],q[6], rotationVector)
    dataPlot[k, 18] = rotationVector.getValue(0)
    dataPlot[k, 19] = rotationVector.getValue(1)
    dataPlot[k, 20] = rotationVector.getValue(2)
    
    
    dataPlot[k, 22] = h* omega[0]
    dataPlot[k, 23] = h* omega[1]
    dataPlot[k, 24] = h* omega[2]
    dataPlot[k, 25] = np.linalg.norm(h*omega)

    

    
    
    k = k + 1
    s.nextStep()




dataPlot=np.resize(dataPlot,(k-2,26))


np.savetxt("result-py.dat", dataPlot)
#
# comparison with the reference file
#
from siconos.kernel import SimpleMatrix, getMatrix

#
# plots
#
from matplotlib.pyplot import subplot, title, plot, grid, show, figure


figure(num='Moreau Jean Siconos', figsize=(12, 12))
subplot(321)
title('angular velocities Omega')
plot(dataPlot[:, 0], dataPlot[:, 11])
plot(dataPlot[:, 0], dataPlot[:, 12])
#plot(dataPlot[:, 0], dataPlot[:, 13])

subplot(322)
title('rotation vector')
plot(dataPlot[:, 0], dataPlot[:, 18])
plot(dataPlot[:, 0], dataPlot[:, 19])
plot(dataPlot[:, 0], dataPlot[:, 20])

subplot(323)
title('Theta (h Omega)')
plot(dataPlot[:, 0], dataPlot[:, 22])
plot(dataPlot[:, 0], dataPlot[:, 23])
plot(dataPlot[:, 0], dataPlot[:, 24])

subplot(325)
title('norm of Theta')
plot(dataPlot[:, 0], dataPlot[:, 25])

subplot(324)
title('angular momentum (pi[0])')
plot(dataPlot[:, 0], dataPlot[:, 14])
#plot(dataPlot[:, 0], dataPlot[:, 15])
#plot(dataPlot[:, 0], dataPlot[:, 16])

subplot(326)
title('norm of angular momentum  pi')
plot(dataPlot[:, 0], dataPlot[:, 17])

grid()
show()
