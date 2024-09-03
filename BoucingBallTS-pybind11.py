#!/usr/bin/env python

# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2021 INRIA.
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

import numpy as np
from numpy.linalg import norm
from siconos.algebra import MatrixWrapper, VectorWrapper, MyVectorClass, MyVectorWrapper, MyEigenClass, MyEigenWrapper
from siconos.modeling import LagrangianLinearTIDS, NewtonImpactNSL,\
    LagrangianLinearTIR, Interaction, NonSmoothDynamicalSystem
from siconos.integrators import MoreauJeanOSI
from siconos.simulation import TimeDiscretisation, TimeStepping
from siconos.nonsmooth_formulations import LCP


t0 = 0.       # start time
T = 10.       # end time
h = 0.005    # time step
r = 0.1      # ball radius
g = 9.81     # gravity
m = 1.        # ball mass
e = 0.9      # restitution coeficient
theta = 0.5  # theta scheme

#
# dynamical system
#
initial_position_np = np.array([1, 0, 0], dtype=np.float64, order='F')
initial_velocity_np = np.array([0, 0, 0], dtype=np.float64, order='F')
mass_np = np.eye(10000, dtype=np.float64, order='F')
mass_np[...] = 2. / 5 * r * r

a = 2

# A = np.array([1., 2., 3.])
# initial_position_vc = MyVectorClass(A)


# initial_position_ec = MyEigenClass()
# A = np.array([4., 5., 6.])
# initial_position_ec.contents = A

# A[0] = 21.

# initial_position_vc.print_v()




# initial_position_ew = MyEigenWrapper(initial_position_ec)
# initial_position_ew.print_vector()
# initial_position_ew.change_vector()
# initial_position_ew.print_vector()
# print(A)
# A[1] = 19.
# print(initial_position_ec.contents)


# initial_position_eigen = VectorWrapper(initial_position_np)
# initial_position_eigen.get_vector()
# initial_velocity_eigen = VectorWrapper(initial_velocity_np)
# mass_eigen = MatrixWrapper(mass_np)

# initial_position_eigen_vec = initial_position_eigen.get_vector()
# initial_position_eigen.print_vector()
# initial_position_eigen.get_vector()


# the dynamical system
# ball = LagrangianLinearTIDS(initial_position_eigen, initial_velocity_eigen, mass_eigen)

ball = LagrangianLinearTIDS(initial_position_np, initial_velocity_np, mass_np)

a = 3
# mass_np = ball.mass_python()
mass_np[2, 2] = 27.
print(ball)
# print(ball.mass_python()[2,2])

# a = 1

# # set external forces
# weight_np = np.array([-m * g, 0, 0], dtype=np.float64, order='F')
# weight_eigen = VectorWrapper(weight_np)
# ball.setFExtPtr(weight_eigen)

# #
# # Interactions
# #

# # ball-floor
# H_np = np.array([[1, 0, 0]], dtype=np.float64, order='F')
# H_eigen = MatrixWrapper(H_np)

# nslaw = NewtonImpactNSL(e)
# relation = LagrangianLinearTIR(H_eigen)
# inter = Interaction(nslaw, relation)
# # print(relation)

# #
# # Model
# #
# bouncingBall = NonSmoothDynamicalSystem(t0, T)

# # add the dynamical system to the non smooth dynamical system
# bouncingBall.insertDynamicalSystem(ball)

# # link the interaction and the dynamical system
# bouncingBall.link(inter, ball)

# # print(bouncingBall)

# # Simulation

# (1) OneStepIntegrators
# print(" =====")

OSI = MoreauJeanOSI(theta)
OSI.tonche_mass(ball)
print(ball)
print(" =====")
# OSI.tonche_mass(ball)
print(mass_np[2, 2])


      
mm = ball.mass2()

print("dsldksldslds", mm[2,2])

a = 4


mm[3,4] = 127

print(ball.mass2()[3,4])


# # (2) Time discretisation --
# t = TimeDiscretisation(t0, h)

# # (3) one step non smooth problem
# osnspb = LCP()

# # (4) Simulation setup with (1) (2) (3)
# s = TimeStepping(bouncingBall,t, OSI, osnspb)

# # print(OSI)

# # end of model definition

# #
# # computation
# #


# # the number of time steps
# N = int((T - t0) / h)

# # Get the values to be plotted
# # ->saved in a matrix dataPlot

# dataPlot = np.zeros((N+1, 5))


# # numpy pointers on dense Siconos vectors

# q = ball.q_wrapper().get_vector()
# v = ball.velocity_wrapper().get_vector()
# p = ball.p_wrapper(1).get_vector()
# lambda_ = inter.lambda_wrapper(1).get_vector()



# # initial data

# dataPlot[0, 0] = t0
# dataPlot[0, 1] = q[0]
# dataPlot[0, 2] = v[0]
# dataPlot[0, 3] = p[0]
# dataPlot[0, 4] = lambda_[0]

# k = 1

# # time loop
# while s.hasNextEvent():
#     s.computeOneStep()

#     # q = ball.q_wrapper().get_vector()
#     # v = ball.velocity_wrapper().get_vector()
#     # p = ball.p_wrapper(1).get_vector()
#     # lambda_ = inter.lambda_wrapper(1).get_vector()

#     # print(id(q))

#     dataPlot[k, 0] = s.nextTime()
#     dataPlot[k, 1] = q[0]
#     dataPlot[k, 2] = v[0]
#     dataPlot[k, 3] = p[0]
#     dataPlot[k, 4] = lambda_[0]

#     k += 1
#     s.nextStep()

# #
# # comparison with the reference file
# #
# # ref = getMatrix(SimpleMatrix("BouncingBallTS.ref"))

# # if (norm(dataPlot - ref) > 1e-12):
# #     print("Warning. The result is rather different from the reference file.")


# #
# # plots
# #
# import matplotlib,os
# havedisplay = "DISPLAY" in os.environ
# if not havedisplay:
#     matplotlib.use('Agg')

# import matplotlib.pyplot as plt
# plt.subplot(411)
# plt.title('position')
# plt.plot(dataPlot[:, 0], dataPlot[:, 1])
# plt.grid()
# plt.subplot(412)
# plt.title('velocity')
# plt.plot(dataPlot[:, 0], dataPlot[:, 2])
# plt.grid()
# plt.subplot(413)
# plt.plot(dataPlot[:, 0], dataPlot[:, 3])
# plt.title('reaction')
# plt.grid()
# plt.subplot(414)
# plt.plot(dataPlot[:, 0], dataPlot[:, 4])
# plt.title('lambda')
# plt.grid()

# if havedisplay:
#     plt.show()
# # else:
#     # plt.savefig("bbts.png")
