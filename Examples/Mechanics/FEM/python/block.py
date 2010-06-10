print """

FEM simulation using getfem++ and siconos.

"""

import Kernel
import getfemtosiconos as gts
import numpy as np
import getfem as gf
from matplotlib.pyplot import *

t0 = 0.0      # start time
T = 5.0      # end time
h = 0.005   # time step
g = 9.81    # gravity
e = 0.9     # restitution coeficient
mu=0.3 # Friction coefficient
theta = 0.5 # theta scheme

with_friction = False

sico = gts.SiconosFem()
fem_model = gts.import_fem(sico)


# =======================================
# Create the siconos Dynamical System
# =======================================
# Initial position and velocity
v0 = np.zeros(sico.nbdof)

block = Kernel.LagrangianLinearTIDS(sico.q0,v0,sico.Mass.full())
block.setFExtPtr(sico.RHS)
block.setKPtr(sico.Stiff.full())

# =======================================
# The interaction
# =======================================
dim0 = sico.H.shape[0]

if(with_friction):
    diminter = dim0
    nslaw = Kernel.NewtonImpactFrictionNSL(e,e,mu,3)
    relation = Kernel.LagrangianLinearTIR(sico.H)
else:
    Hnofric = sico.H[0:dim0-1:3,:]
    diminter = Hnofric.shape[0]
    dist = 3.0
    b = np.repeat([dist], diminter)
    nslaw = Kernel.NewtonImpactNSL(e)
    relation = Kernel.LagrangianLinearTIR(Hnofric,b)

inter = Kernel.Interaction(diminter, nslaw, relation)

# =======================================
# The Model
# =======================================
blockModel = Kernel.Model(t0,T)

# add the dynamical system to the non smooth dynamical system
blockModel.nonSmoothDynamicalSystem().insertDynamicalSystem(block)

# link the interaction and the dynamical system
blockModel.nonSmoothDynamicalSystem().link(inter,block);

# =======================================
# The Simulation
# =======================================

# (1) OneStepIntegrators
OSI = Kernel.Moreau(theta)
OSI.insertDynamicalSystem(block)

# (2) Time discretisation --
t = Kernel.TimeDiscretisation(t0,h)

# (3) one step non smooth problem
if(with_friction):
    osnspb = Kernel.FrictionContact(3)
    osnspb.numericsSolverOptions().iparam[0]=100
    osnspb.numericsSolverOptions().iparam[1]=20
    osnspb.numericsSolverOptions().iparam[4]=2
    osnspb.numericsSolverOptions().dparam[0]=1e-6
    osnspb.numericsSolverOptions().dparam[2]=1e-8
    osnspb.setMaxSize(1000)
    osnspb.setMStorageType(1)
    osnspb.setNumericsVerboseMode(0)
    osnspb.setKeepLambdaAndYState(true)
else:
    osnspb = Kernel.LCP()

# (4) Simulation setup with (1) (2) (3)
s = Kernel.TimeStepping(t)
s.insertIntegrator(OSI)
s.insertNonSmoothProblem(osnspb)

# simulation initialization
blockModel.initialize(s)

# the number of time steps
N = (T-t0)/h

# Get the values to be plotted 
# ->saved in a matrix dataPlot

dataPlot = np.empty((N+1,3))

dataPlot[0, 0] = t0
dataPlot[0, 1] = block.q()[3]
dataPlot[0, 2] = block.velocity()[3]

k = 1

# time loop
while(s.nextTime() < T):
    s.computeOneStep()
    name = 'titi'+str(k)+'.vtk'
    dataPlot[k,0]=s.nextTime()
    dataPlot[k,1]=block.q()[3]    
    dataPlot[k,2]=block.velocity()[3]
    k += 1
    s.nextStep()
    fem_model.to_variables(block.q())
    U = fem_model.variable('u')
    sl = gf.Slice(('boundary',),sico.mfu,1)
    sl.export_to_vtk(name, sico.mfu, U,'Displacement')
    print s.nextTime()
    

subplot(211)
title('position')
plot(dataPlot[:,0], dataPlot[:,1])
grid()
subplot(212)
title('velocity')
plot(dataPlot[:,0], dataPlot[:,2])
grid()
show()
