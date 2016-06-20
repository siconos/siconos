print("""

FEM simulation using getfem++ and siconos.

""")

import siconos.kernel as kernel
import getfemtosiconos as gts
import numpy as np
import getfem as gf
from matplotlib.pyplot import *

t0 = 0.0      # start time
T = 0.1      # end time
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

block = kernel.LagrangianLinearTIDS(sico.pos,v0,sico.Mass.full())
F = sico.RHS + sico.K0
block.setFExtPtr(F)
block.setKPtr(sico.Stiff.full())

# =======================================
# The interaction
# =======================================

dist = 3.0
dimH = sico.H.shape[0]
    
if(with_friction):
    diminter = dimH
    nslaw = kernel.NewtonImpactFrictionNSL(e,e,mu,3)
    relation = kernel.LagrangianLinearTIR(sico.H)
else:
    diminter = dimH/3
    Hnofric = np.zeros((diminter, sico.nbdof))
    b = np.repeat([dist], diminter)
    Hnofric = sico.H[0:dimH-1:3,:]
    nslaw = kernel.NewtonImpactNSL(e)
    relation = kernel.LagrangianLinearTIR(Hnofric,b)

inter = kernel.Interaction(diminter, nslaw, relation)

# =======================================
# The Model
# =======================================
blockModel = kernel.Model(t0,T)

# add the dynamical system to the non smooth dynamical system
blockModel.nonSmoothDynamicalSystem().insertDynamicalSystem(block)

# link the interaction and the dynamical system
blockModel.nonSmoothDynamicalSystem().link(inter,block);

# =======================================
# The Simulation
# =======================================

# (1) OneStepIntegrators
OSI = kernel.Moreau(theta)

# (2) Time discretisation --
t = kernel.TimeDiscretisation(t0,h)

# (3) one step non smooth problem
if(with_friction):
    osnspb = kernel.FrictionContact(3)
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
    osnspb = kernel.LCP()

# (4) Simulation setup with (1) (2) (3)
s = kernel.TimeStepping(t)
s.insertIntegrator(OSI)
s.insertNonSmoothProblem(osnspb)

# simulation initialization
blockModel.initialize()

# the number of time steps
N = (T-t0)/h

# Get the values to be plotted 
# ->saved in a matrix dataPlot

dataPlot = np.empty((N+1,4))

dataPlot[0, 0] = t0
dataPlot[0, 1] = block.q()[2]
dataPlot[0, 2] = block.velocity()[2]
dataPlot[0, 3] = block.q()[5]

k = 1

# time loop
while(s.hasNextEvent()):
    s.computeOneStep()
    name = 'titi'+str(k)+'.vtk'
    dataPlot[k,0]=s.nextTime()
    dataPlot[k,1]=block.q()[2]    
    dataPlot[k,2]=block.velocity()[2]
    dataPlot[k,3]=block.q()[5]    
    k += 1
    s.nextStep()
    fem_model.to_variables(block.q())
    U = fem_model.variable('u')
    sl = gf.Slice(('boundary',),sico.mfu,1)
    sl.export_to_vtk(name, sico.mfu, U,'Displacement')
    print(s.nextTime())
    

subplot(211)
title('position')
plot(dataPlot[:,0], dataPlot[:,1])
plot(dataPlot[:,0], dataPlot[:,3])

grid()
subplot(212)
title('velocity')
plot(dataPlot[:,0], dataPlot[:,2])
grid()
show()
