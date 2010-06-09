print """
Programme test pour le "couplage" entre Siconos et Getfem
 
      req: python_getfem.py et siconos ...    
"""

import Kernel
import getfemtosiconos as gts
import numpy as np

t0 = 0.0      # start time
T = 5.0      # end time
h = 0.005   # time step
r = 0.1     # ball radius
g = 9.81    # gravity
m = 1       # ball mass
e = 0.9     # restitution coeficient
theta = 0.5 # theta scheme


sic1 = gts.SiconosFem()
sic2 = gts.SiconosFem()
m = gts.import_fem2(sic1)
m2 = gts.import_fem(sic2)


# =======================================
# Create the siconos Dynamical System
# =======================================
# Initial position and velocity
v0 = np.zeros(sic2.nbdof)

block = Kernel.LagrangianLinearTIDS(sic2.q0,v0,sic2.Mass.full())
block.setFExtPtr(sic2.RHS)
block.setKPtr(sic2.Stiff.full())

# =======================================
# The interaction
# =======================================
nslaw = Kernel.NewtonImpactNSL(e)
relation = Kernel.LagrangianLinearTIR(sic2.H)
inter = Kernel.Interaction(1, nslaw, relation)

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
osnspb = Kernel.LCP()

# (4) Simulation setup with (1) (2) (3)
s = Kernel.TimeStepping(t)
s.insertIntegrator(OSI)
s.insertNonSmoothProblem(osnspb)

# simulation initialization
blockModel.initialize(s)


