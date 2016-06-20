print("""

FEM simulation using getfem++ and siconos.

""")

import siconos.kernel as kernel
import getfemtosiconos as gts
import numpy as np
import getfem as gf
from matplotlib.pyplot import *

t0 = 0.0      # start time
T = 10.0    # end time
h = 0.005   # time step
g = 9.81    # gravity
e = 0.1    # restitution coeficient
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

block = kernel.LagrangianLinearTIDS(sico.q0,v0,sico.Mass.full())
F = sico.RHS
block.setFExtPtr(F)
block.setKPtr(sico.Stiff.full())

position_init = sico.pos

# =======================================
# The interaction
# =======================================

dist = 0.5
diminter = 1
b = np.repeat([dist], diminter)
nslaw = kernel.NewtonImpactNSL(e)
k=0

pos0 = np.dot(sico.H,sico.pos)

dimH = sico.H.shape[0]
relation=[]
inter=[]
hh = np.zeros((diminter,sico.nbdof))
for i in range(0,dimH,3):
    hh[0,:] = sico.H[i,:]
    b2 = b + pos0[i]
    relation.append(kernel.LagrangianLinearTIR(hh,b2))
    inter.append(kernel.Interaction(diminter, nslaw, relation[k]))
    k+=1


    
nbInter=len(inter)


# =======================================
# The Model
# =======================================
blockModel = kernel.Model(t0,T)

# add the dynamical system to the non smooth dynamical system
blockModel.nonSmoothDynamicalSystem().insertDynamicalSystem(block)

# link the interaction and the dynamical system
for i in range(nbInter):
    blockModel.nonSmoothDynamicalSystem().link(inter[i],block);

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

dataPlot = np.empty((N+1,9))

dataPlot[0, 0] = t0
dataPlot[0, 1] = block.q()[2]
#dataPlot[0, 2] = block.velocity()[2]
dataPlot[0, 2] = block.q()[5]
dataPlot[0, 3] = block.q()[8]
dataPlot[0, 4] = block.q()[11]
dataPlot[0, 5] = block.q()[14]
dataPlot[0, 6] = block.q()[17]
dataPlot[0, 7] = block.q()[20]
dataPlot[0, 8] = block.q()[23]

nbNodes = sico.mesh.pts().shape[1]

k = 1
# time loop
while(s.hasNextEvent()):
    s.computeOneStep()
    name = 'titi'+str(k)+'.vtk'
    dataPlot[k,0]=s.nextTime()
    #dataPlot[k,1]=block.q()[2]    
    #dataPlot[k,2]=block.velocity()[2]
    #   dataPlot[k, 2] = block.q()[5]
    #dataPlot[k, 3] = block.q()[8]
    #dataPlot[k, 4] = block.q()[11]
    #dataPlot[k, 5] = block.q()[14]
    #dataPlot[k, 6] = block.q()[17]
    #dataPlot[k, 7] = block.q()[20]
    #dataPlot[k, 8] = block.q()[23]
    siz = sico.pos.size
    reference_pos = sico.pos
    current_pos = block.q() +  reference_pos
    dataPlot[k,1]=current_pos[2]    
    #dataPlot[k,2]=block.velocity()[2]
    dataPlot[k, 2] = current_pos[5]
    dataPlot[k, 3] = current_pos[8]
    dataPlot[k, 4] = current_pos[11]
    dataPlot[k, 5] = current_pos[14]
    dataPlot[k, 6] = current_pos[17]
    dataPlot[k, 7] = current_pos[20]
    dataPlot[k, 8] = current_pos[23]

    #bottom_pos = current_pos[0:12]
    #new_ref = reference_pos
    #translatX = reference_pos[0] - current_pos[0]
    #translatY = reference_pos[1] - current_pos[1]
    #translatZ = reference_pos[2] - current_pos[2]
    #trans = np.repeat([translatX],nbNodes)
    #new_ref[0:siz:3] = reference_pos[0:siz:3] - trans
    #trans = np.repeat([translatY],nbNodes)
    #new_ref[1:siz:3] = reference_pos[1:siz:3] - trans
    #trans = np.repeat([translatZ],nbNodes)
    #new_ref[2:siz:3] = reference_pos[2:siz:3] - trans
    
    #depl = current_pos - new_ref
    #block.setQPtr(depl)
    #correction = np.dot(sico.H,new_ref)
    #bnew = correction[0:dimH:3]
    #for i in range(nbInter):
    #    bb = b + bnew[i]
    #    relation[i].setEPtr(bb)

    #print "RESUME"
    #relation[0].computeOutput(s.nextTime(),0)
    #inter[0].y(0).display()
    #print current_pos[0:3]
    #print new_ref[0:3]
    #print depl[0:3]
    
    

    fem_model.to_variables(block.q())
    U = fem_model.variable('u')
    #POS = current_pos
    sl = gf.Slice(('boundary',),sico.mfu,1)
    #sl.export_to_vtk(name, sico.mfu, U,'Displacement', sico.mfu, POS, 'Position')
    sl.export_to_vtk(name, sico.mfu, U,'Displacement')
    #print s.nextTime()
    k += 1
    s.nextStep()



subplot(211)
title('position')
plot(dataPlot[:,0], dataPlot[:,1])
plot(dataPlot[:,0], dataPlot[:,2])
plot(dataPlot[:,0], dataPlot[:,3])
plot(dataPlot[:,0], dataPlot[:,4])
plot(dataPlot[:,0], dataPlot[:,5])
plot(dataPlot[:,0], dataPlot[:,6])
plot(dataPlot[:,0], dataPlot[:,7])
plot(dataPlot[:,0], dataPlot[:,8])

grid()
subplot(212)
title('velocity')
plot(dataPlot[:,0], dataPlot[:,2])
grid()
show()
