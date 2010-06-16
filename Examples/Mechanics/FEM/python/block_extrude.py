print """

FEM simulation using getfem++ and siconos.

"""

import Kernel as Kernel
import numpy as np
import getfem as gf
from matplotlib.pyplot import *
import math

class SiconosFem:
    """ The set of matrices required by Siconos, from a Finite Element Model
    
    """
    def __init__(self):
        self.nbdof = 0
        self.Mass = []
        self.Stiff = []
        self.H = []
        self.q0 = []
        self.RHS=[]

def fillH(pid,sic,nbdof,BOTTOM):
    nb_contacts = np.size(pid)
    rowH = 3 * nb_contacts
    sic.H = np.zeros((rowH,nbdof))
    
    dofbottom = mfu.basic_dof_on_region(BOTTOM)
    for i in range(0,rowH,3):
        sic.H[i,dofbottom[i]+2] = 1.0
        sic.H[i+1,dofbottom[i]] = 1.0
        sic.H[i+2,dofbottom[i]+1] = 1.0
              

def fillH2(pid, sic, nbdof, TOP,alpha,coord):
    nb_contacts = np.size(pid)
    rowH = 3 * nb_contacts
    sic.H2 = np.zeros((rowH,nbdof))
    doftop = mfu.basic_dof_on_region(TOP)
    normH2 = 1.0/sqrt(1+alpha*alpha)
    for i in range(0,rowH,3):
        sic.H2[i,doftop[i]] = norm2*alpha
        sic.H2[i,doftop[i]+1] = norm2
        sic.H2[i+1,doftop[i]+2] = 1.0
        sic.H2[i+2,doftop[i]] = -norm2
        sic.H2[i+2,doftop[i]+1] = alpha*norm2

sico = SiconosFem()

# ===============================
# Model Parameters
# ===============================
E = 2.1e11  # Young modulus
Nu = 0.3 # Poisson coef.
# Lame coeff.
Lambda = E*Nu/((1+Nu)*(1-2*Nu))
Mu = E/(2*(1+Nu))
# Density
Rho=7800
Gravity = -9.81
t0 = 0.0      # start time
T = 0.5    # end time
h = 0.001   # time step
e = 0.0    # restitution coeficient
mu=0.3 # Friction coefficient
theta = 0.5 # theta scheme
with_friction = True

# ===============================
# Build FEM using getfem
# ===============================


xB = 5.0
# ==== The geometry and the mesh ==== 
dimX = 10.01 ; dimY = 10.01 ; dimZ = 3.01
stepX = 1.0 ; stepY = 1.0 ; stepZ = 0.75
alpha = math.atan(dimZ/dimX)
x=np.arange(xB,dimX+xB,stepX)
y=np.arange(0,dimY,stepY)
z=np.arange(0,dimZ,stepZ)
m = gf.Mesh('regular simplices', x,y,z)
m.set('optimize_structure')
# Export the mesh to vtk
m.export_to_vtk("BlockMesh.vtk")

# Create MeshFem objects
# (i.e. assign elements onto the mesh for each variable)
mfu = gf.MeshFem(m,3) # displacement
mff = gf.MeshFem(m,1) # for plot von-mises
# assign the FEM
mfu.set_fem(gf.Fem('FEM_PK(3,1)'))
mff.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(3,1,0.01)'))
# mfu.export_to_vtk("BlockMeshDispl.vtk")

# ==== Set the integration method ==== 
mim = gf.MeshIm(m,gf.Integ('IM_TETRAHEDRON(5)'))

# ==== Summary ==== 
print ' ==================================== \n Mesh details: '
print ' Problem dimension:', mfu.qdim(), '\n Number of elements: ', m.nbcvs(), '\n Number of nodes: ', m.nbpts()
print ' Number of dof: ', mfu.nbdof(), '\n Element type: ', mfu.fem()[0].char()
print ' ===================================='

# ==== Boundaries detection ==== 
allPoints = m.pts()
# Bottom points and faces 
cbot = (abs(allPoints[2,:])  < 1e-6)
pidbot = np.compress(cbot,range(0,m.nbpts()))
fbot = m.faces_from_pid(pidbot)
BOTTOM = 1
m.set_region(BOTTOM,fbot)
# Top points and faces
ctop = (abs(allPoints[2,:]) > dimZ-stepZ)
pidtop = np.compress(ctop,range(0,m.nbpts()))
ftop = m.faces_from_pid(pidtop)
TOP = 2
m.set_region(TOP,ftop)
# Left points and faces
cleft = (abs(allPoints[1,:]) < xB + 1e-6)
pidleft = np.compress(cleft,range(0,m.nbpts()))
fleft= m.faces_from_pid(pidleft)
LEFT = 3
m.set_region(LEFT,fleft)
# Right points and faces
cright = (abs(allPoints[1,:]) > xB + dimX-stepX)
pidright = np.compress(cright,range(0,m.nbpts()))
fright= m.faces_from_pid(pidright)
RIGHT = 4
m.set_region(RIGHT,fright)

# ==== Create getfem models ====
# We use two identical models, one to get the stiffness matrix and the rhs
# and the other to get the mass matrix.
# 
md = gf.Model('real')
# The dof (displacements on nodes)
md.add_fem_variable('u',mfu)
# Add model constants
md.add_initialized_data('lambda',Lambda)
md.add_initialized_data('mu',Mu)
md.add_initialized_data('source_term',[0,0,-10])
md.add_initialized_data('push',[-200000,0,0])
md.add_initialized_data('rho',Rho)
md.add_initialized_data('gravity', Gravity)
md.add_initialized_data('weight',[0,0,Rho*Gravity])
# Build model (linear elasticity)
md.add_isotropic_linearized_elasticity_brick(mim,'u','lambda','mu')
# Add volumic/surfacic source terms
#md.add_source_term_brick(mim,'u','source_term',TOP)
md.add_source_term_brick(mim,'u','push',LEFT)
md.add_source_term_brick(mim,'u','weight')
# Add boundary conditions
#md.add_Dirichlet_condition_with_multipliers(mim,'u',mfu,BOTTOM)

# Assembly
md.assembly()
# Get stiffness matrix
sico.Stiff=md.tangent_matrix().full()
#
# Note: getfem returns sparse matrices. .full() means that we translate sparse to dense.
#
# Get right-hand side
sico.RHS = md.rhs()
# Get initial state 
sico.initial_displacement = md.variable('u')

# Second model for the mass matrix
md2 = gf.Model('real')
md2.add_fem_variable('u',mfu)
md2.add_initialized_data('rho',Rho)
md2.add_mass_brick(mim,'u','rho')
md2.assembly()
# Get mass matrix
sico.Mass = md2.tangent_matrix().full()
# number of dof
sico.nbdof = mfu.nbdof()
sico.mfu=mfu
sico.mesh=m

# ===============================
# Here starts the Siconos stuff
# ===============================
#
# From getfem, we have Mass, Stiffness and RHS
# saved in object sico.

# H-Matrix 
fillH(pidbot,sico,mfu.nbdof(),BOTTOM)

# =======================================
# Create the siconos Dynamical System
# 
# Mass.ddot q + Kq = fExt
# 
# q: dof vector (displacements)
# =======================================
# Initial displacement and velocity
v0 = np.zeros(sico.nbdof)
block = Kernel.LagrangianLinearTIDS(sico.initial_displacement,v0,sico.Mass)
# set fExt and K
block.setFExtPtr(sico.RHS)
block.setKPtr(sico.Stiff)

# =======================================
# The interactions
# A contact is defined for each node at
# the bottom of the block
# =======================================
# Create one relation/interaction for each point 
# in the bottom surface
# Each interaction is of size three with a 
# relation between local coordinates at contact and global coordinates given by: 
# y = Hq + b
# y = [ normal component, first tangent component, second tangent component] 
# 
# The friction-contact non-smooth law
if(with_friction):
    nslaw = Kernel.NewtonImpactFrictionNSL(e,e,mu,3)
    diminter = 3
else:
    nslaw = Kernel.NewtonImpactNSL(e)
    diminter = 1

hh = np.zeros((diminter,sico.nbdof))
b = np.zeros(diminter)
b[0] = 0.0
k = 0
relation=[]
inter=[]
hh = np.zeros((diminter,sico.nbdof))
nbInter = pidbot.shape[0]
if(with_friction):
    for i in range(nbInter):
        # hh is a submatrix of sico.H with 3 rows. 
        hh[:,:] = sico.H[k:k+3,:]
        k += 3
        relation.append(Kernel.LagrangianLinearTIR(hh,b))
        inter.append(Kernel.Interaction(diminter, nslaw, relation[i]))
    
else:
    for i in range(nbInter):
        # hh is a submatrix of sico.H with 1 row. 
        hh[:,:]= sico.H[k,:]
        k += 3
        relation.append(Kernel.LagrangianLinearTIR(hh,b))
        inter.append(Kernel.Interaction(diminter, nslaw, relation[i]))


hh = np.zeros((diminter,sico.nbdof))
doftop = mfu.basic_dof_on_region(TOP)
nb_contacts = np.size(pidtop)
qB = np.repeat([0],3)
qB[0] = xB

for i in range(nb_contacts):
    hh = np.zeros((diminter,sico.nbdof))
    hh[0,doftop[i]] = math.sin(alpha)*math.sin(alpha)
    hh[0,doftop[i]+1] = -math.cos(alpha)*math.sin(alpha)
    hh[1,doftop[i]+2] = 1.0*math.sin(alpha)
    hh[2,doftop[i]] = -math.cos(alpha)*math.sin(alpha)
    hh[2,doftop[i]+1] = -math.sin(alpha)*math.sin(alpha)
    b = np.dot(hh[:,i:i+3],qB)
    relation.append(Kernel.LagrangianLinearTIR(hh,b))
    inter.append(Kernel.Interaction(diminter, nslaw, relation[i]))
    


nbInter=len(inter)

# =======================================
# The Model
# =======================================
blockModel = Kernel.Model(t0,T)

# add the dynamical system to the non smooth dynamical system
blockModel.nonSmoothDynamicalSystem().insertDynamicalSystem(block)

# link the interactions and the dynamical system
for i in range(nbInter):
    blockModel.nonSmoothDynamicalSystem().link(inter[i],block);

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
 #   osnspb.numericsSolverOptions().iparam[0]=100
#    osnspb.numericsSolverOptions().iparam[1]=20
#    osnspb.numericsSolverOptions().iparam[4]=2
#    osnspb.numericsSolverOptions().dparam[0]=1e-6
#    osnspb.numericsSolverOptions().dparam[2]=1e-8
#    osnspb.setMaxSize(1000)
#    osnspb.setMStorageType(1)
#    osnspb.setNumericsVerboseMode(0)
#    osnspb.setKeepLambdaAndYState(true)
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
dataPlot[0, 8] = 0.0

nbNodes = sico.mesh.pts().shape[1]

k = 1
# time loop
while(s.nextTime() < T):
    s.computeOneStep()
    name = 'friction'+str(k)+'.vtk'
    dataPlot[k,0]=s.nextTime()
    dataPlot[k,1]=block.q()[2]    
    dataPlot[k,2]=block.velocity()[2]
    dataPlot[k, 7] = block.q()[20]
    
    # Post proc for paraview
    md.to_variables(block.q())
    VM=md.compute_isotropic_linearized_Von_Mises_or_Tresca('u','lambda','mu',mff)
    dataPlot[k, 8] = VM[0]
    
    #U = fem_model.variable('u')
    sl = gf.Slice(('boundary',),sico.mfu,1)
    sl.export_to_vtk(name, sico.mfu, block.q(),'Displacement', mff, VM, 'Von Mises Stress')

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
title('VM')
plot(dataPlot[:,0], dataPlot[:,8])
grid()
show()

