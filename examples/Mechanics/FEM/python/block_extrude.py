print("""

FEM simulation using getfem++ and siconos.

""")

import siconos.kernel as kernel
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
Rho=100.0
Gravity = -9.81
t0 = 0.0      # start time
T = 0.5    # end time
h = 0.0005   # time step
e = 0.0    # restitution coeficient
mu=0.3 # Friction coefficient
theta = 0.5 # theta scheme
with_friction = True

# ===============================
# Build FEM using getfem
# ===============================


# ==== The geometry and the mesh ==== 
dimX = 10.01 ; dimY = 10.01 ; dimZ = 3.01
stepX = 2.0 ; stepY = 2.0 ; stepZ = 1.5
alpha = math.atan(dimZ/30.0)
x=np.arange(0,dimX,stepX)
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
print(' ==================================== \n Mesh details: ')
print(' Problem dimension:', mfu.qdim(), '\n Number of elements: ', m.nbcvs(), '\n Number of nodes: ', m.nbpts())
print(' Number of dof: ', mfu.nbdof(), '\n Element type: ', mfu.fem()[0].char())
print(' ====================================')

# ==== Boundaries detection ==== 
allPoints = m.pts()
# Bottom points and faces 
cbot = (abs(allPoints[2,:])  < 1e-6)
pidbot = np.compress(cbot,list(range(0,m.nbpts())))
fbot = m.faces_from_pid(pidbot)
BOTTOM = 1
m.set_region(BOTTOM,fbot)
# Top points and faces
ctop = (abs(allPoints[2,:]) > dimZ-stepZ)
pidtop = np.compress(ctop,list(range(0,m.nbpts())))
ftop = m.faces_from_pid(pidtop)
TOP = 2
m.set_region(TOP,ftop)
# Left points and faces
cleft = (abs(allPoints[0,:]) <  1e-6)
pidleft = np.compress(cleft,list(range(0,m.nbpts())))
fleft= m.faces_from_pid(pidleft)
LEFT = 3
m.set_region(LEFT,fleft)
# Right points and faces
cright = (abs(allPoints[0,:]) > dimX-stepX)
pidright = np.compress(cright,list(range(0,m.nbpts())))
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
md.add_initialized_data('source_term',[0.0,0.0,-10])
md.add_initialized_data('push',[50000000.0,0.0,0.0])
md.add_initialized_data('rho',Rho)
md.add_initialized_data('gravity', Gravity)
md.add_initialized_data('weight',[0,0,Rho*Gravity])
# Build model (linear elasticity)
md.add_isotropic_linearized_elasticity_brick(mim,'u','lambda','mu')
# Add volumic/surfacic source terms
#md.add_source_term_brick(mim,'u','source_term',TOP)
md.add_source_term_brick(mim,'u','push',LEFT)
md.add_source_term_brick(mim,'u','weight')

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

# ===============================
# Here starts the Siconos stuff
# ===============================
#
# From getfem, we have Mass, Stiffness and RHS
# saved in object sico.

# =======================================
# Create the siconos Dynamical System
# 
# Mass.ddot q + Kq = fExt
# 
# q: dof vector (displacements)
# =======================================
# Initial displacement and velocity
v0 = np.zeros(sico.nbdof)
block = kernel.LagrangianLinearTIDS(sico.initial_displacement,v0,sico.Mass)
# set fExt and K
block.setFExtPtr(sico.RHS)
block.setKPtr(sico.Stiff)

print("Start Siconos process ...")

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

nslaw = kernel.NewtonImpactFrictionNSL(e,e,mu,3)
diminter = 3
b = np.zeros(diminter)
relation=[]
inter=[]
q0 = mfu.basic_dof_nodes() 
# We create an interaction for each point at the bottom of the block
nbInter = pidbot.shape[0]
# The list of dof in bottom face
dofbottom = mfu.basic_dof_on_region(BOTTOM)
for i in range(nbInter):
    index = dofbottom[i]
    hh = np.zeros((diminter,sico.nbdof))
    hh[0,index+2] = 1.0
    hh[1,index] = 1.0
    hh[2,index+1] = 1.0
    qRef =np.repeat([0],3)
    qRef[0] = q0[0,index]
    qRef[1] = q0[1,index]
    qRef[2] = q0[2,index]
    b[0] = np.dot(hh[0,index:index+3],qRef)
    relation.append(kernel.LagrangianLinearTIR(hh,b))
    inter.append(kernel.Interaction(diminter, nslaw, relation[i]))
    

# Now we treat the top face in the same way
relationTop =[]
qRef = np.repeat([0],3)

nbInterTop = pidtop.shape[0]
doftop = mfu.basic_dof_on_region(TOP)
for i in range(nbInterTop):
    index = doftop[i]
    hh = np.zeros((diminter,sico.nbdof))
    hh[0,index] = -math.sin(alpha)
    hh[0,index+2] = -math.cos(alpha)
    hh[1,index] = math.cos(alpha)
    hh[1,index+2] = -math.sin(alpha)
    hh[2,index+1] = -1.0
    qRef = np.repeat([0],3)
    qRef[0] = dimX - q0[0,index]
    b = np.dot(hh[:,index:index+3],qRef)
    relationTop.append(kernel.LagrangianLinearTIR(hh,b))
    inter.append(kernel.Interaction(diminter, nslaw, relationTop[i]))

print(nbInterTop)
nbInter=len(inter)

print(nbInter)

# =======================================
# The Model
# =======================================
blockModel = kernel.Model(t0,T)

# add the dynamical system to the non smooth dynamical system
blockModel.nonSmoothDynamicalSystem().insertDynamicalSystem(block)

# link the interactions and the dynamical system
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
osnspb = kernel.FrictionContact(3)

# (4) Simulation setup with (1) (2) (3)
s = kernel.TimeStepping(t)
s.insertIntegrator(OSI)
s.insertNonSmoothProblem(osnspb)
blockModel.setSimulation(s)

# simulation initialization
blockModel.initialize()

# the number of time steps
N = (T-t0)/h

# Get the values to be plotted 
# ->saved in a matrix dataPlot

dataPlot = np.empty((N+1,3))

dataPlot[0, 0] = t0
dataPlot[0, 1] = block.q()[0]
nbNodes = m.pts().shape[1]

k = 1
# time loop
while(s.hasNextEvent()):
    s.computeOneStep()
    name = 'extr'+str(k)+'.vtk'
    #dataPlot[k,0]=s.nextTime()
    #dataPlot[k,1]=block.q()[0]    
    
    # Post proc for paraview
    md.to_variables(block.q())
    VM=md.compute_isotropic_linearized_Von_Mises_or_Tresca('u','lambda','mu',mff)
    
    #U = fem_model.variable('u')
    sl = gf.Slice(('boundary',),mfu,1)
    sl.export_to_vtk(name, mfu, block.q(),'Displacement', mff, VM, 'Von Mises Stress')

    #print s.nextTime()
    k += 1
    s.nextStep()



#subplot(211)
#title('position')
#plot(dataPlot[:,0], dataPlot[:,1])
#plot(dataPlot[:,0], dataPlot[:,2])
#grid()
#subplot(212)
#title('VM')
#show()

