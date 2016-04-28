print("""

FEM simulation using getfem++ and siconos.

""")

import siconos.kernel as kernel
import numpy as np
import getfem as gf
from matplotlib.pyplot import *

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
              

sico = SiconosFem()

# ===============================
# Model Parameters
# ===============================
E = 1e3  # Young modulus
Nu = 0.3 # Poisson coef.
# Lame coeff.
Lambda = E*Nu/((1+Nu)*(1-2*Nu))
Mu = E/(2*(1+Nu))
# Density
Rho=1.0#7.800
Gravity = -9.81
t0 = 0.0      # start time
T = 10.0    # end time
h = 0.005   # time step
e = 0.0   # restitution coeficient
mu=0.3 # Friction coefficient
theta = 0.5 # theta scheme

# ===============================
# Build FEM using getfem
# ===============================

# ==== The geometry and the mesh ==== 
dimX = 10.01 ; dimY = 10.01 ; dimZ = 3.01
stepX = 2.0 ; stepY = 2.0 ; stepZ = 1.5
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
cleft = (abs(allPoints[1,:]) < 1e-6)
pidleft = np.compress(cleft,list(range(0,m.nbpts())))
fleft= m.faces_from_pid(pidleft)
LEFT = 3
m.set_region(LEFT,fleft)

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
md.add_initialized_data('source_term',[0,0,-100])
md.add_initialized_data('push',[0,100,0])
md.add_initialized_data('rho',Rho)
md.add_initialized_data('gravity', Gravity)
md.add_initialized_data('weight',[0,0,Rho*Gravity])
# Build model (linear elasticity)
md.add_isotropic_linearized_elasticity_brick(mim,'u','lambda','mu')
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
block = kernel.LagrangianLinearTIDS(sico.initial_displacement,v0,sico.Mass)
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
# The non-smooth law
nslaw = kernel.NewtonImpactNSL(e)
diminter = 1

hh = np.zeros((diminter,sico.nbdof))
b = np.zeros(diminter)
b[0] = 1.0
k = 0
relation=[]
inter=[]
hh = np.zeros((diminter,sico.nbdof))
nbInter = pidbot.shape[0]
for i in range(nbInter):
    # hh is a submatrix of sico.H with 1 row. 
    hh[:,:]= sico.H[k,:]
    k += 3
    relation.append(kernel.LagrangianLinearTIR(hh,b))
    inter.append(kernel.Interaction(diminter, nslaw, relation[i]))
    
nbInter=len(inter)

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

osnspb = kernel.LCP()

# (4) Simulation setup with (1) (2) (3)
s = kernel.TimeStepping(t)
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

nbNodes = sico.mesh.pts().shape[1]

k = 1
# time loop
while(s.hasNextEvent()):
    s.computeOneStep()
    name = 'bounce'+str(k)+'.vtk'
    dataPlot[k,0]=s.nextTime()
    dataPlot[k,1]=block.q()[2]    
    dataPlot[k,2]=block.velocity()[2]
    
    # Post proc for paraview
    md.to_variables(block.q())
    VM=md.compute_isotropic_linearized_Von_Mises_or_Tresca('u','lambda','mu',mff)
    #U = fem_model.variable('u')
    sl = gf.Slice(('boundary',),sico.mfu,1)
    sl.export_to_vtk(name, sico.mfu, block.q(),'Displacement', mff, VM, 'Von Mises Stress')
    #print s.nextTime()
    k += 1
    s.nextStep()



#subplot(211)
#title('position')
#plot(dataPlot[:,0], dataPlot[:,1])

#grid()
#subplot(212)
#title('velocity')
#plot(dataPlot[:,0], dataPlot[:,2])
#grid()
#show()

