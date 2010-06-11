# Getfem tools to build the finite element model
import getfem as gf
import numpy as np


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

def import_fem(sico):
    """ Build a mesh object using getfem.
    
    We use getfem++ to build the finite element model and
    to fill in the operators required by siconos:
    - the mass matrix (Mass)
    - the stiffness matrix (Stiff)
    - the matrix that links global coordinates and local coord. at contact points (H)
    
    """
    ############################
    # The geometry and the mesh
    ############################
    dimX = 10.01 ; dimY = 10.01 ; dimZ = 3.01
    stepX = 1.0 ; stepY = 1.0 ; stepZ = 1.0
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
    mfd = gf.MeshFem(m,1) # data
    mff = gf.MeshFem(m,1) # for plot von-mises
    # assign the FEM
    mfu.set_fem(gf.Fem('FEM_PK(3,1)'))
    mfd.set_fem(gf.Fem('FEM_PK(3,0)'))
    mff.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(3,1,0.01)'))
    # mfu.export_to_vtk("BlockMeshDispl.vtk")

    # Set the integration method
    mim = gf.MeshIm(m,gf.Integ('IM_TETRAHEDRON(5)'))

    # Summary
    print ' ==================================== \n Mesh details: '
    print ' Problem dimension:', mfu.qdim(), '\n Number of elements: ', m.nbcvs(), '\n Number of nodes: ', m.nbpts()
    print ' Number of dof: ', mfu.nbdof(), '\n Element type: ', mfu.fem()[0].char()
    print ' ===================================='

    ###########################
    # Set the parameters
    # for the constitutive law
    ###########################
    E = 1e3  # Young modulus
    Nu = 0.3 # Poisson coef.
    # Lame coeff.
    Lambda = E*Nu/((1+Nu)*(1-2*Nu))
    Mu = E/(2*(1+Nu))
    # Density
    Rho=1.0#7.800
    Gravity = -9.81
    ############################
    # Boundaries detection
    ############################
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
    # Top-Left points and faces
    cleft = (abs(allPoints[1,:]) < 1e-6)
    clefttop=cleft*ctop
    pidlefttop = np.compress(clefttop,range(0,m.nbpts()))
    flefttop = m.faces_from_pid(pidlefttop)
    pidleft = np.compress(cleft,range(0,m.nbpts()))
    fleft= m.faces_from_pid(pidleft)
    LEFTTOP = 3
    m.set_region(LEFTTOP,flefttop)
    LEFT = 4
    m.set_region(LEFT,fleft)

    # Create a model
    md = gf.Model('real')
    md.add_fem_variable('u',mfu)
    md.add_initialized_data('lambda',Lambda)
    md.add_initialized_data('mu',Mu)
    md.add_initialized_data('source_term',[0,0,-100])
    md.add_initialized_data('push',[0,100,0])
    md.add_initialized_data('rho',Rho)
    md.add_initialized_data('gravity', Gravity)
#    Weight = np.zeros(mfu.nbdof())
##    Weight = []
    md.add_initialized_data('weight',[0,0,Rho*Gravity])
    md.add_isotropic_linearized_elasticity_brick(mim,'u','lambda','mu')
    #md.add_source_term_brick(mim,'u','source_term',TOP)
    #md.add_source_term_brick(mim,'u','push',LEFT)
    md.add_source_term_brick(mim,'u','weight')
    #md.add_Dirichlet_condition_with_multipliers(mim,'u',mfu,BOTTOM)
    
    md.assembly()
    sico.Stiff=md.tangent_matrix()
    sico.RHS = md.rhs()
    sico.q0 = md.variable('u')
    md2 = gf.Model('real')
    md2.add_fem_variable('u',mfu)
    md2.add_initialized_data('rho',Rho)
    md2.add_mass_brick(mim,'u','rho')
    md2.assembly()
    sico.Mass = md2.tangent_matrix()
    sico.nbdof = mfu.nbdof()
    sico.mfu=mfu
    sico.mesh=m
    sico.pos = np.zeros(sico.nbdof)

    sico.pos[0:sico.nbdof:3] = m.pts()[0,:]
    sico.pos[1:sico.nbdof:3] = m.pts()[1,:]
    sico.pos[2:sico.nbdof:3] = m.pts()[2,:]
    sico.K0 = np.dot(sico.Stiff.full(),sico.pos)
    sico.bot = pidbot
    # running solve...
    #md.solve()
    
    # post-processing
    #VM=md.compute_isotropic_linearized_Von_Mises_or_Tresca('u','lambda','mu',mff)
    # extracted solution
    #U = md.variable('u')
    # export U and VM in a pos file
    #sl = gf.Slice(('boundary',),mfu,1)
    #sl.export_to_vtk('toto.vtk', mfu, U, 'Displacement', mff, VM, 'Von Mises Stress')
    
    # H-Matrix 
    fillH(pidbot,sico,mfu.nbdof())

    return md

def import_fem2(sico):
    """ Build a mesh object using getfem.
    
    We use getfem++ to build the finite element model and
    to fill in the operators required by siconos:
    - the mass matrix (Mass)
    - the stiffness matrix (Stiff)
    - the matrix that links global coordinates and local coord. at contact points (H)
    
    """
    ############################
    # The geometry and the mesh
    ############################
    dimX = 10.01 ; dimY = 10.01 ; dimZ = 10.01
    stepX = 1.0 ; stepY = 1.0 ; stepZ = 1.0
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
    mfd = gf.MeshFem(m,1) # data
    mff = gf.MeshFem(m,1) # for plot von-mises
    # assign the FEM
    mfu.set_fem(gf.Fem('FEM_PK(3,1)'))
    mfd.set_fem(gf.Fem('FEM_PK(3,0)'))
    mff.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(3,1,0.01)'))
    # mfu.export_to_vtk("BlockMeshDispl.vtk")

    # Set the integration method
    mim = gf.MeshIm(m,gf.Integ('IM_TETRAHEDRON(5)'))

    # Summary
    print ' ==================================== \n Mesh details: '
    print ' Problem dimension:', mfu.qdim(), '\n Number of elements: ', m.nbcvs(), '\n Number of nodes: ', m.nbpts()
    print ' Number of dof: ', mfu.nbdof(), '\n Element type: ', mfu.fem()[0].char()
    print ' ===================================='

    ###########################
    # Set the parameters
    # for the constitutive law
    ###########################
    E = 1e3  # Young modulus
    Nu = 0.3 # Poisson coef.
    # Lame coeff.
    Lambda = E*Nu/((1+Nu)*(1-2*Nu))
    Mu = E/(2*(1+Nu))
    # Density
    Rho=7800
    ############################
    # Boundaries detection
    ############################
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
    # Top-Left points and faces
    cleft = (abs(allPoints[1,:]) < 1e-6)
    clefttop=cleft*ctop
    pidlefttop = np.compress(clefttop,range(0,m.nbpts()))
    flefttop = m.faces_from_pid(pidlefttop)
    pidleft = np.compress(cleft,range(0,m.nbpts()))
    fleft= m.faces_from_pid(pidleft)
    LEFTTOP = 3
    m.set_region(LEFTTOP,flefttop)
    LEFT = 4
    m.set_region(LEFT,fleft)
    
    ############################
    # Assembly
    ############################
    nbd =  mfd.nbdof()
    # Stiffness matrix
    sico.Stiff=gf.asm_linear_elasticity(mim, mfu, mfd,np.repeat([Lambda], nbd), np.repeat([Mu], nbd))
    # Mass matrix
    sico.Mass=Rho*gf.asm_mass_matrix(mim, mfu)
    # Right-hand side
    Ftop = gf.asm_boundary_source(TOP, mim, mfu, mfd, np.repeat([[0],[0],[-1]],nbd,1))
    Fleft= gf.asm_boundary_source(LEFT, mim, mfu, mfd, np.repeat([[0],[10],[0]],nbd,1))
    sico.RHS = Ftop + Fleft

    sico.nbdof = mfu.nbdof()
    sico.q0 = mfu.basic_dof_from_cvid()
    
    sico.bot = pidbot
    

    # H-Matrix 
    fillH(pidbot,sico,mfu.nbdof())
    return m

def fillH(pid,sic,nbdof):
    nb_contacts = np.size(pid)
    sic.H = np.zeros((3*nb_contacts,nbdof))
    num = 0 # Current contact number
    for i in range(0,sic.H.shape[0]-1,3):
        sic.H[i,3*pid[num]+2] = 1
        sic.H[i+1,3*pid[num]] = 1
        sic.H[i+2,3*pid[num]+1] = 1
        num = num +1
