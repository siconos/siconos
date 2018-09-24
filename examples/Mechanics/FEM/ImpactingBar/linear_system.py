

class LinearSystem(object):
    """ Description of a mechanical linear sytem with its mass and rigidity matrices and some boundary conditions. """
    
    def __init__(self,Mass,Stiffness,Bc,fc,fext=None):
        """  """
        self.Mass = Mass
        self.Stiffness = Stiffness
        ## Boundary conditions : Bc U = fc
        self.Bc = Bc
        self.fc = fc
        self.fext = None
        if(fext is not None):
            self.fext = fext
        
    def computeLinearModes(self):
        """ Computation of the eigenfrequencies of a system defined as :
        (K - w^2 M)phi = 0
        """
        # Computes Cholesky decompostion for M = LL^t
        L0 = LA.cholesky(Mass)
        L0inv = LA.inv(L0)
        # Singular value decomposition
        U, s, V = LA.svd(Mass)
        L = np.dot(U,np.diag(np.sqrt(s)))
        Linv = LA.inv(L)
        # System in its "standard" form, Ax = lambda x, A = L^-1 Stiffness L^-t
        system = np.dot(Linv,np.dot(Stiffness,Linv.T))
        system0 =  np.dot(L0inv,np.dot(Stiffness,L0inv.T))
        
        # t0 = time.clock()
        # eigenFreq = LA.eigvals(system)
        # print "elapsed with eigvals ...", time.clock()-t0
        
        # t0 = time.clock()
        # eigenFreq2 = LA.eigvalsh(system)
        # print "elapsed with eigvalsh ...", time.clock()-t0
        
        t0 = time.clock()
        omeg0, phi0 = LA.eigh(system0)
        print("elapsed with eig ...", time.clock()-t0)
        t0 = time.clock()
        omeg2,phi = LA.eigh(system)
        print("elapsed with eigh ...", time.clock()-t0)
        
        # freq00 = np.sqrt(eigenFreq)/(2.*pi)
        # freq01 = np.sqrt(eigenFreq2)/(2.*pi)
        # freq02 = np.sqrt(freq2)/(2.*pi)
        freq = np.sqrt(omeg2)/(2.*pi)
        freq0 = np.sqrt(omeg0)/(2.*pi)
        # print "frequencies : ", freq00.min(),freq01.min(),freq02.min(),freq.min()
        
        # Back to the eigenvectors of the initial system
        vec = np.dot(Linv.T,phi)
        vec0 = np.dot(L0inv.T,phi0)
        
        # Check results ...
        # for i in range(vec.shape[1]):
        #    temp = np.dot(Stiffness-omeg2[i]*Mass,vec[:,i])
        #    if(not np.allclose(LA.norm(temp),0.0)):
        #        print "SVD rate ..."
        #    temp = np.dot(Stiffness-omeg0[i]*Mass,vec0[:,i])
        #    if(not np.allclose(LA.norm(temp),0.0)):
        #        print "Chol rate ..."
        
        return freq,vec

    def computeG(self,i,omega):
        """ computes G(i,omega) = K - (i*omega)^2 Mass """
        
        return (self.Stiffness - (i*omega)**2*self.Mass)
    
    

