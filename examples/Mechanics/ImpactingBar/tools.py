import numpy as np
import cmath
pi = cmath.pi


def computeBigH(H,Nfft):

    ## x = cmath.exp(2*1j*pi/Nfft)*np.ones((Nfft),dtype='complex64')
    
    ## sfft = (np.vander(x,Nfft)).T
    ## sfft = np.vander(sfft[:,0],Nfft)

    ## print sfft
    ## print sfft.shape
    bigH = np.zeros(np.asarray(H.shape)*Nfft,dtype='complex64')
    print("shape of bigH:", bigH.shape)

    nc = H.shape[0]
    ndof = H.shape[1]
    
    for i in range(Nfft):
        for j in range(Nfft):
            bigH[i*nc:(i+1)*nc,j*ndof:(j+1)*ndof] = H*cmath.exp(2*1j*pi*i*j/Nfft)
            
    return bigH


def energyNormalisation(mass,U,gamma,Nfft):

    res = -gamma
    for i in range(Nfft):
        res += np.dot((U[i*ndof:(i+1)*ndof]).conjugate.T,np.dot(Mass*U[i*ndof:(i+1)*ndof]))
    
        
