from Siconos.Numerics import *
from scipy.sparse.csr import csr_matrix

M = csr_matrix([1,2,3])
    
MM = SparseMatrix(M)    

printSparse(M)
SBM=newFromFileSBM('SBM1.dat')

A = SBMtoSparse(SBM)

R = sparseToSBM(3,A)
