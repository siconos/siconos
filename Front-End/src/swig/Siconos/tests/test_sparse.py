# sparseToSBM
# scipy csr => 1x1 block
from numpy import finfo, double
eps = finfo(double).eps


def test_from_csr1():

    from Siconos.Numerics import sparseToSBM, getValueSBM
    from scipy.sparse.csr import csr_matrix
    
    M = csr_matrix([[1,2,3],
                    [4,5,6],
                    [7,8,9]])

    print M.indices
    print M.indptr
    print M.data

    r,SBM = sparseToSBM(3,M)

    assert getValueSBM(SBM,0,0) == 1
    assert getValueSBM(SBM,0,1) == 2
    assert getValueSBM(SBM,0,2) == 3
    assert getValueSBM(SBM,1,0) == 4
    assert getValueSBM(SBM,1,1) == 5
    assert getValueSBM(SBM,1,2) == 6    
    assert getValueSBM(SBM,2,0) == 7
    assert getValueSBM(SBM,2,1) == 8
    assert getValueSBM(SBM,2,2) == 9    

# scipy csr 3x3 block
def test_from_csr2():

    from Siconos.Numerics import sparseToSBM, getValueSBM
    from scipy.sparse.csr import csr_matrix
    
    M = csr_matrix([[1,2,3],
                    [4,5,6],
                    [7,8,9]])

    print M.indices
    print M.indptr
    print M.data

    r,SBM = sparseToSBM(1,M)

    assert getValueSBM(SBM,0,0) == 1
    assert getValueSBM(SBM,0,1) == 2
    assert getValueSBM(SBM,0,2) == 3
    assert getValueSBM(SBM,1,0) == 4
    assert getValueSBM(SBM,1,1) == 5
    assert getValueSBM(SBM,1,2) == 6    
    assert getValueSBM(SBM,2,0) == 7
    assert getValueSBM(SBM,2,1) == 8
    assert getValueSBM(SBM,2,2) == 9    


def test_SBMtoSparse1():
    from Siconos.Numerics import getValueSBM, newFromFileSBM, printSBM, SBMtoSparse
    from scipy.sparse.csr import csr_matrix

    SBM=newFromFileSBM('SBM1.dat')

    printSBM(SBM)

    r,A = SBMtoSparse(SBM)

    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            assert abs(A[i,j] - getValueSBM(SBM,i,j)) < eps


def test_sparseToSBM1():
    from Siconos.Numerics import sparseToSBM,getValueSBM, newFromFileSBM, printSBM, SBMtoSparse
    from scipy.sparse import csr_matrix, lil_matrix

    A = lil_matrix((100, 100))
    A.setdiag(range(100))
    A[0, :10] = range(10)
    A[1, 10:20] = A[0, :10]

    M = csr_matrix(A)

    v,SBM=sparseToSBM(2,M)

    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            print i,j, getValueSBM(SBM,i,j), M[i,j]
            assert abs(getValueSBM(SBM,i,j) - M[i,j]) < eps

def test_SBMtoSparseToSBM():

    from Siconos.Numerics import getValueSBM, newFromFileSBM, printSBM, SBMtoSparse, sparseToSBM
    from scipy.sparse.csr import csr_matrix

    SBM1=newFromFileSBM('SBM1.dat')

    printSBM(SBM1)

    r,SPARSE = SBMtoSparse(SBM1)

    v,SBM2 = sparseToSBM(3,SPARSE)

    printSBM(SBM2)

    assert SBM1.nbblocks == SBM2.nbblocks
    assert SBM1.blocknumber0 == SBM2.blocknumber0
    assert SBM1.blocknumber1 == SBM2.blocknumber1

    for i in range(SPARSE.shape[0]):
        for j in range(SPARSE.shape[1]):
            print i,j, getValueSBM(SBM1,i,j), getValueSBM(SBM2,i,j)
            assert (getValueSBM(SBM1,i,j) - getValueSBM(SBM2,i,j)) < eps



