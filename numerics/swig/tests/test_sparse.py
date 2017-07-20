# SBM_from_csparse
# scipy csr => 1x1 block
from numpy import finfo, double
eps = finfo(double).eps
from siconos.tests_setup import working_dir
import os

def test_from_csr1():

    from siconos.numerics import SBM_from_csparse, SBM_get_value
    from scipy.sparse.csr import csr_matrix

    M = csr_matrix([[1,2,3],
                    [4,5,6],
                    [7,8,9]])

    print(M.indices)
    print(M.indptr)
    print(M.data)

    r,SBM = SBM_from_csparse(3,M)

    assert SBM_get_value(SBM,0,0) == 1
    assert SBM_get_value(SBM,0,1) == 2
    assert SBM_get_value(SBM,0,2) == 3
    assert SBM_get_value(SBM,1,0) == 4
    assert SBM_get_value(SBM,1,1) == 5
    assert SBM_get_value(SBM,1,2) == 6
    assert SBM_get_value(SBM,2,0) == 7
    assert SBM_get_value(SBM,2,1) == 8
    assert SBM_get_value(SBM,2,2) == 9

# scipy csr 3x3 block
def test_from_csr2():

    from siconos.numerics import SBM_from_csparse, SBM_get_value
    from scipy.sparse.csr import csr_matrix

    M = csr_matrix([[1,2,3],
                    [4,5,6],
                    [7,8,9]])

    print(M.indices)
    print(M.indptr)
    print(M.data)

    r,SBM = SBM_from_csparse(1,M)

    assert SBM_get_value(SBM,0,0) == 1
    assert SBM_get_value(SBM,0,1) == 2
    assert SBM_get_value(SBM,0,2) == 3
    assert SBM_get_value(SBM,1,0) == 4
    assert SBM_get_value(SBM,1,1) == 5
    assert SBM_get_value(SBM,1,2) == 6
    assert SBM_get_value(SBM,2,0) == 7
    assert SBM_get_value(SBM,2,1) == 8
    assert SBM_get_value(SBM,2,2) == 9


def test_SBM_to_sparse1():
    from siconos.numerics import SBM_get_value, SBM_new_from_file, SBM_print, SBM_to_sparse
    from scipy.sparse.csr import csr_matrix

    SBM = SBM_new_from_file(os.path.join(working_dir, 'data/SBM1.dat'))

    r,A = SBM_to_sparse(SBM)

    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            assert abs(A[i,j] - SBM_get_value(SBM,i,j)) < eps


def test_SBM_from_csparse1():
    from siconos.numerics import SBM_from_csparse,SBM_get_value, SBM_new_from_file, SBM_print, SBM_to_sparse
    from scipy.sparse import csr_matrix, lil_matrix

    A = lil_matrix((100, 100))
    A.setdiag(range(100))
    A[0, :10] = range(10)
    A[1, 10:20] = A[0, :10]

    M = csr_matrix(A)

    v,SBM=SBM_from_csparse(2,M)

    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            assert abs(SBM_get_value(SBM,i,j) - M[i,j]) < eps

def test_SBM_to_SBM_from_csparse():

    from siconos.numerics import SBM_get_value, SBM_new_from_file, SBM_print, SBM_to_sparse, SBM_from_csparse
    from scipy.sparse.csr import csr_matrix

    SBM1 = SBM_new_from_file(os.path.join(working_dir, 'data/SBM1.dat'))

    r,SPARSE = SBM_to_sparse(SBM1)

    v,SBM2 = SBM_from_csparse(3,SPARSE)

    assert SBM1.nbblocks == SBM2.nbblocks
    assert SBM1.blocknumber0 == SBM2.blocknumber0
    assert SBM1.blocknumber1 == SBM2.blocknumber1

    for i in range(SPARSE.shape[0]):
        for j in range(SPARSE.shape[1]):
            assert (SBM_get_value(SBM1,i,j) - SBM_get_value(SBM2,i,j)) < eps



