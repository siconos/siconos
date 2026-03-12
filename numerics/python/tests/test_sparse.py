# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2026 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from numpy import finfo, double

from siconos.numerics import SBM_from_csparse, SBM_get_value
from scipy.sparse.csr import csr_matrix, csc_matrix
from scipy.sparse import lil_matrix
from siconos.numerics import (
    SBM_get_value,
    SBM_new_from_file,
    SBM_to_sparse,
    SBM_from_csparse,
)

eps = finfo(double).eps


def test_from_csr1():

    M = csr_matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    print(M.indices)
    print(M.indptr)
    print(M.data)

    blocksize = 3

    r, SBM = SBM_from_csparse(blocksize, M)

    assert SBM_get_value(SBM, 0, 0) == 1
    assert SBM_get_value(SBM, 0, 1) == 2
    assert SBM_get_value(SBM, 0, 2) == 3
    assert SBM_get_value(SBM, 1, 0) == 4
    assert SBM_get_value(SBM, 1, 1) == 5
    assert SBM_get_value(SBM, 1, 2) == 6
    assert SBM_get_value(SBM, 2, 0) == 7
    assert SBM_get_value(SBM, 2, 1) == 8
    assert SBM_get_value(SBM, 2, 2) == 9


def test_from_csc1():

    M = csc_matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    # print(M.indices)
    # print(M.indptr)
    # print(M.data)

    blocksize = 3

    r, SBM = SBM_from_csparse(blocksize, M)

    assert SBM_get_value(SBM, 0, 0) == 1
    assert SBM_get_value(SBM, 0, 1) == 2
    assert SBM_get_value(SBM, 0, 2) == 3
    assert SBM_get_value(SBM, 1, 0) == 4
    assert SBM_get_value(SBM, 1, 1) == 5
    assert SBM_get_value(SBM, 1, 2) == 6
    assert SBM_get_value(SBM, 2, 0) == 7
    assert SBM_get_value(SBM, 2, 1) == 8
    assert SBM_get_value(SBM, 2, 2) == 9


# scipy csr 3x3 block
def test_from_csr2():

    M = csr_matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    print(M.indices)
    print(M.indptr)
    print(M.data)
    blocksize = 1
    r, SBM = SBM_from_csparse(blocksize, M)

    assert SBM_get_value(SBM, 0, 0) == 1
    assert SBM_get_value(SBM, 0, 1) == 2
    assert SBM_get_value(SBM, 0, 2) == 3
    assert SBM_get_value(SBM, 1, 0) == 4
    assert SBM_get_value(SBM, 1, 1) == 5
    assert SBM_get_value(SBM, 1, 2) == 6
    assert SBM_get_value(SBM, 2, 0) == 7
    assert SBM_get_value(SBM, 2, 1) == 8
    assert SBM_get_value(SBM, 2, 2) == 9


# scipy csr 3x3 block
def test_from_csc162x162(datafile):  # uses datafile pytest fixture

    try:
        from scipy.sparse import load_npz
    except ImportError:
        return 0

    M = load_npz(datafile("csc162x162.npz"))
    # M = load_npz('data/csc162x162.npz')

    blocksize = 9
    r, SBM = SBM_from_csparse(blocksize, M)
    assert SBM_get_value(SBM, 0, 0) == M[0, 0]
    assert SBM_get_value(SBM, 161, 161) == M[161, 161]


#  # scipy csr 3x3 block
# def test_from_csc162x162_to_dense():
#     import siconos.numerics as sn
#     from siconos.numerics import SBM_from_csparse, SBM_get_value, NM_display,
# NM_read_in_filename
#     from scipy.sparse import csr_matrix, load_npz, linalg

#     M = load_npz(os.path.join(working_dir, 'data/csc162x162.npz'))

#     #M = load_npz('data/csc162x162.npz')
#     #print(linalg.eigs(M+M.transpose(),  which='LR')[1])
#     #print(M.indices)
#     #print(M.indptr)
#     #print(M.data)
#     NM_csc=None
#     NM_read_in_filename(NM_csc,'data/NM_csc_162x162.dat')


#     Mdense= M.todense()
#     print(Mdense)


def test_SBM_to_sparse1(datafile):  # uses datafile pytest fixture

    SBM = SBM_new_from_file(datafile("SBM1.dat"))

    r, A = SBM_to_sparse(SBM)

    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            assert abs(A[i, j] - SBM_get_value(SBM, i, j)) < eps


def test_SBM_from_csparse1():

    A = lil_matrix((100, 100))
    A.setdiag(range(100))
    A[0, :10] = range(10)
    A[1, 10:20] = A[0, :10]

    M = csr_matrix(A)

    v, SBM = SBM_from_csparse(2, M)

    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            assert abs(SBM_get_value(SBM, i, j) - M[i, j]) < eps


def test_SBM_to_SBM_from_csparse(datafile):  # uses datafile pytest fixture

    SBM1 = SBM_new_from_file(datafile("SBM1.dat"))

    r, SPARSE = SBM_to_sparse(SBM1)

    v, SBM2 = SBM_from_csparse(3, SPARSE)

    assert SBM1.nbblocks == SBM2.nbblocks
    assert SBM1.blocknumber0 == SBM2.blocknumber0
    assert SBM1.blocknumber1 == SBM2.blocknumber1

    for i in range(SPARSE.shape[0]):
        for j in range(SPARSE.shape[1]):
            assert (SBM_get_value(SBM1, i, j) - SBM_get_value(SBM2, i, j)) < eps


if __name__ == "__main__":
    test_from_csc1()
    test_from_csr1()
    test_from_csr2()
    test_SBM_to_SBM_from_csparse()
    test_from_csc162x162()
    # test_from_csc162x162_to_dense()
    test_SBM_to_sparse1()
    test_SBM_from_csparse1()
