# sparseToSBM
# scipy csr => 1x1 block
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
    from Siconos.Numerics import sparseToSBM, getValueSBM, newFromFileSBM, printSBM, SBMtoSparse
    from scipy.sparse.csr import csr_matrix

    SBM=newFromFileSBM('SBM1.dat')

    printSBM(SBM)

    A = SBMtoSparse(SBM)

    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            assert A[i,j] == getValueSBM(SBM,i,j)




#v,R = sparseToSBM(3,A)

#printSBM(R)

#RA = SBMtoSparse(R)


