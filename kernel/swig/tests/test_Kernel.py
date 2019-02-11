#!/usr/bin/env python
import numpy as np
import siconos.kernel as kernel

def test_autocast():
    dsA = kernel.LagrangianDS([0],[0],[[1]])
    dsB = kernel.FirstOrderLinearDS([0],[[1]])
    nsds = kernel.NonSmoothDynamicalSystem(0, 0)
    nsds.insertDynamicalSystem(dsA)
    nsds.insertDynamicalSystem(dsB)

    assert(type(nsds.dynamicalSystem(dsA.number())) == kernel.LagrangianDS)
    assert(type(nsds.dynamicalSystem(dsB.number())) == kernel.FirstOrderLinearDS)


def test_getVector():
    assert (kernel.getVector([1,2,3]) == np.array([1,2,3])).all()
    v = kernel.SiconosVector(3)
    v.setValue(0,1)
    v.setValue(1,2)
    v.setValue(2,4)

    assert (kernel.getVector(v) != np.array([1,2,3])).any()

    assert (kernel.getVector(v) == np.array([1,2,4])).all()

    v1 = kernel.SiconosVector([1, 2, 3])
    v2 = kernel.SiconosVector(np.asarray([1, 2, 3]))

    assert (kernel.getVector(v1) == kernel.getVector(v2)).all()


def test_castVector():
    i = [1.0,4.0,3.0]
    v = kernel.SiconosVector([1,2,3])
    assert str(v) == '[3](1,2,3)'
    repr(v)
    assert v[0] == 1.0
    try:
        v[5]
        raise Exception("expected IndexError")
    except IndexError:
        pass
    v[1] = 4
    assert v[1] == 4.0
    try:
        v[4] = 5
        raise Exception("expected IndexError")
    except IndexError:
        pass
    for x,y in zip(v,i):
        assert x == y
    for x,y in zip(list(v),i):
        assert x == y
    for x,y in zip(np.array(v),i):
        assert x == y
    assert 3.0 in v
    assert 5.0 not in v


def test_getMatrix():
    assert (kernel.getMatrix([[1,2,3]]) == np.array([[1,2,3]])).all()

    m = kernel.SimpleMatrix(1,3)

    m.setValue(0,0,1)

    m.setValue(0,1,2)

    m.setValue(0,2,3)

    assert (kernel.getMatrix(m) == np.array([[1,2,3]])).all()

    assert (kernel.getMatrix(m) != np.array([[1,0,3]])).any()

    m1 = kernel.SimpleMatrix(((1,2,3), (4,5,6)))
    m2 = kernel.SimpleMatrix(np.array([[1,2,3],[4,5,6]]))
    assert (kernel.getMatrix(m1) == kernel.getMatrix(kernel.SimpleMatrix(m2))).all()

def test_matrix_bracket_operator():
    M = kernel.SimpleMatrix(10,10)
    M.zero()
    def fill_matrix(M):
        for i in range(M.size(0)):
            for j in range(M.size(1)):
                M[i,j] = i+j
         
        return

    fill_matrix(M)
    print(M, type(M))
    for i in range(M.size(0)):
        for j in range(M.size(1)):
            assert(M[i,j]==i+j)

    M[0,1]=266.0
    assert(M[0,1]== 266.0)
    
    try:
        M[0:1]= [0,2]
    except Exception as e:
        print(e)
        pass
    try:
        M[1.0,1]= 4.0
    except Exception as e:
        print(e)
        pass



    
def test_LagrangianDS_setMassPtr():
    class LDS(kernel.LagrangianDS):
        pass

    lds = LDS()

    lds.setMassPtr([[1,2,3],[4,5,6]])

    assert (lds.mass() == np.array([[1,2,3],[4,5,6]])).all()


def test_LagrangianScleronomousR_setJachqPtr():
    class Rel(kernel.LagrangianScleronomousR):
        pass

    r = Rel()
    j = np.array([[1,2,3],[4,5,6]])
    r.setJachqPtr(j)
    # C is transposed()
    r.C()

    assert np.max(r.C() - np.array([[1,2,3],[4,5,6]])) == 0.
    assert np.max(r.C() - np.array([[0,2,3],[4,5,6]])) == 1.

    r.setJachqPtr(r.C())

    r.C()

    assert np.max(r.C() - np.array([[1,2,3],[4,5,6]])) == 0.
    assert np.max(r.C() - np.array([[0,2,3],[4,5,6]])) == 1.


def test_SolverOption():
    lcp = kernel.LCP()

    i0 = lcp.numericsSolverOptions().iparam[0]

    lcp.numericsSolverOptions().iparam[0] = i0+1

    assert lcp.numericsSolverOptions().iparam[0] != i0

    assert lcp.numericsSolverOptions().iparam[0] == i0+1

    d0 = lcp.numericsSolverOptions().dparam[0]

    lcp.numericsSolverOptions().dparam[0] = 0.5 * d0


    assert lcp.numericsSolverOptions().dparam[0] !=  d0

    assert lcp.numericsSolverOptions().dparam[0] ==  0.5 * d0


def test_BoundaryCondition():
    B = kernel.BoundaryCondition([1,2,3])

    print(B)

    print(B.velocityIndices())

    B.velocityIndices()[2]=5

    assert (B.velocityIndices() == [1, 2, 5]).all()
if __name__ == "__main__":
    # execute only if run as a script
    test_matrix_bracket_operator()

