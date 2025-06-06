import numpy as np

import siconos.kernel as sk
import siconos.numerics as sn
import warnings


def test_autocast():
    dsA = sk.LagrangianDS([0], [0], [[1]])
    dsB = sk.FirstOrderLinearDS([0], [[1]])
    nsds = sk.NonSmoothDynamicalSystem(0, 0)
    nsds.insertDynamicalSystem(dsA)
    nsds.insertDynamicalSystem(dsB)

    not_failed = 1
    if not isinstance(nsds.dynamicalSystem(dsA.number()), sk.LagrangianDS):
        not_failed = 0
    if not isinstance(nsds.dynamicalSystem(dsB.number()), sk.FirstOrderLinearDS):
        not_failed = 0

    assert(not_failed)

    #assert(type(nsds.dynamicalSystem(dsA.number())) == sk.LagrangianDS)
    #assert(type(nsds.dynamicalSystem(dsB.number())) == sk.FirstOrderLinearDS)


def test_getVector():
    assert (sk.getVector([1, 2, 3]) == np.array([1, 2, 3])).all()
    v = sk.SiconosVector(3)
    v.setValue(0, 1)
    v.setValue(1, 2)
    v.setValue(2, 4)

    assert (sk.getVector(v) != np.array([1, 2, 3])).any()

    assert (sk.getVector(v) == np.array([1, 2, 4])).all()

    v1 = sk.SiconosVector([1, 2, 3])
    v2 = sk.SiconosVector(np.asarray([1, 2, 3]))

    assert (sk.getVector(v1) == sk.getVector(v2)).all()


def test_castVector():
    i = [1.0, 4.0, 3.0]
    v = sk.SiconosVector([1, 2, 3])
    warnings.simplefilter("always")
    assert str(v) == "[3](1,2,3)"
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
    for x, y in zip(v, i):
        assert x == y
    
    for x, y in zip(list(v), i):
        assert x == y
    assert 3.0 in v
    assert 5.0 not in v


def test_getMatrix():

    assert (sk.getMatrix([[1, 2, 3]]) == np.array([[1, 2, 3]])).all()

    m = sk.SimpleMatrix(1, 3)

    m.setValue(0, 0, 1)

    m.setValue(0, 1, 2)

    m.setValue(0, 2, 3)

    assert (sk.getMatrix(m) == np.array([[1, 2, 3]])).all()

    assert (sk.getMatrix(m) != np.array([[1, 0, 3]])).any()

    m1 = sk.SimpleMatrix(((1, 2, 3), (4, 5, 6)))
    m2 = sk.SimpleMatrix(np.array([[1, 2, 3], [4, 5, 6]]))
    assert (sk.getMatrix(m1) == sk.getMatrix(sk.SimpleMatrix(m2))).all()


def test_matrix_bracket_operator():
    M = sk.SimpleMatrix(10, 10)
    M.zero()
    M_nparray = np.zeros((10, 10))
    print(M_nparray)

    def fill_matrix(M, M_nparray):
        for i in range(M.size(0)):
            for j in range(M.size(1)):
                M[i, j] = float(i + j)
                M_nparray[i, j] = float(i + j)

        return

    fill_matrix(M, M_nparray)

    # print(M, type(M))
    # print(M_nparray, type(M_nparray))
    for i in range(M.size(0)):
        for j in range(M.size(1)):
            if M[i, j] != i + j:
                return False

    # assert((M == M_nparray).all())

    M[0, 1] = 266.0
    if M[0, 1] != 266.0:
        return M[0, 1] == 266.0

    try:
        # Slice indexing for SimpleMatrix is not yet implemented.
        M[0:1] = [0, 2]
    except Exception as e:
        print(e)

    try:
        M[1.0, 1] = 4.0
        # SiconosMatrix must be indexed by integer
    except Exception as e:
        print(e)


def test_LagrangianDS_setMassPtr():
    class LDS(sk.LagrangianDS):

        pass

    lds = LDS()

    lds.setMassPtr([[1, 2, 3], [4, 5, 6]])

    assert (lds.mass() == np.array([[1, 2, 3], [4, 5, 6]])).all()


def test_LagrangianScleronomousR_setJachqPtr():
    class Rel(sk.LagrangianScleronomousR):
        pass

    r = Rel()
    j = np.array([[1, 2, 3], [4, 5, 6]])
    r.setJachqPtr(j)
    # C is transposed()
    r.C()

    assert np.max(r.C() - np.array([[1, 2, 3], [4, 5, 6]])) == 0.0
    assert np.max(r.C() - np.array([[0, 2, 3], [4, 5, 6]])) == 1.0

    r.setJachqPtr(r.C())

    r.C()

    assert np.max(r.C() - np.array([[1, 2, 3], [4, 5, 6]])) == 0.0
    assert np.max(r.C() - np.array([[0, 2, 3], [4, 5, 6]])) == 1.0


def test_SolverOption():

    lcp = sk.LCP()

    # Check default solver
    assert lcp.numericsSolverOptions().solverId == sn.SICONOS_LCP_LEMKE

    lcp.numericsSolverOptions().iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
    lcp.numericsSolverOptions().dparam[sn.SICONOS_DPARAM_TOL] = 1e-12

    iparam = lcp.numericsSolverOptions().iparam
    dparam = lcp.numericsSolverOptions().dparam

    assert iparam[sn.SICONOS_IPARAM_MAX_ITER] == 1000
    assert dparam[sn.SICONOS_DPARAM_TOL] == 1e-12


def test_BoundaryCondition():

    B = sk.BoundaryCondition([1, 2, 3])

    print(B)

    print(B.velocityIndices())

    B.velocityIndices()[2] = 5

    assert (B.velocityIndices() == [1, 2, 5]).all()


if __name__ == "__main__":

    # execute only if run as a script
    test_autocast()
    test_matrix_bracket_operator()
