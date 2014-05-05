try:
    import Siconos.Kernel as SK
    import Siconos.Numerics as SN

except (ImportError):
    print('Could not import Siconos.* module')

import numpy as np
import matplotlib.pyplot as plt

    ## \brief Constructor
    #
    # \param  is a  (optional)
_kappa = .5
_g = 9.81

h = 1e-3

def computeg(x, lpm):
    #print 'call computeg'
    #print(l)
    #print(R)
    r = np.empty((2,))
    l = (lpm[0] - lpm[1], lpm[2] - lpm[3])
    r[0] = -_kappa*l[0]*l[1]*x[1]
    r[1] = _g*l[0]/(1-_kappa*l[0]*l[1])
    #print(R)
    #print('computeg done')
    return r

def computeJachx(x, l):
    #print('call computeJachx')
    #print('x=',x)
    #print(l)
    #print(C)
    C = np.zeros((4, 2))
    C[0, 0] = 1 if x[0] < 0.0 else 0.0
    C[1, 0] = 1 if x[0] > 0.0 else 0.0
    C[2, 1] = 1 if x[1] < 0.0 else 0.0
    C[3, 1] = 1 if x[1] > 0.0 else 0.0
    return C

def computeJacglambda(x, lpm):
    #print('call computeJacglambda')
    #print(B)
    B = np.empty((2, 4))
    l = (lpm[0] - lpm[1], lpm[2] - lpm[3])
    B[0, 0] = -_kappa*l[1]*x[1]
    B[0, 1] = -B[0, 0]
    B[0, 2] = -_kappa*l[0]*x[1]
    B[0, 3] = -B[0, 2]
    B[1, 0] = _g/(1-_kappa*l[0]*l[1])**2
    B[1, 1] = -B[1, 0]
    B[1, 2] = (_g*_kappa*l[0]**2)/(1-_kappa*l[0]*l[1])**2
    B[1, 3] = -B[1, 2]
    return B

def computeJacgx(x, lpm):
    #print('call computeJacgx')

    K = np.zeros((2, 2))
    l = (lpm[0] - lpm[1], lpm[2] - lpm[3])
    K[0, 1] = -_kappa*l[0]*l[1]
    return K

def computeJachlambda(x, l):
    #print('call computeJachlambda')
    #print(D)
    return np.zeros((4,4))


A = np.array(((1, h), (0, 1)))
xk = np.empty((2,))
xkp1 = np.empty((2,))
lmp = np.empty((4,))
F = np.zeros((6,))
JacF = np.zeros((6,6))

def mcp_function(z):
    xkp1[:] = z[0:2]
    lmp[:] = z[2:]
    F[0:2] = xkp1
    F[0:2] -= A.dot(xk)
    F[0:2] -= computeg(xkp1, lmp)
    F[2:5:2] = np.min(xkp1, 0)
    F[3:6:2] = np.max(xkp1, 0)
    res = np.copy(F)
    return res

def mcp_Nablafunction(z):
    xkp1[:] = z[0:2]
    lmp[:] = z[2:]
    JacF[0:2, 0:2] = np.eye(2)
    JacF[0:2, 0:2] -= computeJacgx(xkp1, lmp)
    JacF[0:2, 2:] -= computeJacglambda(xkp1, lmp)
    JacF[2:, 0:2] = computeJachx(xkp1, lmp)
    JacF[2:, 2:] = np.zeros((4,4))
    res = np.copy(JacF)
    return res


if __name__ == '__main__':
    xk[0] = 1.
    xk[1] = 10.

    T = 10.
    t = 0.
    z = np.zeros((6,))
    w = np.empty((6,))
    z[0:2] = xk

    mcp_function(z)
    mcp_Nablafunction(z)
    N = int(T/h + 10)
    print(N)
    mcp = SN.MCP(2,4, mcp_function, mcp_Nablafunction)
    SO=SN.SolverOptions(mcp, SN.SICONOS_MCP_FB)
    SN.mcp_driver_init(mcp, SO)

    lambdaPM = np.empty((N, 4))
    signs = np.empty((N, 2))
    sol = np.empty((N, 2))
    sol[0, :2] = xk

    k = 0

    while t <= T:
        k += 1
        info = SN.mcp_FischerBurmeister(mcp, z, w, SO)
        if info > 0:
            print('MCP solver failed !')

        print('before sol')
        sol[k, :] = z[0:2]
        print('before xk')
        xk[:] = z[0:2]
        print('before lambdaPM')
        lambdaPM[k, :] = z[2:]
        print('before signs')
        signs[k, 0] = z[2] - z[3]
        print('before signs2')
        signs[k, 1] = z[4] - z[5]
        t += h


    np.savetxt("dataZIsol.txt", sol)
    np.savetxt("dataZIlambdaPM.txt", lambdaPM)
    np.savetxt("dataZIsign.txt", signs)

    plt.figure()
    plt.plot(sol[:, 0], sol[:, 1])

