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

B = np.empty((2, 4))
def computeJacglambda(x, lpm):
    #print('call computeJacglambda')
    #print(B)
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

K = np.zeros((2, 4))
def computeJacgx(x, lpm):
    #print('call computeJacgx')

    l = (lpm[0] - lpm[1], lpm[2] - lpm[3])
    K[0, 2] = -_kappa*l[0]*l[1]
    K[0, 3] = -K[0, 2]
    return K

def computeJachlambda(x, l):
    #print('call computeJachlambda')
    #print(D)
    return np.zeros((4,4))


A = np.array(((1, h), (0, 1)))
xk = np.empty((2,))
xkp1 = np.empty((2,))
lmp = np.empty((4,))
JacF = np.zeros((8,8))

def mcp_function(n1, n2, z, F):
#    print("z --------------------------------------")
#    print(z)
    xkp1[0] = z[0] - z[1]
    xkp1[1] = z[2] - z[3]
    lmp[:] = z[4:]
    F[0:2] = xkp1
#    print('xkp1')
#    print(xkp1)
    F[0:2] -= A.dot(xk)
#    print('xk')
#    print(xk)
#    print(F[0:2])
    F[0:2] -= computeg(xkp1, lmp)
    F[2] = lmp[0] + lmp[1] - 1.0
    F[3] = lmp[2] + lmp[3] - 1.0
    F[4] = z[0]
    F[5] = z[1]
    F[6] = z[2]
    F[7] = z[3]
#    print("F --------------------------------------")
#    print(F)
    pass

def mcp_Nablafunction(n1, n2, z, nabla_Fmcp):
    xkp1[0] = z[0] - z[1]
    xkp1[1] = z[2] - z[3]
    lmp[:] = z[4:]
    JacF[0, 0] = 1.0
    JacF[0, 1] = -1.0
    JacF[0, 2:4] = 0.0
    JacF[1, 0:2] = 0.0
    JacF[1, 2] = 1.0
    JacF[1, 3] = -1.0
    JacF[0:2, 0:4] -= computeJacgx(xkp1, lmp)
    JacF[2:4, 0:4] = 0.0
    JacF[0:2, 4:] = -computeJacglambda(xkp1, lmp)
    JacF[2, 4:6] = 1.0
    JacF[2, 6:8] = 0.0
    JacF[3, 4:6] = 0.0
    JacF[3, 6:8] = 1.0
    JacF[4:, 0:4] = np.eye(4)
    JacF[4:, 4:] = 0.0
    #nabla_Fmcp[:] = np.transpose(JacF)
    #print(JacF)
    nabla_Fmcp[:] = JacF


if __name__ == '__main__':
    xk[0] = 1.
    xk[1] = 10.

    T = 10.
    t = 0.
    z = np.zeros((8,))
    w = np.empty((8,))
    z[0] = np.max(xk[0], 0.0)
    z[1] = np.min(xk[0], 0.0)
    z[2] = np.max(xk[1], 0.0)
    z[3] = np.min(xk[1], 0.0)

#    mcp_function(z)
#    mcp_Nablafunction(z)
    N = int(T/h + 10)
    print(N)
    mcp = SN.MixedComplementarityProblem2(4, 4, mcp_function, mcp_Nablafunction)
    SO=SN.SolverOptions(mcp, SN.SICONOS_MCP_NEWTON_FBLSA)

    lambdaPM = np.empty((N, 4))
    signs = np.empty((N, 2))
    sol = np.empty((N, 2))
    sol[0, :] = xk

    k = 0

    while t <= T:
        k += 1
        info = SN.mcp_newton_minFBLSA(mcp, z, w, SO)
        #info = SN.mcp_newton_FBLSA(mcp, z, w, SO)
        print('iter {:} ; solver iter = {:} ; prec = {:}'.format(k, SO.iparam[1], SO.dparam[1]))
        if info > 0:
            print('MCP solver failed ! info = {:}'.format(info))

        sol[k, 0] = z[0] - z[1]
        sol[k, 1] = z[2] - z[3]
        xk[:] = sol[k, :]
        lambdaPM[k, :] = z[4:]
        signs[k, 0] = z[4] - z[5]
        signs[k, 1] = z[6] - z[7]
        t += h


#    np.savetxt("dataZIsol.txt", sol)
#    np.savetxt("dataZIlambdaPM.txt", lambdaPM)
#    np.savetxt("dataZIsign.txt", signs)

    plt.figure()
    plt.plot(sol[:, 0], sol[:, 1])
    plt.figure()
    plt.plot(sol[:, 0])
    plt.plot(sol[:, 1])
    plt.figure()
    plt.plot(signs[:, 0])
    plt.plot(signs[:, 1])
    plt.show()

