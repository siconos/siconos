import Siconos.Kernel as SK
import Siconos.Numerics as SN

import numpy as np
import matplotlib.pyplot as plt

    ## \brief Constructor
    #
    # \param  is a  (optional)
_kappa = .925
_g = 9.81

h = 1e-3

theta = 1.0
gamma = 1.0

r = np.empty((2,))
def computeg(l):
    #print 'call computeg'
    #print(l)
    #print(R)
    r[1] = _g*l[0]/(1.0 - _kappa*l[0]*l[1])
    v_gamma = ((1.0 - gamma)*xk[1] + gamma*(xk[1] + h*r[1]))
    r[0] = -_kappa*l[0]*l[1]*(v_gamma)
    #print(R)
    #print('computeg done')
    return r

B = np.empty((2, 4))
def computeJacg(l):
    #print('call computeJacglambda')
    #print(B)
    rr = computeg(l)
    v_gamma = ((1.0 - gamma)*xk[1] + gamma*(xk[1] + h*rr[1]))
    B[1, 0] = 2*_g/(1-_kappa*l[0]*l[1])**2
    B[1, 1] = 0.0
    B[1, 2] = 2*(_g*_kappa*l[0]**2)/(1-_kappa*l[0]*l[1])**2
    B[1, 3] = 0.0
    B[0, 0] = -2*_kappa*l[1]*(v_gamma) - gamma*h*_kappa*l[0]*l[1]*B[1, 0]
    B[0, 1] = 0.0
    B[0, 2] = -2*_kappa*l[0]*(v_gamma) - gamma*h*_kappa*l[0]*l[1]*B[1, 2]
    B[0, 3] = 0.0
    return B

A = np.array(((1, h), (0, 1)))
xk = np.empty((2,))
xkp1 = np.empty((2,))
Fmcp = np.empty((4,))
l = np.empty((2,))
JacF = np.zeros((4,4))

def mcp_function(n1, n2, z, F):
#    print("z --------------------------------------")
#    print(z)

    l[0] = 2*z[0] - 1.0
    l[1] = 2*z[2] - 1.0
    r = computeg(l)
    v_theta = ((1.0 -theta)*xk[1] + theta*(xk[1] + h*r[1]))
    F[0] = xk[0] + h*v_theta + h*r[0] + z[1]
    F[2] = xk[1] + h*r[1] + z[3]
    F[1] = 1.0 - z[0]
    F[3] = 1.0 - z[2]
#    print("F --------------------------------------")
#    print(F)
    pass

def mcp_Nablafunction(n1, n2, z, nabla_Fmcp):
    mcp_function(n1, n2, z, Fmcp)
    l[0] = 2*z[0] - 1.0
    l[1] = 2*z[2] - 1.0
    Jacg = h*computeJacg(l)
    JacF[0, :] = Jacg[0, :] + theta*h*Jacg[1, :]
    JacF[0, 1] += 1.0
    JacF[1, 0] = -1.0
    JacF[1, 1:] = 0.0
    JacF[2, :] = Jacg[1, :]
    JacF[2, 3] += 1.0
    JacF[3, :] = 0.0
    JacF[3, 2] = -1.0
    nabla_Fmcp[:] = JacF
    pass


if __name__ == '__main__':
    xk[0] = 1.
    xk[1] = 10.

    T = 10.0
    t = 0.0
    z = np.zeros((4,))
    w = np.empty((4,))

#    mcp_function(z)
#    mcp_Nablafunction(z)
    N = int(T/h + 10)
    print(N)
    mcp = SN.MixedComplementarityProblem2(0, 4, mcp_function, mcp_Nablafunction)
    SO=SN.SolverOptions(mcp, SN.SICONOS_MCP_NEWTON_FBLSA)
    SO.dparam[0] = 1e-24
    SO.iparam[0] = 150
    SO.iparam[3] = 2
    SO.iparam[4] = 10

    lambdaPM = np.empty((N, 4))
    signs = np.empty((N, 2))
    sol = np.empty((N, 2))
    sol[0, :] = xk

    k = 0

    while t <= T:
        k += 1
        info = SN.mcp_newton_minFBLSA(mcp, z, w, SO)
        #info = SN.mcp_newton_FBLSA(mcp, z, w, SO)
        if info > 0:
            mcp_function(0, 4, z, w)
            sol[k, 0] = w[0] - z[1]
            sol[k, 1] = w[2] - z[3]
            if sol[k, 0] < -1e-7 and np.abs(z[1]) < 1e-10:
                z[1] = -sol[k, 0]
                z[0] = 1.0
            if xk[1] < -1e-7 and np.abs(z[3]) < 1e-10:
                z[3] = -sol[k, 1]
                z[2] = 1.0
            if z[1] < -1e-7:
                z[1] = 0.0
                z[0] = 0.0
            if z[3] < -1e-7:
                z[3] = 0.0
                z[2] = 0.0
            if z[1] > 1e-7 and z[0] < 1.0 - 1e-7:
                z[0] = 1.0
            if z[3] > 1e-7 and z[2] < 1.0 - 1e-7:
                z[2] = 1.0

            info = SN.mcp_newton_minFBLSA(mcp, z, w, SO)
            print('iter {:} ; solver iter = {:} ; prec = {:}'.format(k, SO.iparam[1], SO.dparam[1]))
            if info >0:
                print('MCP solver failed ! info = {:}'.format(info))
                print(xk)
                print(z)
                print(w)
#        else:
#            print('iter {:} ; solver iter = {:} ; prec = {:}'.format(k, SO.iparam[1], SO.dparam[1]))

        mcp_function(0 ,4, z, w)
        sol[k, 0] = w[0] - z[1]
        sol[k, 1] = w[2] - z[3]
        xk[:] = sol[k, :]
        signs[k, 0] = z[0] - w[1]
        signs[k, 1] = z[2] - w[3]
        t = k*h
        #z[:] = 0.0

#    np.savetxt("dataZIsol.txt", sol)
#    np.savetxt("dataZIlambdaPM.txt", lambdaPM)
#    np.savetxt("dataZIsign.txt", signs)

    plt.figure()
    plt.plot(sol[:, 0], sol[:, 1], 'b-*')
    plt.figure()
    plt.plot(sol[:, 0])
    plt.plot(sol[:, 1])
    plt.figure()
    plt.plot(signs[:, 0])
    plt.plot(signs[:, 1])
    plt.show()

    pos = np.abs(sol[:, 0])
    velocity = (1 - _kappa*np.sign(sol[:, 0]*sol[:, 1]))*sol[:, 1]*np.sign(sol[:, 0])

    plt.subplot(311)
    plt.title('position')
    plt.plot(pos)
    plt.grid()
    plt.subplot(312)
    plt.title('velocity')
    plt.plot(velocity)
    plt.grid()
#    plt.subplot(313)
#    plt.title('control input')
#    plt.plot(dataPlot[:,0], control)
#    plt.grid()
    plt.show()

#    indx = np.nonzero(dataPlot[:, 0]>30)
#    ttt = dataPlot[indx, 0].flatten()
#
#    plt.subplot(311)
#    plt.title('position')
#    plt.plot(ttt, pos[indx])
#    plt.grid()
#    plt.subplot(312)
#    plt.title('velocity')
#    plt.plot(ttt, velocity[indx])
#    plt.grid()
##    plt.subplot(313)
##    plt.title('control input')
##    plt.plot(ttt, control[indx])
#    plt.grid()
#    plt.show()
