import siconos.numerics as SN

import numpy as np
#import ctypes

try:
    from numba import jit
except:
    import sys
    print('numba not found, exiting')
    sys.exit(0)
#import numba

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



    ## \brief Constructor
    #
    # \param  is a  (optional)
hh = 1e-3

withPlot = False

#A = np.array(((1, h), (0, 1)))
xk = np.empty((2,))
xkp1 = np.empty((2,))
Fmcp = np.empty((4,))

#c_float_p = ctypes.POINTER(ctypes.c_double)

#@jit("void(i8,i8,double[:],double[:])",nopython=True)
#@jit("void(i8,i8,double[:],double[:])")
def mcp_function(n, z, F):
#    print("z --------------------------------------")
#    print(z)
    _kappa = .7
    _g = 9.81
    theta = 1.0
    gamma = 1.0
    #xk0 = z[4]
    #xk1 = z[5]
    #h = z[6]
    xk0 = xk[0]
    xk1 = xk[1]
    h = hh

    l0 = 2*z[0] - 1.0
    l1 = 2*z[2] - 1.0
    r1 = _g*l0/(1.0 - _kappa*l0*l1)
    v_gamma = ((1.0 - gamma)*xk1 + gamma*(xk1 + h*r1))
    r0 = -_kappa*l0*l1*(v_gamma)

    v_theta = ((1.0 -theta)*xk1 + theta*(xk1 + h*r1))
    F[0] = xk0 + h*v_theta + h*r0 + z[1]
    F[2] = xk1 + h*r1 + z[3]
    F[1] = 1.0 - z[0]
    F[3] = 1.0 - z[2]
#    print("F --------------------------------------")
#    print(F)
    pass

#@jit("void(i8,i8,double[:],double[:,:])",nopython=True)
#@jit("void(i8,i8,double[:],double[:,:])")
def mcp_Nablafunction(n, z, nabla_Fmcp):
    _kappa = .7
    _g = 9.81
    theta = 1.0
    gamma = 1.0
    l0 = 2*z[0] - 1.0
    l1 = 2*z[2] - 1.0
    #xk1 = z[5]
    #h = z[6]
    xk0 = xk[0]
    xk1 = xk[1]
    h = hh
    r1 = _g*l0/(1.0 - _kappa*l0*l1)
    v_gamma = ((1.0 - gamma)*xk1 + gamma*(xk1 + h*r1))
    nabla_Fmcp[2, 0] = h*2*_g/(1-_kappa*l0*l1)**2
    nabla_Fmcp[2, 1] = 0.0
    nabla_Fmcp[2, 2] = h*2*(_g*_kappa*l0**2)/(1-_kappa*l0*l1)**2
    nabla_Fmcp[2, 3] = 1.0
    nabla_Fmcp[0, 0] = h*(-2*_kappa*l1*(v_gamma) + (theta - gamma*h*_kappa*l0*l1)*nabla_Fmcp[2, 0])
    nabla_Fmcp[0, 1] = 1.0
    nabla_Fmcp[0, 2] = h*(-2*_kappa*l0*(v_gamma) + (theta - gamma*h*_kappa*l0*l1)*nabla_Fmcp[2, 2])
    nabla_Fmcp[0, 3] = 0.0
    nabla_Fmcp[1, 0] = -1.0
    nabla_Fmcp[1, 1] = 0.0
    nabla_Fmcp[1, 2] = 0.0
    nabla_Fmcp[1, 3] = 0.0
    nabla_Fmcp[3, 0] = 0.0
    nabla_Fmcp[3, 1] = 0.0
    nabla_Fmcp[3, 2] = -1.0
    nabla_Fmcp[3, 3] = 0.0
    pass


if __name__ == '__main__':

    print('This is an experimental file !')
    xk[0] = 1.
    xk[1] = 10.

    #numba.jit('void(i8,i8,double[:],double[:,:])', nopython=True)(mcp_Nablafunction).inspect_types()
#    mcp_Nablafunction.inspect_types()
    T = 10.0
    t = 0.0
    z = np.zeros((7,), dtype=np.dtype('d'))
    w = np.zeros((7,), dtype=np.dtype('d'))

#    mcp_function(z)
#    mcp_Nablafunction(z)
    N = int(T/hh + 10)
    print(N)
    mcp = SN.MixedComplementarityProblem2(0, 4, mcp_function, mcp_Nablafunction)
    SO=SN.SolverOptions(mcp, SN.SICONOS_MCP_NEWTON_FBLSA)
    SO.dparam[0] = 1e-13
    SO.iparam[0] = 100
    SO.iparam[3] = 2
    SO.iparam[4] = 10

    lambdaPM = np.empty((N, 4))
    signs = np.empty((N, 2))
    sol = np.empty((N, 2))
    sol[0, :] = xk
    z[4:6] = xk
    z[6] = hh

    k = 0

    while t <= T:
        k += 1
        info = SN.mcp_newton_minFBLSA(mcp, z, w, SO)
        #info = SN.mcp_newton_FBLSA(mcp, z, w, SO)
#        print('iter {:} ; solver iter = {:} ; prec = {:}'.format(k, SO.iparam[1], SO.dparam[1]))
        if info > 0:
            mcp_function(0, 4, z, w)
            sol[k, 0] = w[0] - z[1]
            sol[k, 1] = w[2] - z[3]
            if sol[k, 0] < -1e-7 and np.abs(z[1]) < 1e-10:
                z[1] = -sol[k, 0]
                z[0] = 1.0
            if z[4] < -1e-7 and np.abs(z[3]) < 1e-10:
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
                print(z[4:6])
                print(z)
                print(w)
#        else:
#            print('iter {:} ; solver iter = {:} ; prec = {:}'.format(k, SO.iparam[1], SO.dparam[1]))

        mcp_function(0 ,4, z, w)
        sol[k, 0] = w[0] - z[1]
        sol[k, 1] = w[2] - z[3]
        z[4:6] = sol[k, :]
        xk[:] = sol[k, :]
        signs[k, 0] = z[0] - w[1]
        signs[k, 1] = z[2] - w[3]
        t = k*hh
        #z[:] = 0.0

#    np.savetxt("dataZIsol.txt", sol)
#    np.savetxt("dataZIlambdaPM.txt", lambdaPM)
#    np.savetxt("dataZIsign.txt", signs)

    if withPlot:
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
        _kappa = .7
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
