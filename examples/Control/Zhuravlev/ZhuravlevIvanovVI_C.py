#!/usr/bin/env python

import siconos.numerics as SN

import numpy as np
import matplotlib.pyplot as plt

withPlot = False

try:
    from cffi import FFI
except:
    import sys
    print('no cffi module installed, exiting')
    sys.exit(0)

if __name__ == '__main__':

    h = 1e-3
    T = 10.0
    t = 0.0

    theta = 1.0
    gamma = 1.0
    g = 9.81
    kappa = 0.4

    xk = np.array((1., 10.))
    ffi = FFI()
    ffi.cdef('void set_cstruct(uintptr_t p_env, void* p_struct);')
    ffi.cdef('''typedef struct
             {
             int id;
             double* xk;
             double h;
             double theta;
             double gamma;
             double g;
             double kappa;
             double f_eval;
             double nabla_eval;
              } data;
             ''')

    data_struct = ffi.new('data*')
    data_struct.id = -1  # to avoid freeing the data in the destructor
    data_struct.xk = ffi.cast('double *', xk.ctypes.data)
    data_struct.h = h
    data_struct.theta = theta
    data_struct.gamma = gamma
    data_struct.g = g
    data_struct.kappa = kappa

    vi = SN.VI(2)
    D = ffi.dlopen(SN._numerics.__file__)
    D.set_cstruct(vi.get_env_as_long(), ffi.cast('void*', data_struct))
    vi.set_compute_F_and_nabla_F_as_C_functions('ZhuravlevIvanov.so', 'compute_F', 'compute_nabla_F')

    lambda_ = np.zeros((2,))
    xkp1 = np.zeros((2,))

    SO = SN.SolverOptions(vi, SN.SICONOS_VI_BOX_PATH)
    lb = np.array((-1.0, -1.0))
    ub = np.array((1.0, 1.0))
    vi.set_box_constraints(lb, ub)

    N = int(T/h + 10)
    print(N)
    SO.dparam[0] = 1e-12
    SO.iparam[0] = 50
    SO.iparam[2] = 1
    SO.iparam[3] = 0
    SO.iparam[4] = 10

    signs = np.empty((N, 2))
    sol = np.empty((N, 2))
    sol[0, :] = xk

    k = 0
    #SN.numerics_set_verbose(3)

    while t <= T:
        k += 1
        # info = SN.variationalInequality_box_newton_QiLSA(vi, lambda_, xkp1, SO)
        info = SN.vi_box_path(vi, lambda_, xkp1, SO)
        #print('iter {:} ; solver iter = {:} ; prec = {:}'.format(k, SO.iparam[1], SO.dparam[1]))
        #info = SN.mcp_newton_FBLSA(mcp, z, w, SO)
        if info > 0:
            print(lambda_)
#            vi_function(2, signs[k-1, :], xkp1)
            lambda_[0] = -np.sign(xkp1[0])
            lambda_[1] = -np.sign(xkp1[1])
            if np.abs(xk[0]) < 1e-10:
                lambda_[0] = 0.01
            if np.abs(xk[1]) < 1e-10:
                lambda_[1] = 0.01
                print('ok lambda')
                print(lambda_)
            # info = SN.variationalInequality_box_newton_QiLSA(vi, lambda_, xkp1, SO)
            info = SN.vi_box_path(vi, lambda_, xkp1, SO)
#            info = SN.mcp_newton_minFBLSA(mcp, z, w, SO)
            print('iter {:} ; solver iter = {:} ; prec = {:}'.format(k, SO.iparam[1], SO.dparam[1]))
            if info >0:
                print('VI solver failed ! info = {:}'.format(info))
                print(xk)
                print(lambda_)
                print(xkp1)
                #kaboom()
#        else:
#            print('iter {:} ; solver iter = {:} ; prec = {:}'.format(k, SO.iparam[1], SO.dparam[1]))

 #       vi_function(2, lambda_, xkp1)
        sol[k, 0:2] = xkp1
        np.copyto(xk, xkp1, casting='no')
        signs[k, 0:2] = lambda_
        t = k*h
        #z[:] = 0.0

    print('f_eval', data_struct.f_eval, 'nabla_eval', data_struct.nabla_eval)
#    np.savetxt("dataZIsol.txt", sol)
#    np.savetxt("dataZIlambdaPM.txt", lambdaPM)
#    np.savetxt("dataZIsign.txt", signs)

    if withPlot:
        plt.figure()
        plt.plot(sol[:, 0], sol[:, 1], 'b-*')
        plt.xlabel('s')
        plt.ylabel('v')
        plt.figure()
        plt.plot(sol[:, 0], label=r's')
        plt.plot(sol[:, 1], label=r'v')
        plt.legend(loc='best')
        plt.figure()
        plt.plot(signs[:, 0], label=r'$\lambda_1$')
        plt.plot(signs[:, 1], label=r'$\lambda_2$')
        plt.legend(loc='best')
        plt.show()

        pos = np.abs(sol[:, 0])
        velocity = (1 - kappa*np.sign(sol[:, 0]*sol[:, 1]))*sol[:, 1]*np.sign(sol[:, 0])

        plt.subplot(211)
        plt.title('position')
        plt.plot(pos)
        plt.grid()
        plt.subplot(212)
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
