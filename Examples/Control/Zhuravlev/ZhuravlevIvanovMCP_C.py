import siconos.numerics as SN

import numpy as np
import matplotlib.pyplot as plt

from cffi import FFI

withPlot = True



if __name__ == '__main__':
    xk = np.array((1., 10.))

    T = 10.0
    t = 0.0
    h = 1e-3
    z = np.zeros((4,))
    w = np.empty((4,))

    kappa = 0.9
    g = 9.81
    theta = 1.0
    gamma = 1.0

    mcp = SN.MixedComplementarityProblem2(0, 4)

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
             unsigned int f_eval;
             unsigned int nabla_eval;
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

    D = ffi.dlopen(SN._numerics.__file__)
    D.set_cstruct(mcp.get_env_as_long(), ffi.cast('void*', data_struct))
    mcp.set_compute_F_and_nabla_F_as_C_functions('ZhuravlevIvanov.so', 'compute_Fmcp', 'compute_nabla_Fmcp')






    SO=SN.SolverOptions(mcp, SN.SICONOS_MCP_NEWTON_FBLSA)
    SO.dparam[0] = 1.0e-24
    SO.iparam[0] = 150
    SO.iparam[3] = 2
    SO.iparam[4] = 10

    N = int(T/h + 10)
    print(N)
    lambdaPM = np.empty((N, 4))
    signs = np.empty((N, 2))
    sol = np.empty((N, 2))
    sol[0, :] = xk

    k = 0

    while t <= T:
        k += 1
        info = SN.mcp_newton_minFBLSA(mcp, z, w, SO)
        #info = SN.mcp_newton_FBLSA(mcp, z, w, SO)
        #print('iter {:} ; solver iter = {:} ; prec = {:}'.format(k, SO.iparam[1], SO.dparam[1]))
        if info > 0:
            #zi_syst.compute_Fmcp(0, 4, z, w)
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

        #zi_syst.compute_Fmcp(0 ,4, z, w)
        sol[k, 0] = w[0] - z[1]
        sol[k, 1] = w[2] - z[3]
        xk[:] = sol[k, :]
        signs[k, 0] = z[0] - w[1]
        signs[k, 1] = z[2] - w[3]
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
