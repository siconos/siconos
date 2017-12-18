import siconos.kernel as SK
import siconos.numerics as SN

import numpy as np

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
    ## \brief Constructor
    #
    # \param  is a  (optional)
_kappa = .7
_g = 9.81

h = 1e-3

theta = 1.0
gamma = 1.0

r = np.empty((2,))

B = np.empty((2, 2))

A = np.array(((1, h), (0, 1)))
xk = np.empty((2,))
xkp1 = np.empty((2,))
Fmcp = np.empty((4,))
l = np.empty((2,))
JacF = np.zeros((2,2))

def vi_function(n, l, F):
#    print("z --------------------------------------")
#    print(z)
    # first compute g(lambda)
    xk0 = xk[0]
    xk1 = xk[1]
    l0 = l[0]
    l1 = l[1]
    invR = 1.0/(1.0 - _kappa*l0*l1)
    F[1] = h*_g*l0*invR
    v_gamma = xk1 + gamma*F[1]
    F[0] = -h*_kappa*l0*l1*v_gamma

    v_theta = xk1 + theta*F[1]
    F[0] += xk0 + h*v_theta
    F[1] += xk1
    pass

def vi_nabla_function(n, l, nabla_Fmcp):
    l0 = l[0]
    l1 = l[1]
    xk1 = xk[1]
    invR = 1.0/(1.0 - _kappa*l0*l1)
    invR2 = invR*invR
    rr1 = _g*l0*invR
    v_gamma = xk1 + gamma*h*rr1
    nabla_Fmcp[1, 0] = h*_g*invR2
    nabla_Fmcp[1, 1] = h*(_g*_kappa*l0**2)/invR2
    nabla_Fmcp[0, 0] = -h*(_kappa*l1*v_gamma + (-1.0 + gamma*h*_kappa*l0*l1)*nabla_Fmcp[1, 0])
    nabla_Fmcp[0, 1] = -h*(_kappa*l0*v_gamma + (-1.0 + gamma*h*_kappa*l0*l1)*nabla_Fmcp[1, 1])
    B[:] = nabla_Fmcp
    pass


if __name__ == '__main__':
    xk[0] = 1.
    xk[1] = 10.

    T = 4.0
    t = 0.0
    vi = SN.VI(2, vi_function)
    vi.set_compute_nabla_F(vi_nabla_function)
    lambda_ = np.array((0., 0.))
    xkp1 = np.array((0., 0.))

    SO = SN.SolverOptions(vi, SN.SICONOS_VI_BOX_QI)
    lb = np.array((-1.0, -1.0))
    ub = np.array((1.0, 1.0))
    vi.set_box_constraints(lb, ub)

    N = int(T/h + 10)
    print(N)
    SO.dparam[0] = 1e-16
    SO.iparam[0] = 50
    SO.iparam[3] = 0
    SO.iparam[4] = 10

    signs = np.empty((N, 2))
    sol = np.empty((N, 2))
    sol[0, :] = xk

    k = 0
    #SN.numerics_set_verbose(3)

    while t <= T:
        k += 1
        info = SN.variationalInequality_box_newton_QiLSA(vi, lambda_, xkp1, SO)
        #print('iter {:} ; solver iter = {:} ; prec = {:}'.format(k, SO.iparam[1], SO.dparam[1]))
        #info = SN.mcp_newton_FBLSA(mcp, z, w, SO)
        if info > 0:
#            mcp_function(0, 4, z, w)
#            sol[k, 0] = w[0] - z[1]
#            sol[k, 1] = w[2] - z[3]
#            if sol[k, 0] < -1e-7 and np.abs(z[1]) < 1e-10:
#                z[1] = -sol[k, 0]
#                z[0] = 1.0
#            if xk[1] < -1e-7 and np.abs(z[3]) < 1e-10:
#                z[3] = -sol[k, 1]
#                z[2] = 1.0
#            if z[1] < -1e-7:
#                z[1] = 0.0
#                z[0] = 0.0
#            if z[3] < -1e-7:
#                z[3] = 0.0
#                z[2] = 0.0
#            if z[1] > 1e-7 and z[0] < 1.0 - 1e-7:
#                z[0] = 1.0
#            if z[3] > 1e-7 and z[2] < 1.0 - 1e-7:
#                z[2] = 1.0
            print(lambda_)
            vi_function(2, signs[k-1, :], xkp1)
            #lambda_[0] = -np.sign(xkp1[0])
            #lambda_[1] = *np.sign(xkp1[1])
            if np.abs(xk[0]) < 1e-10:
                lambda_[0] = 0.01
            if np.abs(xk[1]) < 1e-10:
                lambda_[1] = 0.01
                print('ok lambda')
                print(lambda_)
            info = SN.variationalInequality_box_newton_QiLSA(vi, lambda_, xkp1, SO)
#            info = SN.mcp_newton_minFBLSA(mcp, z, w, SO)
            print('iter {:} ; solver iter = {:} ; prec = {:}'.format(k, SO.iparam[1], SO.dparam[1]))
            if info >0:
                print('VI solver failed ! info = {:}'.format(info))
                print(B)
                print(np.linalg.eig(B))
                print(np.linalg.cond(B))
                print(xk)
                print(lambda_)
                print(xkp1)
                # kaboom()
#        else:
#            print('iter {:} ; solver iter = {:} ; prec = {:}'.format(k, SO.iparam[1], SO.dparam[1]))

        vi_function(2, lambda_, xkp1)
        sol[k, :] = xkp1
        xk[:] = xkp1
        signs[k, :] = lambda_
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
    #plt.savefig('Zhuravlev_pv.png')

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
