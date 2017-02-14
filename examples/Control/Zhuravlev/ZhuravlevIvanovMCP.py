import siconos.kernel as SK
import siconos.numerics as SN

import numpy as np
import matplotlib.pyplot as plt

    ## \brief Constructor
    #
    # \param  is a  (optional)

h = 1e-3

withPlot = True


class ZI(object):
    def __init__(self, h, xk, theta, gamma, kappa, g):
        self.xk = xk
        self.h = h
        self.theta = theta
        self.gamma = gamma
        self._kappa = kappa
        self._g = g
        self.Jacg = np.zeros((2, 4))
        self.r = np.zeros((2,))

        self.f_eval = 0
        self.nabla_eval = 0

    def compute_Fmcp(self, n, z, F):
        l0 = 2.0*z[0] - 1.0
        l1 = 2.0*z[2] - 1.0
        r = self.computeg(l0, l1)
        v_theta = ((1.0 - self.theta)*self.xk[1] + self.theta*(self.xk[1] + self.h*self.r[1]))
        F[0] = self.xk[0] + self.h*v_theta + self.h*self.r[0] + z[1]
        F[2] = self.xk[1] + self.h*self.r[1] + z[3]
        F[1] = 1.0 - z[0]
        F[3] = 1.0 - z[2]
        self.f_eval += 1
        pass

    def compute_nabla_Fmcp(self, n, z, nabla_Fmcp):
        l0 = 2.0*z[0] - 1.0
        l1 = 2.0*z[2] - 1.0
        Jacg = self.h*self.computeJacg(l0, l1)
        nabla_Fmcp[0, :] = Jacg[0, :] + self.theta*self.h*Jacg[1, :]
        nabla_Fmcp[0, 1] += 1.0
        nabla_Fmcp[1, 0] = -1.0
        nabla_Fmcp[1, 1:] = 0.0
        nabla_Fmcp[2, :] = Jacg[1, :]
        nabla_Fmcp[2, 3] += 1.0
        nabla_Fmcp[3, :] = 0.0
        nabla_Fmcp[3, 2] = -1.0
        self.nabla_eval += 1
        pass

    def computeg(self, l0, l1):
        self.r[1] = self._g*l0/(1.0 - self._kappa*l0*l1)
        v_gamma = ((1.0 - self.gamma)*self.xk[1] + self.gamma*(self.xk[1] + self.h*self.r[1]))
        self.r[0] = -self._kappa*l0*l1*(v_gamma)
        #print(R)
        #print('computeg done')
        return self.r

    def computeJacg(self, l0, l1):
        #print('call computeJacglambda')
        #print(B)
        rr = self.computeg(l0, l1)
        v_gamma = ((1.0 - self.gamma)*self.xk[1] + self.gamma*(self.xk[1] + self.h*rr[1]))
        B = self.Jacg
        B[1, 0] = 2.0*self._g/(1-self._kappa*l0*l1)**2
        B[1, 2] = 2.0*(self._g*self._kappa*l0**2)/(1.0-self._kappa*l0*l1)**2
        B[0, 0] = -2.0*self._kappa*l1*(v_gamma) - self.gamma*self.h*self._kappa*l0*l1*B[1, 0]
        B[0, 2] = -2.0*self._kappa*l0*(v_gamma) - self.gamma*self.h*self._kappa*l0*l1*B[1, 2]
        return B



if __name__ == '__main__':
    xk = np.array((1., 10.))

    T = 4.0
    t = 0.0
    z = np.zeros((4,))
    w = np.empty((4,))

    kappa = 0.41
    g = 9.81
    theta = 1.0
    gamma = 1.0

    zi_syst = ZI(h, xk, theta, gamma, kappa, g)
    mcp = SN.MixedComplementarityProblem2(0, 4, zi_syst)




    SO=SN.SolverOptions(mcp, SN.SICONOS_MCP_NEWTON_FBLSA)
    SO.dparam[0] = 1e-24
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
            zi_syst.compute_Fmcp(0, 4, z, w)
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

        zi_syst.compute_Fmcp(0 ,4, z, w)
        sol[k, 0] = w[0] - z[1]
        sol[k, 1] = w[2] - z[3]
        xk[:] = sol[k, :]
        signs[k, 0] = z[0] - w[1]
        signs[k, 1] = z[2] - w[3]
        t = k*h
        #z[:] = 0.0

    print('f_eval', zi_syst.f_eval, 'nabla_eval', zi_syst.nabla_eval)

    pos = np.abs(sol[:, 0])
    velocity = (1 - kappa*np.sign(sol[:, 0]*sol[:, 1]))*sol[:, 1]*np.sign(sol[:, 0])

    if True:
        np.savetxt("dataZIsol.txt", sol)
        np.savetxt("dataZIlambdaPM.txt", lambdaPM)
        np.savetxt("dataZIsign.txt", signs)
        np.savetxt("dataZIpos.txt", pos)
        np.savetxt("dataZIvel.txt", velocity)

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
