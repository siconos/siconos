try:
    import siconos.kernel as SK

except (ImportError):
    print('Could not import Siconos.* module')

import numpy


class MyR(SK.FirstOrderNonLinearR):
    # \brief Constructor
    #
    # \param  is a  (optional)
    def __init__(self, C, B):
        SK.FirstOrderNonLinearR.__init__(self)
        self._kappa = .5
        self._g = 9.81
        self.setCPtr(C)
        return

    def computeh(self, time, x, l, z, y):
        # print 'call computeh'
        # print('x=',x)
        # print(y)
        y[:] = self.C().dot(x)
        # print('computeh done')
        pass

    def computeg(self, time, x, l, z, R):
        # print 'call computeg'
        # print(l)
        # print(R)
        R[0] = -self._kappa*l[0]*l[1]*x[1]
        R[1] = self._g*l[0]/(1-self._kappa*l[0]*l[1])
        # print(R)
        # print('computeg done')
        pass

    def computeJachx(self, time, x, l, z, C):
        # print('call computeJachx')
        # print('x=',x)
        # print(l)
        # print(C)
        C[:] = numpy.eye(2)
        pass

    def computeJacglambda(self, time, x, l, z, B):
        # print('call computeJacglambda')
        # print(B)
        B[0, 0] = -self._kappa*l[1]*x[1]
        B[0, 1] = -self._kappa*l[0]*x[1]
        B[1, 0] = self._g/(1-self._kappa*l[0]*l[1])**2
        B[1, 1] = (self._g*self._kappa*l[0]**2)/(1-self._kappa*l[0]*l[1])**2
        pass

    def computeJacgx(self, time, x, l, z, K):
        # print('call computeJacgx')

        K[:]=numpy.zeros(2)
        K[0, 1] = -self._kappa*l[0]*l[1]

        # print(K)
        pass

    def computeJachlambda(self, time, x, l, z, D):
        # print('call computeJachlambda')
        # print(D)
        D[:]=numpy.zeros(2)
        pass
