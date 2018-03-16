try:
    import siconos.kernel as SK

except (ImportError):
    print('Could not import Siconos.* module')

import numpy


class MyNonLinearR(SK.FirstOrderNonLinearR):
    # \brief Constructor
    #
    # \param  is a  (optional)
    def __init__(self, C, B):
        SK.FirstOrderNonLinearR.__init__(self)
        self.setBPtr(B)
        self.setCPtr(C)
        self._verbose=False
        return

    def computeh(self, time, x, l, z, y):
        if self._verbose:
            print('################################## call computeh')
            print(x)
        y[:] = self.C().dot(x)[:]
        if self._verbose:
            print(y)
            print('################################## computeh done')
        pass

    def computeg(self, time, x, l, z, R):
        if self._verbose:
            print('################################## call computeg')
        R[:] = self.B().dot(l)[:]
        if self._verbose:
            print('################################## call computeg done')
        pass

    def computeJachx(self, time, x, l, z, C):
        # since C is given, this should not be called
        if self._verbose:
            print('********************************* call computeJachx')
            print(C)
        C[:] = numpy.eye(2)[:]
        if self._verbose:
            print('********************************* call computeJachx done ')
        pass

    def computeJacglambda(self, time, x, l, z, B):
        # since B is given, this should not be called
        B[:] = self.B()
        pass

    def computeJacgx(self, time, x, l, z, K):
        if self._verbose:
            print('********************************* call computeJacgx')
        K[:] = numpy.zeros(2)[:]
        if self._verbose:
            print(K)
            print('********************************* call computeJacgx done')
        pass

    def computeJachlambda(self, time, x, l, z, D):
        # since D is given, this should not be called
        if self._verbose:
            print('********************************* call computeJachlambda')
        D[:] = numpy.zeros(2)[:]
        if self._verbose:
            print(D)
            print('********************************* call computeJachlambda done')
        pass
