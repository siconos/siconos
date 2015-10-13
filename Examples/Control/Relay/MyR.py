try:
    import siconos.kernel as SK
    import siconos.numerics as N

except (ImportError):
    print 'Could not import Siconos.* module'

import numpy

class MyR(SK.FirstOrderType2R):
    ## \brief Constructor
    #
    # \param  is a  (optional)
    def __init__(self, C, B):
        SK.FirstOrderType2R.__init__(self)
        self.setBPtr(B)
        self.setCPtr(C)
        return


    def computeh(self,time, x, l, y):
        # should be called in the new Olivier Version
        print 'call computeh'
        print(x)
        print(y)
        #numpy.copyto(y, numpy.dot(self.C(), x))
        y[:] = self.C().dot(x)[:]
        print('computeh done')
        pass

    def computeg(self,time, l, R):
        # should be called in the new Olivier Version
        print 'call computeg'
        print(l)
        print(R)
        #numpy.copyto(R, numpy.dot(self.B(), l))
        R[:] = self.B().dot(l)[:]
        print(R)
        print('computeg done')
        pass

    def computeJachx(self,time, x, l, C):
        #self.setJachxPtr(C) not exiting ??
        #numpy.copyto(self._C,self.jachx() ) version 1.7
        print('call computeJachx')
        print(C)
        C[:] = numpy.eye(2)[:]

#        numpy.copyto(SK.getMatrix(C), self.C())
        #print 'self.jachx()', self.jachx()
        pass

    def computeJacglambda(self,time, l, B):
        print('call computeJachlambda')
#        numpy.copyto(SK.getMatrix(B), self.B())
        #print 'self.jacglambda() = ', self.jacglambda()
        #self.setJachglambdaPtr(self._B) not callable in that form ?
        pass

    def computeJachlambda(self, time, x, l, D):
        print('call computeJachlambda')
        print(D)
#        numpy.copyto(SK.getMatrix(B), self.B())
        #print 'self.jacglambda() = ', self.jacglambda()
        #self.setJachglambdaPtr(self._B) not callable in that form ?
        pass




