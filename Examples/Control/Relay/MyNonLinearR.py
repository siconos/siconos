try:
    import Siconos.Kernel
    import Siconos.Numerics as N
    
except (ImportError):
    print 'Could not import Siconos.* module'

import numpy

class MyNonLinearR(Siconos.Kernel.FirstOrderNonLinearR):
    ## \brief Constructor
    #
    # \param  is a  (optional)
    def __init__(self, C, B):
        Siconos.Kernel.FirstOrderNonLinearR.__init__(self)
        self.setBPtr(B)
        self.setCPtr(C)
        return

        
    def computeh(self,time, x, l, y):    
        # should be called in the new Olivier Version
        print 'call computeh'
        y = numpy.dot(self._C,x)
        pass

    def computeg(self,time, l, R):
        # should be called in the new Olivier Version
        print 'call computeg'
        R = numpy.dot(self._C,l)
        pass
 
    def computeJachx(self,time, x, l, C):
        #self.setJachxPtr(C) not exiting ??
        #numpy.copyto(self._C,self.jachx() ) version 1.7
        print 'self.jachx()', self.jachx()
        C = self._C
        #print 'self.jachx()', self.jachx()
        pass

    def computeJacglambda(self,time, l, B):
        print 'self.jacglambda() = ', self.jacglambda()
        B = self._B
        #print 'self.jacglambda() = ', self.jacglambda()
        #self.setJachglambdaPtr(self._B) not callable in that form ?
        pass
 
    
        
        
