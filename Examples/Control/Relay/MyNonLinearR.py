from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_Kernel', [dirname(__file__)])
        except ImportError:
            import _Kernel
            return _Kernel
        if fp is not None:
            try:
                _mod = imp.load_module('_Kernel', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _Kernel = swig_import_helper()
    del swig_import_helper
else:
    import _Kernel


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
        #self.setBPtr(B)
        #self.setCPtr(C)
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
 
    def computeOutput(self, *args):
        """
        computeOutput(FirstOrderNonLinearR self, double time, Interaction inter, VectorOfBlockVectors & DSlink, VectorOfVectors & workV, 
        VectorOfSMatrices & workM, SiconosMatrix osnsM, unsigned int level)
        
        FirstOrderNonLinearR.computeOutput(time,inter,DSlink,workV,workM,osnsM,level)
        default function to compute y
        
        Parameters:
        -----------
        
        time:  the current time
        
        inter:  the interaction using this relation
        
        derivativeNumber:  number of the derivative to compute (optional,
        default = 0) 
        """
        print "computeOutput"
        print args
        for i in range(7):
            print "args[",i,"]", args[i]
            
        args[1].display()

            

        return 0 # _Kernel.FirstOrderNonLinearR_computeOutput(self, *args)
        

    def computeInput(self, *args):
        """
        computeInput(FirstOrderNonLinearR self, double time, Interaction inter, VectorOfBlockVectors & DSlink, VectorOfVectors & workV, 
            VectorOfSMatrices & workM, SiconosMatrix osnsM, unsigned int level)

        FirstOrderNonLinearR.computeInput(time,inter,DSlink,workV,workM,osnsM,level)
        default function to compute r

        Parameters:
        -----------

        time:  the current time

        inter:  the interaction using this relation

        level:  the input "derivative" order of lambda used to compute input

        """
        return 0 # _Kernel.FirstOrderNonLinearR_computeInput(self, *args)

    def initComponents(self, *args):
        """initComponents(FirstOrderNonLinearR self, Interaction inter, VectorOfBlockVectors & DSlink, VectorOfVectors & workV, VectorOfSMatrices & workM)"""
        print "initComponents"
        print args
        for i in range(3):
            print "args[",i,"]", args[i]
            
        return _Kernel.FirstOrderNonLinearR_initComponents(self, *args)
