import numpy 
import scipy
import math

# Explicit simulation 

# \dot x = Ax +a + B \sigma
# \sigma \in Sgn (Cx+D)

def computeOneStepExplicit(x,ti,tf,A,B,C,D,a):
    info=0
    xi=numpy.array(x)
    h=tf-ti
    y=numpy.dot(C,xi)+D
    print("y=", y)
    sigma = numpy.array(y)
    print(numpy.size(x))
    for i in range(numpy.size(y)):
        if (y[i] > 0.0) :
            sigma[i] = +1.0
        elif (y[i] < 0.0) :
            sigma[i] = 0.0
        else :
            sigma[i] = 0.5
    
    xtmp = xi + h*numpy.dot(A,xi) +h* a  + h*numpy.dot(B,sigma) 
    print("xtmp=",xtmp)      
    for i in range(numpy.size(y)):
        x[i]=xtmp[i]
    
        
    print("sigma =" , sigma)

    print("x=",x)
    return info 
