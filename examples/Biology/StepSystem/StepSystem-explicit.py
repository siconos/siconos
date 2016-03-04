import numpy 
import scipy

import sys

from ComputeOneStep import *

#state
n=2
x=numpy.zeros(n, 'd')


# initial condition
x[0]=10
x[1]=5

gamma1=-1.5
gamma2=-4.5

kappa1=40
kappa2=40

theta1=8.0
theta2=8.0


A= numpy.array([[gamma1,0],[0, gamma2]])

a= numpy.array([kappa1,kappa2])

B= numpy.array([[-kappa1,0],[0, -kappa2]])

C= numpy.array([[1.0,0],[0, 1.0]])

D= numpy.array([-theta1,-theta2])


t0=0.0
T=1.0
h=0.001
N= int((T-t0)/h)+1
#N=5000
noutput=14
noutputfine=8
dataPlot = numpy.empty((N+1,noutput))

try :
    k = 0     
    dataPlot[k, 0] = t0

    print("x", x)
    dataPlot[k, 1] = x[0]
    dataPlot[k, 2] = x[1]


 #
    # Time loop
    #
    ti=t0
    tf=ti+h
    while ((tf <= T + 1e-10 ) and k < N) :
        print("=========== Iteration k ===============", k, ti, tf) 
        info = computeOneStepExplicit(x,ti,tf,A,B,C,D,a) 
            
        dataPlot[k+1, 0] = tf
        # print "xd", xd(tf)
        print("x", x)
        dataPlot[k+1, 1] = x[0]
        dataPlot[k+1, 2] = x[1]
        
        ti=ti+h
        tf=tf+h
        k=k+1

    dataPlot.resize(k+1,noutput)
    numpy.savetxt("./result.dat",dataPlot)
    

except :
    print("Unexpected error:", sys.exc_info()[0])
    raise
