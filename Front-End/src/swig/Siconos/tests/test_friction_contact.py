#!/usr/bin/env python

from numpy import *
from Siconos.Numerics import *


NC = 3

M = eye(9)

q = array([-1, 1, 3, -1, 1, 3, -1, 1, 3])

mu = array([0.1, 0.1, 0.1]);

z = array([0,0,0,0,0,0,0,0,0])

setNumericsVerbose(2)

# projection on positive domain
def positive(x):
    return 1.*(x>=0)

def p_positive(x):
    return x*positive(x)


# euclidean ball projection
def p_ball(r,x,y):
    rp = p_positive(r)
    nv = maximum(sqrt(x*x+y*y),finfo(double).eps)
    k0 = nv - rp
    k1 = (k0 <= 0)
    inr = k1*x, k1*y
    k2 = (k0 > 0)
    inoutr = inr[0]+k2*rp*x/nv, inr[1]+k2*rp*y/nv        
    return inoutr

# vectorized ACF
def fac_n(rn,rx,ry,un,ux,uy,rhon,rhox,rhoy,mu):
    return p_positive(rn-rhon*un) - rn
    
def fac_t(rn,rx,ry,un,ux,uy,rhon,rhox,rhoy,mu):
    bx, by = p_ball(mu*rn,rx-rhox*ux,ry-rhoy*uy) 
    return bx - rx, by - ry

def Python_ACF(r,v,rho,mu):
    n = fac_n(r[0],r[1],r[2],v[0],v[1],v[2],rho[0],rho[1],rho[2],mu[0])

    t1,t2 = fac_t(r[0],r[1],r[2],v[0],v[1],v[2],rho[0],rho[1],rho[2],mu[0])
    return array([n, t1, t2])

#def test_nsgs():
#    SO = SolverOptions()
#    frictionContact3D_nsgs_setDefaultSolverOptions(SO)
#    info = frictionContact3D_driver(FrictionContactProblem(3,M,q,mu),q,z,SO, NumericsOptions())
#    print q
#    print z

def test_GlobalAlartCurnierFun():
    reactions = array(([ [ 0., 0., 0. ], 
                         [ 1., 1., 0.], 
                         [ 1., 1., 1.], 
                         [0., 1., 1.] ]))
                       
    velocities = array(([ [ 0., 0., 0. ],
                          [ 1., 3., 4. ],
                          [ 1., 1., 1. ],
                          [ 0., 1., 0. ] ]))

    rho      = ones(12)
    mu = array(([ 0.7, 0.7, 0.7, 0.7 ]))

    n_ac = frictionContact3D_GlobalAlartCurnierFunction(reactions.flat[:], 
                                                        velocities.flat[:], 
                                                        rho.flat[:], mu)
    
    
    p_ac = Python_ACF(reactions.transpose(), velocities.transpose(), rho.transpose(), mu.transpose())
    
    print n_ac 
    print p_ac.transpose().flat[:]
    assert ( (n_ac - p_ac.transpose().flat[:]).all() <= finfo(double).eps )


test_GlobalAlartCurnierFun()
