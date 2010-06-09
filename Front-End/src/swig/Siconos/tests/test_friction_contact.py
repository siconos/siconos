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

# positive(x) = 1 if x>=0 and 0 otherwise
def positive(x):
    return 1.*(x>=0)

# p_positive(x) = x if x>=0 and 0 otherwise
def p_positive(x):
    return x*positive(x)

# disk(r,x,y) = 1 if (x,y) inside disk(radius=r)
def inside_disk(r,x,y):
    return sqrt(x*x+y*y) <= r

assert inside_disk(0,0,0)
assert inside_disk(1,1,0)

# euclidean disk projection
def p_disk(r,x,y):
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
    bx, by = p_disk(mu*rn,rx-rhox*ux,ry-rhoy*uy) 
    return bx - rx, by - ry

def Python_ACF(r,v,rho,mu):
    n = fac_n(r[0],r[1],r[2],v[0],v[1],v[2],rho[0],rho[1],rho[2],mu[0])

    t1,t2 = fac_t(r[0],r[1],r[2],v[0],v[1],v[2],rho[0],rho[1],rho[2],mu[0])
    return array([n, t1, t2])

def gamma(x):
    assert (len(x) == 2)
    p = linalg.norm(x)

    assert p > 0

#    p3 = p*p*p
#    r = zeros([2,2])
#    r[0,0] = 1./p - x[0]*x[0]/p3
#    r[0,1] = -x[0]*x[1]/p3
#    r[1,0] = r[0,1]
#    r[1,1] = 1/p - x[1]*x[1]/p3

    r = eye(2)/p - array([x])*array([x]).transpose()/(p ** 3)

    return r

def dunphi2(rn,rx,ry,un,ux,uy,rhon,rhox,rhoy,mu):
    return positive(rn - rhon*un)*rhon
    
def drnphi2(rn,rx,ry,un,ux,uy,rhon,rhox,rhoy,mu):
    return positive(rhon*un - rn)

# not vectorized
def dutphi3(rn,rx,ry,un,ux,uy,rhon,rhox,rhoy,mu):
    if inside_disk(mu*rn, rx-rhox*ux,ry-rhoy*uy):
        return array([[rhox,0],[0,rhoy]])
    else:
        result = mu*rn*gamma(array([rx-rhox*ux,ry-rhoy*uy]))
        result[0,0] *= rhox
        result[1,1] *= rhoy

        assert not isnan(result[0,0])
        assert not isnan(result[1,1])

        return result


def drnphi3(rn,rx,ry,un,ux,uy,rhon,rhox,rhoy,mu):
    if inside_disk(mu*rn, rx-rhox*ux,ry-rhoy*uy):
        return array([[0.], [0.]])
    else:
        p = linalg.norm(array([rx-rhox*ux,ry-rhoy*uy]))
        return -mu*rn*array([[rx-rhox*ux],[ry-rhoy*uy]])/p

def drtphi3(rn,rx,ry,un,ux,uy,rhon,rhox,rhoy,mu):
    if inside_disk(mu*rn, rx-rhox*ux,ry-rhoy*uy):
        return zeros([2,2])
    else:
        result = mu*rn*gamma(array([rx-rhox*ux,ry-rhoy*uy]))
        result[0,0] *= rhox
        result[1,1] *= rhoy
        return eye(2) - result

#def test_nsgs():
#    SO = SolverOptions()
#    frictionContact3D_nsgs_setDefaultSolverOptions(SO)
#    info = frictionContact3D_driver(FrictionContactProblem(3,M,q,mu),q,z,SO, NumericsOptions())
#    print q
#    print z


def i3x3(ind,i,j):
    return ind*9 + 3*i + j

def nearly_equal(x,y):
    return abs(x-y) < 1e-10


def test_GlobalAlartCurnierFun():

    #
    # test equivalence of two implementations (Numerics/C, python)
    # over random samples
    #
    random.seed(1) # so it can be reproduced
    problem_dimension = 3000

    reactions = random.sample(problem_dimension) 
    velocities = random.sample(problem_dimension)

    rho      = ones(problem_dimension)
    mu = random.sample(problem_dimension/3)

    reactions = reactions.reshape(problem_dimension/3,3)
    velocities = velocities.reshape(problem_dimension/3,3)
    rho = rho.reshape(problem_dimension/3,3)

    # get AC values and elements of the jacobian
    n_ac, A, B = frictionContact3D_GlobalAlartCurnierFunction(reactions.flat[:], 
                                                              velocities.flat[:], 
                                                              rho.flat[:], mu)
    
    p_ac = Python_ACF(reactions.transpose(), velocities.transpose(), rho.transpose(), mu.transpose())
    
    for i in range(0,len(reactions)):
        r = reactions[i]
        v = velocities[i]
        mu_ = mu[i]
        rho_ = rho[i]
        assert nearly_equal(A[i3x3(i,0,0)],dunphi2(r[0],r[1],r[2],v[0],v[1],v[2],rho_[0],rho_[1],rho_[2],mu_))

        assert nearly_equal(B[i3x3(i,0,0)],drnphi2(r[0],r[1],r[2],v[0],v[1],v[2],rho_[0],rho_[1],rho_[2],mu_))
        
        assert nearly_equal(A[i3x3(i,1,0)],0.)
        
        assert nearly_equal(A[i3x3(i,2,0)],0.)
        
        dutphi3_ = dutphi3(r[0],r[1],r[2],v[0],v[1],v[2],rho_[0],rho_[1],rho_[2],mu_)
        
        assert nearly_equal(A[i3x3(i,1,1)], dutphi3_[0,0])
        
        assert nearly_equal(A[i3x3(i,1,2)],dutphi3_[0,1])
        
        assert nearly_equal(A[i3x3(i,2,1)],dutphi3_[1,0])
        
        assert nearly_equal(A[i3x3(i,2,2)],dutphi3_[1,1])

        drnphi3_ = drnphi3(r[0],r[1],r[2],v[0],v[1],v[2],rho_[0],rho_[1],rho_[2],mu_)
        
        assert nearly_equal(B[i3x3(i,1,0)],drnphi3_[0,0] )
        
        assert nearly_equal(B[i3x3(i,2,0)],drnphi3_[1,0])
        
        drtphi3_ = drtphi3(r[0],r[1],r[2],v[0],v[1],v[2],rho_[0],rho_[1],rho_[2],mu_)
        
        assert nearly_equal(B[i3x3(i,1,1)],drtphi3_[0,0])

        assert nearly_equal(B[i3x3(i,1,2)],drtphi3_[0,1])
        
        assert nearly_equal(B[i3x3(i,2,1)],drtphi3_[1,0])

        assert nearly_equal(B[i3x3(i,2,2)],drtphi3_[1,1])



    assert ( (n_ac - p_ac.transpose().flat[:]).all() <= finfo(double).eps )


test_GlobalAlartCurnierFun()
