#!/usr/bin/env python

from numpy import *

# import Siconos.Numerics * fails with py.test!
import Siconos.Numerics as N


NC = 1

M = eye(3*NC)

q = array([-1., 1., 3.])

mu = array([0.1]);

z = array([0.,0.,0.])

reactions = array([0.,0.,0.])

velocities = array([0.,0.,0.])


# gp linesearch
def lineSearch(problem, reaction, rho, delta_reaction, maxiter):

    reaction0=reaction

    alphamin = 0.
    alphamax = inf
    alpha = 1

    m1 = 0.01
    m2 = 0.99

    W = matrix(problem.M)

    velocity = W * reaction + problem.q

    F,A,B = N.frictionContact3D_globalAlartCurnierFunction(reaction, velocity, problem.mu, rho)

    q0 = 0.5 * inner(F.transpose(),F.transpose())

    AWpB = N.computeAWpB(A,W,B)

    dqdt0 = matrix(F.transpose()) * (matrix(AWpB) * delta_reaction)

    iter=0

    while (iter < maxiter):

        iter += 1

        reaction = reaction0 + alpha*delta_reaction 
        
        velocity = W * reaction + problem.q
        

        F,A,B = N.frictionContact3D_globalAlartCurnierFunction(reaction, velocity, problem.mu, rho)

        q = 0.5 * inner(F.transpose(), F.transpose())

        slope = (q-q0) / alpha

        C1 = (slope >= m2*dqdt0)[0,0]
        C2 = (slope <= m1*dqdt0)[0,0]

        if (C1 and C2):
            return alpha

        else:
            if not C1:
                alphamin = alpha
            else:
                alphamax = alpha

            if (alpha < inf):
                alpha = 0.5 *(alphamin + alphamax)

            else:
                alpha = alphamin
            
    return alpha

def globalAlartCurnier(problem, reaction, velocity, options):

    problemSize = 3 * problem.numberOfContacts
    itermax = options.iparam[0]
    erritermax = options.iparam[1]
    rho = matrix(ones(problemSize)).transpose()
    
    tolerance = options.dparam[0]

    W = problem.M
    
    velocity = W*reaction + problem.q

    iter = 0
    error = inf

    while (iter < itermax and error > tolerance):

        F,A,B=N.frictionContact3D_globalAlartCurnierFunction(
            reaction,
            velocity,
            problem.mu,
            rho)

        AWpB = N.computeAWpB(A,W,B)

        AA=matrix(AWpB)

        delta_reaction = linalg.solve(AA,-F)

        alpha = lineSearch(problem, reaction, rho, delta_reaction, 100)

        if alpha is None:
            reaction += delta_reaction
        else:
            reaction += alpha*delta_reaction

        velocity = W*reaction + problem.q

        error = N.FrictionContact3D_compute_error(problem, reaction, velocity,
                                                tolerance, options)

        iter += 1

    return reaction, velocity, iter, error



def test_fc3dnsgs():
    N.setNumericsVerbose(2)
    FCP = N.FrictionContactProblem(3,M,q,mu)
    SO=N.SolverOptions(N.SICONOS_FRICTION_3D_NSGS)
    r=N.frictionContact3D_nsgs(FCP, reactions, velocities, SO)
    assert SO.dparam[1] < 1e-10
    assert not r 


def test_fc3dglobalac_simple_case():
    N.setNumericsVerbose(2)
    FCP = N.FrictionContactProblem(3,M,q,mu)
    SO=N.SolverOptions(N.SICONOS_FRICTION_3D_GLOBALAC)
    
    r=N.frictionContact3D_globalAlartCurnier(FCP, reactions, velocities, SO)
    assert SO.dparam[1] < 1e-10
    assert not r

def test_fc3dglobalac_full():
    N.setNumericsVerbose(0)
    problem = N.frictionContactProblemFromFile("Rover1039.dat")
    print 'problem dimension :',problem.dimension
    print 'number of contacts : ', problem.numberOfContacts

    SO=N.SolverOptions(N.SICONOS_FRICTION_3D_GLOBALAC)

    SO.iparam[1]=1

    reactions = matrix(zeros(problem.dimension*problem.numberOfContacts)).transpose()
    velocities = matrix(zeros(problem.dimension*problem.numberOfContacts)).transpose()

    r,v,iter,err=globalAlartCurnier(problem, reactions, velocities, SO)

    reactions = matrix(zeros(problem.dimension*problem.numberOfContacts)).transpose()
    velocities = matrix(zeros(problem.dimension*problem.numberOfContacts)).transpose()

    i=N.frictionContact3D_globalAlartCurnier(problem,reactions,velocities, SO)

    print 'iter=', iter
    print 'iter=', SO.iparam[7]
    assert (SO.dparam[1] - err) < 1e-7
