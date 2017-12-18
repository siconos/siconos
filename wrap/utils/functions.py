"""Set of functions used in python tests"""

import siconos.kernel as SK
import numpy as np


def compute_dt_matrices(A, B, h, TV=False):
    # variable declaration
    t0 = 0.0  # start time
    T = 1  # end time
    n, m = B.shape

    # Matrix declaration
    x0 = np.random.random(n)
    Csurface = np.random.random((m, n))

    # Declaration of the Dynamical System
    if TV:
        process_ds = SK.FirstOrderLinearDS(x0, A)
    else:
        process_ds = SK.FirstOrderLinearTIDS(x0, A)
    # Model
    process = SK.Model(t0, T)
    process.nonSmoothDynamicalSystem().insertDynamicalSystem(process_ds)
    # time discretisation
    process_time_discretisation = SK.TimeDiscretisation(t0, h)
    # Creation of the Simulation
    process_simu = SK.TimeStepping(process_time_discretisation, 0)
    process_simu.setName("plant simulation")
    # Declaration of the integrator
    process_integrator = SK.ZeroOrderHoldOSI()
    process_simu.insertIntegrator(process_integrator)

    rel = SK.FirstOrderLinearTIR(Csurface, B)
    nslaw = SK.RelayNSL(m)
    inter = SK.Interaction(nslaw, rel)

    #process.nonSmoothDynamicalSystem().insertInteraction(inter, True)
    process.nonSmoothDynamicalSystem().link(inter, process_ds)
    process.nonSmoothDynamicalSystem().setControlProperty(inter, True)
    # Initialization
    process.setSimulation(process_simu)
    process.initialize()

    # Main loop
    process_simu.computeOneStep()
    Ad = SK.getMatrix(process_integrator.Ad(process_ds)).copy()
    Bd = SK.getMatrix(process_integrator.Bd(process_ds)).copy()

    return (Ad, Bd)


def pole_placement(A, B, P):
    """ Compute the column vector K such that the eigenvalues of A - B*F are
    the ones given by the vector P
    A is an nxn matrix, B and P are vectors
    Please note that if you want to specify complex eigenvalues, P is a matrix
    with 2 rows: the first for the real part of the eigenvalues and the second
    for the imaginary part.

    To use this function, you need a python interface to the SEVAS fortran
    routine, also known as algorithm 718 in TOMS. You can freely obtain the
    code at http://calgo.acm.org/718.gz or http://www.netlib.org/toms/718.
    Then extract the file DSEVAS.F from the main file.
    Please take a look at the license of all the algorithm on ACM ToMS
    http://toms.acm.org/AlgPolicy.html

    To generate a python interface, use f2py from numpy. A quick and dirty way
    to do it is to use the following command:

        f2py -c DSEVAS.F -m DSEVAS -l<your blas lib> --noopt

    replace <your blas lib> with your favorite blas library. Please note that
    there are some glitches with this method and then the optimisations have
    to be disabled (--noopt switch).
    """

    import DSEVAS
    n = A.shape[0]
    k = np.int(np.ceil((n**2 - 2*n + 1.0)/4))+1
    w = np.int(np.ceil((n**2 + 3*n - 4.0)/2))+1
    AA = np.array(A, order='F', dtype='f8')
    BB = np.array(B, order='F', dtype='f8')
    ieigal = np.zeros((1), dtype=int)
    rstor = np.zeros((4, k), dtype='f8', order='F')
    istor = np.zeros((k), dtype=int)
    cstor = np.zeros((w), dtype='f8')
    if P.ndim == 1:
        PP = np.zeros((2, n), dtype='f8', order='F')
        PP[0, :] = P
    else:
        PP = np.array(P, dtype='f8', order='F')
    K = np.zeros((n), dtype='f8', order='F')

    DSEVAS.dsevas(AA, BB, PP, rstor, istor, cstor, 0, 0, ieigal, K, n, n)

    errPoles = np.linalg.norm(P-np.linalg.eig(A - B*K)[0], ord=np.inf)
    if errPoles > 1e-10:
        print("Error, the poles are not placed correctly")
        print("The error is: " + str(errPoles))
        print("Desired poles: " + str(P))
        print("Obtained poles: " + str(np.linalg.eig(A - B*K)[0]))

    return (K, errPoles)
