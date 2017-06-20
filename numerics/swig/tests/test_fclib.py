import siconos.numerics as sn
import os.path
import scipy.sparse
import scipy.sparse.linalg
import numpy as np
from siconos.tests_setup import working_dir


# this could be changed at install time
data_dir = working_dir + '/data/'
# solvers = [sn.SICONOS_GLOBAL_FRICTION_3D_NSGS, sn.SICONOS_GLOBAL_FRICTION_3D_NSN_AC]
solvers = (sn.SICONOS_GLOBAL_FRICTION_3D_NSGS,)
solvers_reduced1 = (sn.SICONOS_FRICTION_3D_NSGS, sn.SICONOS_FRICTION_3D_NSN_AC)
solvers_reduced2 = (sn.SICONOS_FRICTION_3D_NSN_AC,)  # sn.SICONOS_FRICTION_3D_NSN_FB)
solvers_reduced3 = (sn.SICONOS_FRICTION_3D_NSGS,)

def condensed_from_global(fcp):
    # spsolve expect the indices to be cint aka 32 bits int
    # Hence, we do some magic to cirumvent those issues.
    Mcoo = scipy.sparse.coo_matrix(fcp.M)
    Mcsc = scipy.sparse.csc_matrix(Mcoo)
    Hcoo = scipy.sparse.coo_matrix(fcp.H)
    Hcsc = scipy.sparse.csc_matrix(Hcoo)
    WW = Hcsc.T.dot(scipy.sparse.linalg.spsolve(Mcsc, Hcsc))
    qprime = fcp.H.T.dot(scipy.sparse.linalg.spsolve(Mcsc, fcp.q))

    fcp_reduced = sn.FrictionContactProblem()
    fcp_reduced.dimension = fcp.dimension
    fcp_reduced.numberOfContacts = fcp.numberOfContacts

    # this is a hack to deal with the inability of fc3d solvers to work with
    # sparse matrices
    _, Wsbm = sn.sparseToSBM(3, WW)
    fcp_reduced.M = Wsbm
    fcp_reduced.mu = fcp.mu
    fcp_reduced.q = fcp.b + qprime
    return fcp_reduced


def solve_reduced(fcp, solver_reduced):
    SO_reduced = sn.SolverOptions(solver_reduced)
    SO_reduced.iparam[0] =100000
    SO_reduced.dparam[0] = np.sqrt(fcp.numberOfContacts)*1e-9
    size_reaction = fcp.numberOfContacts * 3
    reaction_reduced = np.zeros((size_reaction,))
    velocities_reduced = np.zeros((size_reaction,))

    return sn.fc3d_driver(fcp, reaction_reduced, velocities_reduced, SO_reduced)


def solve_global(fcp, solver):
    SO = sn.SolverOptions(solver)
    SO.dparam[0] = np.sqrt(fcp.numberOfContacts)*1e-9

    n, m = fcp.H.shape
    size_reaction = m
    reaction = np.zeros((size_reaction,))
    velocities = np.zeros((size_reaction,))
    global_velocities = np.zeros((n,))
    return sn.gfc3d_driver(fcp, reaction, velocities, global_velocities, SO)


def test_gfc3d():
    data_files = ('CubeH8.hdf5', 'LMGC_GlobalFrictionContactProblem00046.hdf5')
    mark_as_failed = False

    for d in data_files:
        full_path = data_dir + d
        if os.path.isfile(full_path):

            fcp = sn.globalFrictionContact_fclib_read(full_path)
            for s in solvers:
                res = solve_global(fcp, s)
                if res:
                    print('Solver {:} on problem {:} failed with info = {:}'.format(sn.solver_options_id_to_name(s), d, res))
                    mark_as_failed = True
                else:
                    print('Solver {:} on problem {:} is ok'.format(sn.solver_options_id_to_name(s), d))

            fcp_reduced = condensed_from_global(fcp)
            for s in solvers_reduced1:
                res = solve_reduced(fcp_reduced, s)
                if res:
                    print('Solver {:} on problem {:} in reduced form failed with info = {:}'.format(sn.solver_options_id_to_name(s), d, res))
                    mark_as_failed = True
                else:
                    print('Solver {:} on problem {:} is ok'.format(sn.solver_options_id_to_name(s), d))

    assert mark_as_failed is False


def test_fc3d():
    data_files = ['Capsules-i125-1213.hdf5']
    mark_as_failed = False

    for d in data_files:
        full_path = data_dir + d
        if os.path.isfile(full_path):
            sn.numerics_set_verbose(1)
            fcp = sn.frictionContact_fclib_read(full_path)
            for s in solvers_reduced3:
                res = solve_reduced(fcp, s)
                if res:
                    print('Solver {:} on problem {:} failed with info = {:}'.format(sn.solver_options_id_to_name(s), d, res))
                    mark_as_failed = True
                else:
                    print('Solver {:} on problem {:} is ok'.format(sn.solver_options_id_to_name(s), d))

    assert mark_as_failed is False
