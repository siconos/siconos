#!/usr/bin/env python3

# kept alone as we may want to run it with mpirun

def test_nm_mumps():

    # this may be run with mpirun -np <nb processes>
    from siconos import numerics
    from mpi4py import MPI
    import scipy.sparse
    import numpy

    comm = MPI.COMM_WORLD

    # build an empty sparse matrix
    M = numerics.NumericsMatrix(scipy.sparse.eye(0))

    numerics.NM_MPI_set_comm(M, comm)
    numerics.NM_MUMPS_set_control_params(M)

    # start parallel MUMPS, processes with rank > 0 execute NM_MUMPS(M, -1)
    # then all next NM_MUMPS(M, x) with x!=0 then return when process with
    # rank 0 executes NM_MUMPS(M, 0)
    numerics.NM_MUMPS(M, -1)

    if (comm.Get_rank() > 0):
        # when rank 0 executes NM_MUMPS(M, 0)
        exit(0)

    # only on rank 0:

    # fill the matrix
    #     2*x - y = 1
    #     x   + y = 1
    # solution x: 2/3, y: 1/3 */
    numerics.NM_zentry(M, 0, 0, 2.)
    numerics.NM_zentry(M, 0, 1, -1.)
    numerics.NM_zentry(M, 1, 0, 1.)
    numerics.NM_zentry(M, 1, 1, 1.)

    b = numpy.array([1.0, 1.0])
    numerics.NM_MUMPS_set_problem(M, b)

    # analysis, factorization, solve.
    numerics.NM_MUMPS(M, 6)

    # end of the work with this matrix
    numerics.NM_MUMPS(M, -2)

    # send end of listening for ranks > 0
    numerics.NM_MUMPS(M, 0)

    print('solution:', b)

    assert(abs(b[0] - 2./3.) < 1e-7)
    assert(abs(b[1] - 1./3.) < 1e-7)
