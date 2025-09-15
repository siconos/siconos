#include "NM_MUMPS.h"

#include "NM_MPI.h"
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "SiconosNumerics.h"
#define SIZE 3
#include <stdio.h>
#include <stdlib.h>
#ifdef SICONOS_HAS_MPI
#include <mpi.h>
#endif
#include <math.h>
int main(int argc, char* argv[]) {
  int rval = 0;
  double b[SIZE * 2];

#ifdef SICONOS_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  int rank;

#ifdef SICONOS_HAS_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  rank = 0;
#endif
  NumericsMatrix* M;

  /* unsymmetric test */
  M = NM_create(NM_SPARSE, 2, 2);

#ifdef SICONOS_HAS_MPI
  NM_MPI_set_comm(M, MPI_COMM_WORLD);
#endif
  NM_MUMPS_set_control_params(M);
  NM_triplet_alloc(M, 0);
  M->matrix2->origin = NSM_TRIPLET;

  NM_MUMPS(M, -1);

#ifdef SICONOS_HAS_MPI
  if (rank == 0)
#endif
  {
    assert(rank == 0);
    NM_MUMPS_set_verbosity(M, 1);

    /*
      2*x - y = 1
      x   + y = 1
    */
    /* solution x: 2/3, y: 1/3 */

    NM_entry(M, 0, 0, 2.);
    NM_entry(M, 0, 1, -1.);
    NM_entry(M, 1, 0, 1.);
    NM_entry(M, 1, 1, 1.);

    /* solution x: 2/3, y: 1/3 */
    b[0] = 1.;
    b[1] = 1.;

    /* solution x: 1/3, y: -1/3 */
    b[2] = 1.;
    b[3] = 0.;

    NM_MUMPS_set_matrix(M);
    NM_MUMPS_set_dense_rhs(M, 2, b);
    NM_MUMPS(M, 6);
    NM_MUMPS(M, -2);
    NM_MUMPS(M, 0);

    for (unsigned int i = 0; i < 4; ++i) {
      printf("solution b[%u] = %g\n", i, b[i]);
    }

    rval += (fabs(b[0] - 2. / 3.) > 1e-7);
    rval = rval << (fabs(b[1] - 1. / 3.) > 1e-7);

    rval = rval << (fabs(b[2] - 1. / 3.) > 1e-7);
    rval = rval << (fabs(b[3] + 1. / 3.) > 1e-7);
  }

  NM_clear(M);

  /* symmetric test */
  M = NM_create(NM_SPARSE, 3, 3);

#ifdef SICONOS_HAS_MPI
  NM_MPI_set_comm(M, MPI_COMM_WORLD);
#endif
  NM_MUMPS_set_control_params(M);
  NM_triplet_alloc(M, 0);
  M->matrix2->origin = NSM_TRIPLET;
  NM_MUMPS_set_sym(M, 1);

  NM_MUMPS(M, -1);

#ifdef SICONOS_HAS_MPI
  if (rank == 0)
#endif
  {
    assert(rank == 0);
    NM_MUMPS_set_verbosity(M, 1);

    /*
      M=[[1, 0, 0],
         [0, 1, .5],
         [0, .5, 1]]

      b=[1, 1, 1]
    */

    NM_entry(M, 0, 0, 1.);
    NM_entry(M, 1, 1, 1.);
    NM_entry(M, 2, 1, 0.5);
    /*    NM_entry(M, 1, 2, 0.5);*/
    NM_entry(M, 2, 2, 1.);

    /*    NM_display(M);*/

    /* solution: [1, 2/3, 2/3] */
    b[0] = 1.;
    b[1] = 1.;
    b[2] = 1.;

    /* solution: [1, -2/3, 5/3] */
    b[3] = 1.;
    b[4] = 0.;
    b[5] = 1.;

    NM_MUMPS_set_matrix(M);
    NM_MUMPS_set_dense_rhs(M, 2, b);
    NM_MUMPS(M, 6);
    NM_MUMPS(M, -2);
    NM_MUMPS(M, 0);

    for (unsigned int i = 0; i < 6; ++i) {
      printf("solution b[%u] = %g\n", i, b[i]);
    }

    rval = rval << (fabs(b[0] - 1.) > 1e-7);
    rval = rval << (fabs(b[1] - 2. / 3.) > 1e-7);
    rval = rval << (fabs(b[2] - 2. / 3.) > 1e-7);

    rval = rval << (fabs(b[3] - 1) > 1e-7);
    rval = rval << (fabs(b[3] + 2. / 3.) > 1e-7);
    rval = rval << (fabs(b[3] - 5. / 3.) > 1e-7);
  }

  NM_clear(M);
  free(M);

#ifdef SICONOS_HAS_MPI
  MPI_Finalize();
#endif

  return (rval);
}
