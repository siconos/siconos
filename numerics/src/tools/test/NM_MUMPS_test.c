#include "SiconosNumerics.h"
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "NM_MUMPS.h"
#include "NM_MPI.h"
#define SIZE 2
#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include <math.h>
int main(int argc, char *argv[])
{
  double b[SIZE];

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  int rank;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
  rank = 0;
#endif

  NumericsMatrix* M = NM_create(NM_SPARSE, SIZE, SIZE);

#ifdef HAVE_MPI
  NM_MPI_set_comm(M, MPI_COMM_WORLD);
#endif
  NM_MUMPS_set_control_params(M);
  NM_triplet_alloc(M, 0);
  M->matrix2->origin = NSM_TRIPLET;

  NM_MUMPS(M, -1);
#ifdef HAVE_MPI
  if(rank > 0)
  {
    MPI_Finalize();
    NM_free(M);
    exit(0);
  }
#endif

  assert (rank == 0);
  NM_MUMPS_set_verbosity(M, 1);

  /*
     2*x - y = 1
     x   + y = 1
  */
  /* solution x: 2/3, y: 1/3 */

  NM_zentry(M, 0, 0, 2.);
  NM_zentry(M, 0, 1, -1.);
  NM_zentry(M, 1, 0, 1.);
  NM_zentry(M, 1, 1, 1.);

  for(unsigned int i=0; i<SIZE; ++i)
  {
    b[i] = 1.;
  }

  NM_MUMPS_set_problem(M, b);
  NM_MUMPS(M, 6);
  NM_MUMPS(M, -2);
  NM_MUMPS(M, 0);

  for (unsigned int i=0; i<SIZE; ++i)
  {
    printf("solution b[%u] = %g\n", i, b[i]);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  NM_free(M);

  if (fabs(b[0] - 2./3.) > 1e-7) return(1);
  if (fabs(b[1] - 1./3.) > 1e-7) return(1);

  return(0);
}
