#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "LA.h"
#include "NumericsOptions.h"
#include "FrictionContact3D_Solvers.h"
#include "NonSmoothDrivers.h"


int frictionContact3D_LmgcDriver(double *reaction,
                                 double *velocity,
                                 double *q,
                                 double *mu,
                                 double* W,
                                 unsigned int *row,
                                 unsigned int *column,
                                 unsigned int nc,
                                 unsigned int nb,
                                 int solver_id,
                                 double tolerance,
                                 int itermax)
{

  SparseBlockCoordinateMatrix* MC = newSparseBlockCoordinateMatrix3x3fortran(nc, nc, nb, row, column, W);

  SparseBlockStructuredMatrix* M = SBCMToSBM(MC);

  NumericsMatrix* NM = newSparseNumericsMatrix(nc * 3, nc * 3, M);

  FrictionContactProblem* FC = frictionContactProblem_new(3, nc, NM, q, mu);

  frictionContact_display(FC);

  NumericsOptions numerics_options;
  numerics_options.verboseMode = 0; // turn verbose mode to off by default

  //  uncomment to save FrictionContactProblem
  //  FILE * ff =  fopen("LMGC.dat", "w");
  //  frictionContact_printInFile(FC, ff);
  //  fclose(ff);

  SolverOptions numerics_solver_options;

  frictionContact3D_setDefaultSolverOptions(&numerics_solver_options, solver_id);

  numerics_solver_options.dparam[0] = tolerance;
  numerics_solver_options.iparam[0] = itermax;

  int info = frictionContact3D_driver(FC,
                                      reaction , velocity,
                                      &numerics_solver_options, &numerics_options);


  freeSparseBlockCoordinateMatrix3x3fortran(MC);

  free(M->index1_data);
  free(M->index2_data);
  free(M->block);
  free(M);
  free(FC);
  deleteSolverOptions(&numerics_solver_options);
  free(NM);
  free(MC);

  return info;
}
