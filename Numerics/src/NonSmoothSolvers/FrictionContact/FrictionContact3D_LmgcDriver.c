#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "LA.h"
#include "NumericsOptions.h"
#include "FrictionContact3D_Solvers.h"
#include "NonSmoothDrivers.h"
#include "fclib_interface.h"

static int fccounter = 0;

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
                                 int itermax,
                                 int verbose,
                                 int outputFile)
{

  SparseBlockCoordinateMatrix* MC = newSparseBlockCoordinateMatrix3x3fortran(nc, nc, nb, row, column, W);

  SparseBlockStructuredMatrix* M = SBCMToSBM(MC);

  NumericsMatrix* NM = newSparseNumericsMatrix(nc * 3, nc * 3, M);

  FrictionContactProblem* FC = frictionContactProblem_new(3, nc, NM, q, mu);

  /* frictionContact_display(FC); */

  NumericsOptions numerics_options;
  numerics_options.verboseMode = verbose; // turn verbose mode to off by default

  //  uncomment to save FrictionContactProblem
  if (outputFile == 1)
  {
    char fname[256];
    sprintf(fname, "LMGC_FrictionContactProblem%.5d.dat", fccounter++);
    printf("Dump of LMGC_FrictionContactProblem%.5d.dat", fccounter);

    FILE * foutput  =  fopen(fname, "w");
    frictionContact_printInFile(FC, foutput);
    fclose(foutput);
  }
  else if (outputFile == 2)
  {
#ifdef WITH_FCLIB
    char fname[256];
    sprintf(fname, "LMGC_FrictionContactProblem%.5d.hdf5", fccounter++);
    printf("Dump of LMGC_FrictionContactProblem%.5d.hdf5", fccounter);

    FILE * foutput  =  fopen(fname, "w");
    int n = 100;
    char * title = (char *)malloc(n * sizeof(char *));
    strcpy(title, "LMGC dump in hdf5");
    char * description = (char *)malloc(n * sizeof(char *));

    strcat(description, "Rewriting in hdf5 through siconos of  ");
    strcat(description, fname);
    strcat(description, " in FCLIB format");
    char * math_info = (char *)malloc(n * sizeof(char *));
    strcpy(math_info,  "unknown");

    frictionContact_fclib_write(FC,
                                title,
                                description,
                                math_info,
                                fname);



    fclose(foutput);
#else
    printf("Fclib is not available ...\n");
#endif

  }

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
