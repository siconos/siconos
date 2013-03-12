#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "NumericsOptions.h"
#include "FrictionContact3D_Solvers.h"
#include "NonSmoothDrivers.h"
#include "fclib_interface.h"

#ifdef WITH_FCLIB
#include <fclib.h>
#endif

#ifdef WITH_FCLIB
int globalFrictionContact_fclib_write(
  GlobalFrictionContactProblem* problem, 
  char * title, 
  char * description, 
  char * math_info,
  const char *path)
{
  int rinfo = 0;

  struct fclib_global fclib_problem;
  struct fclib_info info;

  fclib_problem.spacedim = problem->dimension;
  fclib_problem.mu =  problem->mu;
  fclib_problem.b =  problem->b;
  fclib_problem.f =  problem->q;

  fclib_problem.info = &info;
  info.title = title;
  info.description = description;
  info.math_info = math_info;

  struct fclib_matrix M, H;

  fclib_problem.M = &M;
  fclib_problem.H = &H;

  /* only coordinates */
  assert(problem->M->matrix2);
  assert(problem->H->matrix2);

  SparseMatrix* _M = problem->M->matrix2;
  SparseMatrix* _H = problem->H->matrix2;

  M.m = _M->m;
  M.n = _M->n;
  M.p = _M->p;
  M.i = _M->i;
  M.x = _M->x;
  M.nzmax = _M->nzmax;
  M.nz = _M->nz;

  H.m = _H->m;
  H.n = _H->n;
  H.p = _H->p;
  H.i = _H->i;
  H.x = _H->x;
  H.nzmax = _H->nzmax;
  H.nz = _H->nz;

  rinfo = fclib_write_global(&fclib_problem, path);

  return rinfo;

}
#endif

static int fccounter = 0;

int globalFrictionContact3D_LmgcDriver(double *reaction,
                                       double *velocity,
                                       double *globalVelocity,
                                       double *q,
                                       double *b,
                                       double *mu,
                                       double *Mdata,
                                       unsigned int nzM,
                                       unsigned int *rowM,
                                       unsigned int *colM,
                                       double* Hdata,
                                       unsigned int nzH,
                                       unsigned int *rowH,
                                       unsigned int *colH,
                                       unsigned int n,
                                       unsigned int nc,
                                       int solver_id,
                                       double tolerance,
                                       int itermax,
                                       int verbose,
                                       int outputFile)
{
  GlobalFrictionContactProblem problem;
  NumericsMatrix M, H;

  SparseMatrix _M, _H;

  M.storageType = 2;
  M.size0 = n;
  M.size1 = n;

  _M.nzmax = nzM;
  _M.nz = nzM;
  _M.m = M.size0;
  _M.n = M.size1;
  _M.p = (int *) rowM;
  _M.i = (int *) colM;
  _M.x = Mdata;

  H.storageType = 2;
  H.size0 = M.size0;
  H.size1 = 3 * nc;

  _H.nzmax = nzH;
  _H.nz = nzH;
  _H.m = H.size0;
  _H.n = H.size1;
  _H.p = (int *) rowH;
  _H.i = (int *) colH;
  _H.x = Hdata;

  M.matrix2 = &_M;
  H.matrix2 = &_H;

  problem.dimension = 3;
  problem.numberOfContacts = nc;

  problem.M = &M;
  problem.H = &H;
  problem.q = q;
  problem.b = b;
  problem.mu = mu;

  NumericsOptions numerics_options; 
  numerics_options.verboseMode = verbose;

  SolverOptions numerics_solver_options;

  frictionContact3D_setDefaultSolverOptions(&numerics_solver_options, solver_id);

  numerics_solver_options.dparam[0] = tolerance;
  numerics_solver_options.iparam[0] = itermax;

  int rinfo = globalFrictionContact3D_driver(&problem,
                                             reaction, 
                                             velocity,
                                             globalVelocity,
                                             &numerics_solver_options,
                                             &numerics_options);

  if(outputFile == 1)
  {
    /* dump in C format */
  }
  else if (outputFile == 2)
  {
    /* dump in Numerics .dat format */
  }
  else if (outputFile == 3)
  {
#ifdef WITH_FCLIB
    char fname[256];
    sprintf(fname, "LMGC_GlobalFrictionContactProblem%.5d.hdf5", fccounter++);
    printf("Dump of LMGC_GlobalFrictionContactProblem%.5d.hdf5", fccounter);

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

    globalFrictionContact_fclib_write(&problem,
                                      title,
                                      description,
                                      math_info,
                                      fname);



    fclose(foutput);
#else
    printf("Fclib is not available ...\n");
#endif

  }

  return rinfo; 
}

