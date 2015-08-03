#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "NumericsOptions.h"
#include "FrictionContact3D_Solvers.h"
#include "NonSmoothDrivers.h"
#include "fclib_interface.h"

//#define DEBUG_MESSAGES 1
#include "debug.h"

#ifdef WITH_FCLIB
#include <fclib.h>
#endif

#ifdef WITH_FCLIB
int globalFrictionContact_fclib_write(
  GlobalFrictionContactProblem* problem,
  char * title,
  char * description,
  char * mathInfo,
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
  info.math_info = mathInfo;

  struct fclib_matrix M, H;

  fclib_problem.M = &M;
  fclib_problem.H = &H;

  /* only coordinates */
  assert(problem->M->matrix2);
  assert(problem->H->matrix2);

  CSparseMatrix* _M = problem->M->matrix2;
  CSparseMatrix* _H = problem->H->matrix2;

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
#ifdef WITH_FCLIB
static int fccounter = 0;
#endif

double * alloc_memory_double(unsigned int size, double *p);
double * alloc_memory_double(unsigned int size, double *p)
{
  double * r = (double *) malloc (size * sizeof(double));
  memcpy(r, p, size * sizeof(double));
  return r;
}

unsigned int * alloc_memory_int(unsigned int size, unsigned int *p);
unsigned int * alloc_memory_int(unsigned int size, unsigned int *p)
{
  unsigned int * r = (unsigned int *) malloc (size * sizeof(int));
  memcpy(r, p, size * sizeof(int));
  return r;
}



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
                                       int isize,
                                       int *iparam,
                                       int dsize,
                                       double *dparam,
                                       int verbose,
                                       int outputFile)
{
  GlobalFrictionContactProblem problem;
  NumericsMatrix M, H;

  CSparseMatrix _M, _H;

  unsigned int * _colM = alloc_memory_int(nzM, colM);
  unsigned int * _rowM = alloc_memory_int(nzM, rowM);

  unsigned int * _colH = alloc_memory_int(nzH, colH);
  unsigned int * _rowH = alloc_memory_int(nzH, rowH);


  M.matrix0 = NULL;
  M.matrix1 = NULL;
  M.matrix2 = NULL;
  M.matrix3 = NULL;
  M.matrix4 = NULL;

  M.storageType = 2; /* csc */
  M.size0 = n;
  M.size1 = n;

  _M.nzmax = nzM;
  _M.nz = nzM;
  _M.m = M.size0;
  _M.n = M.size1;
  _M.p = (int *) _colM;
  _M.i = (int *) _rowM;
  _M.x = Mdata;


  H.matrix0 = NULL;
  H.matrix1 = NULL;
  H.matrix2 = NULL;
  H.matrix3 = NULL;
  H.matrix4 = NULL;

  H.storageType = 2;
  H.size0 = M.size0;
  H.size1 = 3 * nc;

  _H.nzmax = nzH;
  _H.nz = nzH;
  _H.m = H.size0;
  _H.n = H.size1;

  _H.p = (int *) _colH;
  _H.i = (int *) _rowH;
  _H.x = Hdata;

  M.matrix2 = &_M;
  H.matrix2 = &_H;

  for (unsigned int i=0; i< _M.nz; ++i)
  {
    _M.p[i] --;
    _M.i[i] --;
  }

  for (unsigned int i=0; i< _H.nz; ++i)
  {
    _H.p[i] --;

    _H.i[i] --;

    DEBUG_PRINTF("%d -> %d,%d\n", i, _H.p[i], _H.i[i]);

  }


  problem.dimension = 3;
  problem.numberOfContacts = nc;

  problem.M = &M;
  problem.H = &H;
  problem.q = q;
  problem.b = b;
  problem.mu = mu;

  NumericsOptions numerics_options;
  setDefaultNumericsOptions(&numerics_options);
  numerics_options.verboseMode = verbose;

  SolverOptions numerics_solver_options;

  globalFrictionContact3D_setDefaultSolverOptions(&numerics_solver_options, solver_id);

  assert(isize <= numerics_solver_options.iSize);
  assert(dsize <= numerics_solver_options.dSize);

  for (unsigned int i=0; i<isize; ++i)
    numerics_solver_options.iparam[i] = iparam[i];

  for (unsigned int i=0; i<dsize; ++i)
    numerics_solver_options.dparam[i] = dparam[i];

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
    char * mathInfo = (char *)malloc(n * sizeof(char *));
    strcpy(mathInfo,  "unknown");

    globalFrictionContact_fclib_write(&problem,
                                      title,
                                      description,
                                      mathInfo,
                                      fname);


    fclose(foutput);
#else
    printf("Fclib is not available ...\n");
#endif

  }

  if (M.matrix3)
  {
    cs_free(M.matrix3);
    M.matrix3 = NULL;
  }

  if (M.matrix4)
  {
    cs_free(M.matrix4);
    M.matrix4 = NULL;
  }

  if (H.matrix3)
  {
    cs_free(H.matrix3);
    H.matrix3 = NULL;
  }

  if (H.matrix4)
  {
    cs_free(H.matrix4);
    H.matrix4 = NULL;
  }

  free(_colM);
  free(_colH);

  free(_rowM);
  free(_rowH);

  return rinfo;
}

