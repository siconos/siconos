#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "NumericsOptions.h"
#include "FrictionContact3D_Solvers.h"
#include "NonSmoothDrivers.h"
#include "fclib_interface.h"

#define DEBUG_MESSAGES 1
#define DEBUG_STDOUT
#include "debug.h"

static int fccounter =0;
static double * alloc_memory_double(unsigned int size, double *p)
{
  double * r = (double *) malloc (size * sizeof(double));
  memcpy(r, p, size * sizeof(double));
  return r;
}

static csi * alloc_memory_csi(unsigned int size, unsigned int *p)
{
  csi * r = (csi *) malloc (size * sizeof(csi));
  for(unsigned int i=0; i<size; ++i)
  {
    r[i] = (csi) p[i];
  }
  return r;
}

int globalFrictionContact_fclib_write(
  GlobalFrictionContactProblem* problem,
  char * title,
  char * description,
  char * mathInfo,
  const char *path);

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

  /* NumericsMatrix M, H; */
  NumericsMatrix * M =newNumericsMatrix();
  M->storageType = 2; /* csc */
  M->size0 = n;
  M->size1 = n;


  NumericsMatrix * H =newNumericsMatrix();
  H->storageType = 2;
  H->size0 = M->size0;
  H->size1 = 3 * nc;

  NumericsSparseMatrix * SM =newNumericsSparseMatrix();
  M->matrix2 = SM;
  SM->triplet =   (CSparseMatrix * )malloc(sizeof(CSparseMatrix));
  CSparseMatrix * _M = SM->triplet;

  csi * _colM = alloc_memory_csi(nzM, colM);
  csi * _rowM = alloc_memory_csi(nzM, rowM);

  _M->nzmax = nzM;
  _M->nz = nzM;
  _M->m = M->size0;
  _M->n = M->size1;
  _M->p = (csi *) _colM;
  _M->i = (csi *) _rowM;
  double * _Mdata = alloc_memory_double(nzM, Mdata);
  _M->x = _Mdata;

  DEBUG_PRINTF("_M->n=%li\t",_M->n);
  DEBUG_PRINTF("_M->m=%li\n",_M->m);

  NumericsSparseMatrix * SH =newNumericsSparseMatrix();
  H->matrix2 = SH;
  SH->triplet =   (CSparseMatrix * )malloc(sizeof(CSparseMatrix));
  CSparseMatrix * _H = SH->triplet;

  csi * _colH = alloc_memory_csi(nzH, colH);
  csi * _rowH = alloc_memory_csi(nzH, rowH);

  _H->nzmax = nzH;
  _H->nz = nzH;
  _H->m = H->size0;
  _H->n = H->size1;

  _H->p = (csi *) _colH;
  _H->i = (csi *) _rowH;
  double * _Hdata = alloc_memory_double(nzH, Hdata);
  _H->x = _Hdata;

  for (int i=0; i< _M->nz; ++i)
  {
    _M->p[i] --;
    _M->i[i] --;
    /* DEBUG_PRINTF("%d -> %d,%d\n", i, _M->p[i], _M->i[i]); */

  }

  for (int i=0; i< _H->nz; ++i)
  {
    _H->p[i] --;
    _H->i[i] --;
    /* DEBUG_PRINTF("%d -> %d,%d\n", i, _H->p[i], _H->i[i]); */
  }

  GlobalFrictionContactProblem * problem =(GlobalFrictionContactProblem*)malloc(sizeof(GlobalFrictionContactProblem));

  problem->dimension = 3;
  problem->numberOfContacts = nc;

  problem->M = M;
  problem->H = H;
  problem->q = q;
  problem->b = b;
  problem->mu = mu;

  NumericsOptions numerics_options;
  setDefaultNumericsOptions(&numerics_options);
  numerics_options.verboseMode = verbose;

  SolverOptions numerics_solver_options;

  globalFrictionContact3D_setDefaultSolverOptions(&numerics_solver_options, solver_id);

  assert(isize <= numerics_solver_options.iSize);
  assert(dsize <= numerics_solver_options.dSize);
  /* for (int i=0; i<isize; ++i) */
  /*   numerics_solver_options.iparam[i] = iparam[i]; */

  /* for (int i=0; i<dsize; ++i) */
  /*   numerics_solver_options.dparam[i] = dparam[i]; */

  printSolverOptions(&numerics_solver_options);
  FILE * file  =  fopen("toto.dat", "w");
  globalFrictionContact_printInFile(problem, file);
  fclose(file);
  int rinfo =  globalFrictionContact3D_driver(problem,
                                             reaction,
                                             velocity,
                                             globalVelocity,
                                             &numerics_solver_options,
                                             &numerics_options);

  FILE * file1  =  fopen("tutu.dat", "w");
  globalFrictionContact_printInFile(problem, file1);
  fclose(file1);
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

    globalFrictionContact_fclib_write(problem,
                                      title,
                                      description,
                                      mathInfo,
                                      fname);


    fclose(foutput);
#else
    printf("Fclib is not available ...\n");
#endif

  }


  freeNumericsMatrix(M);
  freeNumericsMatrix(H);
  free(problem);

  /* free(_colM); */
  /* free(_colH); */

  /* free(_rowM); */
  /* free(_rowH); */

  return rinfo;
}


/* int globalFrictionContact3D_LmgcDriver_SBM(double *reaction, */
/*                                        double *velocity, */
/*                                        double *globalVelocity, */
/*                                        double *q, */
/*                                        double *b, */
/*                                        double *mu, */
/*                                        double *Mdata, */
/*                                        unsigned int nzM, */
/*                                        unsigned int *rowM, */
/*                                        unsigned int *colM, */
/*                                        double* Hdata, */
/*                                        unsigned int nzH, */
/*                                        unsigned int *rowH, */
/*                                        unsigned int *colH, */
/*                                        unsigned int n, */
/*                                        unsigned int nc, */
/*                                        int solver_id, */
/*                                        int isize, */
/*                                        int *iparam, */
/*                                        int dsize, */
/*                                        double *dparam, */
/*                                        int verbose, */
/*                                        int outputFile) */
/* { */
/*   GlobalFrictionContactProblem problem; */
/*   NumericsMatrix M, H; */

/*   CSparseMatrix _M, _H; */

/*   unsigned int * _colM = alloc_memory_int(nzM, colM); */
/*   unsigned int * _rowM = alloc_memory_int(nzM, rowM); */

/*   unsigned int * _colH = alloc_memory_int(nzH, colH); */
/*   unsigned int * _rowH = alloc_memory_int(nzH, rowH); */


/*   M.matrix0 = NULL; */
/*   M.matrix1 = NULL; */
/*   M.matrix2 = NULL; */
/*   M.internalData = NULL; */

/*   M.storageType = 2; /\* csc *\/ */
/*   M.size0 = n; */
/*   M.size1 = n; */

/*   _M.nzmax = nzM; */
/*   _M.nz = nzM; */
/*   _M.m = M.size0; */
/*   _M.n = M.size1; */
/*   _M.p = (int *) _colM; */
/*   _M.i = (int *) _rowM; */
/*   _M.x = Mdata; */


/*   H.matrix0 = NULL; */
/*   H.matrix1 = NULL; */
/*   H.matrix2 = NULL; */
/*   H.internalData = NULL; */

/*   H.storageType = 2; */
/*   H.size0 = M.size0; */
/*   H.size1 = 3 * nc; */

/*   _H.nzmax = nzH; */
/*   _H.nz = nzH; */
/*   _H.m = H.size0; */
/*   _H.n = H.size1; */

/*   _H.p = (int *) _colH; */
/*   _H.i = (int *) _rowH; */
/*   _H.x = Hdata; */

/*   M.matrix2->triplet = &_M; */
/*   H.matrix2->triplet = &_H; */

/*   for (int i=0; i< _M.nz; ++i) */
/*   { */
/*     _M.p[i] --; */
/*     _M.i[i] --; */
/*   } */

/*   for (int i=0; i< _H.nz; ++i) */
/*   { */
/*     _H.p[i] --; */

/*     _H.i[i] --; */

/*     DEBUG_PRINTF("%d -> %d,%d\n", i, _H.p[i], _H.i[i]); */

/*   } */


/*   problem.dimension = 3; */
/*   problem.numberOfContacts = nc; */

/*   problem.M = &M; */
/*   problem.H = &H; */
/*   problem.q = q; */
/*   problem.b = b; */
/*   problem.mu = mu; */

/*   NumericsOptions numerics_options; */
/*   setDefaultNumericsOptions(&numerics_options); */
/*   numerics_options.verboseMode = verbose; */

/*   SolverOptions numerics_solver_options; */

/*   globalFrictionContact3D_setDefaultSolverOptions(&numerics_solver_options, solver_id); */

/*   assert(isize <= numerics_solver_options.iSize); */
/*   assert(dsize <= numerics_solver_options.dSize); */

/*   for (int i=0; i<isize; ++i) */
/*     numerics_solver_options.iparam[i] = iparam[i]; */

/*   for (int i=0; i<dsize; ++i) */
/*     numerics_solver_options.dparam[i] = dparam[i]; */

/*   int rinfo = globalFrictionContact3D_driver(&problem, */
/*                                              reaction, */
/*                                              velocity, */
/*                                              globalVelocity, */
/*                                              &numerics_solver_options, */
/*                                              &numerics_options); */

/*   if(outputFile == 1) */
/*   { */
/*     /\* dump in C format *\/ */
/*   } */
/*   else if (outputFile == 2) */
/*   { */
/*     /\* dump in Numerics .dat format *\/ */
/*   } */
/*   else if (outputFile == 3) */
/*   { */
/* #ifdef WITH_FCLIB */
/*     char fname[256]; */
/*     sprintf(fname, "LMGC_GlobalFrictionContactProblem%.5d.hdf5", fccounter++); */
/*     printf("Dump of LMGC_GlobalFrictionContactProblem%.5d.hdf5", fccounter); */

/*     FILE * foutput  =  fopen(fname, "w"); */
/*     int n = 100; */
/*     char * title = (char *)malloc(n * sizeof(char *)); */
/*     strcpy(title, "LMGC dump in hdf5"); */
/*     char * description = (char *)malloc(n * sizeof(char *)); */

/*     strcat(description, "Rewriting in hdf5 through siconos of  "); */
/*     strcat(description, fname); */
/*     strcat(description, " in FCLIB format"); */
/*     char * mathInfo = (char *)malloc(n * sizeof(char *)); */
/*     strcpy(mathInfo,  "unknown"); */

/*     globalFrictionContact_fclib_write(&problem, */
/*                                       title, */
/*                                       description, */
/*                                       mathInfo, */
/*                                       fname); */


/*     fclose(foutput); */
/* #else */
/*     printf("Fclib is not available ...\n"); */
/* #endif */

/*   } */


/*   freeNumericsMatrix(&M); */
/*   freeNumericsMatrix(&H); */

/*   free(_colM); */
/*   free(_colH); */

/*   free(_rowM); */
/*   free(_rowH); */

/*   return rinfo; */
/* } */



/* int globalFrictionContact3D_LmgcDriver(double *reaction, */
/*                                        double *velocity, */
/*                                        double *globalVelocity, */
/*                                        double *q, */
/*                                        double *b, */
/*                                        double *mu, */
/*                                        double *Mdata, */
/*                                        unsigned int nzM, */
/*                                        unsigned int *rowM, */
/*                                        unsigned int *colM, */
/*                                        double* Hdata, */
/*                                        unsigned int nzH, */
/*                                        unsigned int *rowH, */
/*                                        unsigned int *colH, */
/*                                        unsigned int n, */
/*                                        unsigned int nc, */
/*                                        int solver_id, */
/*                                        int isize, */
/*                                        int *iparam, */
/*                                        int dsize, */
/*                                        double *dparam, */
/*                                        int verbose, */
/*                                        int outputFile) */
/* { */
/*   return globalFrictionContact3D_LmgcDriver_CSC(reaction, */
/*                                                 velocity, */
/*                                                 globalVelocity, */
/*                                                 q, */
/*                                                 b, */
/*                                                 mu, */
/*                                                 Mdata, */
/*                                                 nzM, */
/*                                                 rowM, */
/*                                                 colM, */
/*                                                 Hdata, */
/*                                                 nzH, */
/*                                                 rowH, */
/*                                                 colH, */
/*                                                 n, */
/*                                                 nc, */
/*                                                 solver_id, */
/*                                                 isize, */
/*                                                 iparam, */
/*                                                 dsize, */
/*                                                 dparam, */
/*                                                 verbose, */
/*                                                 outputFile); */

/* } */
