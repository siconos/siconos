#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "csparse.h"
#include "fc3d_Solvers.h"
#include "NonSmoothDrivers.h"
#include "fclib_interface.h"
#include "GlobalFrictionContactProblem.h"
#include "gfc3d_Solvers.h"
#include "NumericsSparseMatrix.h"
#include "numerics_verbose.h"
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "debug.h"

#ifdef WITH_FCLIB
static int gfccounter =-1;

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#endif

//#define USE_NM_DENSE


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

int gfc3d_LmgcDriver(double *reaction,
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
                     int verbose_in,
                     int outputFile,
                     int freq_output)
{


  verbose = verbose_in;
  
  /* NumericsMatrix M, H; */
  NumericsMatrix * M =NM_new();
  M->storageType = 2; /* sparse */
  M->size0 = n;
  M->size1 = n;


  NumericsMatrix * H =NM_new();
  H->storageType = 2;
  H->size0 = M->size0;
  H->size1 = 3 * nc;

  NumericsSparseMatrix * SM =newNumericsSparseMatrix();
  M->matrix2 = SM;
  SM->triplet =   (CSparseMatrix * )malloc(sizeof(CSparseMatrix));
  CSparseMatrix * _M = SM->triplet;
  SM->origin = NS_TRIPLET;

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

  DEBUG_PRINTF("_M->n=%lli\t",_M->n);
  DEBUG_PRINTF("_M->m=%lli\n",_M->m);

  NumericsSparseMatrix * SH =newNumericsSparseMatrix();
  H->matrix2 = SH;
  SH->triplet =   (CSparseMatrix * )malloc(sizeof(CSparseMatrix));
  CSparseMatrix * _H = SH->triplet;
  SH->origin = NS_TRIPLET;

  csi * _colH = alloc_memory_csi(nzH, colH);
  csi * _rowH = alloc_memory_csi(nzH, rowH);

  _H->nzmax = nzH;
  _H->nz = nzH;
  _H->m = H->size0;
  _H->n = H->size1;

  _H->p = _colH;
  _H->i =  _rowH;
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



#ifdef USE_NM_DENSE
  assert(M);
  assert(H);

  NumericsMatrix *MMtmp = NM_new();
  NumericsMatrix *HHtmp = NM_new();
  
  NM_copy(M,MMtmp);
  NM_copy(H,HHtmp);

  NM_clearSparse(M);
  NM_clearSparse(H);

  M = NM_create(NM_DENSE, H->size0, H->size0);
  H = NM_create(NM_DENSE, H->size0, H->size1);
 
  NM_to_dense(MMtmp,M);
  NM_to_dense(HHtmp,H);

  /* NM_display(M); */
  /* NM_display(H); */

#endif


  
  GlobalFrictionContactProblem * problem =(GlobalFrictionContactProblem*)malloc(sizeof(GlobalFrictionContactProblem));

  problem->dimension = 3;
  problem->numberOfContacts = nc;
  problem->env = NULL;
  problem->workspace = NULL;

  problem->M = M;
  problem->H = H;
  problem->q = q;
  problem->b = b;
  problem->mu = mu;

  SolverOptions numerics_solver_options;
  
  int infi = gfc3d_setDefaultSolverOptions(&numerics_solver_options, solver_id);
  assert(!infi);
  int iSize_min = isize < numerics_solver_options.iSize ? isize : numerics_solver_options.iSize;
  DEBUG_PRINTF("iSize_min = %i", iSize_min);
  for (int i = 0; i < iSize_min; ++i) 
    numerics_solver_options.iparam[i] = iparam[i];

  int dSize_min = dsize <  numerics_solver_options.dSize ? dsize : numerics_solver_options.dSize;
  for (int i=0; i < dSize_min; ++i)
    numerics_solver_options.dparam[i] = dparam[i];

  /* solver_options_print(&numerics_solver_options); */
  /* FILE * file  =  fopen("toto.dat", "w"); */
  /* globalFrictionContact_printInFile(problem, file); */
  /* fclose(file); */
  int rinfo =  gfc3d_driver(problem,
			    reaction,
			    velocity,
			    globalVelocity,
			    &numerics_solver_options);


  iparam[SICONOS_IPARAM_ITER_DONE] = numerics_solver_options.iparam[SICONOS_IPARAM_ITER_DONE];
  dparam[SICONOS_DPARAM_TOL] = numerics_solver_options.dparam[SICONOS_DPARAM_TOL];

  /* FILE * file1  =  fopen("tutu.dat", "w"); */
  /* globalFrictionContact_printInFile(problem, file1); */
  /* fclose(file1); */
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
    gfccounter++;
    struct stat st = {};
    if (stat("./fclib-hdf5/", &st) == -1) {
      mkdir("./fclib-hdf5/", 0700);
    }
    printf("################################## gfcccounter = %i\n", gfccounter);
    if (gfccounter % freq_output == 0)
    {
      char fname[256];
      snprintf(fname, sizeof(fname), "./fclib-hdf5/LMGC_GFC3D-i%.5d-%i-%.5d.hdf5", numerics_solver_options.iparam[SICONOS_IPARAM_ITER_DONE], nc, gfccounter);
      printf("Dump ./fclib-hdf5/LMGC_GFC3D-i%.5d-%i-%.5d.hdf5.\n", numerics_solver_options.iparam[SICONOS_IPARAM_ITER_DONE], nc, gfccounter);
      /* printf("ndof = %i.\n", ndof); */

      FILE * foutput  =  fopen(fname, "w");
      int n = 100;
      char * title = (char *)malloc(n * sizeof(char *));
      strncpy(title, "LMGC dump in hdf5", n);
      char * description = (char *)malloc(n * sizeof(char *));

      snprintf(description, n, "Rewriting in hdf5 through siconos of %s in FCLIB format", fname);
      char * mathInfo = (char *)malloc(n * sizeof(char *));
      strncpy(mathInfo, "unknown", n);

      globalFrictionContact_fclib_write(problem,
                                        title,
                                        description,
                                        mathInfo,
                                        fname);


      fclose(foutput);
    }
#else
    printf("Fclib is not available ...\n");
#endif

  }


  NM_free(M);
  NM_free(H);
  free(M);
  free(H);
  free(problem);

  /* free(_colM); */
  /* free(_colH); */

  /* free(_rowM); */
  /* free(_rowH); */

  return rinfo;
}


/* int gfc3d_LmgcDriver_SBM(double *reaction, */
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

/*   SolverOptions numerics_solver_options; */

/*   gfc3d_setDefaultSolverOptions(&numerics_solver_options, solver_id); */

/*   assert(isize <= numerics_solver_options.iSize); */
/*   assert(dsize <= numerics_solver_options.dSize); */

/*   for (int i=0; i<isize; ++i) */
/*     numerics_solver_options.iparam[i] = iparam[i]; */

/*   for (int i=0; i<dsize; ++i) */
/*     numerics_solver_options.dparam[i] = dparam[i]; */

/*   int rinfo = gfc3d_driver(&problem, */
/*                                              reaction, */
/*                                              velocity, */
/*                                              globalVelocity, */
/*                                              &numerics_solver_options); */

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
/*     sprintf(fname, "LMGC_GlobalFrictionContactProblem%.5d.hdf5", gfccounter++); */
/*     printf("Dump of LMGC_GlobalFrictionContactProblem%.5d.hdf5", gfccounter); */

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


/*   NM_free(&M); */
/*   NM_free(&H); */

/*   free(_colM); */
/*   free(_colH); */

/*   free(_rowM); */
/*   free(_rowH); */

/*   return rinfo; */
/* } */



/* int gfc3d_LmgcDriver(double *reaction, */
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
/*   return gfc3d_LmgcDriver_CSC(reaction, */
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
