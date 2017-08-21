/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <assert.h>

#include "SiconosLapack.h"
#include "gfc3d_Solvers.h"
#include "NonSmoothDrivers.h"
#include "fc3d_Solvers.h"
#include "cond.h"
#include "pinv.h"
#include <string.h>

#include "NumericsSparseMatrix.h"
#include "NumericsVector.h"

#include "sanitizer.h"
#include "numerics_verbose.h"
//#define TEST_COND
/* #define OUTPUT_DEBUG */


/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "debug.h"

//#define USE_LAPACK_DGETRS

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

/* Global Variable for the reformulation of the problem */

int gfc3d_reformulation_local_problem(GlobalFrictionContactProblem* problem, FrictionContactProblem* localproblem)
{
  int info = -1;

  NumericsMatrix *M = problem->M;
  NumericsMatrix *H = problem->H;

  int n = M->size0;
  int m = H->size1;

  localproblem->numberOfContacts = problem->numberOfContacts;
  localproblem->dimension =  problem->dimension;
  localproblem->mu =  problem->mu;

  /* assert(M); */
  /* assert(H); */
  /* NM_display(M); */
  /* NM_display(H); */
  /* NumericsMatrix *MMtmp = NM_new(); */
  /* NumericsMatrix *HHtmp = NM_new(); */

  /* NM_copy(M,MMtmp); */
  /* NM_copy(H,HHtmp); */

  /* NM_clearSparse(M); */
  /* NM_clearSparse(H); */

  /* M = NM_create(NM_DENSE, n, n); */
  /* H = NM_create(NM_DENSE, n, m); */

  /* NM_to_dense(MMtmp,M); */
  /* NM_to_dense(HHtmp,H); */

  /* NM_display(M); */
  /* NM_display(H); */

  if (H->storageType != M->storageType)
  {
    //     if(verbose==1)
    printf(" ->storageType != M->storageType :This case is not taken into account\n");
    return info;
  }
#ifdef OUTPUT_DEBUG
  FILE * fileout;
#endif
  if (M->storageType == NM_DENSE)
  {



    int nm = n * m;

    double *Htmp = (double*)malloc(nm * sizeof(double));
    // compute W = H^T M^-1 H
    //Copy Htmp <- H
    cblas_dcopy_msan(nm,  H->matrix0 , 1, Htmp, 1);

    //Compute Htmp   <- M^-1 Htmp
#ifdef USE_LAPACK_DGETRS
    lapack_int* ipiv = (lapack_int*)NM_iWork(M, M->size0, sizeof(lapack_int));
    lapack_int infoDGETRF;
    lapack_int infoDGETRS;
    DGETRF(n, n, M->matrix0, n, ipiv, &infoDGETRF);
    assert(!infoDGETRF);
    NM_internalData(M)->isLUfactorized = true;
    DGETRS(LA_NOTRANS, n, m,  M->matrix0, n, ipiv, Htmp, n, &infoDGETRS);
#else
    NM_gesv_expert_multiple_rhs(M,Htmp,m,NM_KEEP_FACTORS);
#endif

    /* assert(!infoDGETRS); */
    /*      DGESV(n, m, M->matrix0, n, ipiv, Htmp, n, infoDGESV); */

    localproblem->M = NM_new();
    NumericsMatrix *Wnum = localproblem->M;
    Wnum->storageType = 0;
    Wnum-> size0 = m;
    Wnum-> size1 = m;
    Wnum->matrix0 = (double*)calloc(m * m , sizeof(double));
    Wnum->matrix1 = NULL;
    Wnum->matrix2 = NULL;
    Wnum->internalData = NULL;
    // Compute W <-  H^T M^1 H

    assert(H->matrix0);
    assert(Htmp);
    assert(Wnum->matrix0);

    cblas_dgemm(CblasColMajor,CblasTrans, CblasNoTrans, m, m, n, 1.0, H->matrix0, n, Htmp, n, 0.0, Wnum->matrix0, m);
    /*     DGEMM(CblasTrans,CblasNoTrans,m,m,n,1.0,H->matrix0,n,Htmp,n,0.0,Wnum->matrix0,m); */

    // compute localq = H^T M^(-1) q +b

    //Copy localq <- b
    localproblem->q = (double*)malloc(m * sizeof(double));
    cblas_dcopy_msan(m, problem->b , 1, localproblem->q, 1);

    double* qtmp = (double*)malloc(n * sizeof(double));
    cblas_dcopy_msan(n,  problem->q, 1, qtmp, 1);

    // compute H^T M^(-1) q + b
#ifdef USE_LAPACK_DGETRS
    DGETRS(LA_NOTRANS, n, 1,  M->matrix0, n, ipiv, qtmp, n, &infoDGETRS);
#else
    NM_gesv_expert(M,qtmp,NM_KEEP_FACTORS);
#endif

    cblas_dgemv(CblasColMajor,CblasTrans, n, m, 1.0, H->matrix0 , n, qtmp, 1, 1.0, localproblem->q, 1);
    // Copy mu
    localproblem->mu = problem->mu;

    frictionContact_display(localproblem);

    free(Htmp);
    free(qtmp);


  }

  else if (M->storageType == NM_SPARSE_BLOCK)
  {
    int n = M->size0;
    int m = H->size1;

    // compute W = H^T M^-1 H


    NumericsMatrix * MinvH= NM_create(NM_SPARSE_BLOCK, m, m);

    //Compute MinvH   <- M^-1 H


    int infoMInv = NM_inverse_diagonal_block_matrix_in_place(M);
    assert(!infoMInv);

    DEBUG_PRINT("M inverse :");
    DEBUG_EXPR(NM_display(M));

    double alpha = 1.0, beta = 1.0;
    NM_gemm(alpha, M, H, beta, MinvH);

    NumericsMatrix * Htrans= NM_create(NM_SPARSE_BLOCK, H->size1, H->size0);

    SBM_transpose(H->matrix1, Htrans->matrix1);

    localproblem->M = NM_create(NM_SPARSE_BLOCK, m, m );

    NumericsMatrix *W = localproblem->M;
    NM_gemm(alpha, Htrans, MinvH, beta, W);



#ifdef OUTPUT_DEBUG
    FILE * fileout;
    fileout = fopen("dataW.sci", "w");
    NM_write_in_file_scilab(W, fileout);
    fclose(fileout);
#endif

#ifdef TEST_COND
    NumericsMatrix *WnumInverse = NM_new();
    WnumInverse->storageType = 0;
    WnumInverse-> size0 = m;
    WnumInverse-> size1 = m;
    WnumInverse->matrix1 = NULL;
    WnumInverse->matrix2 = NULL;
    WnumInverse->internalData = NULL;
    WnumInverse->matrix0 = (double*)malloc(m * m * sizeof(double));
    double * WInverse = WnumInverse->matrix0;
    SBM_to_dense(W, WnumInverse->matrix0);

    FILE * file1 = fopen("dataW.dat", "w");
    NM_write_in_file_scilab(WnumInverse, file1);
    fclose(file1);

    double * WInversetmp = (double*)malloc(m * m * sizeof(double));
    memcpy(WInversetmp, WInverse, m * m * sizeof(double));
    double  condW;
    condW = cond(WInverse, m, m);

    lapack_int* ipiv = (lapack_int *)malloc(m * sizeof(lapack_int));
    int infoDGETRF = 0;
    DGETRF(m, m, WInverse, m, ipiv, &infoDGETRF);
    assert(!infoDGETRF);
    int infoDGETRI = 0;
    DGETRI(m, WInverse, m, ipiv, &infoDGETRI);


    free(ipiv);
    assert(!infoDGETRI);


    double  condWInverse;
    condWInverse = cond(WInverse, m, m);




    FILE * file2 = fopen("dataWInverse.dat", "w");
    NM_write_in_file_scilab(WnumInverse, file2);
    fclose(file2);

    double tol = 1e-24;
    pinv(WInversetmp, m, m, tol);
    NumericsMatrix *WnumInversetmp = NM_new();
    WnumInversetmp->storageType = 0;
    WnumInversetmp-> size0 = m;
    WnumInversetmp-> size1 = m;
    WnumInversetmp->matrix1 = NULL;
    WnumInversetmp->matrix2 = NULL;
    WnumInversetmp->internalData = NULL;
    WnumInversetmp->matrix0 = WInversetmp ;

    FILE * file3 = fopen("dataWPseudoInverse.dat", "w");
    NM_write_in_file_scilab(WnumInversetmp, file3);
    fclose(file3);


    free(WInverse);
    free(WInversetmp);
    free(WnumInverse);
#endif

    localproblem->q = (double*)malloc(m * sizeof(double));

    //Copy localq<- b
    cblas_dcopy_msan(m, problem->b  , 1, localproblem->q, 1);

    // compute H^T M^-1 q+ b
    double* qtmp = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) qtmp[i] = 0.0;

    double beta2 = 0.0;
    NM_gemv(alpha, M, problem->q, beta2, qtmp);
    NM_gemv(alpha, Htrans, qtmp, beta, localproblem->q);

    localproblem->mu = problem->mu;

    NM_free(MinvH);
    NM_free(Htrans);
    free(MinvH);
    free(Htrans);
    free(qtmp);
  }
  else if (M->storageType == NM_SPARSE)
  {

#ifdef OUTPUT_DEBUG
    fileout = fopen("dataM.py", "w");
    NM_write_in_file_python(M, fileout);
    fclose(fileout);
    fileout = fopen("dataH.py", "w");
    NM_write_in_file_python(H, fileout);
    fclose(fileout);
    fileout = fopen("dataq.py", "w");
    NV_write_in_file_python(problem->q, M->size0  , fileout);
    fclose(fileout);

#endif




    // Product M^-1 H
    DEBUG_EXPR(NM_display(H););

    NumericsMatrix * Minv  = NM_new();
    Minv->size0 = n;
    Minv->size1 = n;
    Minv->storageType = NM_SPARSE;
    NM_inv(M, Minv);
    DEBUG_EXPR(NM_display(Minv););


    NumericsMatrix* MinvH = NM_create(NM_SPARSE,n,m);
    NM_triplet_alloc(MinvH, n);
    MinvH->matrix2->origin = NS_TRIPLET;
    DEBUG_EXPR(NM_display(MinvH););
    NM_gemm(1.0, Minv, H, 0.0, MinvH);
    DEBUG_EXPR(NM_display(MinvH););

    // Product H^T M^-1 H
    NM_csc_trans(H);
    NumericsMatrix* Htrans = NM_new();
    Htrans->storageType = NM_SPARSE;
    Htrans-> size0 = m;
    Htrans-> size1 = n;
    NM_csc_alloc(Htrans, 0);
    Htrans->matrix2->origin = NS_CSC;
    Htrans->matrix2->csc = NM_csc_trans(H);
    DEBUG_EXPR(NM_display(Htrans););

    localproblem->M = NM_create(NM_SPARSE, m, m);

    NumericsMatrix *W = localproblem->M;
    int nzmax= m*m;
    NM_csc_empty_alloc(W, nzmax);
    W->matrix2->origin = NS_CSC;

    NM_gemm(1.0, Htrans, MinvH, 0.0, W);
    DEBUG_EXPR(NM_display(W););

    // compute localq = H^T M^(-1) q +b

    //Copy localq <- b
    //DEBUG_PRINT("Compute locaproblem q\n");
    localproblem->q = (double*)malloc(m * sizeof(double));
    cblas_dcopy_msan(m, problem->b , 1, localproblem->q, 1);

    double* qtmp = (double*)malloc(n * sizeof(double));
    //cblas_dcopy_msan(n,  problem->q, 1, qtmp, 1);

    // compute H^T M^(-1) q + b
    NM_gemv(1.0, Minv, problem->q, 0.0, qtmp);
    NM_gemv(1.0, Htrans, qtmp, 1.0, localproblem->q);

    // Copy mu
    localproblem->mu = problem->mu;
    DEBUG_EXPR(frictionContact_display(localproblem););
  }
  else
  {
    printf("gfc3d_reformulation_local_problem :: unknown matrix storage");
    exit(EXIT_FAILURE);
  }
  return info;
}
int computeGlobalVelocity(GlobalFrictionContactProblem* problem, double * reaction, double * globalVelocity)
{
  int info = -1;

  if (problem->M->storageType == NM_DENSE)
  {
    int n = problem->M->size0;
    int m = problem->H->size1;


    /* Compute globalVelocity   <- H reaction + q*/

    /* globalVelocity <- problem->q */
    cblas_dcopy(n,  problem->q , 1, globalVelocity, 1);

    // We compute only if the local problem has contacts
    if (m>0)
    {
      /* globalVelocity <-  H*reaction + globalVelocity*/
      cblas_dgemv(CblasColMajor,CblasNoTrans, n, m, 1.0, problem->H->matrix0 , n, reaction , 1, 1.0, globalVelocity, 1);
    }

    /* Compute globalVelocity <- M^(-1) globalVelocity*/

    assert(NM_internalData(problem->M)->isLUfactorized);
#ifdef USE_LAPACK_DGETRS
    lapack_int infoDGETRS = 0;
    lapack_int* ipiv = (lapack_int*)NM_iWork(problem->M, problem->M->size0, sizeof(lapack_int));
    DGETRS(LA_NOTRANS, n, 1,   problem->M->matrix0, n, ipiv, globalVelocity , n, &infoDGETRS);
    assert(!infoDGETRS);
#else
    NM_gesv_expert(problem->M,globalVelocity,NM_KEEP_FACTORS);
#endif

  }
  else if (problem->M->storageType == NM_SPARSE_BLOCK)
  {
    int n = problem->M->size0;
    int m = problem->H->size1;

    /* Compute qtmp   <- H reaction + q*/

    double* qtmp = (double*)malloc(n * sizeof(double));
    double alpha = 1.0;
    double beta = 1.0;

    cblas_dcopy_msan(n,  problem->q , 1, qtmp, 1);
    SBM_gemv(m, n, alpha, problem->H->matrix1, reaction, beta, qtmp);
    /* Compute global velocity = M^(-1) qtmp*/


    /*      SBM_inverse_diagonal_block_matrix(M->matrix1); We assume that M->matrix1 is already inverse*/
    assert(NM_internalData(problem->M)->isInversed);

    double beta2 = 0.0;
    SBM_gemv(n, n, alpha,  problem->M->matrix1, qtmp, beta2, globalVelocity);

    free(qtmp);

  }
  else if (problem->M->storageType == NM_SPARSE)
  {
    int n = problem->M->size0;
    int m = problem->H->size1;


    /* Compute globalVelocity   <- H reaction + q*/

    /* globalVelocity <- problem->q */
    cblas_dcopy(n,  problem->q , 1, globalVelocity, 1);
    // We compute only if the local problem has contacts
    if (m>0)
    {
      /* globalVelocity <-  H*reaction + globalVelocity*/
      NM_gemv(1.0, problem->H, reaction, 1.0, globalVelocity);
      DEBUG_EXPR(NM_vector_display(reaction, m));
    }
    /* Compute globalVelocity <- M^(-1) globalVelocity*/
    info = NM_gesv_expert(problem->M, globalVelocity, NM_PRESERVE);
    DEBUG_EXPR(NM_vector_display(globalVelocity, n));
  }
  else
  {
    printf("gfc3d_reformulation_local_problem :: unknown matrix storage");
    exit(EXIT_FAILURE);
  }

  return info;
}

int freeLocalProblem(FrictionContactProblem* localproblem)
{
  int info = -1;

  /*    if (!localproblem->M->storageType) */
  /*  { */
  if (localproblem->M->matrix0)
    free(localproblem->M->matrix0);
  /*  } */
  /*     else */
  /*  { */
  if (localproblem->M->matrix1)
  {
    SBM_free(localproblem->M->matrix1);
    free(localproblem->M->matrix1);
  }
  /*  } */
  free(localproblem->M);
  free(localproblem->q);
  free(localproblem);
  return info;
}



void  gfc3d_nsgs_wr(GlobalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  DEBUG_BEGIN("gfc3d_nsgs_wr\n");
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if (H->size1 > 0)
  {
    // Reformulation
    FrictionContactProblem* localproblem = (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
    if (verbose)
    {
      printf("Reformulation info a reduced problem onto local variables ...\n");
    }
    gfc3d_reformulation_local_problem(problem, localproblem);
    DEBUG_EXPR(frictionContact_display(localproblem););
    if (verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_nsgs(localproblem, reaction , velocity , info , options->internalSolvers);

    options->iparam[1] =  options->internalSolvers->iparam[1];
    options->dparam[1] =  options->internalSolvers->dparam[1];
    computeGlobalVelocity(problem, reaction, globalVelocity);

    freeLocalProblem(localproblem);
  }
  else
  {
    computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
  DEBUG_END("gfc3d_nsgs_wr\n");
}
int gfc3d_nsgs_wr_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the NSGS_WR Solver\n");
  }
  fc3d_nsgs_setDefaultSolverOptions(options);
  options->solverId = SICONOS_GLOBAL_FRICTION_3D_NSGS_WR;
  options->numberOfInternalSolvers = 1;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  fc3d_nsgs_setDefaultSolverOptions(options->internalSolvers);
  return 0;
}



void  gfc3d_nonsmooth_Newton_AlartCurnier_wr(GlobalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  DEBUG_BEGIN("gfc3d_nonsmooth_Newton_AlartCurnier_wr(...)\n");
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if (H->size1 > 0)
  {
    // Reformulation
    FrictionContactProblem* localproblem = (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
    if (verbose)
    {
      printf("Reformulation info a reduced problem onto local variables ...\n");
    }
    gfc3d_reformulation_local_problem(problem, localproblem);
    DEBUG_EXPR(frictionContact_display(localproblem););
    if (verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_nonsmooth_Newton_AlartCurnier(localproblem, reaction , velocity , info , options->internalSolvers);

    options->iparam[1] =  options->internalSolvers->iparam[1];
    options->dparam[1] =  options->internalSolvers->dparam[1];
    computeGlobalVelocity(problem, reaction, globalVelocity);

    freeLocalProblem(localproblem);
  }
  else
  {
    computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }

  DEBUG_END("gfc3d_nonsmooth_Newton_AlartCurnier_wr(...)\n")


}
int gfc3d_nonsmooth_Newton_AlartCurnier_wr_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the NSN_AC_WR Solver\n");
  }
  gfc3d_nonsmooth_Newton_AlartCurnier_setDefaultSolverOptions(options);
  options->solverId = SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR;
  options->numberOfInternalSolvers = 1;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  gfc3d_nonsmooth_Newton_AlartCurnier_setDefaultSolverOptions(options->internalSolvers);
  return 0;
}

void  gfc3d_nsgs_velocity_wr(GlobalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if (H->size1 > 0)
  {
    // Reformulation
    FrictionContactProblem* localproblem = (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
    if (verbose)
    {
      printf("Reformulation info a reduced problem onto local variables ...\n");
    }
    gfc3d_reformulation_local_problem(problem, localproblem);
    DEBUG_EXPR(frictionContact_display(localproblem););
    if (verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_nsgs_velocity(localproblem, reaction , velocity , info , options->internalSolvers);

    options->iparam[1] =  options->internalSolvers->iparam[1];
    options->dparam[1] =  options->internalSolvers->dparam[1];
    computeGlobalVelocity(problem, reaction, globalVelocity);

    freeLocalProblem(localproblem);
  }
  else
  {
    computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
}
int gfc3d_nsgs_velocity_wr_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the NSGSV_WR Solver\n");
  }
  fc3d_nsgs_velocity_setDefaultSolverOptions(options);
  options->solverId = SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR;
  options->numberOfInternalSolvers = 1;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  fc3d_nsgs_velocity_setDefaultSolverOptions(options->internalSolvers);
  return 0;

}

void  gfc3d_proximal_wr(GlobalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if (H->size1 > 0)
  {
    // Reformulation
    FrictionContactProblem* localproblem = (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
    if (verbose)
    {
      printf("Reformulation info a reduced problem onto local variables ...\n");
    }
    gfc3d_reformulation_local_problem(problem, localproblem);
    DEBUG_EXPR(frictionContact_display(localproblem););
    if (verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_proximal(localproblem, reaction , velocity , info , options->internalSolvers);

    options->iparam[1] =  options->internalSolvers->iparam[1];
    options->dparam[1] =  options->internalSolvers->dparam[1];
    computeGlobalVelocity(problem, reaction, globalVelocity);

    freeLocalProblem(localproblem);
  }
  else
  {
    computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
}
int gfc3d_proximal_wr_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the PROX_WR Solver\n");
  }
  fc3d_proximal_setDefaultSolverOptions(options);
  options->solverId = SICONOS_GLOBAL_FRICTION_3D_PROX_WR;
  options->numberOfInternalSolvers = 1;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  fc3d_proximal_setDefaultSolverOptions(options->internalSolvers);
  return 0;
}
void  gfc3d_DeSaxceFixedPoint_wr(GlobalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if (H->size1 > 0)
  {
    // Reformulation
    FrictionContactProblem* localproblem = (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
    if (verbose)
    {
      printf("Reformulation info a reduced problem onto local variables ...\n");
    }
    gfc3d_reformulation_local_problem(problem, localproblem);
    DEBUG_EXPR(frictionContact_display(localproblem););
    if (verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_DeSaxceFixedPoint(localproblem, reaction , velocity , info , options->internalSolvers);
    options->iparam[1] =  options->internalSolvers->iparam[1];
    options->dparam[1] =  options->internalSolvers->dparam[1];
    computeGlobalVelocity(problem, reaction, globalVelocity);

    freeLocalProblem(localproblem);
  }
  else
  {
    computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
}
int gfc3d_DeSaxceFixedPoint_setDefaultSolverOptions(SolverOptions* options)
{


  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the DSFP_WR Solver\n");
  }
  fc3d_DeSaxceFixedPoint_setDefaultSolverOptions(options);
  options->solverId = SICONOS_GLOBAL_FRICTION_3D_DSFP_WR;
  options->numberOfInternalSolvers = 1;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  fc3d_DeSaxceFixedPoint_setDefaultSolverOptions(options->internalSolvers);
  return 0;
}

void  gfc3d_TrescaFixedPoint_wr(GlobalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if (H->size1 > 0)
  {
    // Reformulation
    FrictionContactProblem* localproblem = (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
    if (verbose)
    {
      printf("Reformulation info a reduced problem onto local variables ...\n");
    }
    gfc3d_reformulation_local_problem(problem, localproblem);
    DEBUG_EXPR(frictionContact_display(localproblem););
    if (verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_TrescaFixedPoint(localproblem, reaction , velocity , info , options->internalSolvers);
    options->iparam[1] =  options->internalSolvers->iparam[1];
    options->dparam[1] =  options->internalSolvers->dparam[1];
    computeGlobalVelocity(problem, reaction, globalVelocity);

    freeLocalProblem(localproblem);
  }
  else
  {
    computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }



}
int gfc3d_TrescaFixedPoint_setDefaultSolverOptions(SolverOptions* options)
{


  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the DSFP_WR Solver\n");
  }
  fc3d_TrescaFixedPoint_setDefaultSolverOptions(options);
  options->solverId = SICONOS_GLOBAL_FRICTION_3D_TFP_WR;
  options->numberOfInternalSolvers = 1;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  fc3d_TrescaFixedPoint_setDefaultSolverOptions(options->internalSolvers);
  return 0;
}
