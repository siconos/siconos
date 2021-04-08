/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include "GlobalFrictionContactProblem.h"
#include "FrictionContactProblem.h"
#include <assert.h>            // for assert
#include <stdlib.h>            // for free, malloc, exit, EXIT_FAILURE
#include <sys/errno.h>         // for errno
#include <string.h>            // for memcpy
#include "SiconosBlas.h"         // for cblas_dscal, cblas_dcopy
#include "NumericsMatrix.h"    // for NumericsMatrix, NM_display, NM_clear
#include "NumericsSparseMatrix.h"    // for NumericsMatrix, NM_display, NM_clear
#include "SparseBlockMatrix.h"                   // for SBM_gemv, SBM_free
#include "sanitizer.h"                           // for cblas_dcopy_msan
#include "numerics_verbose.h"  // for CHECK_IO
#include "io_tools.h"

#include "siconos_debug.h"
#if defined(WITH_FCLIB)
#include "fclib_interface.h"
#endif

//#define OUTPUT_DEBUG 1
#ifdef OUTPUT_DEBUG
#include "NumericsVector.h"
#endif

GlobalFrictionContactProblem* globalFrictionContactProblem_new(void)
{
  GlobalFrictionContactProblem* problem = malloc(sizeof(GlobalFrictionContactProblem));
  problem->M = NULL;
  problem->H = NULL;
  problem->q = NULL;
  problem->b = NULL;
  problem->mu = NULL;
  problem->env = NULL;
  problem->numberOfContacts = 0;
  problem->dimension = 0;
  return problem;
}

int globalFrictionContact_printInFile(GlobalFrictionContactProblem*  problem, FILE* file)
{
  if(! problem)
  {
    fprintf(stderr, "Numerics, GlobalFrictionContactProblem printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i;

  int d  = problem->dimension;
  fprintf(file, "%d\n", d);
  int nc = problem->numberOfContacts;
  fprintf(file, "%d\n", nc);
  NM_write_in_file(problem->M, file);
  NM_write_in_file(problem->H, file);
  for(i = 0; i < problem->M->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->q[i]);
  }
  fprintf(file, "\n");
  for(i = 0; i < problem->H->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->b[i]);
  }
  fprintf(file, "\n");
  for(i = 0; i < nc; i++)
  {
    fprintf(file, "%32.24e ", problem->mu[i]);
  }
  fprintf(file, "\n");
  return 0;
}

GlobalFrictionContactProblem* globalFrictionContact_newFromFile(FILE* file)
{
  GlobalFrictionContactProblem* problem = globalFrictionContactProblem_new();
  int nc = 0, d = 0;
  int info = 0;
  CHECK_IO(fscanf(file, "%d\n", &d), &info);
  problem->dimension = d;
  CHECK_IO(fscanf(file, "%d\n", &nc), &info);
  problem->numberOfContacts = nc;
  problem->M = NM_new_from_file(file);

  problem->H =  NM_new_from_file(file);

  problem->q = (double *) malloc(problem->M->size1 * sizeof(double));
  for(int i = 0; i < problem->M->size1; ++i)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->q[i])), &info);
  }
  problem->b = (double *) malloc(problem->H->size1 * sizeof(double));
  for(int i = 0; i < problem->H->size1; ++i)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->b[i])), &info);
  }

  problem->mu = (double *) malloc(nc * sizeof(double));
  for(int i = 0; i < nc; ++i)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->mu[i])), &info);
  }

  if(info)
  {
    problem = NULL;
  }
  return problem;
}

GlobalFrictionContactProblem* globalFrictionContact_new_from_filename(const char* filename)
{
  GlobalFrictionContactProblem* problem = NULL;
  int is_hdf5 = check_hdf5_file(filename);
  if(is_hdf5)
  {
#if defined(WITH_FCLIB)
    problem = globalFrictionContact_fclib_read(filename);
#else
    numerics_error("GlobalFrictionContactProblem",
                   "Try to read an hdf5 file, while fclib interface is not active. Recompile Siconos with fclib.",
                   filename);
#endif
  }
  else
  {
    FILE * file = fopen(filename, "r");
    if(!file)
      numerics_error("GlobalFrictionContactProblem", "Can not open file ", filename);

    problem = globalFrictionContact_newFromFile(file);
    fclose(file);
  }
  return problem;

}

int globalFrictionContact_printInFileName(GlobalFrictionContactProblem* problem, const char* filename)
{
  int info = 0;
  FILE * file = fopen(filename, "w");

  if(!file)
  {
    return errno;
  }

  info = globalFrictionContact_printInFile(problem, file);

  fclose(file);
  return info;
}

void globalFrictionContact_free(GlobalFrictionContactProblem* problem)
{

  if(problem->M)
  {
    NM_clear(problem->M);
    free(problem->M);
  }
  problem->M = NULL;

  if(problem->H)
  {
    NM_clear(problem->H);
    free(problem->H);
  }
  problem->H = NULL;

  if(problem->mu)
  {
    free(problem->mu);
  }
  problem->mu = NULL;

  if(problem->q)
  {
    free(problem->q);
  }
  problem->q = NULL;

  if(problem->b)
  {
    free(problem->b);
  }
  problem->b = NULL;

  if(problem->env) assert(0 && "globalFrictionContact_free :: problem->env != NULL, don't know what to do");

  free(problem);
}

void globalFrictionContact_display(GlobalFrictionContactProblem* problem)
{

  assert(problem);
  int i, n = problem->dimension * problem->numberOfContacts;
  printf("GlobalFrictionContact Display :\n-------------\n");
  printf("dimension :%d \n", problem->dimension);
  printf("numberOfContacts:%d \n", problem->numberOfContacts);
  int m = problem->M->size0;
  if(problem->M)
  {
    printf("M matrix:\n");
    NM_display(problem->M);
  }
  else
    printf("No M matrix:\n");
  if(problem->H)
  {
    printf("H matrix:\n");
    NM_display(problem->H);
  }
  else
    printf("No H matrix:\n");

  if(problem->q)
  {
    printf("q vector:\n");
    for(i = 0; i < m; i++) printf("q[ %i ] = %12.8e\n", i, problem->q[i]);
  }
  else
    printf("No q vector:\n");

  if(problem->b)
  {
    printf("b vector:\n");
    for(i = 0; i < n; i++) printf("b[ %i ] = %12.8e\n", i, problem->b[i]);
  }
  else
    printf("No q vector:\n");

  if(problem->mu)
  {
    printf("mu vector:\n");
    for(i = 0; i < problem->numberOfContacts; i++) printf("mu[ %i ] = %12.8e\n", i, problem->mu[i]);
  }
  else
    printf("No mu vector:\n");

}

GlobalFrictionContactProblem* globalFrictionContact_copy(GlobalFrictionContactProblem* problem)
{
  assert(problem);

  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  int m = 3 * nc;

  GlobalFrictionContactProblem* new = (GlobalFrictionContactProblem*) malloc(sizeof(GlobalFrictionContactProblem));
  new->dimension = problem->dimension;
  new->numberOfContacts = problem->numberOfContacts;
  new->M = NM_new();
  NM_copy(problem->M, new->M);
  new->H = NM_new();
  NM_copy(problem->H, new->H);
  new->q = (double*)malloc(n*sizeof(double));
  memcpy(new->q, problem->q, n*sizeof(double));
  new->b = (double*)malloc(m*sizeof(double));
  memcpy(new->b, problem->b, m*sizeof(double));
  new->mu = (double*)malloc(nc*sizeof(double));
  memcpy(new->mu, problem->mu, nc*sizeof(double));
  new->env = NULL;
  return new;
}

int globalFrictionContact_computeGlobalVelocity(
  GlobalFrictionContactProblem* problem,
  double * reaction,
  double * globalVelocity)
{
  int info = -1;

  int n = problem->M->size0;
  int m = problem->H->size1;

  /* globalVelocity <- problem->q */
  cblas_dcopy(n,  problem->q, 1, globalVelocity, 1);

  // We add the reaction part only if the problem has contacts
  if(m>0)
  {
    /* globalVelocity <-  H*reaction + globalVelocity*/
    NM_gemv(1.0, problem->H, reaction, 1.0, globalVelocity);
    DEBUG_EXPR(NM_vector_display(reaction, m));
  }

  /* Compute globalVelocity <- M^(-1) globalVelocity*/
  //info = NM_gesv_expert(problem->M, globalVelocity, NM_PRESERVE);
  info = NM_LU_solve(problem->M, globalVelocity, 1);
  DEBUG_EXPR(NM_vector_display(globalVelocity, n));

  return info;
}

FrictionContactProblem * globalFrictionContact_reformulation_FrictionContact(GlobalFrictionContactProblem* problem)
{
  int info = -1;

  NumericsMatrix *M = problem->M;
  NumericsMatrix *H = problem->H;

  int n = M->size0;
  int m = H->size1;

  FrictionContactProblem* localproblem = frictionContactProblem_new();

  localproblem->numberOfContacts = problem->numberOfContacts;
  localproblem->dimension =  problem->dimension;
  localproblem->mu = (double*)calloc(problem->numberOfContacts,sizeof(double));
  cblas_dcopy_msan(problem->numberOfContacts, problem->mu, 1, localproblem->mu, 1);

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

  if(H->storageType != M->storageType)
  {
    //     if(verbose==1)
    printf(" ->storageType != M->storageType :This case is not taken into account\n");
    return NULL;
  }
#ifdef OUTPUT_DEBUG
  FILE * fileout;
#endif
  if(M->storageType == NM_DENSE)
  {



    int nm = n * m;

    double *Htmp = (double*)calloc(nm, sizeof(double));
    // compute W = H^T M^-1 H
    //Copy Htmp <- H
    cblas_dcopy_msan(nm,  H->matrix0, 1, Htmp, 1);

    //Compute Htmp   <- M^-1 Htmp
#ifdef USE_LAPACK_DGETRS
    lapack_int* ipiv = (lapack_int*)NM_iWork(M, M->size0, sizeof(lapack_int));
    lapack_int infoDGETRF;
    lapack_int infoDGETRS;
    DGETRF(n, n, M->matrix0, n, ipiv, &infoDGETRF);
    assert(!infoDGETRF);
    NM_set_LU_factorized(M, true);
    DGETRS(LA_NOTRANS, n, m,  M->matrix0, n, ipiv, Htmp, n, &infoDGETRS);
#else
    // NM_gesv_expert_multiple_rhs(M,Htmp,m,NM_KEEP_FACTORS);
    NM_LU_solve(M, Htmp, m);
#endif

    /* assert(!infoDGETRS); */
    /*      DGESV(n, m, M->matrix0, n, ipiv, Htmp, n, infoDGESV); */

    localproblem->M = NM_new();
    NumericsMatrix *Wnum = localproblem->M;
    Wnum->storageType = 0;
    Wnum-> size0 = m;
    Wnum-> size1 = m;
    Wnum->matrix0 = (double*)calloc(m * m, sizeof(double));
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
    localproblem->q = (double*)calloc(m, sizeof(double));
    cblas_dcopy_msan(m, problem->b, 1, localproblem->q, 1);

    double* qtmp = (double*)calloc(n , sizeof(double));
    cblas_dcopy_msan(n,  problem->q, 1, qtmp, 1);

    // compute H^T M^(-1) q + b
#ifdef USE_LAPACK_DGETRS
    DGETRS(LA_NOTRANS, n, 1,  M->matrix0, n, ipiv, qtmp, n, &infoDGETRS);
#else
    // NM_gesv_expert(M,qtmp,NM_KEEP_FACTORS);
    NM_LU_solve(M, qtmp, 1);
#endif

    cblas_dgemv(CblasColMajor,CblasTrans, n, m, 1.0, H->matrix0, n, qtmp, 1, 1.0, localproblem->q, 1);

    frictionContact_display(localproblem);

    free(Htmp);
    free(qtmp);


  }

  else if(M->storageType == NM_SPARSE_BLOCK)
  {
    int n = M->size0;
    int m = H->size1;

    // compute W = H^T M^-1 H
    // compute MinvH   <- M^-1 H
    /* int infoMInv = 0; */
    /* infoMInv = NM_inverse_diagonal_block_matrix_in_place(M); */
    assert(!NM_inverse_diagonal_block_matrix_in_place(M));

    DEBUG_PRINT("M inverse :");
    DEBUG_EXPR(NM_display(M));


    NumericsMatrix * MinvH = NM_multiply(M,H);

    /* NumericsMatrix * MinvH= NM_create(NM_SPARSE_BLOCK, m, m); */
    /* double alpha = 1.0, beta = 0.0; */
    /* NM_gemm(alpha, M, H, beta, MinvH); */

    NumericsMatrix * Htrans= NM_create(NM_SPARSE_BLOCK, H->size1, H->size0);
    SBM_transpose(H->matrix1, Htrans->matrix1);

    /* localproblem->M = NM_create(NM_SPARSE_BLOCK, m, m ); */
    /* NumericsMatrix *W = localproblem->M; */
    /* NM_gemm(alpha, Htrans, MinvH, beta, W); */

    localproblem->M =  NM_multiply(Htrans,MinvH);


#ifdef OUTPUT_DEBUG
    FILE * fileout;
    fileout = fopen("dataW.sci", "w");
    NM_write_in_file_scilab(localproblem->M, fileout);
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
    WnumInverse->matrix0 = (double*)calloc(m * m , sizeof(double));
    double * WInverse = WnumInverse->matrix0;
    SBM_to_dense(W, WnumInverse->matrix0);

    FILE * file1 = fopen("dataW.dat", "w");
    NM_write_in_file_scilab(WnumInverse, file1);
    fclose(file1);

    double * WInversetmp = (double*)calloc(m * m,  sizeof(double));
    memcpy(WInversetmp, WInverse, m * m * sizeof(double));
    double  condW;
    condW = cond(WInverse, m, m);

    lapack_int* ipiv = (lapack_int *)calloc(m , sizeof(lapack_int));
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

    localproblem->q = (double*)calloc(m, sizeof(double));

    //Copy localq<- b
    cblas_dcopy_msan(m, problem->b, 1, localproblem->q, 1);

    // compute H^T M^-1 q+ b
    double* qtmp = (double*)calloc(n,  sizeof(double));
    double alpha = 1.0, beta = 1.0;
    double beta2 = 0.0;
    NM_gemv(alpha, M, problem->q, beta2, qtmp); /* Warning M contains Minv */
    NM_gemv(alpha, Htrans, qtmp, beta, localproblem->q);


    NM_clear(MinvH);
    NM_clear(Htrans);
    free(MinvH);
    free(Htrans);
    free(qtmp);
  }
  else if(M->storageType == NM_SPARSE)
  {

#ifdef OUTPUT_DEBUG
    fileout = fopen("dataM.py", "w");
    NM_write_in_file_python(M, fileout);
    fclose(fileout);
    fileout = fopen("dataH.py", "w");
    NM_write_in_file_python(H, fileout);
    fclose(fileout);
    fileout = fopen("dataq.py", "w");
    NV_write_in_file_python(problem->q, M->size0, fileout);
    fclose(fileout);
    fileout = fopen("datab.py", "w");
    NV_write_in_file_python(problem->b, H->size1, fileout);
    fclose(fileout);
#endif




    // Product M^-1 H
    DEBUG_EXPR(NM_display(H););
    numerics_printf_verbose(1,"inversion of the matrix M ...");
    NumericsMatrix * Minv  = NM_LU_inv(M);
    DEBUG_EXPR(NM_display(Minv););


    /* NumericsMatrix* MinvH = NM_create(NM_SPARSE,n,m); */
    /* NM_triplet_alloc(MinvH, n); */
    /* MinvH->matrix2->origin = NSM_TRIPLET; */
    /* DEBUG_EXPR(NM_display(MinvH);); */
    /* NM_gemm(1.0, Minv, H, 0.0, MinvH); */
    numerics_printf_verbose(1,"multiplication  H^T M^{-1} H ...");
    NumericsMatrix* MinvH = NM_multiply(Minv,H);
    DEBUG_EXPR(NM_display(MinvH););

    // Product H^T M^-1 H
    NM_csc_trans(H);


    NumericsMatrix* Htrans = NM_new();
    Htrans->storageType = NM_SPARSE;
    Htrans-> size0 = m;
    Htrans-> size1 = n;
    NM_csc_alloc(Htrans, 0);
    Htrans->matrix2->origin = NSM_CSC;
    Htrans->matrix2->csc = NM_csc_trans(H);
    DEBUG_EXPR(NM_display(Htrans););

    /* localproblem->M = NM_create(NM_SPARSE, m, m); */
    /* NumericsMatrix *W = localproblem->M; */
    /* int nzmax= m*m; */
    /* NM_csc_empty_alloc(W, nzmax); */
    /* W->matrix2->origin = NSM_CSC; */
    /* NM_gemm(1.0, Htrans, MinvH, 0.0, W); */

    localproblem->M = NM_multiply(Htrans,MinvH);
    DEBUG_EXPR(NM_display(localproblem->M););

    NSM_fix_csc(NM_csc(localproblem->M));

#ifdef OUTPUT_DEBUG
    fileout = fopen("dataW.py", "w");
    NM_write_in_file_python(localproblem->M, fileout);
    fclose(fileout);
#endif

    /* compute localq <- H^T M^(-1) q + b */
    numerics_printf_verbose(1,"Compute localq = H^T M^(-1) q +b  ...");
    double* qtmp = (double*)calloc(n, sizeof(double));
    NM_gemv(1.0, Minv, problem->q, 0.0, qtmp);   /* qtmp <- M^(-1) q  */
    DEBUG_EXPR(NV_display(qtmp,n););
    DEBUG_EXPR(NV_display(problem->q,n););
    localproblem->q = (double*)calloc(m, sizeof(double));
    cblas_dcopy_msan(m, problem->b, 1, localproblem->q, 1);     /*Copy localq <- b */
    NM_gemv(1.0, Htrans, qtmp, 1.0, localproblem->q); /* localq <- H^T qtmp + localq   */

#ifdef OUTPUT_DEBUG
    fileout = fopen("dataqloc.py", "w");
    NV_write_in_file_python(localproblem->q, H->size1, fileout);
    fclose(fileout);
#endif


    DEBUG_EXPR(frictionContact_display(localproblem););
    //getchar();
  }
  else
  {
    printf("globalFrictionContact_reformulation_FrictionContact :: unknown matrix storage");
    exit(EXIT_FAILURE);
  }
  return localproblem;
}
