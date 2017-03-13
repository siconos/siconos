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

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "avi_caoferris.h"

#include "AVI_Solvers.h"
#include "AVI_cst.h"
#include "pivot-utils.h"
#include "LinearComplementarityProblem.h"
#include "vertex_extraction.h"
#include "numerics_verbose.h"

#include "SiconosLapack.h"

#include "sanitizer.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"


int avi_caoferris(AffineVariationalInequalities* problem, double *z, double *w, SolverOptions* options)
{
  assert(problem);
  assert(problem->M);
  assert(problem->q);
  assert(problem->poly.set);
  if (problem->poly.set->id != SICONOS_SET_POLYHEDRON)
  {
    numerics_error_nonfatal("avi_caoferris", "unsupported set type %d", problem->poly.set->id);
    return -1;
  }
  unsigned n = problem->size;
  assert(n > 0);
  unsigned nrows = problem->poly.split->size_ineq;
  assert(nrows - n > 0);
  unsigned n_I = nrows - n; /* Number of inactive constraints */

  /* Create the data  problem */
  LinearComplementarityProblem lcplike_pb;
  lcplike_pb.size = nrows;
  NumericsMatrix num_mat;
  fillNumericsMatrix(&num_mat, NM_DENSE, nrows, nrows, calloc(nrows*nrows, sizeof(double)));

  lcplike_pb.M = &num_mat;

  lcplike_pb.q = (double *)calloc(nrows, sizeof(double));
  double* a_bar = (double *)malloc(nrows*sizeof(double));

  double* B_A_T = (double*)malloc(n*n*sizeof(double));
  double* copyA = (double*)malloc(n*n*sizeof(double));
  double* B_I_T = (double*)malloc(n*(n_I)*sizeof(double));
  double* d_vec = (double *)malloc(nrows*sizeof(double));
  lapack_int* basis = (lapack_int *)malloc((2*nrows+1)*sizeof(lapack_int));

  siconos_find_vertex(problem->poly.split, n, basis);
  DEBUG_PRINT_VEC_INT(basis, nrows+1);
  const double* H = problem->poly.split->H->matrix0;
  assert(H);
  const double* K = problem->poly.split->K;
  /* Set of active constraints */
  unsigned* A = (unsigned*)malloc(n*sizeof(unsigned));
  lapack_int* active_constraints = &basis[nrows+1];

  /* set active_constraints to 1 at the beginning */
  memset(active_constraints, -1, nrows*sizeof(lapack_int));
  DEBUG_PRINT_VEC_INT(active_constraints, nrows);
  unsigned indx_B_I_T = 0;
  for (unsigned i = 1; i <= nrows; ++i)
  {
    assert((unsigned)abs(basis[i]) > nrows); /* we don't want slack variable here */
    int indx = abs(basis[i]) - nrows - 1 - n;
    if (indx >= 0)
    {
      /* this is an inactive constraint */
      assert(indx_B_I_T < n_I);
      assert((unsigned)indx < nrows); 
      cblas_dcopy(n, &H[indx], nrows, &B_I_T[indx_B_I_T*n], 1); /* form B_I_T */
      active_constraints[indx] = 0; /* desactivate the constraint */
      lcplike_pb.q[n+indx_B_I_T] = -K[indx]; /* partial construction of q[n:nrows] as -K_I  */
      indx_B_I_T++;
    }
  }
  DEBUG_PRINT_VEC_INT(active_constraints, nrows);

  unsigned indx_B_A_T = 0;
  for (unsigned i = 0; i < nrows; ++i)
  {
    if (active_constraints[i] == -1)
    {
      assert(indx_B_A_T < n);
      A[indx_B_A_T] = i+1; /* note which constraints is active */
      cblas_dcopy(n, &H[i], nrows, &B_A_T[indx_B_A_T*n], 1); /* form B_A_T */
      d_vec[indx_B_A_T] = K[i]; /* save K_A */
      indx_B_A_T++;
    }
  }
  assert(indx_B_A_T == n && "there were not enough active constraints");
  DEBUG_PRINT_VEC_STR("K_A", d_vec, n);
  cblas_dcopy(n*n, problem->M->matrix0, 1, copyA, 1);

  DEBUG_PRINT_MAT(B_A_T, n, n);
  DEBUG_PRINT_MAT(B_I_T, n, n_I);

  /* get LU for B_A_T */
  lapack_int* ipiv = basis;
  lapack_int infoLAPACK = 0;

  /* LU factorisation of B_A_T  */
  DGETRF(n, n, B_A_T, n, ipiv, &infoLAPACK);
  assert(infoLAPACK <= 0 && "avi_caoferris :: info from DGETRF > 0, this should not append !\n");

  /* compute B_A_T^{-1}B_I_T  */
  DGETRS(LA_NOTRANS, n, n_I, B_A_T, n, ipiv, B_I_T, n, &infoLAPACK);
  assert(infoLAPACK == 0 && "avi_caoferris :: info from DGETRS for solving B_A_T X = B_I_T is not zero!\n");

  DEBUG_PRINT("B_A_T^{-1}B_I_T\n");
  DEBUG_PRINT_MAT(B_I_T, n, n_I);

  /* Compute B_A_T^{-1} A */
  DGETRS(LA_NOTRANS, n, n, B_A_T, n, ipiv, copyA, n, &infoLAPACK);
  assert(infoLAPACK == 0 && "avi_caoferris :: info from DGETRS for solving B_A_T X = A is not zero!\n");

  DEBUG_PRINT("B_A_T^{-1}A\n");
  DEBUG_PRINT_MAT(copyA, n, n);

  /* do some precomputation for \bar{q}: B_A_T^{-1}q_{AVI} */
  cblas_dcopy_msan(n, problem->q, 1, a_bar, 1);
  DGETRS(LA_NOTRANS, n, 1, B_A_T, n, ipiv, a_bar, n, &infoLAPACK);
  assert(infoLAPACK == 0  && "avi_caoferris :: info from DGETRS for solving B_A_T X = a_bar is not zero!\n");
  DEBUG_PRINT_VEC_STR("B_A_T{-1}q_{AVI}", a_bar, n);

  /* Do the transpose of B_A_T^{-1} A */
  double* basepointer = &num_mat.matrix0[nrows*nrows - n*n];
  for (unsigned i = 0; i < n; ++i) cblas_dcopy(n, &copyA[i*n], 1, &basepointer[i], n);

  /* Compute B_A_T^{-1}(B_A_T^{-1}M)_T */
  DGETRS(LA_NOTRANS, n, n, B_A_T, n, ipiv, basepointer, n, &infoLAPACK);
  assert(infoLAPACK == 0  && "avi_caoferris :: info from DGETRS for solving B_A_T X = (B_A_T^{-1}M)_T is not zero!\n");

  DEBUG_PRINT("B_A_T^{-1}(B_A_T^{-1}M)_T\n");
  DEBUG_PRINT_MAT(basepointer, n, n);

  for (unsigned i = 0; i < n; ++i) cblas_dcopy(n, &basepointer[n*i], 1, &copyA[i], n);

  DEBUG_PRINT_VEC_STR("b_I =: q[n:nrows]", (&lcplike_pb.q[n]), n_I);
  /* partial construction of q: q[n:nrows] += (B_A_T^{-1}*B_I_T)_T K_A */
  cblas_dgemv(CblasColMajor, CblasTrans, n_I, n, 1.0, B_I_T, n_I, d_vec, 1, 1.0, &lcplike_pb.q[n], 1);
  DEBUG_PRINT_VEC_STR("final q[n:nrows] as b_I + B_I B_A^{-1}b_A", (&lcplike_pb.q[n]), n_I);

  /* Compute B_A_T^{-1} M B_A^{-1} K_A 
   * We have to set CblasTrans since we still have a transpose */
  /* XXX It looks like we could have 2 here, but not it does not work w/ it. Investigate why -- xhub  */
  cblas_dgemv(CblasColMajor, CblasTrans, n, n, 1.0, basepointer, n, d_vec, 1, 0.0, lcplike_pb.q, 1);
  DEBUG_PRINT_VEC_STR("B_A_T^{-1} M B_A^{-1} K_A =: q[0:n]", lcplike_pb.q, n);

  /* q[0:n] = 2 B_A_T^{-1} A B_A^{-1}b_A  + B_A_T{-1} q_{AVI} */
  /*  XXX about the + or -: we do not follow the convention of Cao & Ferris */
  cblas_daxpy(n, 1.0, a_bar, 1, lcplike_pb.q, 1);
  DEBUG_PRINT("final q\n");
  DEBUG_PRINT_VEC(lcplike_pb.q, nrows);

  /* q is now ready, let's deal with M */

  /* set some pointers to sub-matrices */
  double* upper_left_mat = num_mat.matrix0;
  double* upper_right_mat = &num_mat.matrix0[n*nrows];
  double* lower_left_mat = &num_mat.matrix0[n];
  double* lower_right_mat = &upper_right_mat[n];


  /* copy the B_A_T^{-1} B_I_T (twice) and set the lower-left part to 0*/
  for (unsigned i = 0, j = 0, k = 0; i < n_I; ++i, j += n_I, k += nrows)
  {
    cblas_dcopy(n, &copyA[n*i], 1, &upper_right_mat[k], 1);/* copy into the right location B_A_T^{-1} M B_A^{-1} */
    cblas_dcopy(n_I, &B_I_T[j], 1, &upper_left_mat[k], 1); /* copy B_A_T^{-1}*B_I_T to the upper-right block */
    cblas_dscal(n, -1.0, &upper_left_mat[k], 1); /* take the opposite of the matrix */
    cblas_dcopy(n_I, &B_I_T[j], 1, &lower_right_mat[i], nrows); /*  copy B_IB_A^{-1} to the lower-left block */
    memset(&lower_left_mat[k], 0, sizeof(double)*(n_I)); /* set the lower-left block to 0 */
  }

  DEBUG_PRINT_MAT(num_mat.matrix0, nrows, nrows);


  /* Matrix M is now ready */

  /* Save K_A */
  double* K_A = a_bar;
  cblas_dcopy(n, d_vec, 1, K_A, 1);
  DEBUG_PRINT_VEC(K_A, n);
  /* We put -1 because we directly copy it in stage 3 */
  for (unsigned int i = 0; i < n; ++i) d_vec[i] =  -1.0;
  memset(&d_vec[n], 0, n_I*sizeof(double));

  DEBUG_PRINT_VEC_INT_STR("Active set", A, n);
  double* u_vec = (double *)calloc(nrows, sizeof(double));
  double* s_vec = (double *)calloc(nrows, sizeof(double));
  /* Call directly the 3rd stage 
   * Here w is used as u and z as s in the AVI */
  int info = avi_caoferris_stage3(&lcplike_pb, u_vec, s_vec, d_vec, n, A, options);

  /* Update z  */
  /* XXX why no w ?  */
  DEBUG_PRINT_VEC_INT(A, n);
  for (unsigned i = 0; i < n; ++i) z[i] = s_vec[A[i]-1] + K_A[i];
  DEBUG_PRINT_VEC_STR("s_A + K_A", z, n);
  DGETRS(LA_TRANS, n, 1, B_A_T, n, ipiv, z, n, &infoLAPACK);
  assert(infoLAPACK == 0  && "avi_caoferris :: info from DGETRS for solving B_A X = s_A + K_A is not zero!\n");
  DEBUG_PRINT_VEC_STR("solution z", z, n);

  /* free allocated stuff */
  free(u_vec);
  free(s_vec);
  free(A);
  free(basis);
  free(d_vec);
  free(B_I_T);
  free(copyA);
  free(B_A_T);
  freeNumericsMatrix(lcplike_pb.M);
  free(lcplike_pb.q);
  free(a_bar);

  return info;
}

int avi_caoferris_stage3(LinearComplementarityProblem* problem, double* restrict u , double* restrict s, double* restrict d, unsigned size_x, unsigned* restrict A, SolverOptions* options)
{

  /* returned value */
  int info = 0;
  assert(size_x > 0);
  /* matrix M of the avi */
  assert(problem);
  assert(problem->M);
  double * M = problem->M->matrix0;
  assert(M);
  /* size of the AVI */

  unsigned int dim = problem->size;
  assert(dim>0);
  unsigned int dim2 = 2 * (dim + 1);

  unsigned int drive = dim+1;
  unsigned int drive_number = dim+1;
  int block = -1;
  unsigned int has_sol = 0;
  unsigned int nb_iter = 0;
  unsigned int leaving;
  unsigned int itermax = options->iparam[0];
  unsigned aux_indx = 0;

  double pivot, pivot_inv;
  double tmp;
  unsigned int* basis;
  unsigned int* u_indx;
  unsigned int* s_indx;
  double* mat;

  /*output*/

  options->iparam[1] = 0;

  /* Allocation */
  basis = (unsigned int *)malloc(dim * sizeof(unsigned int));
  mat = (double *)malloc(dim * dim2 * sizeof(double));

  u_indx = (unsigned int *)malloc(dim * sizeof(unsigned int));
  s_indx = (unsigned int *)malloc(dim * sizeof(unsigned int));


  /* construction of mat matrix such that
   * mat = [ q | Id | -d | -M ] with d_i = i in A ? 1 : 0
   */

  /* We need to init only the part corresponding to Id */
  for (unsigned int i = 0 ; i < dim; ++i)
    for (unsigned int j = dim ; j <= dim*dim; j += dim)
      mat[i + j] = 0.0;

  /*  Copy M but mat[dim+2:, :] = -M */
  for (unsigned int i = 0 ; i < dim; ++i)
    for (unsigned int j = 0 ; j < dim*dim; j += dim)
      mat[i + j + dim*(dim + 2)] = -M[j + i]; // Siconos is in column major

  assert(problem->q);

  for (unsigned int i = 0, j = dim ; i < dim; ++i, j += dim)
  {
    mat[i] = problem->q[i];
    mat[i + j] =  1.0;
  }

  /** Add covering vector */
  assert(d);
  for (unsigned int i = 0; i < size_x  ; ++i) mat[i + dim*(dim + 1)] = d[i];
  for (unsigned int i = size_x; i < dim; ++i) mat[i + dim*(dim + 1)] = 0.0;


  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("% 2.2e ", mat[i + j*dim]) }
      DEBUG_PRINT("\n")});
  /* End of construction of mat */


  unsigned int val_A = A[0];

  /** Contruct the basis and the index maps for u and s
   * basis = (u_A, s_I) and nonbasic variables are (u_I, s_A)*/
  for (unsigned int i = 0, indx_A = 0, indx_I = 0; i < dim; ++i)
  {
    if (i == val_A-1) // i is in A, u_i is basic and s_i is nonbasic
    {
      basis[indx_A] = val_A;
      u_indx[i] = indx_A + 1;
      s_indx[i] = dim2 - size_x + indx_A;
      if (++indx_A < size_x)
        val_A = A[indx_A];
    }
    else // i is not in A (therefore in I), s_i is basic, u_i is nonbasic
    {
      basis[size_x+indx_I] = dim + 2 + i;
      u_indx[i] = dim + 2 + indx_I;
      s_indx[i] = size_x + indx_I + 1;
      ++indx_I;
    }
  }
  DEBUG_PRINT("basis u_indx s_indx\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i %i %i\n", basis[i], u_indx[i], s_indx[i]) });

  /* Start research of argmin lexico
   * lexicographic order is simple: we just look for the min of index in case
   * of tie */
  /* With this first step the covering vector enter in the basis */

  block = pivot_init_lemke(mat, size_x);

  /* Stop research of argmin lexico */


  if (block == -1)
  {
    /** exit, the solution is at hand with the current basis */
    DEBUG_PRINT("Trivial solution\n");
    has_sol = 1;
    goto exit_caoferris;
  }
  /* save the position of the auxiliary variable */
  aux_indx = drive;

  /* Pivot < mu , driver > */

  DEBUG_PRINTF("Pivoting %i and %i\n", block, drive);
  pivot = mat[block + drive*dim];
  pivot_inv = 1.0/pivot;

  /* Update column mat[block, :] */
  mat[block + drive*dim] = 1;
  for (unsigned int i = 0        ; i < drive*dim; i += dim) mat[block + i] *= pivot_inv;
  for (unsigned int i = dim*(drive + 1); i < dim2*dim ; i += dim) mat[block + i] *= pivot_inv;

  /* Update other columns*/
  for (int i = 0; i < block; ++i)
  {
    tmp = mat[i + drive*dim];
    for (unsigned int j = 0; j < dim2*dim; j += dim) mat[i + j] -= tmp*mat[block + j];
  }
  for (unsigned int i = block + 1; i < dim; ++i)
  {
    tmp = mat[i + drive*dim];
    for (unsigned int j = 0; j < dim2*dim; j += dim) mat[i + j] -= tmp*mat[block + j];
  }

  /** one basic u is leaving and mu enters the basis */
  leaving = basis[block];
  DEBUG_PRINTF("leaving variable: %d\n", leaving);
  basis[block] = drive;

  DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("% 2.2e ", mat[i + j*dim]) }
      DEBUG_PRINT("\n")});

  while (nb_iter < itermax && !has_sol)
  {

    ++nb_iter;

    if (leaving < dim + 1)
    {
      drive_number = leaving + dim + 1;
      drive = s_indx[leaving-1];
    }
    else if (leaving > dim + 1)
    {
      drive_number = leaving - (dim + 1);
      drive = u_indx[leaving - (dim + 2)];
    }

    DEBUG_PRINTF("driving variable %i \n", drive_number);
    assert(drive_number < dim2);

    block = pivot_selection_lemke(mat, dim, drive, aux_indx);

    if (block == -1) break;

    if (basis[block] == dim + 1) has_sol = 1;

    /* Pivot < block , drive > */
    DEBUG_PRINTF("Pivoting %i and %i\n", block, drive);

    pivot = mat[block + drive*dim];
    pivot_inv = 1.0/pivot;

    /* Update column mat[block, :] */
    mat[block + drive*dim] = 1;
    for (unsigned int i = 0        ; i < drive*dim; i += dim) mat[block + i] *= pivot_inv;
    for (unsigned int i = dim*(drive + 1); i < dim2*dim ; i += dim) mat[block + i] *= pivot_inv;

    /* Update other columns*/
    for (int i = 0; i < block; ++i)
    {
      tmp = mat[i + drive*dim];
      for (unsigned int j = 0; j < dim2*dim; j += dim) mat[i + j] -= tmp*mat[block + j];
    }
    for (unsigned int i = block + 1; i < dim; ++i)
    {
      tmp = mat[i + drive*dim];
      for (unsigned int j = 0; j < dim2*dim; j += dim) mat[i + j] -= tmp*mat[block + j];
    }

    /** one basic variable is leaving and driver enters the basis */
    leaving = basis[block];
    basis[block] = drive_number;

  DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("% 2.2e ", mat[i + j*dim]) }
      DEBUG_PRINT("\n")});
  } /* end while*/

  DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("% 2.2e ", mat[i + j*dim]) }
      DEBUG_PRINT("\n")});

exit_caoferris:

  DEBUG_PRINT("final basis\n");
  DEBUG_PRINT_VEC_INT(basis, dim);
  /* Recover solution */
  for (unsigned int i = 0 ; i < dim; ++i)
  {
    drive = basis[i];
    if (drive < dim + 1)
    {
      s[drive - 1] = 0.0;
      u[drive - 1] = mat[i];
    }
    else if (drive > dim + 1)
    {
      s[drive - dim - 2] = mat[i];
      u[drive - dim - 2] = 0.0;
    }
  }

  DEBUG_PRINT("u s\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%e %e\n", u[i], s[i]) });

  options->iparam[1] = nb_iter;

  if (has_sol) info = 0;
  else info = 1;

  free(basis);
  free(u_indx);
  free(s_indx);

  free(mat);

  return info;
}
