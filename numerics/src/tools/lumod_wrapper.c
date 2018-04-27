/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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


#include "lumod_wrapper.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "SiconosBlas.h"
#include "SiconosLapack.h"
#include "lumod_dense.h"


//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"


#define TOL_BLU 1e-30
#define BASIS_OFFSET 1
#define DEEP_DEBUG_Ck
#include <fenv.h>
#include "pivot-utils.h"

/* TODO :
 *  - Yk should be row-major (be careful with gemv and lda)
 *  - in solve, we should try to better handle the resolution with C
 *  - the solution to Hx = col should be saved only when col correspond to a
 *  driving variable not in the factorized basis*/ 

#ifdef DEBUG_MESSAGES
inline static void lumod_full_check(SN_lumod_dense_data* lumod_data)
{
#ifdef DEEP_DEBUG_Ck
  unsigned n = lumod_data->n;
  unsigned k = lumod_data->k;
  unsigned maxmod = lumod_data->maxmod;
  assert(k > 0 && "lumod_full_check k is 0, nothing to do!");

  double* Ctmp = (double*)calloc(k*k, sizeof(double));
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, k, k, n, 1.0, lumod_data->Uk, n, lumod_data->Yk, n, 0., Ctmp, k);
  for (unsigned i = 0; i < k; ++i)
  {
    memset(lumod_data->y, 0, k*sizeof(double));
    lumod_data->y[i] = 1.;
    Lprod_dense(1, maxmod, k, lumod_data->L_C, lumod_data->y, lumod_data->z);
    /*          solve U y = z */
    Usolve_dense(1, maxmod, k, lumod_data->U_C-1, lumod_data->z-1);
    /* real check */
    cblas_dgemv(CblasColMajor, CblasNoTrans, k, k, 1.0, Ctmp, k, lumod_data->z, 1, 0., lumod_data->y, 1);
#define TOL_CHECK 1e-10
    unsigned failcount = 0;
    for (unsigned j = 0; j < k; ++j)
    {
      if ((j == i) && (fabs(lumod_data->y[j] - 1.) > TOL_CHECK)) { ++failcount; }
      else if ((j != i) && (fabs(lumod_data->y[j]) > TOL_CHECK)) { ++failcount; }
    }
    if (failcount > 0)
    {
      DEBUG_PRINT_MAT_STR("Ck", Ctmp, k, k);
      DEBUG_PRINT_VEC_STR("ep", lumod_data->y, k);
      DEBUG_PRINT_VEC_STR("solution", lumod_data->z, k);
    }
    assert(failcount == 0);
  }

  free(Ctmp);
#endif
}
#endif
/* Wrapper on dense matrix */
SN_lumod_dense_data* SN_lumod_dense_allocate(unsigned n, unsigned maxmod)
{
  SN_lumod_dense_data* lumod_data = (SN_lumod_dense_data*)malloc(sizeof(SN_lumod_dense_data));
  lumod_data->maxmod = maxmod;
  lumod_data->n = n;
  lumod_data->k = 0;

  /* Perform only one big allocation
   * Formula: size = H + Uk + Yk + L_C + U_C + y + z + w*/
  unsigned size_H = n*n;
  unsigned size_Uk = n*maxmod;
  unsigned size_Yk = n*maxmod;
  unsigned size_L_C = maxmod*maxmod;
  /* the original formula is ceil(maxmod*(maxmod + 1)/2) */
  unsigned size_U_C = maxmod*(maxmod + 1)/2 + 1;
  unsigned size_y = maxmod;
  unsigned size_z = maxmod;
  unsigned size_w = maxmod;
  unsigned size = size_H + size_Uk + size_Yk + size_L_C + size_U_C + size_y + size_z + size_w;
  double* data = (double*)calloc(size, sizeof(double));

  unsigned current_pointer = 0;
  /* H matrix */
  lumod_data->LU_H = &data[current_pointer];
  current_pointer += size_H;
  lumod_data->ipiv_LU_H = (lapack_int*)malloc(n*sizeof(lapack_int));
  lumod_data->factorized_basis = (unsigned*)malloc((4*n+2)*sizeof(unsigned));
  lumod_data->row_col_indx = (int*)malloc((2*n+1)*sizeof(int));

  /* matrices for the BLU updates */
  lumod_data->Uk = &data[current_pointer];
  current_pointer += size_Uk;
  lumod_data->Yk = &data[current_pointer];
  current_pointer += size_Yk;

  /* C matrix allocation. Sizes are given in lumod*/
  lumod_data->L_C = &data[current_pointer];
  current_pointer += size_L_C;
  lumod_data->U_C = &data[current_pointer];
  current_pointer += size_U_C;
  /* needed for LUMOD */
  lumod_data->y = &data[current_pointer];
  current_pointer += size_y;
  lumod_data->z = &data[current_pointer];
  current_pointer += size_z;
  lumod_data->w = &data[current_pointer];
  current_pointer += size_w;

  return lumod_data;
}

void SM_lumod_dense_free(SN_lumod_dense_data* lumod_data)
{
  assert(lumod_data->LU_H);
  free(lumod_data->LU_H);
  assert(lumod_data->ipiv_LU_H);
  free(lumod_data->ipiv_LU_H);
  assert(lumod_data->factorized_basis);
  free(lumod_data->factorized_basis);
  assert(lumod_data->row_col_indx);
  free(lumod_data->row_col_indx);
  /* Let's do things by the book */
  lumod_data->LU_H = NULL;
  lumod_data->ipiv_LU_H = NULL;
  lumod_data->factorized_basis = NULL;
  lumod_data->row_col_indx = NULL;
  lumod_data->Uk = NULL;
  lumod_data->Yk = NULL;
  lumod_data->L_C = NULL;
  lumod_data->U_C = NULL;
  lumod_data->y = NULL;
  lumod_data->z = NULL;
  lumod_data->w = NULL;
  free(lumod_data);
}

int SN_lumod_dense_solve(SN_lumod_dense_data* restrict lumod_data, double* restrict x, double* restrict col_tilde)
{
  /* Steps to solve H_k x = b:
   * 1. Solve H_0 x1 = b
   * 2. Solve C x2 = Uk^T x1
   * 3. Compute x3 = x1 - Yk x2
   * 4. Compute x = x3 + Uk x2
   */
  unsigned n = lumod_data->n;
  unsigned k = lumod_data->k;
  unsigned maxmod = lumod_data->maxmod;
  lapack_int infoLAPACK = 0;

  /*  Step 1. */
  DEBUG_PRINT_VEC_STR("col", x, n);
  DGETRS(LA_NOTRANS, n, 1, lumod_data->LU_H, n, lumod_data->ipiv_LU_H, x, n, &infoLAPACK);
  assert(infoLAPACK == 0  && "SN_lumod_solve :: info from DGETRS for solving H_0 X = b is not zero!\n");
  DEBUG_PRINT_VEC_STR("x1 sol to H x = col", x, n);

  /* Save H col_tilde = col for a possible BLU
   * Note that in practice, we need this only when the driving variable is not
   * in the basis used for factorisation*/
  if (col_tilde)
  {
    cblas_dcopy(n, x, 1, col_tilde, 1);
  }

  if (k > 0)
  {

    feclearexcept(FE_ALL_EXCEPT);

    DEBUG_EXPR_WE(double* Ctmp = (double*)calloc(k*k, sizeof(double));
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, k, k, n, 1.0, lumod_data->Uk, n, lumod_data->Yk, n, 0., Ctmp, k);
        DEBUG_PRINT_MAT_STR("Ck", Ctmp, k, k);
        free(Ctmp););
    /* Step 2. */
    /* Step 2.a Compute Uk^T x1  */
    /** this is just a permutation, could replace accordingly*/
    cblas_dgemv(CblasColMajor, CblasTrans, n, k, 1.0, lumod_data->Uk, n, x, 1, 0.0, lumod_data->y, 1);
    DEBUG_PRINT_VEC_STR("Uk^T x1", lumod_data->y, k);

    /* XXX: hack */
    if (cblas_ddot(k, lumod_data->y, 1, lumod_data->y, 1) < k*TOL_BLU)
    {
      goto exit_SN_lumod_dense_solve;
    }
    /* Step 2.b Solve C x2 = Uk^T x1 using the LU factors of C*/
    /*          compute z = Ly */
    Lprod_dense(1, maxmod, k, lumod_data->L_C, lumod_data->y, lumod_data->z);
    DEBUG_PRINT_VEC_STR("L Uk^T x1", lumod_data->z, k);
    /*          solve U y = z */
    Usolve_dense(1, maxmod, k, lumod_data->U_C-1, lumod_data->z-1);
    DEBUG_PRINT_VEC_STR("x2 sol to C x = Uk^T x1", lumod_data->z, k);

    if (fetestexcept(FE_ALL_EXCEPT & ~FE_INEXACT & ~FE_UNDERFLOW))
    {
      printf("solution of the small system is spurious!\n");
      return SN_LUMOD_NEED_REFACTORIZATION;
    }
    /* Step 3. Compute x3 = x1 - Yk x2 */
    /* Row-major looks broken on ATLAS ... */
    cblas_dgemv(CblasColMajor, CblasNoTrans, n, k, -1., lumod_data->Yk, n, lumod_data->z, 1, 1., x, 1);
    DEBUG_PRINT_VEC_STR("x3 = x1 - Yk x2", x, n);

    /* Step 4. Compute x = x3 + Uk x2  */
    cblas_dgemv(CblasColMajor, CblasNoTrans, n, k, 1.0, lumod_data->Uk, n, lumod_data->z, 1, 1.0, x, 1);
    DEBUG_PRINT_VEC_STR("x = x1 - Yk x2 + Uk x2", x, n);
  }
exit_SN_lumod_dense_solve:
  DEBUG_PRINT_VEC_STR("SN_lumod_dense_solve H x = col final", x, n);
  return infoLAPACK;
}

int SN_lumod_factorize(SN_lumod_dense_data* restrict lumod_data, unsigned* restrict basis, NumericsMatrix* restrict M, double* covering_vector)
{
  /* Construct the basis matrix  */
  unsigned n = lumod_data->n;
  assert(n > 0);
  double* H = lumod_data->LU_H;
  unsigned* factorized_basis = lumod_data->factorized_basis;
  double* Mlcp =  M->matrix0;
  assert(Mlcp);

  /* Reset the factorized_basis */
  memset(factorized_basis, 0, (2*n+1)*sizeof(unsigned));
  memset(lumod_data->row_col_indx, 0, (2*n+1)*sizeof(int));
  DEBUG_PRINT("Variables in factorized basis\n")

  for (unsigned i = 0, j = 0; i < n; ++i, j += n)
  {
    unsigned var = basis[i] - BASIS_OFFSET;
    DEBUG_PRINTF("%s%d ", basis_to_name(basis[i], n), basis_to_number(basis[i], n))
    factorized_basis[var] = i + BASIS_OFFSET;
    if (var > n) /* z var */
    {
      unsigned z_indx = var - n - 1;
      assert(var >=  n + 1 );
      assert(var - n - 1 < n);
      cblas_dcopy(n, &Mlcp[z_indx*n], 1, &H[j], 1);
    }
    else if (var < n)
    {
      unsigned w_indx = var;
      memset(&H[j], 0, n*sizeof(double));
      H[j + w_indx] = -1.;
    }
    else /* we have the auxiliary variable  */
    {
      cblas_dcopy(n, covering_vector, 1, &H[j], 1);
    }
  }
  DEBUG_PRINT("\n")

  DEBUG_PRINT_MAT_STR("basis for factorization", H, n, n);

  /*  Compute LU factorisation of basis */
  lapack_int infoDGETRF = 0;
  DGETRF(n, n, H, n, lumod_data->ipiv_LU_H, &infoDGETRF);
  if (infoDGETRF > 0)
  {
    printf("crash_pivot_basis :: the crash basis is singular, cannot inverse the matrix.\n\
        The (first) diagonal element of U to be 0.0 is %d\n\
        A remedy remains to be implemented! Do not hesitate to report to the developers!\n", infoDGETRF);
    return infoDGETRF;
  }
  else if (infoDGETRF < 0)
  {
    printf("SN_lumod_factorize :: wrong call to DGETRF: %d parameter has an illegal value", infoDGETRF);
    return infoDGETRF;
  }

  /* Reset C */
  lumod_data->k = 0;
  memset(lumod_data->Uk, 0, lumod_data->maxmod*n*sizeof(double));

  return 0;
}

void SN_lumod_add_row_col(SN_lumod_dense_data* restrict lumod_data, unsigned leaving_indx_in_H, double* restrict col)
{
  /* We have a few things to update now: Uk, Vk and the LU factorization of Ck.
   * - Uk+1 = (Uk, e_{leaving_indx})
   * - Vk+1 = (Vk, col)
   * - Ck+1 has to be update using LUMOD
   */

  unsigned n = lumod_data->n;
  unsigned k = lumod_data->k;
  unsigned maxmod = lumod_data->maxmod;

  /*  We first update Ck: we have to prepare
   *  y = ((Yk)_{leaving} col_{leaving})
   *  z = Uk^T col*/
  if (k > 0)
  {
    DEBUG_PRINT_MAT_STR("Uk before expansion", lumod_data->Uk, n, (k+1));
    DEBUG_PRINT_MAT_STR("Yk before expansion", lumod_data->Yk, n, (k+1));
    cblas_dcopy(k, &lumod_data->Yk[leaving_indx_in_H], n, lumod_data->y, 1);
    cblas_dgemv(CblasColMajor, CblasTrans, n, k, 1.0, lumod_data->Uk, n, col, 1, 0., lumod_data->z, 1);
    DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS_STR("L_C before expansion", lumod_data->L_C, k, maxmod, k);
  }

  lumod_data->y[k] = col[leaving_indx_in_H];

  /* remember that the dimension increases from k-1 to k */
  LUmod_dense(1, maxmod, k+1, 0, 0, lumod_data->L_C-1, lumod_data->U_C-1, lumod_data->y-1, lumod_data->z-1, lumod_data->w-1);

  /* Now update Yk and Uk */
  DEBUG_PRINT_VEC_STR("new col", col, n);
  cblas_dcopy(n, col, 1, &lumod_data->Yk[k*n], 1);
  lumod_data->Uk[leaving_indx_in_H + n*k] = 1.;

  DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS_STR("L_C after expansion", lumod_data->L_C, k+1, maxmod, k+1);
  DEBUG_PRINT_VEC_STR("U_C", lumod_data->U_C, 10);

  DEBUG_PRINT_MAT(lumod_data->Uk, n, (k+1));
  DEBUG_PRINT_MAT(lumod_data->Yk, n, (k+1));
  ++(lumod_data->k);

  DEBUG_EXPR_WE(lumod_full_check(lumod_data););

  /* costly check */
  DEBUG_EXPR_WE(for (unsigned i = 0; i < k*n; i += n) assert(cblas_ddot(n, &lumod_data->Uk[i], 1, &lumod_data->Uk[i], 1) == 1));
}

void SN_lumod_replace_col(SN_lumod_dense_data* restrict lumod_data,  unsigned index_col, double* restrict col)
{
  /*  We have to update Yk and Ck
   *  Yk: col(Yk, index_col) <- col
   *  Ck: col(Ck, index_col) <- Uk^T col
   */

  unsigned n = lumod_data->n;
  unsigned k = lumod_data->k;
  unsigned maxmod = lumod_data->maxmod;

  assert(k > 0);
  assert(index_col < k);

  /* Update a column of Yk */
  DEBUG_PRINT_MAT_STR("Yk before replacement", lumod_data->Yk, n, k);
  cblas_dcopy(n, col, 1, &lumod_data->Yk[index_col*n], 1);
  DEBUG_PRINT_MAT_STR("Yk after replacement", lumod_data->Yk, n, k);

  /* Update of Ck: the new column has to be in z = Uk^T col */
  cblas_dgemv(CblasColMajor, CblasTrans, n, k, 1.0, lumod_data->Uk, n, col, 1, 0., lumod_data->z, 1);

  DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS_STR("L_C before col replace", lumod_data->L_C, k, maxmod, k);
  LUmod_dense(2, maxmod, k, 0, index_col+1, lumod_data->L_C-1, lumod_data->U_C-1, lumod_data->y-1, lumod_data->z-1, lumod_data->w-1);
  DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS_STR("L_C after col replace", lumod_data->L_C, k, maxmod, k);

  DEBUG_EXPR_WE(lumod_full_check(lumod_data););
}

void SN_lumod_replace_row(SN_lumod_dense_data* restrict lumod_data,  unsigned index_row, unsigned leaving_indx_in_H)
{
  /*  We have to update Uk and Ck
   *  Uk: col(Uk, index_row) <- e_{leaving}
   *  Ck: row(Ck, index_row) <- row(Yk, leaving_indx)
   */

  unsigned n = lumod_data->n;
  unsigned k = lumod_data->k;
  unsigned maxmod = lumod_data->maxmod;
  assert(k > 0);
  assert(index_row < k);
  assert(leaving_indx_in_H < n);

  /* Update a column of Uk  */
  DEBUG_PRINTF("SN_lumod_replace_row :: replacing row %d in Uk", index_row);
  DEBUG_PRINT_MAT_STR("Uk before replacement", lumod_data->Uk, n, k);
  memset(&lumod_data->Uk[index_row*n], 0, n*sizeof(double));
  lumod_data->Uk[index_row*n + leaving_indx_in_H] = 1.;
  DEBUG_PRINT_MAT_STR("Uk after replacement", lumod_data->Uk, n, k);

  /* Update of Ck: the new row has to be in z = row(Yk, leaving_indx) */
  cblas_dcopy(k, &lumod_data->Yk[leaving_indx_in_H], n, lumod_data->y, 1);

  DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS_STR("L_C before row replace", lumod_data->L_C, k, maxmod, k);
  LUmod_dense(3, maxmod, k, index_row+1, 0, lumod_data->L_C-1, lumod_data->U_C-1, lumod_data->y-1, lumod_data->z-1, lumod_data->w-1);
  DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS_STR("L_C after row replace", lumod_data->L_C, k, maxmod, k);

  DEBUG_EXPR_WE(lumod_full_check(lumod_data););
}

void SN_lumod_delete_row_col(SN_lumod_dense_data* restrict lumod_data, unsigned index_row, unsigned index_col)
{
  /* Delete a column and a row from Ck. Delete also a column in Yk and Uk. */

  unsigned n = lumod_data->n;
  unsigned k = lumod_data->k;
  unsigned maxmod = lumod_data->maxmod;
  assert(k > 0);

  DEBUG_PRINTF("Deleting row %d and col %d\n", index_row, index_col);

  DEBUG_PRINT_MAT_STR("Uk before deletion", lumod_data->Uk, n, k);
  if ((index_row + 1) != k)
  {
    assert(index_row < k);
    cblas_dcopy(n, &lumod_data->Uk[(k-1)*n], 1, &lumod_data->Uk[(index_row)*n], 1);
  }
  memset(&lumod_data->Uk[(k-1)*n], 0, n*sizeof(double));
  DEBUG_PRINT_MAT_STR("Uk after deletion", lumod_data->Uk, n, k-1);

  DEBUG_PRINT_MAT_STR("Yk before deletion", lumod_data->Yk, n, k);
  /* Row major logic */
/*  if (index_col + 1 != k)
  {
    assert(index_col < k);
    unsigned len = k - 1;
    for (unsigned i = 0, indx_dest = index_col, indx_src = index_col+1; i < n; ++i, indx_dest += len, indx_src += k)
    {
      memmove(&lumod_data->Yk[indx_dest], &lumod_data->Yk[indx_src], len);
    }
  }*/
  if ((index_col + 1) != k)
  {
    assert(index_col < k);
    cblas_dcopy(n, &lumod_data->Yk[(k-1)*n], 1, &lumod_data->Yk[index_col*n], 1);
  }
  DEBUG_PRINT_MAT_STR("Yk after deletion", lumod_data->Yk, n, k-1);

  DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS_STR("L_C before suppr", lumod_data->L_C, k, maxmod, k);
  LUmod_dense(4, maxmod, k, index_row+1, index_col+1, lumod_data->L_C-1, lumod_data->U_C-1, lumod_data->y-1, lumod_data->z-1, lumod_data->w-1);
  DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS_STR("L_C after suppr", lumod_data->L_C, k-1, maxmod, k-1);
  --(lumod_data->k);
  DEBUG_EXPR_WE(if (k > 0) lumod_full_check(lumod_data););
}
