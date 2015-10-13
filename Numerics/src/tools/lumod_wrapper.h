
#ifndef LUMOD_WRAPPER_H
#define LUMOD_WRAPPER_H

#include "NumericsMatrix.h"
#include "assert.h"

#define SN_LUMOD_NEED_REFACTORIZATION 1

typedef struct {
  unsigned n; /**< size of the matrix H*/
  unsigned maxmod; /**< maximum number of changes */
  unsigned k; /**< number of rows (or columns) of C */
  double* LU_H; /**< LU factors of the initial matrix H */
  int* ipiv_LU_H; /**< pivot for the LU factorization of H*/
  unsigned* factorized_basis; /**< basis when H was factorized and storing for the info where the columns of the non basic variables are in U and when a basic variable exited */
  int* row_col_indx; /**< Store the information to which column or row the variable correspond */
  double* Uk; /**< matrix which keeps track of modified columns in Hk*/
  double* Yk; /**< matrix of changed columns in Hk (row major)*/
  double* L_C; /**< L matrix of the LU factorization for Ck*/
  double* U_C; /**< U matrix of the LU factorization  for Ck*/
  double* y; /**< y vector for LUMOD */
  double* z; /**< z vector for LUMOD */
  double* w; /**< work vector for LUMOD */
} SN_lumod_dense_data;

static inline unsigned SN_lumod_need_refactorization(int info) { return (info == SN_LUMOD_NEED_REFACTORIZATION ? 1 : 0); }

static inline unsigned SN_lumod_find_arg_var(int* array, int p, unsigned n)
{
  for (unsigned i = 0; i < 2*n + 1; ++i)
  {
    if (array[i] == p) return i;
  }
  assert(0 && "We should not be here");
  return 0;
}

SN_lumod_dense_data* SN_lumod_dense_allocate(unsigned n, unsigned maxmod);
void SM_lumod_dense_free(SN_lumod_dense_data* lumod_data);
int SN_lumod_dense_solve(SN_lumod_dense_data* lumod_data, double* x, double* col_tilde);
void SN_lumod_add_row_col(SN_lumod_dense_data* lumod_data, unsigned leaving_indx, double* col);
void SN_lumod_replace_col(SN_lumod_dense_data* lumod_data, unsigned index_col, double* col);
void SN_lumod_replace_row(SN_lumod_dense_data* lumod_data, unsigned index_row, unsigned leaving_indx);
void SN_lumod_delete_row_col(SN_lumod_dense_data* lumod_data, unsigned index_row, unsigned index_col);
int SN_lumod_factorize(SN_lumod_dense_data* lumod_data, unsigned* basis, NumericsMatrix* M, double* covering_vector);

#endif
