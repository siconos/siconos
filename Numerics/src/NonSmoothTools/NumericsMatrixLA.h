

#ifndef NM_LA_H
#define NM_LA_H

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif


static inline double* NMD_row_rmajor(double* restrict mat, unsigned ncols, unsigned rindx)
{ 
  return &mat[rindx*ncols];
}

static inline void NMD_copycol_rmajor(unsigned nrows, double* col, double* restrict mat, unsigned ncols, unsigned cindx)
{
  cblas_dcopy(nrows, col, 1, &mat[cindx], ncols);
}

static inline void NMD_dense_gemv(unsigned nrows, unsigned ncols, double alpha, double* restrict mat, double* restrict y, double beta, double* restrict x)
{
cblas_dgemv(CblasColMajor, CblasTrans, ncols, nrows, alpha, mat, ncols, y, 1, beta, x, 1);
}
#endif
