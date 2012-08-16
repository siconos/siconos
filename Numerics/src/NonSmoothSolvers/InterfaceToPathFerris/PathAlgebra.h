#ifndef PATHALGEBRA_H
#define PATHALGEBRA_H

/*!\file PathAlgebra.h
  \brief functions and tools used to convert matrices between Numerics and Path format
  \author Franck Perignon , 26/05/2008
*/



#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** Numerics dense matrix to Path Sparse matrix
      \param size0 number of rows of input matrix
      \param size1 number of columns of input matrix
      \param matIn matrix to convert (col. major)
      \param col_start
      \param col_len
      \param row
      \param data
  */
  void convertToPathSparse(int size0, int size1, double* matIn,
                           int* col_start, int* col_len, int* row, double* data);
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif
#endif
