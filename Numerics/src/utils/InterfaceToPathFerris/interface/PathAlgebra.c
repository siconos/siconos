#include "PathAlgebra.h"
#include "NumericsOptions.h"
#include "math.h"

/** tolerance value for zero */
static double zeroTol = 1e-15;

void convertToPathSparse(int size0, int size1, double* matIn, int* col_start, int* col_len, int* row, double* data)
{

  int pos = 0;
  col_start[0] = 1;
  for (int j = 0; j < size1 ; ++j)
  {
    if (j > 0)
    {
      col_start[j] = col_start[j - 1] + col_len[j - 1];
      if (col_start[j] == col_start[j - 1])
        numericsError("PathAlgebra::convertToPathSparse()", "Null column in input matrix");
    }
    col_len[j] = 0;
    for (int i = 0; i < size0 ; ++i)
    {
      if (fabs(matIn[i + j * size0]) > zeroTol)
      {
        data[pos] = matIn[i + j * size0];
        col_len[j] ++;
        row[pos] = i + 1; // Warning: indices start from 1 for row
        pos++;
      }
    }
  }
}
