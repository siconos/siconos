#include <stdio.h>
#include <stdlib.h>
#include "SparseBlockMatrix.h"
#include "LA.h"


void prod(int size, const SparseBlockStructuredMatrix* const A, const double* const x, double* y, int init)
{
  /* Product SparseMat - vector, y = A*x (init = 1 = true) or y += A*x (init = 0 = false) */

  if (A == NULL || x == NULL || y == NULL)
  {
    fprintf(stderr, "Numerics, SparseBlockMatrix, product sparse matrix - vector prod(A,x,y) failed, A, x or y is a null pointer\n");
    exit(EXIT_FAILURE);
  }
  /* Checks sizes */
  if (size != A->blocksize[A->size - 1])
  {
    fprintf(stderr, "Numerics, SparseBlockMatrix, product sparse matrix - vector prod(A,x,y) failed. Unconsistent sizes between A and x or y\n");
    exit(EXIT_FAILURE);
  }

  /* Number of non-null blocks in the matrix */
  int nbblocks = A->nbblocks;
  /* Row (block) position of the current block */
  int currentRowNumber = 0;
  /* Column (block) position of the current block*/
  int colNumber = 0;
  /* Number of rows/columns of the current block */
  int nbRows, nbColumns;
  /* Position of the sub-block of x multiplied by the sub-block of A */
  int posInX = 0;
  /* Position of the sub-block of y, result of the product */
  int posInY = 0;

  int incx = 1, incy = 1;
  /* Loop over all non-null blocks
     Works whatever the ordering order of the block is, in A->block
     But it requires a set to 0 of all y components
  */

  /* Set y to 0 */
  if (init == 1)
    DSCAL(size, 0.0, y, incy);

  for (int blockNum = 0; blockNum < nbblocks; ++blockNum)
  {
    /* Get row/column position of the current block */
    currentRowNumber = A->RowIndex[blockNum];
    colNumber = A->ColumnIndex[blockNum];
    /* Get dim. of the current block */
    nbRows = A->blocksize[currentRowNumber];
    if (currentRowNumber != 0)
      nbRows -= A->blocksize[currentRowNumber - 1];

    nbColumns = A->blocksize[colNumber];
    if (colNumber != 0)
      nbColumns -= A->blocksize[colNumber - 1];

    /* Get position in x of the sub-block multiplied by A sub-block */
    posInX = 0;
    if (colNumber != 0)
      posInX += A->blocksize[colNumber - 1];
    /* Get position in y for the ouput sub-block, result of the product */
    posInY = 0;
    if (currentRowNumber != 0)
      posInY += A->blocksize[currentRowNumber - 1];
    /* Computes y[] += currentBlock*x[] */
    DGEMV(LA_NOTRANS, nbRows, nbColumns, 1.0, A->block[blockNum], nbRows, &x[posInX], incx, 1.0, &y[posInY], incy);
  }

}

void subRowProd(int sizeX, int sizeY, int currentRowNumber, const SparseBlockStructuredMatrix* const A, const double* const x, double* y, int init)
{
  /*
     Product (Row of blocks of a SparseMat) - vector, y = rowA*x (init = 1 = true) or y += rowA*x (init = 0 = false)
  */

  if (A == NULL || x == NULL || y == NULL)
  {
    fprintf(stderr, "Numerics, SparseBlockMatrix, product sparse matrix - vector prod(A,x,y) failed, A, x or y is a null pointer\n");
    exit(EXIT_FAILURE);
  }

  /* Checks sizes */
  if (sizeX != A->blocksize[A->size - 1])
  {
    fprintf(stderr, "Numerics, SparseBlockMatrix, product sparse matrix - vector prod(A,x,y) failed. Unconsistent sizes between A and x\n");
    exit(EXIT_FAILURE);
  }

  /* Number of non-null blocks in the matrix */
  int nbblocks = A->nbblocks;
  /* Column (block) position of the current block*/
  int colNumber = 0;
  /* Number of rows/columns of the current block */
  int nbRows, nbColumns;
  /* Position of the sub-block of x multiplied by the sub-block of A */
  int posInX = 0;

  int incx = 1, incy = 1;

  /* Check if currentRowNumber is lower or equal to the number of rows in A */
  if (currentRowNumber > A->size)
  {
    fprintf(stderr, "Numerics, SparseBlockMatrix, product (Row of a sparse matrix)-vector subRowProd(rowPos,A,x,y) failed, rowPos is out of range.\n");
    exit(EXIT_FAILURE);
  }

  /* Look for the first element of the wanted row */
  int blockNum = 0;
  while (A->RowIndex[blockNum] != currentRowNumber)
    blockNum++;

  /* Get dim (rows) of the current block */
  nbRows = A->blocksize[currentRowNumber];
  if (currentRowNumber != 0)
    nbRows -= A->blocksize[currentRowNumber - 1];

  if (nbRows != sizeY)
  {
    fprintf(stderr, "Numerics, SparseBlockMatrix, product (Row of a sparse matrix)-vector subRowProd(rowPos,A,x,y) failed, Unconsistent sizes between rowA and y\n");
    exit(EXIT_FAILURE);
  }

  /* Set y to 0, if required */
  if (init == 1)
    DSCAL(sizeY, 0.0, y, incy);

  /* Loop over all non-null blocks
     Works whatever the ordering order of the block is, in A->block
     But it requires a set to 0 of all y components
  */

  int nextRowNumber = currentRowNumber;
  while (nextRowNumber == currentRowNumber && blockNum < nbblocks)
  {
    /* Get row/column position of the current block */
    colNumber = A->ColumnIndex[blockNum];

    /* Get dim(columns) of the current block */
    nbColumns = A->blocksize[colNumber];
    if (colNumber != 0)
      nbColumns -= A->blocksize[colNumber - 1];

    /* Get position in x of the sub-block multiplied by A sub-block */
    posInX = 0;
    if (colNumber != 0)
      posInX += A->blocksize[colNumber - 1];
    /* Computes y[] += currentBlock*x[] */
    DGEMV(LA_NOTRANS, nbRows, nbColumns, 1.0, A->block[blockNum], nbRows, &x[posInX], incx, 1.0, y, incy);

    /* Step to next block*/
    blockNum ++;
    nextRowNumber = A->RowIndex[blockNum];
  }
}

void freeSpBlMat(SparseBlockStructuredMatrix *blmat)
{
  /* Free memory for SparseBlockStructuredMatrix */
  /* Warning: nothing is done to check if memory has really been allocated for each component
   or if it was only pointer links.
   Note that when used from the Kernel, memory is not directly allocated for such structures and
   this function must not be call. See in Kernel/src/simulationsTools/SparseBlockMatrix.cpp for details on
   the way the structure is filled in.
  */

  free(blmat->RowIndex);
  free(blmat->ColumnIndex);
  free(blmat->blocksize);
  /*    for (int i = 0 ; i < blmat->nbblocks ; i++) free(blmat->block[i]);  */
  free(blmat->block);
}

