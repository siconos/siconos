#include <stdio.h>
#include <stdlib.h>
#include "SparseBlockMatrix.h"
#include "LA.h"

void prodSBM(int size, double alpha, const SparseBlockStructuredMatrix* const A, const double* const x, double beta, double* y)
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
  */

  DSCAL(size, beta, y, incy);

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
    DGEMV(LA_NOTRANS, nbRows, nbColumns, alpha, A->block[blockNum], nbRows, &x[posInX], incx, 1.0, &y[posInY], incy);
  }

}

void subRowProdSBM(int sizeX, int sizeY, int currentRowNumber, const SparseBlockStructuredMatrix* const A, const double* const x, double* y, int init)
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

  /* Check if currentRowNumber fits with A dimensions */
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

void rowProdNoDiagSBM(int sizeX, int sizeY, int currentRowNumber, const SparseBlockStructuredMatrix* const A, const double* const x, double* y, int init)
{
  /*
     If: A is a SparseBlockStructuredMatrix matrix, Aij a block at row i and column j
     (Warning: i and j are indices of block position, not scalar component positions)

     Then rowProdNoDiagSBM computes
     y = sum for i not equal to j of Aij.xj over a row of blocks
     (or += if init = false)

     currentRowNumber represents the position (block number) of the required line of blocks in the matrix A.

  */

  if (A == NULL || x == NULL || y == NULL)
  {
    fprintf(stderr, "Numerics, SparseBlockMatrix, product sparse matrix - vector rowProdNoDiagSBM(A,x,y) failed, A, x or y is a null pointer\n");
    exit(EXIT_FAILURE);
  }

  int sizeOfA = A->blocksize[A->size - 1];
  /* Checks sizes */
  if (sizeX != sizeOfA)
  {
    fprintf(stderr, "Numerics, SparseBlockMatrix, product sparse matrix - vector rowProdNoDiagSBM(A,x,y) failed. Unconsistent sizes.\n");
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

  /* Check if currentRowNumber fits with A dimensions */
  if (currentRowNumber > A->size)
  {
    fprintf(stderr, "Numerics, SparseBlockMatrix, product (Row of a sparse matrix)-vector rowProdNoDiagSBM(rowPos,A,x,y) failed, rowPos is out of range.\n");
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
    fprintf(stderr, "Numerics, SparseBlockMatrix, product (Row of a sparse matrix)-vector rowProdNoDiagSBM(rowPos,A,x,y) failed, Unconsistent sizes between rowA and y.\n");
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

    /* Computes product only for extra diagonal blocks */
    if (colNumber != currentRowNumber)
    {
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
    }
    /* Steps to next block*/
    blockNum ++;
    nextRowNumber = A->RowIndex[blockNum];
  }
}

void freeSBM(SparseBlockStructuredMatrix *blmat)
{
  /* Free memory for SparseBlockStructuredMatrix */
  /* Warning: nothing is done to check if memory has really been allocated for each component
   or if it was only pointer links.
   Note that when used from the Kernel, memory is not directly allocated for such structures and
   this function must not be called. See in Kernel/src/simulationsTools/SparseBlockMatrix.cpp for details on
   the way the structure is filled in.
  */
  if (blmat->RowIndex != NULL)
    free(blmat->RowIndex);
  if (blmat->ColumnIndex != NULL)
    free(blmat->ColumnIndex);
  if (blmat->blocksize != NULL)
    free(blmat->blocksize);
  for (int i = 0 ; i < blmat->nbblocks ; i++)
  {
    if (blmat->block[i] != NULL)
      free(blmat->block[i]);
  }
  if (blmat->block != NULL)
    free(blmat->block);
}

void printSBM(const SparseBlockStructuredMatrix* const m)
{
  if (m == NULL)
  {
    fprintf(stderr, "Numerics, SparseBlockStructuredMatrix display failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  if (m->size == 0)
  {
    printf("Numerics, SparseBlockStructuredMatrix display: matrix dim = 0.");
    return;
  }
  int size = m->blocksize[m->size - 1];
  printf("Sparse-Block structured matrix of size %dX%d, with %d blocks in a row(or column)\n", size, size, m->size);
  printf("and %d non null blocks\n", m->nbblocks);
  printf("Diagonal blocks sizes = [");
  for (int i = 0; i < m->size; i++)
  {
    size = m->blocksize[i];
    if (i != 0) size -= m->blocksize[i - 1];
    printf("%d ", size);
  }
  printf("]\n");
  printf("Indices of blocks row-positions = [");
  for (int i = 0; i < m->nbblocks; i++)
    printf("%d ", m->RowIndex[i]);
  printf("]\n");
  printf("Indices of blocks column-positions = [");
  for (int i = 0; i < m->nbblocks; i++)
    printf("%d ", m->ColumnIndex[i]);
  printf("]\n");
  int dimblock, row, col, posR, posC;
  printf("Non null blocks (warning: column-major storage)\n");
  for (int i = 0; i < m->nbblocks; i++)
  {
    posR = m->RowIndex[i];
    posC = m->ColumnIndex[i];
    row = m->blocksize[posR];
    if (posR != 0) row -= m->blocksize[posR - 1];
    col = m->blocksize[posC];
    if (posC != 0) col -= m->blocksize[posC - 1];
    dimblock = row * col;
    printf("block %d,%d (dim:%dX%d)= [", posR, posC, row, col);
    for (int j = 0; j < dimblock; ++j)
      printf("%lf ", m->block[i][j]);
    printf("]\n");
  }
}

void freeSpBlMatPred(SparseBlockStructuredMatrixPred *blmatpred)
{

  int i;

  for (i = 0 ; i < blmatpred->nbbldiag ; i++)
  {
    free(blmatpred->indic[i]);
    free(blmatpred->indicop[i]);
    free(blmatpred->submatlcp[i]);
    free(blmatpred->submatlcpop[i]);
    free(blmatpred->ipiv[i]);
    free(blmatpred->subq[i]);
    free(blmatpred->bufz[i]);
    free(blmatpred->newz[i]);
    free(blmatpred->workspace[i]);
  }
  free(blmatpred->indic);
  free(blmatpred->indicop);
  free(blmatpred->submatlcp);
  free(blmatpred->submatlcpop);
  free(blmatpred->ipiv);
  free(blmatpred->sizesublcp);
  free(blmatpred->sizesublcpop);
  free(blmatpred->subq);
  free(blmatpred->bufz);
  free(blmatpred->newz);
  free(blmatpred->workspace);

}

int getDiagonalBlockPos(const SparseBlockStructuredMatrix* const M, int num)
{
  // Todo: save positions of diagonal blocks in the structure?
  int pos = 0;
  int blockNum = 0;
  /* Look for the first block of row number num */
  while (M->RowIndex[blockNum] != num) blockNum++;
  pos = blockNum;
  while (M->ColumnIndex[pos] != M->RowIndex[blockNum]) pos++;
  return pos;
}
