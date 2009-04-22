#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "SparseBlockMatrix.h"
#include "LA.h"
#include <math.h>
void prodSBM(int sizeX, int sizeY, double alpha, const SparseBlockStructuredMatrix* const A, const double* const x, double beta, double* y)
{
  /* Product SparseMat - vector, y = A*x (init = 1 = true) or y += A*x (init = 0 = false) */

  assert(A);
  assert(x);
  assert(y);
  assert(A->blocksize0);
  assert(A->blocksize1);
  assert(A->index1_data);
  assert(A->index2_data);

  /* Checks sizes */
  assert(sizeX == A->blocksize1[A->blocknumber1 - 1]);
  assert(sizeY == A->blocksize0[A->blocknumber0 - 1]);

  /* Row (block) position of the current block */
  int currentRowNumber ;
  /* Column (block) position of the current block*/
  int colNumber;
  /* Number of rows/columns of the current block */
  int nbRows, nbColumns;
  /* Position of the sub-block of x multiplied by the sub-block of A */
  int posInX = 0;
  /* Position of the sub-block of y, result of the product */
  int posInY = 0;

  /* Loop over all non-null blocks
     Works whatever the ordering order of the block is, in A->block
  */
  DSCAL(sizeY, beta, y, 1);

  for (currentRowNumber = 0 ; currentRowNumber < A->filled1 - 1; ++currentRowNumber)
  {
    for (unsigned int blockNum = A->index1_data[currentRowNumber];
         blockNum < A->index1_data[currentRowNumber + 1]; ++blockNum)
    {
      assert(blockNum < A->filled2);

      colNumber = A->index2_data[blockNum];

      assert(colNumber < sizeX);

      /* Get dim. of the current block */
      nbRows = A->blocksize0[currentRowNumber];


      if (currentRowNumber != 0)
        nbRows -= A->blocksize0[currentRowNumber - 1];

      assert(nbRows <= sizeY & nbRows >= 0);

      nbColumns = A->blocksize1[colNumber];
      if (colNumber != 0)
        nbColumns -= A->blocksize1[colNumber - 1];

      assert(nbColumns <= sizeX & nbColumns >= 0);

      /* Get position in x of the sub-block multiplied by A sub-block */
      posInX = 0;
      if (colNumber != 0)
        posInX += A->blocksize1[colNumber - 1];
      /* Get position in y for the ouput sub-block, result of the product */
      posInY = 0;
      if (currentRowNumber != 0)
        posInY += A->blocksize0[currentRowNumber - 1];
      /* Computes y[] += currentBlock*x[] */
      DGEMV(LA_NOTRANS, nbRows, nbColumns, alpha, A->block[blockNum],
            nbRows, &x[posInX], 1, 1.0, &y[posInY], 1);
    }
  }
}
void prodSBMSBM(double alpha, const SparseBlockStructuredMatrix* const A, const SparseBlockStructuredMatrix* const B,  double beta, SparseBlockStructuredMatrix*  C)
{
  fprintf(stderr, "Numerics, SparseBlockStructuredMatrix, product matrix - matrix prodSBMSBM(alpha,A,B,beta,C) not yet implemented.\n");
  exit(EXIT_FAILURE);

}
void allocateMemoryForProdSBMSBM(const SparseBlockStructuredMatrix* const A, const SparseBlockStructuredMatrix* const B, SparseBlockStructuredMatrix*  C)
{
  fprintf(stderr, "Numerics, allocate memory for SparseBlockStructuredMatrix, product matrix - matrix  AllocateMemoryForProdSBMSBM(alpha,A,B,beta,C) not yet implemented.\n");
  exit(EXIT_FAILURE);

}

void subRowProdSBM(int sizeX, int sizeY, int currentRowNumber,
                   const SparseBlockStructuredMatrix* const A,
                   const double* const x, double* y, int init)
{
  /*
     Product (Row of blocks of a SparseMat) - vector, y = rowA*x (init
     = 1 = true) or y += rowA*x (init = 0 = false)
  */


  assert(A);
  assert(x);
  assert(y);

  /* Checks sizes */
  assert(sizeX == A->blocksize1[A->blocknumber1 - 1]);

  /* Number of non-null blocks in the matrix */
  int nbblocks = A->nbblocks;
  /* Column (block) position of the current block*/
  int colNumber = 0;
  /* Number of rows/columns of the current block */
  int nbRows, nbColumns;
  /* Position of the sub-block of x multiplied by the sub-block of A */
  int posInX = 0;

  /* Check if currentRowNumber fits with A dimensions */
  assert(currentRowNumber <= A->blocknumber0);

  /* Look for the first element of the wanted row */
  int blockNum = A->index1_data[currentRowNumber];

  /* Get dim (rows) of the current block */
  nbRows = sizeY;

  assert(
  {
    nbRows = A->blocksize0[currentRowNumber];
    if (currentRowNumber != 0)
      nbRows -= A->blocksize0[currentRowNumber - 1];
    nbRows == sizeY;
  });

  /* Set y to 0, if required */
  if (init == 1)
    DSCAL(sizeY, 0.0, y, 1);

  /* Loop over all non-null blocks
     Works whatever the ordering order of the block is, in A->block
     But it requires a set to 0 of all y components
  */

  for (unsigned int blockNum = A->index1_data[currentRowNumber];
       blockNum < A->index1_data[currentRowNumber + 1];
       ++blockNum)
  {
    /* Get row/column position of the current block */
    colNumber = A->index2_data[blockNum];

    /* Get dim(columns) of the current block */
    nbColumns = A->blocksize1[colNumber];
    if (colNumber != 0)
      nbColumns -= A->blocksize1[colNumber - 1];

    /* Get position in x of the sub-block multiplied by A sub-block */
    posInX = 0;
    if (colNumber != 0)
      posInX += A->blocksize0[colNumber - 1];
    /* Computes y[] += currentBlock*x[] */
    DGEMV(LA_NOTRANS, nbRows, nbColumns, 1.0, A->block[blockNum], nbRows, &x[posInX], 1, 1.0, y, 1);

  }
}

void rowProdNoDiagSBM(int sizeX, int sizeY, int currentRowNumber, const SparseBlockStructuredMatrix* const A, const double* const x, double* y, int init)
{
  /*
     If: A is a SparseBlockStructuredMatrix matrix, Aij a block at row
     i and column j (Warning: i and j are indices of block position,
     not scalar component positions)

     Then rowProdNoDiagSBM computes y = sum for i not equal to j of
     Aij.xj over a row of blocks (or += if init = false)

     currentRowNumber represents the position (block number) of the
     required line of blocks in the matrix A.

  */

  /* Number of non-null blocks in the matrix */
  int nbblocks = A->nbblocks;

  /* Column (block) position of the current block*/
  int colNumber = 0;

  /* Number of rows/columns of the current block */
  int nbRows, nbColumns;

  /* Position of the sub-block of x multiplied by the sub-block of A */
  int posInX = 0;

  /* Look for the first element of the wanted row */
  int blockNum = 0;

  /* Assertions */
  assert(A);
  assert(x);
  assert(y);
  assert(sizeX == A->blocksize1[A->blocknumber1 - 1]);
  assert(currentRowNumber <= A->blocknumber0);

  /* Get current block number */
  blockNum = A->index1_data[currentRowNumber];

  /* Get dim (rows) of the current block */
  nbRows = sizeY;

  assert(
  {
    nbRows = A->blocksize0[currentRowNumber];
    if (currentRowNumber != 0)
      nbRows -= A->blocksize0[currentRowNumber - 1];
    nbRows == sizeY ;
  });

  /* Set y to 0, if required */
  if (init == 1)
    DSCAL(sizeY, 0.0, y, 1);

  /* Loop over all non-null blocks. Works whatever the ordering order
     of the block is, in A->block, but it requires a set to 0 of all y
     components
  */
  for (unsigned int blockNum = A->index1_data[currentRowNumber];
       blockNum < A->index1_data[currentRowNumber + 1];
       ++blockNum)
  {
    /* Get row/column position of the current block */
    colNumber = A->index2_data[blockNum];

    /* Computes product only for extra diagonal blocks */
    if (colNumber != currentRowNumber)
    {
      /* Get dim(columns) of the current block */
      nbColumns = A->blocksize1[colNumber];
      if (colNumber != 0)
        nbColumns -= A->blocksize1[colNumber - 1];

      /* Get position in x of the sub-block multiplied by A sub-block */
      posInX = 0;
      if (colNumber != 0)
        posInX += A->blocksize0[colNumber - 1];
      /* Computes y[] += currentBlock*x[] */
      DGEMV(LA_NOTRANS, nbRows, nbColumns, 1.0, A->block[blockNum], nbRows, &x[posInX], 1, 1.0, y, 1);
    }
  }
}

void freeSBM(SparseBlockStructuredMatrix *blmat)
{
  /* Free memory for SparseBlockStructuredMatrix */
  /* Warning: nothing is done to check if memory has really been
   allocated for each component or if it was only pointer links.  Note
   that when used from the Kernel, memory is not directly allocated
   for such structures and this function must not be called. See in
   Kernel/src/simulationsTools/SparseBlockMatrix.cpp for details on
   the way the structure is filled in.
  */
  if (blmat->blocksize0)
  {

    free(blmat->blocksize0);
    if (blmat->blocksize0 == blmat->blocksize1)
    {
      blmat->blocksize1 = NULL ;
    }
  }
  if (blmat->blocksize1)
  {
    free(blmat->blocksize1);
  }
  for (int i = 0 ; i < blmat->nbblocks ; i++)
  {
    if (blmat->block[i])
    {
      free(blmat->block[i]);
    }
  }
  if (blmat->block)
    free(blmat->block);
}

void printSBM(const SparseBlockStructuredMatrix* const m)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, SparseBlockStructuredMatrix display failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  if (m->blocknumber0 == 0)
  {
    printf("Numerics, SparseBlockStructuredMatrix display: matrix dim = 0.");
    return;
  }
  assert(m->blocksize0);
  assert(m->blocksize1);

  int size0 = m->blocksize0[m->blocknumber0 - 1];
  int size1 = m->blocksize1[m->blocknumber1 - 1];
  printf("Sparse-Block structured matrix of size %dX%d elements, and  %dX%d blocks\n", size0, size1, m->blocknumber0, m->blocknumber1);
  printf("and %d non null blocks\n", m->nbblocks);
  printf("Diagonal blocks sizes = [ ");
  int diagonalblocknumber  = m->blocknumber1 + ((m->blocknumber0 - m->blocknumber1) & -(m->blocknumber0 < m->blocknumber1)); // min(m->blocknumber0,m->blocknumber1);
  for (int i = 0; i < diagonalblocknumber; i++)
  {
    size0 = m->blocksize0[i];
    if (i != 0) size0 -= m->blocksize0[i - 1];
    size1 = m->blocksize1[i];
    if (i != 0) size1 -= m->blocksize1[i - 1];
    printf("%dX%d ", size0, size1);
  }
  printf("]\n");
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
  /* Look for the first block of row number num */
  int pos;
  int firstBlockOfRow = M->index1_data[num];

  /* Look at the diagonal block */
  for (pos = firstBlockOfRow; M->index2_data[pos] != num; ++pos, assert(pos < M->filled2));

  return pos;
}
