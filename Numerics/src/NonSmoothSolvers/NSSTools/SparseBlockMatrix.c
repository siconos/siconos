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

void allocateMemoryForProdSBMSBM(const SparseBlockStructuredMatrix* const A, const SparseBlockStructuredMatrix* const B, SparseBlockStructuredMatrix*  C)
{

  assert(A);
  assert(B);
  assert(C);

  int i, j, k;
  /*     Check the compatibility of size of matrices */
  assert(A->blocksize0);
  assert(A->blocksize1);
  assert(B->blocksize0);
  assert(B->blocksize1);

  assert(A->blocksize1[A->blocknumber1 - 1] == B->blocksize0[B->blocknumber0 - 1]);

  /*     Check the compatibility of the number and the sizes of blocks */
  int compat = 1;

  if (A->blocknumber1 != B->blocknumber0) compat = 1; /* Compatibility of number of blocks */
  for (i = 0; i < A->blocknumber1; i++)
  {
    if (A->blocksize1[i] != B->blocksize0[i]) compat = 0; /* Compatibility of sizes of blocks */
  }
  if (!compat)
  {
    fprintf(stderr, "Numerics, allocate memory for SparseBlockStructuredMatrix, product matrix - matrix  AllocateMemoryForProdSBMSBM(alpha,A,B,beta,C) not implemented for non compatible blosk sizes.\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    /*      compute blocknumber and block sizes of C */
    C->blocknumber0 = A->blocknumber0;
    C->blocknumber1 = B->blocknumber1;
    C->blocksize0  = (int *)malloc(C->blocknumber0 * sizeof(int));
    C->blocksize1  = (int *)malloc(C->blocknumber1 * sizeof(int));
    for (i = 0; i < C->blocknumber0; i++) C->blocksize0[i] = A->blocksize0[i];
    for (j = 0; j < C->blocknumber1; j++) C->blocksize1[j] = B->blocksize1[j];
  }

  int currentRowNumberofA, currentRowNumberofB;
  int currentcolNumberofA, currentColNumberofB;
  int colNumberofA, rowNumberofB;
  int colNumberofB, rowNumberofA;

  /* Search  the non null blocks of C */

  /*     Creation of filled3,filled4 and index3_data,index4_data for B (indexation by column) */

  int Bfilled3 = B->blocknumber1 + 1;
  int * Bindex3_data = (int *)malloc(Bfilled3 * sizeof(int));
  int Bfilled4 = B->nbblocks;
  int * Bindex4_data = (int *)malloc(Bfilled4 * sizeof(int));
  int blockNumB = -1;
  Bindex3_data[0] = 0;
  int * blockMap  = (int *)malloc(Bfilled4 * sizeof(int));
  for (currentColNumberofB = 0 ; currentColNumberofB < Bfilled3 - 1; ++currentColNumberofB)
  {
    Bindex3_data[currentColNumberofB + 1] = Bindex3_data[currentColNumberofB];
    for (currentRowNumberofB = 0 ; currentRowNumberofB < B->filled1 - 1; ++currentRowNumberofB)
    {
      for (unsigned int blockNum = B->index1_data[currentRowNumberofB];
           blockNum < B->index1_data[currentRowNumberofB + 1]; ++blockNum)
      {
        assert(blockNum < B->filled2);
        colNumberofB = B->index2_data[blockNum];
        if (colNumberofB == currentColNumberofB)
        {
          blockNumB++;
          assert(blockNumB < B->nbblocks);
          Bindex3_data[currentColNumberofB + 1]++;
          Bindex4_data[blockNumB] = currentRowNumberofB;
          blockMap[blockNumB] = blockNum;

        }


      }

    }
  }


  /*    printf("\n"); */
  /*     for (i = 0 ; i< Bfilled3; i++) printf("Bindex3_data[%i]= %i\t", i,Bindex3_data[i] );printf("\n"); */
  /*     for (i = 0 ; i< Bfilled4; i++) printf("Bindex4_data[%i]= %i\t", i,Bindex4_data[i] );printf("\n"); */
  /*     for (i = 0 ; i< Bfilled4; i++) printf("blockMap[%i]= %i\t", i,blockMap[i] );printf("\n"); */
  /*     printf("\n"); */
  int colNumberAA, rowNumberBB;
  C->nbblocks = -1;
  C->filled2 = -1;
  /*     \warning The implementation is chosen to optimize cpu effort rather than memory. Otherwise a two loops are needed */
  int nbblocksmax = A->blocknumber0 * B->blocknumber1;

  double **Cblocktmp = (double**)malloc(nbblocksmax * sizeof(double*));
  size_t *Cindex2_datatmp = (size_t*)malloc(nbblocksmax * sizeof(size_t*));
  C->filled1 = C->blocknumber0 + 1;
  C->index1_data = (size_t*)malloc(C->filled1 * sizeof(size_t));
  C->index1_data[0] = 0;
  for (currentRowNumberofA = 0 ; currentRowNumberofA < A->filled1 - 1; ++currentRowNumberofA)
  {
    C->index1_data[currentRowNumberofA + 1] = C->index1_data[currentRowNumberofA];
    for (currentColNumberofB = 0 ; currentColNumberofB < Bfilled3 - 1; ++currentColNumberofB)
    {
      int BlockCexists = 0;

      for (unsigned int blockNumAA = A->index1_data[currentRowNumberofA];
           blockNumAA < A->index1_data[currentRowNumberofA + 1]; ++blockNumAA)
      {
        assert(blockNumAA < A->filled2);
        colNumberAA = A->index2_data[blockNumAA];

        for (unsigned int blockNumBB = Bindex3_data[currentColNumberofB];
             blockNumBB < Bindex3_data[currentColNumberofB + 1]; ++blockNumBB)
        {
          rowNumberBB = Bindex4_data[blockNumBB];
          if (rowNumberBB == colNumberAA)
          {
            BlockCexists = 1;
            C->nbblocks++;
            /*               printf("C block number %i exists for %i %i ",C->nbblocks, currentRowNumberofA, currentColNumberofB ); */

            C->filled2++;
            int Cblosksize0 = A->blocksize0[currentRowNumberofA];
            if (currentRowNumberofA != 0)
              Cblosksize0  -= A->blocksize0[currentRowNumberofA - 1];
            int Cblocksize1 = B->blocksize1[currentColNumberofB];
            if (currentColNumberofB != 0)
              Cblocksize1 -= B->blocksize1[currentColNumberofB - 1];
            /*               printf("of size %dX%d\n",Cblosksize0,Cblocksize1  ); */
            Cblocktmp[C->nbblocks] = (double*)malloc(Cblosksize0 * Cblocksize1 * sizeof(double));
            for (i = 0; i < Cblosksize0 * Cblocksize1; i++) Cblocktmp[C->nbblocks][i] = 0.0;

            C->index1_data[currentRowNumberofA + 1]++ ;
            Cindex2_datatmp[C->nbblocks] = currentColNumberofB ;
            break;

          };
        }
        if (BlockCexists) break;
      }

    }
  }
  C->nbblocks++;
  C->filled2++;
  assert(C->nbblocks ==  C->filled2);
  C->block = (double**)malloc(C->nbblocks * sizeof(double*));
  C->index2_data = (size_t*)malloc(C->nbblocks * sizeof(size_t*));

  for (i = 0 ; i < C->nbblocks; i++)  C->block[i] = Cblocktmp[i];
  for (i = 0 ; i < C->nbblocks; i++)  C->index2_data[i] = Cindex2_datatmp[i];
  free(Cblocktmp);
  free(Cindex2_datatmp);

  free(Bindex3_data);
  free(Bindex4_data);
  free(blockMap);

  /*   printSBM(C); */


  /*   fprintf(stderr,"Numerics, allocate memory for SparseBlockStructuredMatrix, product matrix - matrix  AllocateMemoryForProdSBMSBM(alpha,A,B,beta,C) not yet implemented.\n"); */
  /*   exit(EXIT_FAILURE); */

}
void prodSBMSBM(double alpha, const SparseBlockStructuredMatrix* const A, const SparseBlockStructuredMatrix* const B,  double beta, SparseBlockStructuredMatrix*  C)
{

  assert(A);
  assert(B);
  assert(C);

  int i, j, k;
  /*     Check the compatibility of size of matrices */
  assert(A->blocksize0);
  assert(A->blocksize1);
  assert(B->blocksize0);
  assert(B->blocksize1);

  assert(A->blocksize1[A->blocknumber1 - 1] == B->blocksize0[B->blocknumber0 - 1]);

  /*     Check the compatibility of the number and the sizes of blocks */
  int compat = 1;

  if (A->blocknumber1 != B->blocknumber0) compat = 1; /* Compatibility of number of blocks */
  for (i = 0; i < A->blocknumber1; i++)
  {
    if (A->blocksize1[i] != B->blocksize0[i]) compat = 0; /* Compatibility of sizes of blocks */
  }
  if (!compat)
  {
    fprintf(stderr, "Numerics, allocate memory for SparseBlockStructuredMatrix, product matrix - matrix  prodSBMSBM(alpha,A,B,beta,C) not implemented for non compatible blosk sizes.\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    assert(C->blocknumber0 == A->blocknumber0);
    assert(C->blocknumber1 == B->blocknumber1);

    assert( {int compat = 1;
             for (i = 0; i < C->blocknumber0; i++)
  {
    if (C->blocksize0[i] != A->blocksize0[i]) compat = 0;
    }
    compat == 1;
            });
    assert( {int compat = 1;
             for (j = 0; j < C->blocknumber1; j++)
  {
    if (C->blocksize1[j] != B->blocksize1[j]) compat = 0;
    }
    compat == 1;
            });
  }

  int currentRowNumberofA, currentRowNumberofB;
  int currentcolNumberofA, currentColNumberofB;
  int colNumberofA, rowNumberofB;
  int colNumberofB, rowNumberofA;

  /* Search  the non null blocks of C */

  /*     Creation of filled3,filled4 and index3_data,index4_data for B (indexation by column) */

  int Bfilled3 = B->blocknumber1 + 1;
  int * Bindex3_data = (int *)malloc(Bfilled3 * sizeof(int));
  int Bfilled4 = B->nbblocks;
  int * Bindex4_data = (int *)malloc(Bfilled4 * sizeof(int));
  int blockNumB = -1;
  Bindex3_data[0] = 0;
  int * blockMap  = (int *)malloc(Bfilled4 * sizeof(int));

  for (currentColNumberofB = 0 ; currentColNumberofB < Bfilled3 - 1; ++currentColNumberofB)
  {
    Bindex3_data[currentColNumberofB + 1] = Bindex3_data[currentColNumberofB];
    for (currentRowNumberofB = 0 ; currentRowNumberofB < B->filled1 - 1; ++currentRowNumberofB)
    {
      for (unsigned int blockNum = B->index1_data[currentRowNumberofB];
           blockNum < B->index1_data[currentRowNumberofB + 1]; ++blockNum)
      {
        assert(blockNum < B->filled2);
        colNumberofB = B->index2_data[blockNum];
        if (colNumberofB == currentColNumberofB)
        {
          blockNumB++;
          assert(blockNumB < B->nbblocks);
          Bindex3_data[currentColNumberofB + 1]++;
          Bindex4_data[blockNumB] = currentRowNumberofB;
          blockMap[blockNumB] = blockNum;
        }


      }

    }
  }
  int colNumberAA, rowNumberBB;

  int Cnbblocks = -1;

  for (currentRowNumberofA = 0 ; currentRowNumberofA < A->filled1 - 1; ++currentRowNumberofA)
  {

    for (currentColNumberofB = 0 ; currentColNumberofB < Bfilled3 - 1; ++currentColNumberofB)
    {

      int CblockPassed = 1;
      for (unsigned int blockNumAA = A->index1_data[currentRowNumberofA];
           blockNumAA < A->index1_data[currentRowNumberofA + 1]; ++blockNumAA)
      {
        assert(blockNumAA < A->filled2);
        colNumberAA = A->index2_data[blockNumAA];
        /*            printf("blockNumAA = %i, colNumberAA = %i\n",blockNumAA,colNumberAA  ); */

        for (unsigned int blockNumBB = Bindex3_data[currentColNumberofB];
             blockNumBB < Bindex3_data[currentColNumberofB + 1]; ++blockNumBB)
        {
          rowNumberBB = Bindex4_data[blockNumBB];
          /*             printf("blockNumBB = %i, rowNumberBB = %i\n",blockNumBB,rowNumberBB  ); */
          /*             printf("blocMap[blockNumBB] = %i, rowNumberBB = %i\n",blockMap[blockNumBB],rowNumberBB  ); */

          if (rowNumberBB == colNumberAA)
          {
            if (CblockPassed)
            {
              Cnbblocks++; /* Find the right C block number*/
              CblockPassed = 0;
            }
            assert(Cnbblocks < C->nbblocks);
            /*               printf("Compute C block number %i for %i %i ", Cnbblocks,currentRowNumberofA, currentColNumberofB ); */

            int Ablocksize0 = A->blocksize0[currentRowNumberofA];
            if (currentRowNumberofA != 0)
              Ablocksize0  -= A->blocksize0[currentRowNumberofA - 1];

            int Bblocksize1 = B->blocksize1[currentColNumberofB];
            if (currentColNumberofB != 0)
              Bblocksize1 -= B->blocksize1[currentColNumberofB - 1];
            /*               printf("of size %dX%d\n",Ablocksize0,Bblocksize1  ); */

            int Ablocksize1 = A->blocksize1[colNumberAA];
            if (colNumberAA != 0)
              Ablocksize1  -= A->blocksize1[colNumberAA - 1];

            int Bblocksize0 = B->blocksize0[rowNumberBB];
            if (rowNumberBB != 0)
              Bblocksize0 -= B->blocksize0[rowNumberBB - 1];


            /*               printf("Contribution of the product of blocks matrices A(%i,%i) and B(%i,%i) of  sizes %dX%d by %dX%d\n",currentRowNumberofA,colNumberAA,rowNumberBB, currentColNumberofB,   Ablocksize0,Ablocksize1,Bblocksize0,Bblocksize1   ); */

            assert(Ablocksize1 == Bblocksize0);
            /* for (i=0;i<Ablocksize0;i++) */
            /*             { */
            /*                 for (j=0;j<Ablocksize1;j++)  { */
            /*               printf("A->block[%i](%i,%i) = %f\n",blockNumAA,i,j, A->block[blockNumAA][i+j*Ablocksize0]); */
            /*                 } */
            /*             } */

            /*                for (i=0;i<Bblocksize0;i++) */
            /*             { */
            /*                 for (j=0;j<Bblocksize1;j++)  { */
            /*               printf("B->block[%i](%i,%i) = %f\n",blockNumBB,i,j, B->block[blockMap[blockNumBB]][i+j*Bblocksize0]); */
            /*                 } */
            /*             } */

            /*                printf("DGEMM call\n"); */
            DGEMM(LA_NOTRANS, LA_NOTRANS, Ablocksize0, Bblocksize1, Ablocksize1, alpha, A->block[blockNumAA], Ablocksize0, B->block[blockMap[blockNumBB]], Bblocksize0, beta, C->block[Cnbblocks], Ablocksize0);


            /*               for (i=0;i<Ablocksize0;i++) */
            /*             { */
            /*                 for (j=0;j<Bblocksize1;j++)  { */
            /*               printf("C->block[%i](%i,%i) = %f\n",Cnbblocks,i,j, C->block[Cnbblocks][i+j*Ablocksize0]); */
            /*                 } */
            /*             } */

          }; /* printf("\n"); */

        } /*  printf("\n"); */

      }

    }
  }

  assert((Cnbblocks + 1) == C->nbblocks);



  free(Bindex3_data);
  free(Bindex4_data);
  free(blockMap);

  /* printSBM(C); */



  /*   fprintf(stderr,"Numerics, SparseBlockStructuredMatrix, product matrix - matrix prodSBMSBM(alpha,A,B,beta,C) not yet implemented.\n"); */
  /*   exit(EXIT_FAILURE); */

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
  if (blmat->index1_data)
    free(blmat->index1_data);
  if (blmat->index2_data)
    free(blmat->index2_data);
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
  assert(m->index1_data);
  assert(m->index2_data);

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
  printf("index1_data of size %i= {", m->filled1);
  for (int i = 0 ; i < m->filled1 - 1; i++) printf("%d,  ", m->index1_data[i]);
  printf("%d}\n", m->index1_data[m->filled1 - 1]);
  printf("index2_data od size %i= {", m->filled2);
  for (int i = 0 ; i < m->filled2 - 1; i++) printf("%d,  ", m->index2_data[i]);
  printf("%d}\n", m->index1_data[m->filled2 - 1]);
  int sizemax = 10;
  int currentRowNumber ;
  int colNumber;
  int nbRows, nbColumns;
  for (currentRowNumber = 0 ; currentRowNumber < m->filled1 - 1; ++currentRowNumber)
  {
    for (unsigned int blockNum = m->index1_data[currentRowNumber];
         blockNum < m->index1_data[currentRowNumber + 1]; ++blockNum)
    {
      assert(blockNum < m->filled2);
      colNumber = m->index2_data[blockNum];
      /* Get dim. of the current block */
      nbRows = m->blocksize0[currentRowNumber];

      if (currentRowNumber != 0)
        nbRows -= m->blocksize0[currentRowNumber - 1];

      nbColumns = m->blocksize1[colNumber];
      if (colNumber != 0)
        nbColumns -= m->blocksize1[colNumber - 1];
      printf("block[%i] of size %dX%d\n", blockNum, nbRows, nbColumns);
      if ((nbRows <= sizemax) & (nbColumns <= sizemax))
      {
        for (int i = 0; i < nbRows; i++)
        {
          for (int j = 0; j < nbColumns; j++)
          {
            printf("block[%i](%i,%i) = %f\n", blockNum, i, j, m->block[blockNum][i + j * nbRows]);
          }
        }
      }
      else
      {
        printf("Block[%i] is too lao=rge to be displayed\n");
      }

    }
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
  /* Look for the first block of row number num */
  int pos;
  int firstBlockOfRow = M->index1_data[num];

  /* Look at the diagonal block */
  for (pos = firstBlockOfRow; M->index2_data[pos] != num; ++pos, assert(pos < M->filled2));

  return pos;
}
double getValueSBM(const SparseBlockStructuredMatrix* const M, int row, int col)
{
  /*      Search the row of blocks and the column of blocks */
  int i, j;
  int rowblocknumber = M->blocknumber0;
  int colblocknumber = M->blocknumber1;
  int rowpos = -1;
  int colpos = -1;
  int blockNumber = -1;
  int colNumber;
  assert(row < M->blocksize0[M->blocknumber0 - 1]);
  assert(col < M->blocksize1[M->blocknumber1 - 1]);


  for (i = 0; i < M->blocknumber0; i++)
  {
    if ((row < M->blocksize0[i]))
    {
      rowblocknumber = i;
      assert(rowblocknumber < M->blocknumber0);
      rowpos = row;
      if (i != 0)
        rowpos  -= M->blocksize0[i - 1];
      break;
    }

  }
  for (j = 0; j < M->blocknumber1; j++)
  {
    if ((col < M->blocksize1[j]))
    {
      colblocknumber = j;
      assert(colblocknumber < M->blocknumber1);
      colpos = col;
      if (j != 0)
        colpos -= M->blocksize1[j - 1];
      blockNumber = -1;
      for (unsigned int blockNum = M->index1_data[rowblocknumber];
           blockNum < M->index1_data[rowblocknumber + 1]; ++blockNum)
      {
        assert(blockNum < M->filled2);

        colNumber = M->index2_data[blockNum];
        if (colblocknumber == colNumber)
        {
          blockNumber = blockNum;
          break;
        }

      }
      break;
    }
  }

  if (blockNumber == -1) return 0.0;

  /* Number of rows/columns of the current block */
  int nbRows, nbColumns;
  nbRows = M->blocksize0[rowblocknumber];
  if (rowblocknumber != 0)
    nbRows -= M->blocksize0[rowblocknumber - 1];
  nbColumns = M->blocksize1[colblocknumber];
  if (colblocknumber != 0)
    nbColumns -= M->blocksize1[colblocknumber - 1];
  assert(rowpos < nbRows);
  assert(colpos < nbColumns);

  return M->block[blockNumber][rowpos + colpos * nbRows];



}

int copySBM(const SparseBlockStructuredMatrix* const A, SparseBlockStructuredMatrix*  B)
{
  assert(A);
  B->nbblocks = A->nbblocks;
  B->blocknumber0 = A->blocknumber0;
  B->blocknumber1 = A->blocknumber1;
  B->blocksize0 = (int*)malloc(B->blocknumber0 * sizeof(int));
  for (int i = 0; i < B->blocknumber0; i++) B->blocksize0[i] = A->blocksize0[i];
  B->blocksize1 = (int*)malloc(B->blocknumber1 * sizeof(int));
  for (int i = 0; i < B->blocknumber1; i++) B->blocksize1[i] = A->blocksize1[i];
  B->filled1 = A->filled1;
  B->filled2 = A->filled2;
  B->index1_data = (size_t*)malloc(B->filled1 * sizeof(size_t));
  for (int i = 0; i < B->filled1; i++) B->index1_data[i] = A->index1_data[i];
  B->index2_data = (size_t*)malloc(B->filled2 * sizeof(size_t));
  for (int i = 0; i < B->filled2; i++) B->index2_data[i] = A->index2_data[i];
  B->block = (double **)malloc(B->nbblocks * sizeof(double*));
  int currentRowNumber ;
  int colNumber;
  int nbRows, nbColumns;
  for (currentRowNumber = 0 ; currentRowNumber < B->filled1 - 1; ++currentRowNumber)
  {
    for (unsigned int blockNum = B->index1_data[currentRowNumber];
         blockNum < B->index1_data[currentRowNumber + 1]; ++blockNum)
    {
      assert(blockNum < B->filled2);
      colNumber = B->index2_data[blockNum];
      /* Get dim. of the current block */
      nbRows = B->blocksize0[currentRowNumber];
      if (currentRowNumber != 0)
        nbRows -= B->blocksize0[currentRowNumber - 1];
      nbColumns = B->blocksize1[colNumber];
      if (colNumber != 0)
        nbColumns -= B->blocksize1[colNumber - 1];
      B->block[blockNum] = (double*)malloc(nbRows * nbColumns * sizeof(double));
      for (int i = 0; i < nbRows * nbColumns; i++)
      {
        B->block[blockNum] [i] = A->block[blockNum] [i] ;
      }


    }
  }


}

int transposeSBM(const SparseBlockStructuredMatrix* const A, SparseBlockStructuredMatrix*  B)
{
  assert(A);
  B->nbblocks = A->nbblocks;
  B->blocknumber0 = A->blocknumber1;
  B->blocknumber1 = A->blocknumber0;
  B->blocksize0 = (int*)malloc(B->blocknumber0 * sizeof(int));
  for (int i = 0; i < B->blocknumber0; i++) B->blocksize0[i] = A->blocksize1[i];
  B->blocksize1 = (int*)malloc(B->blocknumber1 * sizeof(int));
  for (int i = 0; i < B->blocknumber1; i++) B->blocksize1[i] = A->blocksize0[i];

  B->filled1 = A->blocknumber1 + 1;;
  B->filled2 = A->filled2;
  B->index1_data = (size_t *)malloc(B->filled1 * sizeof(size_t));
  B->index2_data = (size_t *)malloc(B->filled2 * sizeof(size_t));


  int currentRowNumberofA;
  int currentColNumberofA ;
  int colNumberofA;

  int * blockMap  = (int *)malloc(B->filled2 * sizeof(int));



  int blockNumB = -1;
  B->index1_data[0] = 0;
  for (currentColNumberofA = 0 ; currentColNumberofA < A->filled2 - 1; ++currentColNumberofA)
  {
    B->index1_data[currentColNumberofA + 1] = B->index1_data[currentColNumberofA];
    for (currentRowNumberofA = 0 ; currentRowNumberofA < A->filled1 - 1; ++currentRowNumberofA)
    {
      for (unsigned int blockNum = A->index1_data[currentRowNumberofA];
           blockNum < A->index1_data[currentRowNumberofA + 1]; ++blockNum)
      {
        assert(blockNum < B->filled2);
        colNumberofA = A->index2_data[blockNum];
        if (colNumberofA == currentColNumberofA)
        {
          blockNumB++;
          assert(blockNumB < B->nbblocks);
          B->index1_data[currentColNumberofA + 1]++;
          B->index2_data[blockNumB] = currentRowNumberofA;
          blockMap[blockNum] = blockNumB;

        }


      }
    }
  }



  B->block = (double **)malloc(B->nbblocks * sizeof(double*));
  int currentRowNumber ;
  int colNumber;
  int nbRows, nbColumns;
  for (currentRowNumber = 0 ; currentRowNumber < B->filled1 - 1; ++currentRowNumber)
  {
    for (unsigned int blockNum = B->index1_data[currentRowNumber];
         blockNum < B->index1_data[currentRowNumber + 1]; ++blockNum)
    {
      assert(blockNum < B->filled2);
      colNumber = B->index2_data[blockNum];
      /* Get dim. of the current block */
      nbRows = B->blocksize0[currentRowNumber];
      if (currentRowNumber != 0)
        nbRows -= B->blocksize0[currentRowNumber - 1];
      nbColumns = B->blocksize1[colNumber];
      if (colNumber != 0)
        nbColumns -= B->blocksize1[colNumber - 1];
      B->block[blockNum] = (double*)malloc(nbRows * nbColumns * sizeof(double));
      for (int i = 0; i < nbRows; i++)
      {
        for (int j = 0; j < nbColumns; j++)
        {
          B->block[blockNum] [i + j * nbRows ] =  A->block[blockMap[blockNum]] [j + i * nbColumns] ;
        }
      }

    }
  }


}
int inverseDiagSBM(const SparseBlockStructuredMatrix*  M)
{

  for (int i = 0; i < M->filled1 - 1; i++)
  {
    int numberofblockperrow = M->index1_data[i + 1] - M->index1_data[i];
    if (numberofblockperrow != 1)
    {
      fprintf(stderr, "SparseBlockMatrix : inverseDiagSBM: Not a diagonal blocks matrix\n");
      exit(EXIT_FAILURE);
    }
  }
  for (int i = 0; i < M->filled2; i++)
  {
    if (M->index2_data[i] != i)
    {
      fprintf(stderr, "SparseBlockMatrix : inverseDiagSBM: Not a diagonal blocks matrix\n");
      exit(EXIT_FAILURE);
    }
  }

  int currentRowNumber ;
  int colNumber;
  int nbRows, nbColumns;
  for (currentRowNumber = 0 ; currentRowNumber < M->filled1 - 1; ++currentRowNumber)
  {
    for (unsigned int blockNum = M->index1_data[currentRowNumber];
         blockNum < M->index1_data[currentRowNumber + 1]; ++blockNum)
    {
      assert(blockNum < M->filled2);
      colNumber = M->index2_data[blockNum];
      /* Get dim. of the current block */
      nbRows = M->blocksize0[currentRowNumber];
      if (currentRowNumber != 0)
        nbRows -= M->blocksize0[currentRowNumber - 1];
      nbColumns = M->blocksize1[colNumber];
      if (colNumber != 0)
        nbColumns -= M->blocksize1[colNumber - 1];

      assert(nbRows == nbColumns);

      int infoDGETRF, infoDGETRI ;
      int* ipiv = (int *)malloc(nbRows * sizeof(*ipiv));

      DGETRF(nbRows, nbColumns, M->block[blockNum], nbRows, ipiv, infoDGETRF);
      DGETRI(nbRows, M->block[blockNum], nbRows, ipiv, infoDGETRF);
      free(ipiv);


    }
  }




}
