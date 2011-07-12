#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "SparseBlockMatrix.h"
#include "LA.h"
#include <math.h>
#include "misc.h"


unsigned int NUMERICS_SBM_FREE_BLOCK = 1 << 2;
unsigned int NUMERICS_SBM_FREE_SBM = 1 << 3;
//#define VERBOSE_DEBUG

void prodSBM(unsigned int sizeX, unsigned int sizeY, double alpha, const SparseBlockStructuredMatrix* const A, const double* const x, double beta, double* y)
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

  /* Column (block) position of the current block*/
  unsigned int colNumber;
  /* Number of rows/columns of the current block */
  unsigned int nbRows, nbColumns;
  /* Position of the sub-block of x multiplied by the sub-block of A */
  unsigned int posInX = 0;
  /* Position of the sub-block of y, result of the product */
  unsigned int posInY = 0;

  /* Loop over all non-null blocks
     Works whatever the ordering order of the block is, in A->block
  */
  DSCAL(sizeY, beta, y, 1);

  for (unsigned int currentRowNumber = 0 ; currentRowNumber < A->filled1 - 1; ++currentRowNumber)
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

      assert((nbRows <= sizeY));

      nbColumns = A->blocksize1[colNumber];
      if (colNumber != 0)
        nbColumns -= A->blocksize1[colNumber - 1];

      assert((nbColumns <= sizeX));

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

  /*     Check the compatibility of size of matrices */
  assert(A->blocksize0);
  assert(A->blocksize1);
  assert(B->blocksize0);
  assert(B->blocksize1);

  assert(A->blocksize1[A->blocknumber1 - 1] == B->blocksize0[B->blocknumber0 - 1]);

  /*     Check the compatibility of the number and the sizes of blocks */
  int compat = 1;

  if (A->blocknumber1 != B->blocknumber0) compat = 1; /* Compatibility of number of blocks */
  for (unsigned int i = 0; i < A->blocknumber1; i++)
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
    C->blocksize0  = (unsigned int *)malloc(C->blocknumber0 * sizeof(int));
    C->blocksize1  = (unsigned int *)malloc(C->blocknumber1 * sizeof(int));
    for (unsigned int i = 0; i < C->blocknumber0; i++) C->blocksize0[i] = A->blocksize0[i];
    for (unsigned int j = 0; j < C->blocknumber1; j++) C->blocksize1[j] = B->blocksize1[j];
  }

  unsigned int currentRowNumberofA, currentRowNumberofB;
  unsigned int colNumberofB;

  /* Search  the non null blocks of C */

  /*     Creation of filled3,filled4 and index3_data,index4_data for B (indexation by column) */

  unsigned int Bfilled3 = B->blocknumber1 + 1;
  unsigned int * Bindex3_data = (unsigned int *)malloc(Bfilled3 * sizeof(int));
  unsigned int Bfilled4 = B->nbblocks;
  unsigned int * Bindex4_data = (unsigned int *)malloc(Bfilled4 * sizeof(int));
  unsigned int blockNumB = -1;
  Bindex3_data[0] = 0;
  int * blockMap  = (int *)malloc(Bfilled4 * sizeof(int));
  for (unsigned int currentColNumberofB = 0 ; currentColNumberofB < Bfilled3 - 1; ++currentColNumberofB)
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
    for (unsigned int currentColNumberofB = 0 ; currentColNumberofB < Bfilled3 - 1; ++currentColNumberofB)
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
            /*           printf("C block number %i exists for %i %i ",C->nbblocks, currentRowNumberofA, currentColNumberofB ); */

            C->filled2++;
            unsigned int Cblosksize0 = A->blocksize0[currentRowNumberofA];
            if (currentRowNumberofA != 0)
              Cblosksize0  -= A->blocksize0[currentRowNumberofA - 1];
            unsigned int Cblocksize1 = B->blocksize1[currentColNumberofB];
            if (currentColNumberofB != 0)
              Cblocksize1 -= B->blocksize1[currentColNumberofB - 1];
            /*           printf("of size %dX%d\n",Cblosksize0,Cblocksize1  ); */
            Cblocktmp[C->nbblocks] = (double*)malloc(Cblosksize0 * Cblocksize1 * sizeof(double));
            for (unsigned int i = 0; i < Cblosksize0 * Cblocksize1; i++) Cblocktmp[C->nbblocks][i] = 0.0;

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

  for (unsigned int i = 0 ; i < C->nbblocks; i++)  C->block[i] = Cblocktmp[i];
  /*   for (i =0 ; i <C->nbblocks; i++)   */
  /*       { */
  /*    C->block[i] =  (double*)malloc(C->blocksize0[i]*C->blocksize1[i]*sizeof(double)); */
  /*    for (j = 0; j<C->blocksize0[i]*C->blocksize1[i]; j++) C->block[i][j]=0.0; */
  /*       } */


  for (unsigned int i = 0 ; i < C->nbblocks; i++)  C->index2_data[i] = Cindex2_datatmp[i];
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

  /*     Check the compatibility of size of matrices */
  assert(A->blocksize0);
  assert(A->blocksize1);
  assert(B->blocksize0);
  assert(B->blocksize1);

  assert(A->blocksize1[A->blocknumber1 - 1] == B->blocksize0[B->blocknumber0 - 1]);

  /*     Check the compatibility of the number and the sizes of blocks */
  int compat = 1;

  if (A->blocknumber1 != B->blocknumber0) compat = 1; /* Compatibility of number of blocks */
  for (unsigned int i = 0; i < A->blocknumber1; i++)
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
             for (unsigned int i = 0; i < C->blocknumber0; i++)
  {
    if (C->blocksize0[i] != A->blocksize0[i]) compat = 0;
    }
    compat == 1;
            });
    assert( {int compat = 1;
             for (unsigned int j = 0; j < C->blocknumber1; j++)
  {
    if (C->blocksize1[j] != B->blocksize1[j]) compat = 0;
    }
    compat == 1;
            });
  }

  unsigned int currentRowNumberofA, currentRowNumberofB;
  unsigned int currentColNumberofB;
  unsigned int colNumberofB ;

  /* Search  the non null blocks of C */

  /*     Creation of filled3,filled4 and index3_data,index4_data for B (indexation by column) */

  unsigned int Bfilled3 = B->blocknumber1 + 1;
  unsigned int * Bindex3_data = (unsigned int *)malloc(Bfilled3 * sizeof(int));
  unsigned int Bfilled4 = B->nbblocks;
  unsigned int * Bindex4_data = (unsigned int *)malloc(Bfilled4 * sizeof(int));
  unsigned int blockNumB = -1;
  Bindex3_data[0] = 0;
  unsigned int * blockMap  = (unsigned int *)malloc(Bfilled4 * sizeof(int));

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
        /*          printf("blockNumAA = %i, colNumberAA = %i\n",blockNumAA,colNumberAA  ); */

        for (unsigned int blockNumBB = Bindex3_data[currentColNumberofB];
             blockNumBB < Bindex3_data[currentColNumberofB + 1]; ++blockNumBB)
        {
          rowNumberBB = Bindex4_data[blockNumBB];
          /*          printf("blockNumBB = %i, rowNumberBB = %i\n",blockNumBB,rowNumberBB  ); */
          /*          printf("blocMap[blockNumBB] = %i, rowNumberBB = %i\n",blockMap[blockNumBB],rowNumberBB  ); */

          if (rowNumberBB == colNumberAA)
          {
            if (CblockPassed)
            {
              Cnbblocks++; /* Find the right C block number*/
              CblockPassed = 0;
            }
            assert(Cnbblocks < (int) C->nbblocks);
            /*           printf("Compute C block number %i for %i %i ", Cnbblocks,currentRowNumberofA, currentColNumberofB ); */

            int Ablocksize0 = A->blocksize0[currentRowNumberofA];
            if (currentRowNumberofA != 0)
              Ablocksize0  -= A->blocksize0[currentRowNumberofA - 1];

            int Bblocksize1 = B->blocksize1[currentColNumberofB];
            if (currentColNumberofB != 0)
              Bblocksize1 -= B->blocksize1[currentColNumberofB - 1];
            /*           printf("of size %dX%d\n",Ablocksize0,Bblocksize1  ); */

            int Ablocksize1 = A->blocksize1[colNumberAA];
            if (colNumberAA != 0)
              Ablocksize1  -= A->blocksize1[colNumberAA - 1];

            int Bblocksize0 = B->blocksize0[rowNumberBB];
            if (rowNumberBB != 0)
              Bblocksize0 -= B->blocksize0[rowNumberBB - 1];


            /*           printf("Contribution of the product of blocks matrices A(%i,%i) and B(%i,%i) of  sizes %dX%d by %dX%d\n",currentRowNumberofA,colNumberAA,rowNumberBB, currentColNumberofB,   Ablocksize0,Ablocksize1,Bblocksize0,Bblocksize1   ); */

            assert(Ablocksize1 == Bblocksize0);
            /* for (i=0;i<Ablocksize0;i++) */
            /*        { */
            /*            for (j=0;j<Ablocksize1;j++)  { */
            /*         printf("A->block[%i](%i,%i) = %f\n",blockNumAA,i,j, A->block[blockNumAA][i+j*Ablocksize0]); */
            /*            } */
            /*        } */

            /*            for (i=0;i<Bblocksize0;i++) */
            /*        { */
            /*            for (j=0;j<Bblocksize1;j++)  { */
            /*         printf("B->block[%i](%i,%i) = %f\n",blockNumBB,i,j, B->block[blockMap[blockNumBB]][i+j*Bblocksize0]); */
            /*            } */
            /*        } */

            /*            printf("DGEMM call\n"); */
            DGEMM(LA_NOTRANS, LA_NOTRANS, Ablocksize0, Bblocksize1, Ablocksize1, alpha, A->block[blockNumAA], Ablocksize0, B->block[blockMap[blockNumBB]], Bblocksize0, beta, C->block[Cnbblocks], Ablocksize0);


            /*           for (i=0;i<Ablocksize0;i++) */
            /*        { */
            /*            for (j=0;j<Bblocksize1;j++)  { */
            /*         printf("C->block[%i](%i,%i) = %f\n",Cnbblocks,i,j, C->block[Cnbblocks][i+j*Ablocksize0]); */
            /*            } */
            /*        } */

          }
          ; /* printf("\n"); */

        } /*  printf("\n"); */

      }

    }
  }

  assert((Cnbblocks + 1) == (int) C->nbblocks);



  free(Bindex3_data);
  free(Bindex4_data);
  free(blockMap);

  /* printSBM(C); */



  /*   fprintf(stderr,"Numerics, SparseBlockStructuredMatrix, product matrix - matrix prodSBMSBM(alpha,A,B,beta,C) not yet implemented.\n"); */
  /*   exit(EXIT_FAILURE); */

}
void subRowProdSBM(unsigned int sizeX, unsigned int sizeY, unsigned int currentRowNumber,
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

  /* Column (block) position of the current block*/
  unsigned int colNumber = 0;
  /* Number of rows/columns of the current block */
  unsigned int nbRows, nbColumns;
  /* Position of the sub-block of x multiplied by the sub-block of A */
  unsigned int posInX = 0;

  /* Check if currentRowNumber fits with A dimensions */
  assert(currentRowNumber <= A->blocknumber0);

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

void rowProdNoDiagSBM(unsigned int sizeX, unsigned int sizeY, unsigned int currentRowNumber, const SparseBlockStructuredMatrix* const A, const double* const x, double* y, int init)
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


  /* Column (block) position of the current block*/
  unsigned int colNumber = 0;

  /* Number of rows/columns of the current block */
  unsigned int nbRows, nbColumns;

  /* Position of the sub-block of x multiplied by the sub-block of
   * A */
  unsigned int posInX = 0;

  /* Look for the first element of the wanted row */

  /* Assertions */
  assert(A);
  assert(x);
  assert(y);
  assert(sizeX == A->blocksize1[A->blocknumber1 - 1]);
  assert(currentRowNumber <= A->blocknumber0);

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
  assert(blmat);

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

  for (unsigned int i = 0 ; i < blmat->nbblocks ; i++)
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
  printf("index1_data of size %li= {", (long int)m->filled1);
  for (unsigned int i = 0 ; i < m->filled1 - 1; i++) printf("%li,  ", (long int)m->index1_data[i]);
  printf("%li}\n", (long int)m->index1_data[m->filled1 - 1]);

  printf("index2_data of size %li= {", (long int)m->filled2);
  for (unsigned int i = 0 ; i < m->filled2 - 1; i++) printf("%li,  ", (long int)m->index2_data[i]);
  printf("%li}\n", (long int)m->index2_data[m->filled2 - 1]);

  printf("blocksize0 of size %li= {", (long int)m->blocknumber0);
  for (unsigned int i = 0 ; i < m->blocknumber0 - 1; i++) printf("%li,  ", (long int)m->blocksize0[i]);
  printf("%li}\n", (long int)m->blocksize0[m->blocknumber0 - 1]);

  printf("blocksize1 of size %li= {", (long int)m->blocknumber1);
  for (unsigned int i = 0 ; i < m->blocknumber1 - 1; i++) printf("%li,  ", (long int)m->blocksize1[i]);
  printf("%li}\n", (long int)m->blocksize1[m->blocknumber1 - 1]);



  unsigned int sizemax = 10;
  unsigned int currentRowNumber ;
  unsigned int colNumber;
  unsigned int nbRows, nbColumns;
  for (currentRowNumber = 0 ; currentRowNumber < m->filled1 - 1; ++currentRowNumber)
  {
    for (unsigned int blockNum = m->index1_data[currentRowNumber];
         blockNum < m->index1_data[currentRowNumber + 1]; ++blockNum)
    {
      assert(blockNum < m->filled2);
      colNumber = m->index2_data[blockNum];
      assert(colNumber < m->blocknumber1);
      /* Get dim. of the current block */
      nbRows = m->blocksize0[currentRowNumber];

      if (currentRowNumber != 0)
        nbRows -= m->blocksize0[currentRowNumber - 1];
      assert(nbRows);
      nbColumns = m->blocksize1[colNumber];
      if (colNumber != 0)
        nbColumns -= m->blocksize1[colNumber - 1];
      assert(nbColumns);

      printf("block[%i] of size %dX%d\n", blockNum, nbRows, nbColumns);
      if ((nbRows <= sizemax) & (nbColumns <= sizemax))
      {
        for (unsigned int i = 0; i < nbRows; i++)
        {
          for (unsigned int j = 0; j < nbColumns; j++)
          {
            printf("block[%i](%i,%i) = %12.8e\n", blockNum, i, j, m->block[blockNum][i + j * nbRows]);
          }
        }
      }
      else
      {
        printf("Block[%i] is too large to be displayed\n", blockNum);
      }

    }
  }

}
void printInFileSBM(const SparseBlockStructuredMatrix* const m, FILE * file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, SparseBlockStructuredMatrix printInFileSBM failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(file, "%i\n", m->nbblocks);
  if (m->nbblocks == 0)  return;

  fprintf(file, "%i\n", m->blocknumber0);
  fprintf(file, "%i\n", m->blocknumber1);

  assert(m->blocksize0);
  assert(m->blocksize1);
  assert(m->index1_data);
  assert(m->index2_data);

  for (unsigned int i = 0; i < m->blocknumber0; i++)
  {
    fprintf(file, "%i\t", m->blocksize0[i]);
  }
  fprintf(file, "\n");
  for (unsigned int i = 0; i < m->blocknumber1; i++)
  {
    fprintf(file, "%i\t", m->blocksize1[i]);
  }
  fprintf(file, "\n");
  fprintf(file, "%li\n", (long int)m->filled1);
  fprintf(file, "%li\n", (long int)m->filled2);
  for (unsigned int i = 0; i < m->filled1; i++)
  {
    fprintf(file, "%li\t", (long int)m->index1_data[i]);
  }
  fprintf(file, "\n");
  for (unsigned int i = 0; i < m->filled2; i++)
  {
    fprintf(file, "%li\t", (long int)m->index2_data[i]);
  }
  fprintf(file, "\n");

  unsigned int currentRowNumber ;
  unsigned int colNumber;
  unsigned int nbRows, nbColumns;
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
      //fprintf(file,"block[%i] of size %dX%d\n", blockNum, nbRows,nbColumns);
      fprintf(file, "%i\n", blockNum);
      for (unsigned int i = 0; i < nbRows * nbColumns; i++)
      {
        fprintf(file, "%32le\n", m->block[blockNum][i]);
      }

    }
  }

}
void printInFileSBMForScilab(const SparseBlockStructuredMatrix* const m, FILE * file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, SparseBlockStructuredMatrix printInFileSBM failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(file, "nbblock = %i;\n", m->nbblocks);
  if (m->nbblocks == 0)  return;

  fprintf(file, "blocknumber0 = %i;\n", m->blocknumber0);
  fprintf(file, "blocknumber1 = %i; \n", m->blocknumber1);

  assert(m->blocksize0);
  assert(m->blocksize1);
  assert(m->index1_data);
  assert(m->index2_data);

  fprintf(file, "blocksize0 = [ \t");
  for (unsigned int i = 0; i < m->blocknumber0; i++)
  {
    fprintf(file, "%i\t", m->blocksize0[i]);
  }
  fprintf(file, "];\n");
  fprintf(file, "blocksize1 = [ \t");
  for (unsigned int i = 0; i < m->blocknumber1; i++)
  {
    fprintf(file, "%i\t", m->blocksize1[i]);
  }
  fprintf(file, "];\n");
  fprintf(file, "filled1 = %li;\n", (long int)m->filled1);
  fprintf(file, "filled2 = %li;\n", (long int)m->filled2);

  fprintf(file, "index1_data = [ \t");
  for (unsigned int i = 0; i < m->filled1; i++)
  {
    fprintf(file, "%li\t", (long int)m->index1_data[i]);
  }
  fprintf(file, "];\n");
  fprintf(file, "index2_data = [ \t");
  for (unsigned int i = 0; i < m->filled2; i++)
  {
    fprintf(file, "%li\t", (long int)m->index2_data[i]);
  }
  fprintf(file, "];\n");

  unsigned int currentRowNumber ;
  unsigned int colNumber;
  unsigned int nbRows, nbColumns;
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
      //fprintf(file,"block[%i] of size %dX%d\n", blockNum, nbRows,nbColumns);
      fprintf(file, "block%i = [ \n", blockNum);

      for (unsigned int i = 0; i < nbRows; i++)
      {
        fprintf(file, "[");
        for (unsigned int j = 0; j < nbColumns; j++)
        {
          fprintf(file, "%32.24e\t ", m->block[blockNum][i + j * nbRows]);
        }
        fprintf(file, "];\n");
      }
      fprintf(file, "];\n");
      /*       for (int i=0;i<nbRows*nbColumns;i++) */
      /*       { */
      /*         fprintf(file,"%32le\n",m->block[blockNum][i]); */
      /*       } */
    }
  }

  fprintf(file, "// Dense version\n");
  int size0 = m->blocksize0[m->blocknumber0 - 1];
  int size1 = m->blocksize1[m->blocknumber1 - 1];
  double * denseMat = (double*)malloc(size0 * size1 * sizeof(double));
  SBMtoDense(m, denseMat);

  fprintf(file, "data= [");
  for (int i = 0; i < size0; i++)
  {
    fprintf(file, "[");
    for (int j = 0; j < size1; j++)
    {
      fprintf(file, "%32.24e,\t ", denseMat[i + j * size0]);
    }
    fprintf(file, "];\n");
  }
  fprintf(file, "]");
  free(denseMat);



}


void printInFileNameSBM(const SparseBlockStructuredMatrix* const m, const char *filename)
{

}
void newFromFileSBM(SparseBlockStructuredMatrix* const m, FILE *file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, SparseBlockStructuredMatrix readInFileSBM failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  CHECK_IO(fscanf(file, "%d", &(m->nbblocks)));

  if (m->nbblocks == 0)  return;

  CHECK_IO(fscanf(file, "%d", &(m->blocknumber0)));
  CHECK_IO(fscanf(file, "%d", &(m->blocknumber1)));

  m->blocksize0 = (unsigned int *)malloc(m->blocknumber0 * sizeof(int));
  m->blocksize1 = (unsigned int *)malloc(m->blocknumber1 * sizeof(int));


  for (unsigned int i = 0; i < m->blocknumber0; i++)
  {
    CHECK_IO(fscanf(file, "%d", &(m->blocksize0[i])));
  }
  for (unsigned int i = 0; i < m->blocknumber1; i++)
  {
    CHECK_IO(fscanf(file, "%d", &(m->blocksize1[i])));
  }

  unsigned int filled1 = 0, filled2 = 0;
  CHECK_IO(fscanf(file, "%d", &(filled1)));
  m->filled1 = filled1;
  CHECK_IO(fscanf(file, "%d", &(filled2)));
  m->filled2 = filled2;
  m->index1_data = (size_t *)malloc(m->filled1 * sizeof(size_t));
  m->index2_data = (size_t *)malloc(m->filled2 * sizeof(size_t));

  int index1_dataCurrent = 0;
  for (unsigned int i = 0; i < m->filled1; i++)
  {
    CHECK_IO(fscanf(file, "%d", &(index1_dataCurrent)));
    m->index1_data[i] = index1_dataCurrent;
  }
  int index2_dataCurrent = 0;
  for (unsigned int i = 0; i < m->filled2; i++)
  {
    CHECK_IO(fscanf(file, "%d", &(index2_dataCurrent)));
    m->index2_data[i] = index2_dataCurrent;
  }
  m->block = (double**)malloc(m->nbblocks * sizeof(double*));
  unsigned int currentRowNumber ;
  unsigned int colNumber;
  unsigned int nbRows, nbColumns;
  unsigned int blockk;
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
      //fprintf(file,"block[%i] of size %dX%d\n", blockNum, nbRows,nbColumns);

      CHECK_IO(fscanf(file, "%d", &(blockk)));
      if (blockk != blockNum)
      {
        printf("Numerics, SparseBlockStructuredMatrix readInFileSBM failed, problem in block numbering. \n");
      }
      m->block[blockNum] = (double*)malloc(nbRows * nbColumns * sizeof(double));
      for (unsigned int i = 0; i < nbRows * nbColumns; i++)
      {
        CHECK_IO(fscanf(file, "%32le\n", &(m->block[blockNum][i])));
      }

    }
  }
}

void readInFileSBM(SparseBlockStructuredMatrix* const m, FILE *file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, SparseBlockStructuredMatrix readInFileSBM failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  CHECK_IO(fscanf(file, "%d", &(m->nbblocks)));

  if (m->nbblocks == 0)  return;

  CHECK_IO(fscanf(file, "%d", &(m->blocknumber0)));
  CHECK_IO(fscanf(file, "%d", &(m->blocknumber1)));

  for (unsigned int i = 0; i < m->blocknumber0; i++)
  {
    CHECK_IO(fscanf(file, "%d", &(m->blocksize0[i])));
  }
  for (unsigned int i = 0; i < m->blocknumber1; i++)
  {
    CHECK_IO(fscanf(file, "%d", &(m->blocksize1[i])));
  }

  unsigned int filled1 = 0, filled2 = 0;
  CHECK_IO(fscanf(file, "%d", &(filled1)));
  m->filled1 = filled1;
  CHECK_IO(fscanf(file, "%d", &(filled2)));
  m->filled2 = filled2;

  unsigned int index1_dataCurrent = 0;
  for (unsigned int i = 0; i < m->filled1; i++)
  {
    CHECK_IO(fscanf(file, "%d", &(index1_dataCurrent)));
    m->index1_data[i] = index1_dataCurrent;
  }
  unsigned int index2_dataCurrent = 0;
  for (unsigned int i = 0; i < m->filled2; i++)
  {
    CHECK_IO(fscanf(file, "%d", &(index2_dataCurrent)));
    m->index2_data[i] = index2_dataCurrent;
  }
  unsigned int currentRowNumber ;
  unsigned int colNumber;
  unsigned int nbRows, nbColumns;
  unsigned int blockk;
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
      //fprintf(file,"block[%i] of size %dX%d\n", blockNum, nbRows,nbColumns);

      CHECK_IO(fscanf(file, "%d", &(blockk)));
      if (blockk != blockNum)
      {
        printf("Numerics, SparseBlockStructuredMatrix readInFileSBM failed, problem in block numbering. \n");
      }

      for (unsigned int i = 0; i < nbRows * nbColumns; i++)
      {
        CHECK_IO(fscanf(file, "%32le\n", &(m->block[blockNum][i])));
      }

    }
  }
}
void readInFileNameSBM(SparseBlockStructuredMatrix* const m, const char *filename)
{

}
void freeSpBlMatPred(SparseBlockStructuredMatrixPred *blmatpred)
{

  for (int i = 0 ; i < blmatpred->nbbldiag ; i++)
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

unsigned int getDiagonalBlockPos(const SparseBlockStructuredMatrix* const M, unsigned int num)
{
  /* Look for the first block of row number num */
  unsigned int pos;
  unsigned int firstBlockOfRow = M->index1_data[num];

  /* Look at the diagonal block */
  for (pos = firstBlockOfRow; M->index2_data[pos] != num; ++pos, assert(pos < M->filled2));

  return pos;
}
double getValueSBM(const SparseBlockStructuredMatrix* const M, unsigned int row, unsigned int col)
{
  /*      Search the row of blocks and the column of blocks */
  unsigned int rowblocknumber = M->blocknumber0;
  unsigned int colblocknumber = M->blocknumber1;
  int rowpos = -1;
  int colpos = -1;
  int blockNumber = -1;
  unsigned int colNumber;
  assert(row < M->blocksize0[M->blocknumber0 - 1]);
  assert(col < M->blocksize1[M->blocknumber1 - 1]);


  for (unsigned int i = 0; i < M->blocknumber0; i++)
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
  for (unsigned int j = 0; j < M->blocknumber1; j++)
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

int copySBM(const SparseBlockStructuredMatrix* const A, SparseBlockStructuredMatrix*  B, unsigned int copyBlock)
{
  assert(A);
  B->nbblocks = A->nbblocks;
  B->blocknumber0 = A->blocknumber0;
  B->blocknumber1 = A->blocknumber1;
  B->blocksize0 = (unsigned int*)malloc(B->blocknumber0 * sizeof(int));
  for (unsigned int i = 0; i < B->blocknumber0; i++) B->blocksize0[i] = A->blocksize0[i];
  B->blocksize1 = (unsigned int*)malloc(B->blocknumber1 * sizeof(int));
  for (unsigned int i = 0; i < B->blocknumber1; i++) B->blocksize1[i] = A->blocksize1[i];
  B->filled1 = A->filled1;
  B->filled2 = A->filled2;
  B->index1_data = (size_t*)malloc(B->filled1 * sizeof(size_t));
  for (unsigned int i = 0; i < B->filled1; i++) B->index1_data[i] = A->index1_data[i];
  B->index2_data = (size_t*)malloc(B->filled2 * sizeof(size_t));
  for (unsigned int i = 0; i < B->filled2; i++) B->index2_data[i] = A->index2_data[i];
  B->block = (double **)malloc(B->nbblocks * sizeof(double*));
  if (copyBlock)
  {
    unsigned int currentRowNumber ;
    unsigned int colNumber;
    unsigned int nbRows, nbColumns;
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
        for (unsigned int i = 0; i < nbRows * nbColumns; i++)
        {
          B->block[blockNum] [i] = A->block[blockNum] [i] ;
        }
      }
    }
  }
  else
  {
    for (unsigned int n = 0; n < B->nbblocks; n++)
      B->block[n] = A->block[n];
  }

  return 0;
}

int transposeSBM(const SparseBlockStructuredMatrix* const A, SparseBlockStructuredMatrix*  B)
{
  assert(A);
  assert(B);
  B->nbblocks = A->nbblocks;
  B->blocknumber0 = A->blocknumber1;
  B->blocknumber1 = A->blocknumber0;
  B->blocksize0 = (unsigned int*)malloc(B->blocknumber0 * sizeof(int));
  for (unsigned int i = 0; i < B->blocknumber0; i++) B->blocksize0[i] = A->blocksize1[i];
  B->blocksize1 = (unsigned int*)malloc(B->blocknumber1 * sizeof(int));
  for (unsigned int i = 0; i < B->blocknumber1; i++) B->blocksize1[i] = A->blocksize0[i];

  B->filled1 = A->blocknumber1 + 1;;
  B->filled2 = A->filled2;
  B->index1_data = (size_t *)malloc(B->filled1 * sizeof(size_t));
  B->index2_data = (size_t *)malloc(B->filled2 * sizeof(size_t));


  unsigned int currentRowNumberofA;
  unsigned int currentColNumberofA ;
  unsigned int colNumberofA;

  unsigned int * blockMap  = (unsigned int *)malloc(B->filled2 * sizeof(unsigned int));



  int blockNumB = -1;
  B->index1_data[0] = 0;
  for (currentColNumberofA = 0 ; currentColNumberofA < A->blocknumber1; ++currentColNumberofA)
  {
    assert(currentColNumberofA + 1 < B->filled1);
    B->index1_data[currentColNumberofA + 1] = B->index1_data[currentColNumberofA];
    for (currentRowNumberofA = 0 ; currentRowNumberofA < A->blocknumber0; ++currentRowNumberofA)
    {
      for (unsigned int blockNum = A->index1_data[currentRowNumberofA];
           blockNum < A->index1_data[currentRowNumberofA + 1]; ++blockNum)
      {
        assert(blockNum < B->filled2);
        colNumberofA = A->index2_data[blockNum];
        if (colNumberofA == currentColNumberofA)
        {
          blockNumB++;
          assert(blockNumB < (int) B->nbblocks && blockNumB >= 0);
          B->index1_data[currentColNumberofA + 1]++;
          B->index2_data[blockNumB] = currentRowNumberofA;
          blockMap[blockNumB] = blockNum;

        }
      }
    }
  }

#ifdef VERBOSE_DEBUG
  printf("----------------- blockMap ---------------\n");
  for (int i = 0; i < B->filled2; i++)
  {
    printf("blockMap[%i] = %i\n", i, blockMap[i]);
  }
  printf("----------------- blockMap ---------------\n");
#endif



  B->block = (double **)malloc(B->nbblocks * sizeof(double*));
  unsigned int currentRowNumber ;
  unsigned int colNumber;
  unsigned int nbRows, nbColumns;
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
      int lengthblock = nbRows * nbColumns;
      B->block[blockNum] = (double*)malloc(lengthblock * sizeof(double));




      for (unsigned int i = 0; i < nbRows; i++)
      {
        for (unsigned int j = 0; j < nbColumns; j++)
        {
          assert(i + j * nbRows < nbRows * nbColumns);
          assert(j + i * nbColumns < nbRows * nbColumns);

          B->block[blockNum] [i + j * nbRows ] = A->block[blockMap[blockNum]] [j + i * nbColumns] ;
        }
      }

    }
  }
  free(blockMap);


  return 0;
}
int inverseDiagSBM(const SparseBlockStructuredMatrix*  M)
{

  for (unsigned int i = 0; i < M->filled1 - 1; i++)
  {
    unsigned int numberofblockperrow = M->index1_data[i + 1] - M->index1_data[i];
    if (numberofblockperrow != 1)
    {
      fprintf(stderr, "SparseBlockMatrix : inverseDiagSBM: Not a diagonal blocks matrix\n");
      exit(EXIT_FAILURE);
    }
  }
  for (unsigned int i = 0; i < M->filled2; i++)
  {
    if (M->index2_data[i] != i)
    {
      fprintf(stderr, "SparseBlockMatrix : inverseDiagSBM: Not a diagonal blocks matrix\n");
      exit(EXIT_FAILURE);
    }
  }

  unsigned int currentRowNumber ;
  unsigned int colNumber;
  unsigned int nbRows, nbColumns;
  int infoDGETRF, infoDGETRI ;
  int info = 0;
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


      int* ipiv = (int *)malloc(nbRows * sizeof(*ipiv));

      DGETRF(nbRows, nbColumns, M->block[blockNum], nbRows, ipiv, infoDGETRF);
      assert(!infoDGETRF);

      DGETRI(nbRows, M->block[blockNum], nbRows, ipiv, infoDGETRI);
      assert(!infoDGETRI);
      free(ipiv);


    }
  }
  if ((!infoDGETRF) || (!infoDGETRF)) info = 0;

  return info;


}

void SBMtoDense(const SparseBlockStructuredMatrix* const A, double *denseMat)
{
  assert(A);
  assert(A->blocksize0);
  assert(A->blocksize1);
  int n = A->blocksize0[A->blocknumber0 - 1];
  int m = A->blocksize1[A->blocknumber1 - 1];

  /*    denseMat = (double*)malloc(n*m*sizeof(double));  */
  for (int i = 0; i < n ; i++)
  {
    for (int j = 0; j < m; j++)
    {
      denseMat[i + j * n] = getValueSBM(A, i, j);
    }
  }
}

int  SBMtoSparseInitMemory(const SparseBlockStructuredMatrix* const A, SparseMatrix *sparseMat)
{
  assert(A);
  assert(A->blocksize0);
  assert(A->blocksize1);
  int n = A->blocksize0[A->blocknumber0 - 1];
  int m = A->blocksize1[A->blocknumber1 - 1];

  sparseMat->m = n;
  sparseMat->n = m;

  sparseMat->nz = -2;
  sparseMat->nzmax = 0;
  sparseMat->p = (int*)malloc((sparseMat->m + 1) * sizeof(int));

  /* Row (block) position of the current block */
  unsigned int currentRowNumber ;
  /* Column (block) position of the current block*/
  unsigned int colNumber;
  /* Number of rows/columns of the current block */
  int nbRows, nbColumns;

  for (currentRowNumber = 0 ; currentRowNumber < A->filled1 - 1; ++currentRowNumber)
  {
    for (unsigned int blockNum = A->index1_data[currentRowNumber];
         blockNum < A->index1_data[currentRowNumber + 1]; ++blockNum)
    {
      assert(blockNum < A->filled2);
      colNumber = A->index2_data[blockNum];

      /* Get dim. of the current block */
      nbRows = A->blocksize0[currentRowNumber];

      if (currentRowNumber != 0)
        nbRows -= A->blocksize0[currentRowNumber - 1];
      assert((nbRows >= 0));
      nbColumns = A->blocksize1[colNumber];
      if (colNumber != 0)
        nbColumns -= A->blocksize1[colNumber - 1];
      assert((nbColumns >= 0));

      sparseMat->nzmax += nbColumns * nbRows;
    }
  }
  sparseMat->i = (int*)malloc((sparseMat->nzmax) * sizeof(int));
  sparseMat->x = (double*)malloc((sparseMat->nzmax) * sizeof(double));

  return 0;
}



int  SBMtoSparse(const SparseBlockStructuredMatrix* const A, SparseMatrix *sparseMat)
{
  assert(A);
  assert(A->blocksize0);
  assert(A->blocksize1);

  assert(sparseMat);
  assert(sparseMat->p);
  assert(sparseMat->i);
  assert(sparseMat->x);

  /* Row (block) position of the current block */
  unsigned int currentRowNumber ;
  /* Column (block) position of the current block*/
  unsigned int colNumber;
  /* Number of rows/columns of the current block */
  int nbRows, nbColumns;

  int nnz = 0;
  int isparserowend, isparserowstart, isparserow;
  int isparsecolumnend, isparsecolumnstart, isparsecolumn;
  sparseMat->p[0] = 0; /* We assume that the first row is non empty */
  for (currentRowNumber = 0 ; currentRowNumber < A->filled1 - 1; ++currentRowNumber)
  {

    /* Get row dim. of the current block line*/
    nbRows = A->blocksize0[currentRowNumber];

    isparserowend = A->blocksize0[currentRowNumber];

    if (currentRowNumber != 0)
    {
      nbRows -= A->blocksize0[currentRowNumber - 1];
      isparserowstart = A->blocksize0[currentRowNumber - 1];
    }
    else
    {
      isparserowstart = 0;

    }
    assert((nbRows >= 0));
#ifdef VERBOSE_DEBUG
    printf("isparserowstart = %i\t", isparserowstart);
    printf("isparserowend = %i\n", isparserowend);

#endif

    for (isparserow = isparserowstart; isparserow < isparserowend ; isparserow++)
    {


      sparseMat->p[isparserow + 1] =  sparseMat->p[isparserow];

      for (unsigned int blockNum = A->index1_data[currentRowNumber];
           blockNum < A->index1_data[currentRowNumber + 1]; ++blockNum)
      {
        colNumber = A->index2_data[blockNum];

        nbColumns = A->blocksize1[colNumber];
        isparsecolumnend = A->blocksize1[colNumber];

        if (colNumber != 0)
        {
          nbColumns -= A->blocksize1[colNumber - 1];
          isparsecolumnstart =  A->blocksize1[colNumber - 1];
        }
        else
          isparsecolumnstart = 0;

        assert((nbColumns >= 0));

        sparseMat->p[isparserow + 1] += nbColumns  ;

#ifdef VERBOSE_DEBUG
        printf("isparsecolumnstart = %i\t", isparsecolumnstart);
        printf("isparsecolumnend = %i\n", isparsecolumnend);

#endif

        for (isparsecolumn = isparsecolumnstart; isparsecolumn < isparsecolumnend; isparsecolumn++)
        {

#ifdef VERBOSE_DEBUG
          /*      printf("isparsecolumn = %i\t", isparsecolumn);*/
          printf("nnz = %i \t", nnz);
          printf("isparserow = %i \t", isparserow);
          printf("isparsecolumn = %i \t", isparsecolumn);
#endif
          sparseMat->i[nnz] = isparsecolumn;

          int rowintheblock = isparserow;
          if (currentRowNumber != 0)
          {
            rowintheblock -= A->blocksize0[currentRowNumber - 1];
          }



          int colintheblock = isparsecolumn;
          if (colNumber != 0)
          {
            colintheblock -= A->blocksize1[colNumber - 1];
          }
          assert(rowintheblock < nbRows);
          assert(colintheblock < nbColumns);
#ifdef VERBOSE_DEBUG
          printf(" rowintheblock= %i \t", rowintheblock);
          printf("colintheblock = %i \n", colintheblock);
#endif

          sparseMat->x[nnz] = A->block[blockNum][colintheblock * nbRows + rowintheblock];
          nnz++;
        }
#ifdef VERBOSE_DEBUG
        //      printf("isparsecolumn = %i\t", isparsecolumn);
        printf("\n");
#endif
      }


    }
  }
  return 0;
}


void SBMfree(SparseBlockStructuredMatrix* A, unsigned int level)
{

  if (level & NUMERICS_SBM_FREE_BLOCK)
  {
    for (unsigned int i = 0; i < A->nbblocks; i++)
      free(A->block[i]);
  }
  free(A->block);
  free(A->blocksize0);
  free(A->blocksize1);
  free(A->index1_data);
  free(A->index2_data);
  if (level & NUMERICS_SBM_FREE_SBM)
    free(A);
}

//#define SBM_DEBUG_SBMRowToDense
void SBMRowToDense(const SparseBlockStructuredMatrix* const A, int row, double *denseMat, int rowPos, int rowNb)
{
  assert(A);
  int BlockRowNb = 0;
  int ColNb = A->blocksize1[A->blocknumber1 - 1];
  if (row)
    BlockRowNb = A->blocksize0[row] - A->blocksize0[row - 1];
  else
    BlockRowNb = A->blocksize0[row];
#ifdef SBM_DEBUG_SBMRowToDense
  printf("SBMRowToDense : copi block row %i, containing %i row and %i col.\n", row, BlockRowNb, ColNb);
#endif

  //zero memory
  for (int numRow = rowPos; numRow < rowPos + BlockRowNb; numRow++)
    for (int numCol = 0; numCol < ColNb; numCol++)
      denseMat[numRow + numCol * rowNb] = 0.0;

  //index1_data[rowNumber]<= blockNumber <index1_data[rowNumber+1]
  for (unsigned int numBlock = A->index1_data[row]; numBlock < A->index1_data[row + 1]; numBlock++)
  {
    double * beginBlock = A->block[numBlock];
    int colNumber = A->index2_data[numBlock];
    int indexColBegin = 0;
    if (colNumber)
      indexColBegin = A->blocksize1[colNumber - 1];
    int BlockColNb = A->blocksize1[colNumber] - indexColBegin;
    for (int i = 0; i < BlockRowNb; i++)
    {
      for (int j = 0; j < BlockColNb; j++)
      {
        denseMat[rowPos + i + (indexColBegin + j)*rowNb] = beginBlock[i + j * BlockRowNb];
      }
    }
  }
#ifdef SBM_DEBUG_SBMRowToDense
  printf("SBMRowToDense : res in file SBMRowToDense.txt.");
  FILE * titi  = fopen("SBMRowToDense.txt", "w");
  printInFileSBMForScilab(A, titi);
  fprintf(titi, "\n//Dense matrix of row block %i:\n", row);
  fprintf(titi, "denseRow = [ \t");
  for (int i = 0; i < BlockRowNb; i++)
  {
    fprintf(titi, "[");
    for (int j = 0; j < ColNb; j++)
    {
      fprintf(titi, "%32.24e\t  ", denseMat[rowPos + i + j * rowNb]);
    }
    fprintf(titi, "];\n");
  }
  fprintf(titi, "];\n");
  fclose(titi);
#endif
}
//#define SBM_DEBUG_SBM_ROW_PERM
void RowPermutationSBM(unsigned int *rowIndex, SparseBlockStructuredMatrix* A, SparseBlockStructuredMatrix*  C)
{
#ifdef SBM_DEBUG_SBM_ROW_PERM
  FILE * titi  = fopen("RowPermutationSBM_input.txt", "w");
  printInFileSBMForScilab(A, titi);
  fclose(titi);
#endif
  int nbRow = A->blocknumber0;
  int nbCol = A->blocknumber1;
  C->nbblocks = A->nbblocks;
  C->block = (double**)malloc(A->nbblocks * sizeof(double*));
  C->blocknumber0 = A->blocknumber0;
  C->blocknumber1 = A->blocknumber1;
  C->blocksize0 = (unsigned int*)malloc(nbRow * sizeof(unsigned int));
  C->blocksize1 = (unsigned int*)malloc(nbCol * sizeof(unsigned int));
  C->filled1 = A->filled1;
  C->filled2 = A->filled2;
  C->index1_data = (size_t*)malloc(C->filled1 * sizeof(size_t));
  C->index2_data = (size_t*)malloc(C->filled2 * sizeof(size_t));
  /*Row permutation ==> same col size*/
  for (int i = 0; i < nbCol; i++)
  {
    C->blocksize1[i] = A->blocksize1[i];
  }
  int curNbBlockC = 0;
  C->index1_data[0] = 0;
  int nbRowInBlock;
  for (int rowC = 0; rowC < nbCol; rowC++)
  {
    int rowA = rowIndex[rowC];
    /*C->blocksize0[rowC+1]=C->blocksize0[rowC]+ number of row of the current block*/
    nbRowInBlock = A->blocksize0[rowA];
    if (rowA)
      nbRowInBlock -= A->blocksize0[rowA - 1];
#ifdef SBM_DEBUG_SBM_ROW_PERM
    printf("RowPermutationSBM rowA=%i, rowC=%i\n", rowA, rowC);
#endif
    if (rowC)
      C->blocksize0[rowC] = C->blocksize0[rowC - 1] + nbRowInBlock;
    else
      C->blocksize0[rowC] = nbRowInBlock;
    size_t NbBlockCurRow = A->index1_data[rowA + 1] - A->index1_data[rowA];

    C->index1_data[rowC + 1] = C->index1_data[rowC] + NbBlockCurRow;
    for (size_t numBlockInRowA = A->index1_data[rowA]; numBlockInRowA < A->index1_data[rowA + 1]; numBlockInRowA++)
    {
      C->index2_data[curNbBlockC] = A->index2_data[numBlockInRowA];
      C->block[curNbBlockC] = A->block[numBlockInRowA];
      curNbBlockC++;
    }
  }
#ifdef SBM_DEBUG_SBM_ROW_PERM
  titi  = fopen("RowPermutationSBM_output.txt", "w");
  printInFileSBMForScilab(C, titi);
  fclose(titi);
#endif
}
//#define SBM_DEBUG_SBM_COL_PERM
void ColPermutationSBM(unsigned int *colIndex, SparseBlockStructuredMatrix* A, SparseBlockStructuredMatrix*  C)
{
#ifdef SBM_DEBUG_SBM_COL_PERM
  FILE * titi  = fopen("ColPermutationSBM_input.txt", "w");
  printInFileSBMForScilab(A, titi);
  fclose(titi);
#endif
  copySBM(A, C, 0);
  for (unsigned int n = 0; n < C->nbblocks; n++)
  {
    C->index2_data[n] = colIndex[C->index2_data[n]];
  }
  //int curColnb=0;
  int nbBlockCol = A->blocknumber1;
  for (int numCol = 0; numCol < nbBlockCol; numCol++)
  {
#ifdef SBM_DEBUG_SBM_COL_PERM
    printf("ColPermutationSBM colA=%i, colC=%i\n", numCol, colIndex[numCol]);
#endif
    int colInA = colIndex[numCol];
    int nbCol = A->blocksize1[colInA];
    if (colInA)
      nbCol -= A->blocksize1[colInA - 1];
    if (numCol)
      C->blocksize1[numCol] = C->blocksize1[numCol - 1] + nbCol;
    else
      C->blocksize1[numCol] = nbCol;
  }
#ifdef SBM_DEBUG_SBM_COL_PERM
  titi  = fopen("ColPermutationSBM_output.txt", "w");
  printInFileSBMForScilab(C, titi);
  fclose(titi);
#endif
}
