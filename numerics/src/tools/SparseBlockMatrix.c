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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "CSparseMatrix_internal.h"
#include "SparseBlockMatrix.h"
#include "SiconosLapack.h"
#include <math.h>
#include <string.h>
#include "numerics_verbose.h"
#include "op3x3.h"
#include "SiconosCompat.h"
#include "NumericsArrays.h"
/* #define DEBUG_NOCOLOR 1 */
/* #define DEBUG_STDOUT 1 */
/* #define DEBUG_MESSAGES 1 */

#include "debug.h"
#include "CSparseMatrix.h"

#ifdef DEBUG_MESSAGES
#include "NumericsMatrix.h"
#endif

//#define VERBOSE_DEBUG

/* Clang 3.6 seems to be at odd with the C99 and C11 std on the conversion from
 * double to _Bool, cf paragraph "6.3.1.2 Boolean type" in the latest C99 or C11 draft  */
#ifdef __clang__
#if (__clang_major__ == 3 && __clang_minor__ == 6)
#pragma clang diagnostic ignored "-Wfloat-conversion"
#endif
#endif
void SBM_null(SparseBlockStructuredMatrix* sbm)
{
  sbm->nbblocks = 0;
  sbm->block = NULL;
  sbm->blocknumber0 = 0;
  sbm->blocknumber1 = 0;
  sbm->blocksize0 = NULL;
  sbm->blocksize1 = NULL;
  sbm->filled1 = 0;
  sbm->filled2 = 0;
  sbm->index1_data = NULL;
  sbm->index2_data = NULL;
  sbm->diagonal_blocks = NULL;
}

SparseBlockStructuredMatrix* SBM_new(void)
{
  SparseBlockStructuredMatrix* sbm = (SparseBlockStructuredMatrix*)
    malloc(sizeof(SparseBlockStructuredMatrix));

  SBM_null(sbm);

  return sbm;
}

void SBM_free(SparseBlockStructuredMatrix *sbm)
{
  /* Free memory for SparseBlockStructuredMatrix */
  /* Warning: nothing is done to check if memory has really been
   allocated for each component or if it was only pointer links.  Note
   that when used from the Kernel, memory is not directly allocated
   for such structures and this function must not be called. See in
   kernel/src/simulationsTools/SparseBlockMatrix.cpp for details on
   the way the structure is filled in.
  */
  assert(sbm);

  if (sbm->blocksize0)
  {

    free(sbm->blocksize0);
    if (sbm->blocksize0 == sbm->blocksize1)
    {
      sbm->blocksize1 = NULL ;
    }
    sbm->blocksize0 = NULL;
  }
  if (sbm->blocksize1)
  {
    free(sbm->blocksize1);
    sbm->blocksize1 = NULL;
  }

  for (unsigned int i = 0 ; i < sbm->nbblocks ; i++)
  {
    if (sbm->block[i])
    {
      free(sbm->block[i]);
      sbm->block[i] = NULL;
    }
  }

  if (sbm->block)
  {
    free(sbm->block);
    sbm->block = NULL;
  }


  if (sbm->index1_data)
  {
    free(sbm->index1_data);
    sbm->index1_data = NULL;
  }

  if (sbm->index2_data)
  {
    free(sbm->index2_data);
    sbm->index2_data = NULL;
  }

  if (sbm->diagonal_blocks)
  {
    free(sbm->diagonal_blocks);
    sbm->diagonal_blocks = NULL;
  }
  sbm->filled1 = 0;
  sbm->filled2 = 0;
  sbm->blocknumber0 = 0;
  sbm->blocknumber1 = 0;
  sbm->nbblocks = 0;
}


void SBM_print(const SparseBlockStructuredMatrix* const m)
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
  if (m->index1_data)
  {
    for (unsigned int i = 0 ; i < m->filled1 - 1; i++) printf("%li,  ", (long int)m->index1_data[i]);
    printf("%li}\n", (long int)m->index1_data[m->filled1 - 1]);
  }
  else
    printf("m->index1_data --> NULL}\n");

  printf("index2_data of size %li= {", (long int)m->filled2);
  if (m->index2_data)
  {
    for (unsigned int i = 0 ; i < m->filled2 - 1; i++) printf("%li,  ", (long int)m->index2_data[i]);
    printf("%li}\n", (long int)m->index2_data[m->filled2 - 1]);
  }
  else
    printf("m->index2_data --> NULL}\n");

  printf("blocksize0 of size %li= {", (long int)m->blocknumber0);
  if (m->blocksize0)
  {
    for (unsigned int i = 0 ; i < m->blocknumber0 - 1; i++) printf("%li,  ", (long int)m->blocksize0[i]);
    printf("%li}\n", (long int)m->blocksize0[m->blocknumber0 - 1]);
  }
  else
    printf("m->blocksize0 --> NULL}\n");

  printf("blocksize1 of size %li= {", (long int)m->blocknumber1);
   if (m->blocksize1)
   {
     for (unsigned int i = 0 ; i < m->blocknumber1 - 1; i++) printf("%li,  ", (long int)m->blocksize1[i]);
     printf("%li}\n", (long int)m->blocksize1[m->blocknumber1 - 1]);
   }
   else
     printf("m->blocksize1 --> NULL}\n");


  unsigned int sizemax = 10;
  unsigned int currentRowNumber ;
  size_t colNumber;
  unsigned int nbRows, nbColumns;
  if (m->block)
  {
    for (currentRowNumber = 0 ; currentRowNumber < m->filled1 - 1; ++currentRowNumber)
    {
      /* Get dim. of the current block */
      nbRows = m->blocksize0[currentRowNumber];
      if (currentRowNumber != 0)
        nbRows -= m->blocksize0[currentRowNumber - 1];
      assert(nbRows);

      for (size_t blockNum = m->index1_data[currentRowNumber];
           blockNum < m->index1_data[currentRowNumber + 1]; ++blockNum)
      {
        assert(blockNum < m->filled2);
        colNumber = m->index2_data[blockNum];
        assert(colNumber < m->blocknumber1);
        

        
        nbColumns = m->blocksize1[colNumber];
        if (colNumber != 0)
          nbColumns -= m->blocksize1[colNumber - 1];
        assert(nbColumns);

        printf("block[" SN_SIZE_T_F "] of size %dX%d\n", blockNum, nbRows, nbColumns);
        if(m->block[blockNum])
        {
          if ((nbRows <= sizemax) & (nbColumns <= sizemax))
          {
            for (unsigned int i = 0; i < nbRows; i++)
            {
              for (unsigned int j = 0; j < nbColumns; j++)
              {
                printf("block[" SN_SIZE_T_F "](%i,%i) = %12.8e\n", blockNum, i, j, m->block[blockNum][i + j * nbRows]);
              }
            }
          }
          else
          {
            printf("Block[" SN_SIZE_T_F "] is too large to be displayed\n", blockNum);
          }
        }
        else
          printf("Block[" SN_SIZE_T_F "] --> NULL \n", blockNum);

      }
    }
  }
  else
    printf("m->block --> NULL");

  printf("m-diagonal_blocks  ");
  if (m->diagonal_blocks)
  {
    printf("[ ");
    for (currentRowNumber = 0 ; currentRowNumber < m->filled1 - 1; ++currentRowNumber)
    {
      printf(" %i",m->diagonal_blocks[currentRowNumber] );
    }
    printf("]\n");
  }
  else
    printf("--> NULL");
  

}


/* a basic iterator scheme for different kind of sparse
 * matrices (csc, csr, triplet) */
typedef struct sparse_matrix_iterator
{
  CS_INT counter1;
  CS_INT counter2;
  CS_INT first;
  CS_INT second;
  double third;
  const CSparseMatrix* mat;
} sparse_matrix_iterator;

static sparse_matrix_iterator sparseMatrixBegin(const CSparseMatrix* const sparseMat);
static int sparseMatrixNext(sparse_matrix_iterator* it);


void SBM_gemv(unsigned int sizeX, unsigned int sizeY, double alpha, const SparseBlockStructuredMatrix* const restrict A, const double* restrict x, double beta, double* restrict y)
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
  size_t colNumber;
  /* Number of rows/columns of the current block */
  unsigned int nbRows, nbColumns;
  /* Position of the sub-block of x multiplied by the sub-block of A */
  unsigned int posInX = 0;
  /* Position of the sub-block of y, result of the product */
  unsigned int posInY = 0;

  /* Loop over all non-null blocks
     Works whatever the ordering order of the block is, in A->block
  */
  cblas_dscal(sizeY, beta, y, 1);

  for (unsigned int currentRowNumber = 0 ; currentRowNumber < A->filled1 - 1; ++currentRowNumber)
  {
    /* Get dim. of the current block */
    nbRows = A->blocksize0[currentRowNumber];
    if (currentRowNumber != 0)
      nbRows -= A->blocksize0[currentRowNumber - 1];
    assert((nbRows <= sizeY));
    for (size_t blockNum = A->index1_data[currentRowNumber];
         blockNum < A->index1_data[currentRowNumber + 1]; ++blockNum)
    {
      assert(blockNum < A->filled2);

      colNumber = A->index2_data[blockNum];

      assert(colNumber < sizeX);

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
      if (nbRows == 3 && nbColumns == 3)
      {
        mvp_alpha3x3(alpha, A->block[blockNum], &x[posInX], &y[posInY]);
      }
      else
      {
        cblas_dgemv(CblasColMajor, CblasNoTrans, nbRows, nbColumns, alpha, A->block[blockNum],
                  nbRows, &x[posInX], 1, 1.0, &y[posInY], 1);
      }
    }
  }
}
void SBM_gemv_3x3(unsigned int sizeX, unsigned int sizeY, const SparseBlockStructuredMatrix* const restrict A,  double* const restrict x, double* restrict y)
{
  /* Product SparseMat - vector, y = vector product y += alpha*A*x  for block of size 3x3 */

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
  size_t colNumber;
  /* Number of rows/columns of the current block */
  unsigned int nbRows, nbColumns;
  /* Position of the sub-block of x multiplied by the sub-block of A */
  unsigned int posInX = 0;
  /* Position of the sub-block of y, result of the product */
  unsigned int posInY = 0;

  /* Loop over all non-null blocks
     Works whatever the ordering order of the block is, in A->block
  */

  for (unsigned int currentRowNumber = 0 ; currentRowNumber < A->filled1 - 1; ++currentRowNumber)
  {
    /* Get dim. of the current block */
    nbRows = A->blocksize0[currentRowNumber];
    if (currentRowNumber != 0)
      nbRows -= A->blocksize0[currentRowNumber - 1];
      
    assert((nbRows <= sizeY));
    for (size_t blockNum = A->index1_data[currentRowNumber];
         blockNum < A->index1_data[currentRowNumber + 1]; ++blockNum)
    {
      assert(blockNum < A->filled2);

      colNumber = A->index2_data[blockNum];

      assert(colNumber < sizeX);



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

      /* cblas_dgemv(CblasColMajor, CblasNoTrans, nbRows, nbColumns, alpha, A->block[blockNum], */
      /*             nbRows, &x[posInX], 1, 1.0, &y[posInY], 1); */
      assert((nbColumns == 3));
      assert((nbRows == 3));
      mvp3x3(A->block[blockNum], &x[posInX], &y[posInY]);

    }
  }
}
void SBM_extract_component_3x3(const SparseBlockStructuredMatrix* const restrict A, SparseBlockStructuredMatrix*  B,
                               unsigned int *row_components, unsigned int row_components_size,
                               unsigned int *col_components, unsigned int col_components_size)
{


  assert(A);
  assert(B);
  assert(A->blocksize0);
  assert(A->blocksize1);
  assert(A->index1_data);
  assert(A->index2_data);

  /* allocation of data */
  B->nbblocks= A->nbblocks;
  B->blocknumber0= A->blocknumber0;
  B->blocknumber1= A->blocknumber1;
  B->filled1= A->filled1;
  B->filled2= A->filled2;

  B->index1_data = (size_t *) malloc(B->filled1*sizeof(size_t));
  memcpy(B->index1_data, A->index1_data , B->filled1*sizeof(size_t));
  B->index2_data = (size_t *) malloc(B->filled2*sizeof(size_t));
  memcpy(B->index2_data, A->index2_data , B->filled2*sizeof(size_t));

  B->blocksize0  = (unsigned int *)malloc(B->blocknumber0*sizeof(unsigned int));
  B->blocksize1  = (unsigned int *)malloc(B->blocknumber1*sizeof(unsigned int));


  int sum =0;

  for (unsigned int row =0; row < B->blocknumber0; row++)
  {
    sum +=row_components_size;
    B->blocksize0[row] = sum ;
  }
  sum =0;
  for (unsigned int col =0; col < B->blocknumber1; col++)
  {
    sum +=col_components_size;
    B->blocksize1[col] = sum ;
  }


  B->block= (double **)malloc(B->nbblocks*sizeof(double*));



  /* Number of rows of the current block */
  unsigned int nbRows;

  /* Loop over all non-null blocks
     Works whatever the ordering order of the block is, in A->block
  */

  for (unsigned int currentRowNumber = 0 ; currentRowNumber < A->filled1 - 1; ++currentRowNumber)
  {
    /* Get dim. of the current block */
    nbRows = A->blocksize0[currentRowNumber];
    if (currentRowNumber != 0)
      nbRows -= A->blocksize0[currentRowNumber - 1];
    
    for (size_t blockNum = A->index1_data[currentRowNumber];
         blockNum < A->index1_data[currentRowNumber + 1]; ++blockNum)
    {
      assert(blockNum < A->filled2);
      
      B->block[blockNum] = (double*) malloc(row_components_size*col_components_size*sizeof(double));

      for (unsigned int i = 0; i < row_components_size; i++ )
      {
        for (unsigned int j = 0; j < col_components_size; j++ )
        {
        B->block[blockNum][i +  row_components_size*j] = A->block[blockNum][row_components[i] + col_components[j] * nbRows];
        }
      }
    }
  }

  DEBUG_EXPR(SBM_print(B));
}

static int SBM_check_compatibility_for_add(const SparseBlockStructuredMatrix* const A, const SparseBlockStructuredMatrix* const B)
{
  DEBUG_BEGIN(" SBM_check_compatibility_for_add(...)\n");
  assert(A);
  assert(B);

  /*     Check the compatibility of size of matrices */
  assert(A->blocksize0);
  assert(A->blocksize1);
  assert(B->blocksize0);
  assert(B->blocksize1);


  assert(A->blocksize1[A->blocknumber1 - 1] == B->blocksize1[B->blocknumber1 - 1]);
  assert(A->blocksize0[A->blocknumber0 - 1] == B->blocksize0[B->blocknumber0 - 1]);

  /*     Check the compatibility of the number and the sizes of blocks */
  int compat = 1;

  if (A->blocknumber0 != B->blocknumber0) compat = 0; /* Compatibility of number of blocks */
  if (A->blocknumber1 != B->blocknumber1) compat = 0; /* Compatibility of number of blocks */


  for (unsigned int i = 0; i < A->blocknumber0; i++)
  {
    if (A->blocksize0[i] != B->blocksize0[i]) compat = 0; /* Compatibility of sizes of blocks */
  }
  for (unsigned int i = 0; i < A->blocknumber1; i++)
  {
    if (A->blocksize1[i] != B->blocksize1[i]) compat = 0; /* Compatibility of sizes of blocks */
  }
  DEBUG_PRINTF("compat = %i \n", compat);
  DEBUG_END("SBM_check_compatibility_for_add(...)\n");
  return compat;

}



static SparseBlockStructuredMatrix * SBM_calloc_for_add(SparseBlockStructuredMatrix * A, SparseBlockStructuredMatrix * B)
{
  DEBUG_BEGIN("SBM_calloc_for_add(...)\n");
  size_t max_number_of_blocks_in_row = 0;
  size_t number_of_blocks_in_row = 0;
  size_t nb_blocks_in_rowA  =0,  nb_blocks_in_rowB =0;

  for (unsigned int currentRowNumber = 0 ; currentRowNumber < A->filled1 - 1; ++currentRowNumber)
  {

    nb_blocks_in_rowA  = A->index1_data[currentRowNumber + 1] - A->index1_data[currentRowNumber];
    nb_blocks_in_rowB  = B->index1_data[currentRowNumber + 1] - B->index1_data[currentRowNumber];
    max_number_of_blocks_in_row = max(max_number_of_blocks_in_row, nb_blocks_in_rowA +nb_blocks_in_rowB);
  }

  size_t  * common_block = (size_t *)malloc(max_number_of_blocks_in_row * sizeof(size_t));

  SparseBlockStructuredMatrix * C = SBM_new();

  C->blocknumber0 = A->blocknumber0;
  C->blocknumber1 = A->blocknumber1;
  C->blocksize0 = (unsigned int*)malloc(C->blocknumber0*sizeof(unsigned int));
  C->blocksize1 = (unsigned int*)malloc(C->blocknumber1*sizeof(unsigned int));

  for (unsigned int i = 0; i < A->blocknumber0; i++)
  {
    C->blocksize0[i] = A->blocksize0[i];
  }
  for (unsigned int i = 0; i < A->blocknumber1; i++)
  {
    C->blocksize1[i] = A->blocksize1[i];
  }

  C->index1_data = (size_t *)malloc(A->filled1 * sizeof(size_t));

  C->index2_data = (size_t *)malloc(max_number_of_blocks_in_row * A->blocknumber0* sizeof(size_t)); /* We oversize a bit */
  C->block = (double **)malloc(max_number_of_blocks_in_row * A->blocknumber0* sizeof(double *)); /* We oversize a bit */

  C->index1_data[0]=0;
  C->filled1=0;
  C->filled2=0;

  for (unsigned int currentRowNumber = 0 ; currentRowNumber < A->filled1 - 1; ++currentRowNumber)
  {

    nb_blocks_in_rowA  = A->index1_data[currentRowNumber + 1] - A->index1_data[currentRowNumber];
    nb_blocks_in_rowB  = B->index1_data[currentRowNumber + 1] - B->index1_data[currentRowNumber];
    DEBUG_PRINTF("nb_blocks_in_rowA = %zu\n",nb_blocks_in_rowA );
    DEBUG_PRINTF("nb_blocks_in_rowB = %zu\n",nb_blocks_in_rowB );

    number_of_blocks_in_row = NA_merge_and_sort_sorted_arrays(
      &(A->index2_data[A->index1_data[currentRowNumber]]),
      &(B->index2_data[B->index1_data[currentRowNumber]]),nb_blocks_in_rowA,nb_blocks_in_rowB,
      common_block);

    DEBUG_EXPR(NA_display(common_block,number_of_blocks_in_row ););

    DEBUG_PRINTF("number_of_blocks_in_row = %zu\n",number_of_blocks_in_row);
    C->index1_data[currentRowNumber+1] = C->index1_data[currentRowNumber] + number_of_blocks_in_row ;
    C->filled1++;

    int blocksize0 = C->blocksize0[currentRowNumber];
    if (currentRowNumber != 0)
      blocksize0  -= C->blocksize0[currentRowNumber - 1];

    int currentColNumber=0;
    int idx=0;
    for (size_t blockNum = C->index1_data[currentRowNumber];
         blockNum < C->index1_data[currentRowNumber + 1]; ++blockNum)
    {

       C->filled2++;
       currentColNumber = common_block[idx];
       C->index2_data[blockNum]=currentColNumber;
       idx++;

       int blocksize1 = C->blocksize1[currentColNumber];
       if (currentColNumber != 0)
          blocksize1 -= C->blocksize1[currentColNumber - 1];

       C->block[blockNum] = (double *)calloc(blocksize0*blocksize1,sizeof(double));
    }

  }
  C->filled1++;
  C->nbblocks = C->filled2;
  DEBUG_END("SBM_calloc_for_add(...)\n");
  return C;
}

void SBM_add_without_allocation(SparseBlockStructuredMatrix * A, SparseBlockStructuredMatrix * B,
                                double alpha, double beta,
                                SparseBlockStructuredMatrix * C, double gamma)
{
  DEBUG_BEGIN("SBM_add_without_allocation(...)\n");
  assert(SBM_check_compatibility_for_add(A,B) && "Non compatible matrices or blocks sizes.\n");
  int nm =0;
  size_t colNumber;
  for (unsigned int currentRowNumber = 0 ; currentRowNumber < C->filled1 - 1; ++currentRowNumber)
  {
    DEBUG_PRINT("\n");
    DEBUG_PRINTF("Computation of the blocks [%i, %i ] in row %i\n", C->index1_data[currentRowNumber], C->index1_data[currentRowNumber+1]-1,  currentRowNumber);

    unsigned int blocksize0 = C->blocksize0[currentRowNumber];
    if (currentRowNumber != 0)
      blocksize0  -= C->blocksize0[currentRowNumber - 1];

    for (size_t blockNum = C->index1_data[currentRowNumber];
         blockNum < C->index1_data[currentRowNumber + 1]; ++blockNum)
    {
      colNumber=C->index2_data[blockNum];
      unsigned int blocksize1 = C->blocksize1[colNumber];
      if (colNumber != 0)
        blocksize1 -= C->blocksize1[colNumber - 1];

      nm = blocksize0 * blocksize1;
      DEBUG_PRINTF("gamma*C for block %zu of size %ix%i\n",blockNum, blocksize0, blocksize1 );
      cblas_dscal(nm, gamma,  C->block[blockNum], 1 );
    }

    size_t jC = C->index1_data[currentRowNumber], jA = A->index1_data[currentRowNumber];

    while(jC < C->index1_data[currentRowNumber+1] && jA < A->index1_data[currentRowNumber+1] )
    {
      /* DEBUG_PRINTF("block number jC = %zu\t, blocknumber  jA = %zu\n", jC, jA); */
      if (C->index2_data[jC] < A->index2_data[jA])
      {
        /* DEBUG_PRINT("no block of A to add Block increment jA\n"); */
        jC++;
      }
      else
      {
        DEBUG_PRINTF("column number C->index2_data[jC] = %zu\t, column number  jA = %zu\n", C->index2_data[jC], A->index2_data[jA]);
        DEBUG_PRINTF("add a block number %zu of A to blocknumber %zu of C\n",jA,jC );

        unsigned int blocksize1 = C->blocksize1[C->index2_data[jC]];
        if (C->index2_data[jC] != 0)
          blocksize1 -= C->blocksize1[C->index2_data[jC] - 1];

        nm = blocksize0*blocksize1;
        cblas_daxpy(nm, alpha, A->block[jA], 1, C->block[jC], 1 );
        jA++;

      }
    }
    jC = C->index1_data[currentRowNumber];
    size_t jB = B->index1_data[currentRowNumber];

    while(jC < C->index1_data[currentRowNumber+1] && jB < B->index1_data[currentRowNumber+1] )
    {

      if (C->index2_data[jC] < B->index2_data[jB])
      {
        /* DEBUG_PRINT("no block of B to add Block increment jB\n"); */
        jC++;
      }
      else
      {
        DEBUG_PRINTF("add a block number %zu of B to blocknumber %zu of C\n",jB,jC );
        DEBUG_PRINTF("column number C->index2_data[jC] = %zu\t, column number  jB = %zu\n", C->index2_data[jC], B->index2_data[jB]);
        unsigned int blocksize1 = C->blocksize1[C->index2_data[jC]];
        if (C->index2_data[jC] != 0)
          blocksize1 -= C->blocksize1[C->index2_data[jC] - 1];

        nm = blocksize0*blocksize1;
        cblas_daxpy(nm, beta, B->block[jB], 1, C->block[jC], 1 );

        jB++;

      }
    }



  }
  DEBUG_END("SBM_add_without_allocation(...)\n");

}

SparseBlockStructuredMatrix * SBM_add(SparseBlockStructuredMatrix * A, SparseBlockStructuredMatrix * B, double alpha, double beta)
{

  /*     Check the compatibility of the number and the sizes of blocks */
  int compat =  SBM_check_compatibility_for_add(A,B);
  if (!compat)
  {
    numerics_error("SBM_add", "Non compatible matrices or blocks sizes.\n");
    return NULL;
  }
  SparseBlockStructuredMatrix * C = SBM_calloc_for_add(A,B);
  SBM_add_without_allocation(A, B, alpha, beta, C, 0.0);

  return C;

}





static struct SBM_index_by_column* SBM_index_by_column_new()
{
  struct SBM_index_by_column* p = (struct SBM_index_by_column*) malloc(sizeof(struct SBM_index_by_column));
  p->filled3=0;
  p->filled4=0;
  p->index3_data=NULL;
  p->index4_data=NULL;
  p->blockMap = NULL;
  return p;
}
static  struct SBM_index_by_column*  SBM_index_by_column_free(struct SBM_index_by_column* p)
{
  if (p->index3_data)
  {
    free(p->index3_data);
    p->index3_data=NULL;
  }
  if (p->index4_data)
  {
    free(p->index4_data);
    p->index4_data = NULL;
  }
  if (p->blockMap)
  {
    free(p->blockMap);
    p->blockMap=NULL;
  }
  free(p);
  p = NULL;
  return p;
}

static void SBM_index_by_column_compute(const SparseBlockStructuredMatrix* const M, struct SBM_index_by_column * SBM_index_by_column_M )
{

  SBM_index_by_column_M->filled3 = M->blocknumber1 + 1;
  SBM_index_by_column_M->index3_data = (size_t *)malloc(SBM_index_by_column_M->filled3 * sizeof(size_t));;

  SBM_index_by_column_M->filled4 =  M->nbblocks;
  SBM_index_by_column_M->index4_data = (size_t *)malloc(SBM_index_by_column_M->filled4 * sizeof(size_t));;

  size_t blockNumM = -1;
  SBM_index_by_column_M->index3_data[0] = 0;
  SBM_index_by_column_M->blockMap  = (size_t *)malloc(SBM_index_by_column_M->filled4 * sizeof(size_t));
  unsigned int currentRowNumberofM;
  unsigned int currentColNumberofM;
  size_t colNumberofM ;
  for (currentColNumberofM = 0 ; currentColNumberofM < SBM_index_by_column_M->filled3 - 1; ++currentColNumberofM)
  {
    SBM_index_by_column_M->index3_data[currentColNumberofM + 1] = SBM_index_by_column_M->index3_data[currentColNumberofM];

    for (currentRowNumberofM = 0 ; currentRowNumberofM < M->filled1 - 1; ++currentRowNumberofM)
    {
      for (size_t blockNum = M->index1_data[currentRowNumberofM];
           blockNum < M->index1_data[currentRowNumberofM + 1]; ++blockNum)
      {
        assert(blockNum < M->filled2);
        colNumberofM = M->index2_data[blockNum];
        if (colNumberofM == currentColNumberofM)
        {
          blockNumM++;
          assert(blockNumM < M->nbblocks);
          SBM_index_by_column_M->index3_data[currentColNumberofM + 1]++;
          SBM_index_by_column_M->index4_data[blockNumM] = currentRowNumberofM;
          SBM_index_by_column_M->blockMap[blockNumM] = blockNum;
        }
      }
    }
  }
}

static int SBM_check_compatibility_for_multiply(const SparseBlockStructuredMatrix* const A, const SparseBlockStructuredMatrix* const B)
{
  DEBUG_BEGIN(" SBM_check_compatibility_for_gemm(...)\n");
  assert(A);
  assert(B);

  /*     Check the compatibility of size of matrices */
  assert(A->blocksize0);
  assert(A->blocksize1);
  assert(B->blocksize0);
  assert(B->blocksize1);

  assert(A->blocksize1[A->blocknumber1 - 1] == B->blocksize0[B->blocknumber0 - 1]);

  /*     Check the compatibility of the number and the sizes of blocks */
  int compat = 1;

  if (A->blocknumber1 != B->blocknumber0) compat = 0; /* Compatibility of number of blocks */


  for (unsigned int i = 0; i < A->blocknumber1; i++)
  {
    if (A->blocksize1[i] != B->blocksize0[i]) compat = 0; /* Compatibility of sizes of blocks */
  }

  DEBUG_PRINTF("compat = %i \n", compat);
  DEBUG_END("SBM_check_compatibility_for_gemm(...)\n");
  return compat;

}


static SparseBlockStructuredMatrix*  SBM_calloc_multiply(const SparseBlockStructuredMatrix* const A,
                                                         const SparseBlockStructuredMatrix* const B,
                                                         struct SBM_index_by_column * SBM_index_by_column_B)
{

  SparseBlockStructuredMatrix* C = SBM_new();

  C->blocknumber0 = A->blocknumber0;
  C->blocknumber1 = B->blocknumber1;
  C->blocksize0  = (unsigned int *)malloc(C->blocknumber0 * sizeof(unsigned int));
  C->blocksize1  = (unsigned int *)malloc(C->blocknumber1 * sizeof(unsigned int));
  for (unsigned int i = 0; i < C->blocknumber0; i++) C->blocksize0[i] = A->blocksize0[i];
  for (unsigned int j = 0; j < C->blocknumber1; j++) C->blocksize1[j] = B->blocksize1[j];

  unsigned int currentRowNumberofA;

  size_t colNumberAA;
  size_t rowNumberBB;
  C->nbblocks = -1;
  C->filled2 = -1;
  /*     \warning The implementation is chosen to optimize cpu effort rather than memory. Otherwise a two loops are needed */
  int nbblocksmax = A->blocknumber0 * B->blocknumber1;

  double **Cblocktmp = (double**)malloc(nbblocksmax * sizeof(double*));
  size_t *Cindex2_datatmp = (size_t*)malloc(nbblocksmax * sizeof(size_t));
  C->filled1 = C->blocknumber0 + 1;
  C->index1_data = (size_t*)malloc(C->filled1 * sizeof(size_t));
  C->index1_data[0] = 0;
  for (currentRowNumberofA = 0 ; currentRowNumberofA < A->filled1 - 1; ++currentRowNumberofA)
  {
    C->index1_data[currentRowNumberofA + 1] = C->index1_data[currentRowNumberofA];
    unsigned int Cblosksize0 = A->blocksize0[currentRowNumberofA];
    if (currentRowNumberofA != 0)
      Cblosksize0  -= A->blocksize0[currentRowNumberofA - 1];
    for (unsigned int currentColNumberofB = 0 ; currentColNumberofB < SBM_index_by_column_B->filled3 - 1; ++currentColNumberofB)
    {
      int BlockCexists = 0;

      for (size_t blockNumAA = A->index1_data[currentRowNumberofA];
           blockNumAA < A->index1_data[currentRowNumberofA + 1]; ++blockNumAA)
      {
        assert(blockNumAA < A->filled2);
        colNumberAA = A->index2_data[blockNumAA];

        for (size_t blockNumBB = SBM_index_by_column_B->index3_data[currentColNumberofB];
             blockNumBB < SBM_index_by_column_B->index3_data[currentColNumberofB + 1]; ++blockNumBB)
        {
          rowNumberBB = SBM_index_by_column_B->index4_data[blockNumBB];
          if (rowNumberBB == colNumberAA)
          {
            BlockCexists = 1;
            C->nbblocks++;
            /*           printf("C block number %i exists for %i %i ",C->nbblocks, currentRowNumberofA, currentColNumberofB ); */

            C->filled2++;

            unsigned int Cblocksize1 = B->blocksize1[currentColNumberofB];
            if (currentColNumberofB != 0)
              Cblocksize1 -= B->blocksize1[currentColNumberofB - 1];
            /*           printf("of size %dX%d\n",Cblosksize0,Cblocksize1  ); */
            Cblocktmp[C->nbblocks] = (double*)calloc(Cblosksize0 * Cblocksize1, sizeof(double));
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
  C->index2_data = (size_t*)malloc(C->nbblocks * sizeof(size_t));

  for (unsigned int i = 0 ; i < C->nbblocks; i++)  C->block[i] = Cblocktmp[i];
  /*   for (i =0 ; i <C->nbblocks; i++)   */
  /*       { */
  /*    C->block[i] =  (double*)malloc(C->blocksize0[i]*C->blocksize1[i]*sizeof(double)); */
  /*    for (j = 0; j<C->blocksize0[i]*C->blocksize1[i]; j++) C->block[i][j]=0.0; */
  /*       } */


  for (unsigned int i = 0 ; i < C->nbblocks; i++)  C->index2_data[i] = Cindex2_datatmp[i];
  free(Cblocktmp);
  free(Cindex2_datatmp);

  return C;
}
SparseBlockStructuredMatrix*  SBM_zero_matrix_for_multiply(const SparseBlockStructuredMatrix* const A, const SparseBlockStructuredMatrix* const B)
{
  DEBUG_BEGIN("SBM_zero_matrix_for_multiply\n");
  assert(A);
  assert(B);

  /*     Check the compatibility of the number and the sizes of blocks */
  int compat =  SBM_check_compatibility_for_multiply(A,B);
  if (!compat)
  {
    numerics_error("SBM_zero_matrix_for_multiply", "Non compatible matrices or blocks sizes.\n");
    return NULL;
  }

  /* indexation of B by column */
  struct SBM_index_by_column* B_index_by_column = SBM_index_by_column_new();
  SBM_index_by_column_compute(B, B_index_by_column);

  /* allocation of C */
  SparseBlockStructuredMatrix*  C = SBM_calloc_multiply(A, B, B_index_by_column);

  SBM_index_by_column_free(B_index_by_column);

  DEBUG_END("SBM_zero_matrix_for_multiply\n");

  return C;
}



SparseBlockStructuredMatrix*  SBM_multiply(const SparseBlockStructuredMatrix* const A, const SparseBlockStructuredMatrix* const B)
{
  DEBUG_BEGIN("SBM_multiply\(...)\n");
  assert(A);
  assert(B);


  /*     Check the compatibility of the number and the sizes of blocks */
  int compat =  SBM_check_compatibility_for_multiply(A,B);
  if (!compat)
  {
    numerics_error("SBM_multiply", "Non compatible matrices or blocks sizes.\n");
    return NULL;
  }

  /* indexation of B by column */
  struct SBM_index_by_column* B_index_by_column = SBM_index_by_column_new();
  SBM_index_by_column_compute(B, B_index_by_column);

  /* allocation of C */
  SparseBlockStructuredMatrix*  C = SBM_calloc_multiply(A, B, B_index_by_column);

  size_t Bfilled3 =  B_index_by_column->filled3;
  size_t * Bindex3_data = B_index_by_column->index3_data;
  /* unsigned int Bfilled4 =  B_index_by_column->filled4; */
  size_t * Bindex4_data = B_index_by_column->index4_data;
  size_t * blockMap  = B_index_by_column->blockMap;

  double alpha = 1.0;
  double beta =  1.0; /* We assume that C is zero */

  size_t colNumberAA;
  size_t rowNumberBB;

  int Cnbblocks = -1;

  for (size_t currentRowNumberofA = 0 ; currentRowNumberofA < A->filled1 - 1; ++currentRowNumberofA)
  {
    int Ablocksize0 = A->blocksize0[currentRowNumberofA];
    if (currentRowNumberofA != 0)
      Ablocksize0  -= A->blocksize0[currentRowNumberofA - 1];
    for (size_t currentColNumberofB = 0 ; currentColNumberofB < Bfilled3 - 1; ++currentColNumberofB)
    {
      DEBUG_PRINTF("\n computation of the block of C(%zu,%zu)\n",currentRowNumberofA,currentColNumberofB );
      int Bblocksize1 = B->blocksize1[currentColNumberofB];
      if (currentColNumberofB != 0)
        Bblocksize1 -= B->blocksize1[currentColNumberofB - 1];
      /*           printf("of size %dX%d\n",Ablocksize0,Bblocksize1  ); */
      int CblockPassed = 1;
      for (size_t blockNumAA = A->index1_data[currentRowNumberofA];
           blockNumAA < A->index1_data[currentRowNumberofA + 1]; ++blockNumAA)
      {
        assert(blockNumAA < A->filled2);
        colNumberAA = A->index2_data[blockNumAA];
        DEBUG_PRINTF("blockNumAA = %zu, colNumberAA = %zu\n",blockNumAA,colNumberAA  );

        for (size_t blockNumBB = Bindex3_data[currentColNumberofB];
             blockNumBB < Bindex3_data[currentColNumberofB + 1]; ++blockNumBB)
        {
          rowNumberBB = Bindex4_data[blockNumBB];
          DEBUG_PRINTF("blockNumBB = %zu, rowNumberBB = %zu\n",blockNumBB,rowNumberBB  );
          DEBUG_PRINTF("blocMap[blockNumBB] = %zu, rowNumberBB = %zu\n",blockMap[blockNumBB],rowNumberBB  );

          if (rowNumberBB == colNumberAA)
          {
            if (CblockPassed)
            {
              Cnbblocks++; /* Find the right C block number*/
              CblockPassed = 0;
            }
            assert(Cnbblocks < (int) C->nbblocks);
            /*           printf("Compute C block number %i for %i %i ", Cnbblocks,currentRowNumberofA, currentColNumberofB ); */

            int Ablocksize1 = A->blocksize1[colNumberAA];
            if (colNumberAA != 0)
              Ablocksize1  -= A->blocksize1[colNumberAA - 1];

            int Bblocksize0 = B->blocksize0[rowNumberBB];
            if (rowNumberBB != 0)
              Bblocksize0 -= B->blocksize0[rowNumberBB - 1];


            DEBUG_PRINTF("Contribution of the product of blocks matrices A(%zu,%zu) (block number %zu) and B(%zu,%zu) (block number %zu) of  sizes %dX%d by %dX%d\n",
                         currentRowNumberofA,colNumberAA, blockNumAA, rowNumberBB, currentColNumberofB, blockMap[blockNumBB],
                         Ablocksize0,Ablocksize1,Bblocksize0,Bblocksize1   );


            assert(Ablocksize1 == Bblocksize0);

            /*            printf("DGEMM call\n"); */
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        Ablocksize0, Bblocksize1, Ablocksize1, alpha, A->block[blockNumAA],
                        Ablocksize0, B->block[blockMap[blockNumBB]], Bblocksize0,
                        beta, C->block[Cnbblocks], Ablocksize0);

            DEBUG_EXPR(NM_dense_display(A->block[blockNumAA], Ablocksize0, Ablocksize1, 0););
            DEBUG_EXPR(NM_dense_display(B->block[blockMap[blockNumBB]], Bblocksize0, Bblocksize1, 0););
            DEBUG_EXPR(NM_dense_display(C->block[Cnbblocks], Ablocksize0, Bblocksize1, 0););

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

  SBM_index_by_column_free(B_index_by_column);

  /* SBM_print(C); */


  DEBUG_END("SBM_multiply\(...)\n");
  return C;
}

static int SBM_check_C_for_gemm(const SparseBlockStructuredMatrix* const A,
                                const SparseBlockStructuredMatrix* const B,
                                SparseBlockStructuredMatrix*  C,
                                struct SBM_index_by_column * SBM_index_by_column_B)
{
  DEBUG_BEGIN("SBM_check_C_for_gemm(...)\n");

  assert(A);
  assert(B);
  assert(C);

  /*     Check the compatibility of size of matrices */
  assert(A->blocksize0);
  assert(A->blocksize1);
  assert(B->blocksize0);
  assert(B->blocksize1);
  assert(C->blocksize0);
  assert(C->blocksize1);

  assert(A->blocksize1[A->blocknumber1 - 1] == B->blocksize0[B->blocknumber0 - 1]);

  /*     Check the compatibility of the number and the sizes of blocks */

  if (!C)
  {
    return 0;
  }
  else
  {
    if (C->blocknumber0 != A->blocknumber0)
    {
      DEBUG_PRINT("blocknumber0 are not consistent\n");
      return 0;
    }
    if (C->blocknumber1 != B->blocknumber1)
    {
       DEBUG_PRINT("blocknumber1 are not consistent\n");
       return 0;
    }
  }

  size_t currentRowNumberofA;

  size_t colNumberAA;
  size_t rowNumberBB;

  for (currentRowNumberofA = 0 ; currentRowNumberofA < A->filled1 - 1; ++currentRowNumberofA)
  {

    if(C->blocksize0[currentRowNumberofA]!= A->blocksize0[currentRowNumberofA])
    {
      DEBUG_PRINTF("C->blocksize0[%zu] is not consistent with A->blocksize0[%zu]\n",currentRowNumberofA,currentRowNumberofA);
      return 0;
    }
    unsigned int Cblocksize0 = A->blocksize0[currentRowNumberofA];
    if (currentRowNumberofA != 0)
      Cblocksize0  -= A->blocksize0[currentRowNumberofA - 1];

    for (size_t currentColNumberofB = 0 ; currentColNumberofB < SBM_index_by_column_B->filled3 - 1; ++currentColNumberofB)
    {
      DEBUG_PRINT("\n");
      DEBUG_PRINTF("check for the block of C(%zu,%zu)\n",currentRowNumberofA,currentColNumberofB );
      int BlockCexists = 0;

      if(C->blocksize1[currentColNumberofB]!= B->blocksize1[currentColNumberofB])
      {
        DEBUG_PRINTF("C->blocksize1[%zu] is not consistent with B->blocksize1[%zu]\n",currentRowNumberofA,currentRowNumberofA);
        return 0;
      }
      for (size_t blockNumAA = A->index1_data[currentRowNumberofA];
           blockNumAA < A->index1_data[currentRowNumberofA + 1]; ++blockNumAA)
      {
        assert(blockNumAA < A->filled2);
        colNumberAA = A->index2_data[blockNumAA];

        for (size_t blockNumBB = SBM_index_by_column_B->index3_data[currentColNumberofB];
             blockNumBB < SBM_index_by_column_B->index3_data[currentColNumberofB + 1]; ++blockNumBB)
        {
          rowNumberBB = SBM_index_by_column_B->index4_data[blockNumBB];

          if (rowNumberBB == colNumberAA)
          {
            BlockCexists = 1;

            unsigned int Cblocksize1 = B->blocksize1[currentColNumberofB];
            if (currentColNumberofB != 0)
              Cblocksize1 -= B->blocksize1[currentColNumberofB - 1];
            DEBUG_PRINTF("C block number is needed for %zu %zu of size %dX%d\n",
                         currentRowNumberofA, currentColNumberofB,Cblocksize0,Cblocksize1);

            /* search for the block */
            size_t C_blockNum = 0;
            int block_C_found = 0;
            /* find the correct blockNum is C */
            for (size_t blockNumCC= C->index1_data[currentRowNumberofA];
                 blockNumCC < C->index1_data[currentRowNumberofA + 1]; ++blockNumCC)
            {
              if (C->index2_data[blockNumCC] == currentColNumberofB)
              {
                C_blockNum = blockNumCC;
                block_C_found = 1;
                DEBUG_PRINTF("Block in C (number %zu ) is found with corresponding row %zu and col %zu\n",
                             C_blockNum, currentRowNumberofA, currentColNumberofB);
              }
            }
            if(!block_C_found)
            {
              DEBUG_PRINTF("Block in C (number %zu ) is found with corresponding row %zu and col %zu\n",
                           C_blockNum, currentRowNumberofA, currentColNumberofB);
              return 0;
            }

            if (!C->block[C_blockNum])
            {
              DEBUG_PRINTF("Block in C (number %zu ) is found but not allocated\n",
                             C_blockNum);
              return 0;
            }

            break;
          };
        }
        if (BlockCexists) break;
      }
    }
  }

  DEBUG_END("SBM_check_C_for_gemm(...)\n");
  return 1;
}

void SBM_gemm_without_allocation(double alpha, const SparseBlockStructuredMatrix* const A,
                                 const SparseBlockStructuredMatrix* const B,
                                 double beta, SparseBlockStructuredMatrix*  C)
{
  DEBUG_BEGIN("SBM_gemm_without_allocation\(...)\n");

  /*  Check the compatibility of the number and the sizes of blocks */
  assert(SBM_check_compatibility_for_multiply(A,B));

  struct SBM_index_by_column* B_index_by_column = SBM_index_by_column_new();
  SBM_index_by_column_compute(B, B_index_by_column);


  /*  Check the memory allocation in C */
  assert(SBM_check_C_for_gemm(A, B, C, B_index_by_column));

  size_t Bfilled3 =  B_index_by_column->filled3;
  size_t  * Bindex3_data = B_index_by_column->index3_data;
  /* unsigned int Bfilled4 =  B_index_by_column->filled4; */
  size_t * Bindex4_data = B_index_by_column->index4_data;
  size_t * blockMap  = B_index_by_column->blockMap;

  size_t colNumberAA;
  size_t rowNumberBB;

  double local_beta=0.0;
  long int C_blockNum = -1;
  for (size_t currentRowNumberofA = 0 ; currentRowNumberofA < A->filled1 - 1; ++currentRowNumberofA)
  {
    int Ablocksize0 = A->blocksize0[currentRowNumberofA];
    if (currentRowNumberofA != 0)
      Ablocksize0  -= A->blocksize0[currentRowNumberofA - 1];
    for (size_t currentColNumberofB = 0 ; currentColNumberofB < Bfilled3 - 1; ++currentColNumberofB)
    {
      DEBUG_PRINTF("\n computation of the block of C(%zu,%zu)\n",currentRowNumberofA,currentColNumberofB );

      int Bblocksize1 = B->blocksize1[currentColNumberofB];
      if (currentColNumberofB != 0)
        Bblocksize1 -= B->blocksize1[currentColNumberofB - 1];
      int CblockPassed = 1;
      local_beta=beta;
      for (size_t blockNumAA = A->index1_data[currentRowNumberofA];
           blockNumAA < A->index1_data[currentRowNumberofA + 1]; ++blockNumAA)
      {
        assert(blockNumAA < A->filled2);
        colNumberAA = A->index2_data[blockNumAA];
        DEBUG_PRINTF("blockNumAA = %zu, colNumberAA = %zu\n",blockNumAA,colNumberAA  );


        for (size_t blockNumBB = Bindex3_data[currentColNumberofB];
             blockNumBB < Bindex3_data[currentColNumberofB + 1]; ++blockNumBB)
        {
          rowNumberBB = Bindex4_data[blockNumBB];
          DEBUG_PRINTF("blockNumBB = %zu, rowNumberBB = %zu\n",blockNumBB,rowNumberBB  );
          DEBUG_PRINTF("blocMap[blockNumBB] = %zu, rowNumberBB = %zu\n",blockMap[blockNumBB],rowNumberBB  );


          if (rowNumberBB == colNumberAA)
          {
            /* search for the block */
            if (CblockPassed)
            {
              CblockPassed = 0;
              C_blockNum = -1;
              /* find the correct blockNum is C */
              for (size_t blockNumCC= C->index1_data[currentRowNumberofA];
                   blockNumCC < C->index1_data[currentRowNumberofA + 1]; ++blockNumCC)
              {
                if (C->index2_data[blockNumCC] == currentColNumberofB)
                {
                  C_blockNum = blockNumCC;
                  DEBUG_PRINTF("Block in C (number %zu ) is found with corresponding row %zu and col %zu\n",
                               C_blockNum, currentRowNumberofA, currentColNumberofB);
                }
              }
            }
            assert(C_blockNum >=0);



            /*           printf("of size %dX%d\n",Ablocksize0,Bblocksize1  ); */

            int Ablocksize1 = A->blocksize1[colNumberAA];
            if (colNumberAA != 0)
              Ablocksize1  -= A->blocksize1[colNumberAA - 1];

            int Bblocksize0 = B->blocksize0[rowNumberBB];
            if (rowNumberBB != 0)
              Bblocksize0 -= B->blocksize0[rowNumberBB - 1];


            DEBUG_PRINTF("Contribution of the product of blocks matrices A(%zu,%zu) (block number %zu) and B(%zu,%zu) (block number %zu) of  sizes %dX%d by %dX%d\n",
                         currentRowNumberofA,colNumberAA, blockNumAA, rowNumberBB, currentColNumberofB, blockMap[blockNumBB],
                         Ablocksize0,Ablocksize1,Bblocksize0,Bblocksize1   );

            assert(Ablocksize1 == Bblocksize0);

            /*            printf("DGEMM call\n"); */
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        Ablocksize0, Bblocksize1, Ablocksize1, alpha, A->block[blockNumAA],
                        Ablocksize0, B->block[blockMap[blockNumBB]], Bblocksize0,
                        local_beta, C->block[C_blockNum], Ablocksize0);
            local_beta=1.0;
            DEBUG_EXPR(NM_dense_display(A->block[blockNumAA], Ablocksize0, Ablocksize1, 0););
            DEBUG_EXPR(NM_dense_display(B->block[blockMap[blockNumBB]], Bblocksize0, Bblocksize1, 0););
            DEBUG_EXPR(NM_dense_display(C->block[C_blockNum], Ablocksize0, Bblocksize1, 0););
          }
        }
      }
    }
  }

  SBM_index_by_column_free(B_index_by_column);

  DEBUG_END("SBM_gemm_without_allocation\(...)\n")
}
void SBM_row_prod(unsigned int sizeX, unsigned int sizeY, unsigned int currentRowNumber,
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
  size_t colNumber = 0;
  /* Number of rows/columns of the current block */
  unsigned int nbRows, nbColumns;
  /* Position of the sub-block of x multiplied by the sub-block of A */
  unsigned int posInX = 0;

  /* Check if currentRowNumber fits with A dimensions */
  assert(currentRowNumber <= A->blocknumber0);

  /* Get dim (rows) of the current block */
  nbRows = sizeY;

  /* if this is important, move it into a function --xhub */
  /*
    assert(
    {
      nbRows = A->blocksize0[currentRowNumber];
      if(currentRowNumber!=0)
        nbRows -= A->blocksize0[currentRowNumber-1];
      nbRows ==sizeY;
    }); */

  /* Set y to 0, if required */
  if (init == 1)
    cblas_dscal(sizeY, 0.0, y, 1);

  /* Loop over all non-null blocks
     Works whatever the ordering order of the block is, in A->block
     But it requires a set to 0 of all y components
  */

  for (size_t blockNum = A->index1_data[currentRowNumber];
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
    cblas_dgemv(CblasColMajor, CblasNoTrans, nbRows, nbColumns, 1.0, A->block[blockNum], nbRows, &x[posInX], 1, 1.0, y, 1);

  }
}
void SBM_row_prod_no_diag(unsigned int sizeX, unsigned int sizeY, unsigned int currentRowNumber, const SparseBlockStructuredMatrix* const A, const double* const x, double* y, int init)
{
  /*
     If: A is a SparseBlockStructuredMatrix matrix, Aij a block at row
     i and column j (Warning: i and j are indices of block position,
     not scalar component positions)

     Then SBM_row_prod_no_diag computes y = sum for i not equal to j of
     Aij.xj over a row of blocks (or += if init = false)

     currentRowNumber represents the position (block number) of the
     required line of blocks in the matrix A.

  */


  /* Column (block) position of the current block*/
  size_t colNumber = 0;

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

  /*  if this is important, move it into a function --xhub */
  /*   assert(
    {
      nbRows = A->blocksize0[currentRowNumber];
      if(currentRowNumber!=0)
        nbRows -= A->blocksize0[currentRowNumber-1];
      nbRows == sizeY ;
    });*/

  /* Set y to 0, if required */
  if (init == 1)
    cblas_dscal(sizeY, 0.0, y, 1);

  /* Loop over all non-null blocks. Works whatever the ordering order
     of the block is, in A->block, but it requires a set to 0 of all y
     components
  */
  for (size_t blockNum = A->index1_data[currentRowNumber];
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
      cblas_dgemv(CblasColMajor,CblasNoTrans, nbRows, nbColumns, 1.0, A->block[blockNum], nbRows, &x[posInX], 1, 1.0, y, 1);
    }
  }
}
void SBM_row_prod_no_diag_3x3(unsigned int sizeX, unsigned int sizeY, unsigned int currentRowNumber, const SparseBlockStructuredMatrix* const A, double* const x, double* y)
{
  /*
     If: A is a SparseBlockStructuredMatrix matrix, Aij a block at row
     i and column j (Warning: i and j are indices of block position,
     not scalar component positions)

     Then SBM_row_prod_no_diag computes y = sum for i not equal to j of
     Aij.xj over a row of blocks (or += if init = false)

     currentRowNumber represents the position (block number) of the
     required line of blocks in the matrix A.

  */


  /* Column (block) position of the current block*/
  size_t colNumber = 0;

  /* Number of columns of the current block */
  unsigned int nbColumns;

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

  /* Loop over all non-null blocks. Works whatever the ordering order
     of the block is, in A->block, but it requires a set to 0 of all y
     components
  */
  for (size_t blockNum = A->index1_data[currentRowNumber];
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
      /* cblas_dgemv(CblasColMajor,CblasNoTrans, nbRows, nbColumns, 1.0, A->block[blockNum], nbRows, &x[posInX], 1, 1.0, y, 1); */
      assert((nbColumns == 3));
      mvp3x3(A->block[blockNum], &x[posInX], y);
    }
  }
}
void SBM_row_prod_no_diag_1x1(unsigned int sizeX, unsigned int sizeY, unsigned int currentRowNumber, const SparseBlockStructuredMatrix* const A, double* const x, double* y)
{
  /*
     If: A is a SparseBlockStructuredMatrix matrix, Aij a block at row
     i and column j (Warning: i and j are indices of block position,
     not scalar component positions)

     Then SBM_row_prod_no_diag computes y = sum for i not equal to j of
     Aij.xj over a row of blocks (or += if init = false)

     currentRowNumber represents the position (block number) of the
     required line of blocks in the matrix A.

  */


  /* Column (block) position of the current block*/
  size_t colNumber = 0;

  /* Number of columns of the current block */
  unsigned int nbColumns;

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

  /* Loop over all non-null blocks. Works whatever the ordering order
     of the block is, in A->block, but it requires a set to 0 of all y
     components
  */
  for (size_t blockNum = A->index1_data[currentRowNumber];
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
      /* cblas_dgemv(CblasColMajor,CblasNoTrans, nbRows, nbColumns, 1.0, A->block[blockNum], nbRows, &x[posInX], 1, 1.0, y, 1); */
      assert((nbColumns == 1));
      y[0] += A->block[blockNum][0] * x[posInX] ;
    }
  }
}


void SBM_write_in_file(const SparseBlockStructuredMatrix* const m, FILE * file)
{
  DEBUG_PRINT("printInFileSBM\n");
  if (! m)
  {
    fprintf(stderr, "Numerics, SparseBlockStructuredMatrix printInFileSBM failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  assert(file);
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
  size_t colNumber;
  unsigned int nbRows, nbColumns;
  for (currentRowNumber = 0 ; currentRowNumber < m->filled1 - 1; ++currentRowNumber)
  {
    DEBUG_PRINTF("currentRowNumber = %i\n", currentRowNumber);
    DEBUG_PRINTF(" m->index1_data[currentRowNumber] = %zu\n",  m->index1_data[currentRowNumber]);
    DEBUG_PRINTF(" m->index1_data[currentRowNumber+1] = %zu\n",  m->index1_data[currentRowNumber+1]);

    for (size_t blockNum = m->index1_data[currentRowNumber];
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
      fprintf(file, SN_SIZE_T_F "\n", blockNum);
      DEBUG_PRINTF("nbRows * nbColumns = %i\n", nbRows * nbColumns);
      assert(m->block[blockNum]);
      for (unsigned int i = 0; i < nbRows * nbColumns; i++)
      {
        DEBUG_PRINTF("i = %i, blockNum = %zu,  m->block[%zu][%i] = %g\n ", i, blockNum, blockNum, i, m->block[blockNum][i] );
        fprintf(file, "%32.24e\n", m->block[blockNum][i]);
      }

    }
  }

}
void SBM_write_in_fileForScilab(const SparseBlockStructuredMatrix* const m, FILE * file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, SparseBlockStructuredMatrix SBM_write_in_file failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  assert(file);
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
  size_t colNumber;
  unsigned int nbRows, nbColumns;
  for (currentRowNumber = 0 ; currentRowNumber < m->filled1 - 1; ++currentRowNumber)
  {
    for (size_t blockNum = m->index1_data[currentRowNumber];
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
      fprintf(file, "block" SN_SIZE_T_F " = [ \n", blockNum);

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
  SBM_to_dense(m, denseMat);

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


void SBM_write_in_filename(const SparseBlockStructuredMatrix* const m, const char *filename)
{

}
void SBM_new_from_file(SparseBlockStructuredMatrix* const m, FILE *file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, SparseBlockStructuredMatrix SBM_read_in_file failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  assert(file);
  CHECK_IO(fscanf(file, "%d", &(m->nbblocks)));

  if (m->nbblocks == 0)  return;

  CHECK_IO(fscanf(file, "%d", &(m->blocknumber0)));
  CHECK_IO(fscanf(file, "%d", &(m->blocknumber1)));

  m->blocksize0 = (unsigned int *)malloc(m->blocknumber0 * sizeof(unsigned int));
  m->blocksize1 = (unsigned int *)malloc(m->blocknumber1 * sizeof(unsigned int));


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
  size_t colNumber;
  unsigned int nbRows, nbColumns;
  unsigned int blockk;
  for (currentRowNumber = 0 ; currentRowNumber < m->filled1 - 1; ++currentRowNumber)
  {
    for (size_t blockNum = m->index1_data[currentRowNumber];
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
        printf("Numerics, SparseBlockStructuredMatrix SBM_read_in_file failed, problem in block numbering. \n");
      }
      m->block[blockNum] = (double*)malloc(nbRows * nbColumns * sizeof(double));
      for (unsigned int i = 0; i < nbRows * nbColumns; i++)
      {
        CHECK_IO(fscanf(file, "%32le\n", &(m->block[blockNum][i])));
      }

    }
  }
}

void SBM_read_in_file(SparseBlockStructuredMatrix* const m, FILE *file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, SparseBlockStructuredMatrix SBM_read_in_file failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  assert(file);
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
  size_t colNumber;
  unsigned int nbRows, nbColumns;
  unsigned int blockk;
  for (currentRowNumber = 0 ; currentRowNumber < m->filled1 - 1; ++currentRowNumber)
  {
    for (size_t blockNum = m->index1_data[currentRowNumber];
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
        printf("Numerics, SparseBlockStructuredMatrix SBM_read_in_file failed, problem in block numbering. \n");
      }

      for (unsigned int i = 0; i < nbRows * nbColumns; i++)
      {
        CHECK_IO(fscanf(file, "%32le\n", &(m->block[blockNum][i])));
      }

    }
  }
}
void SBM_read_in_filename(SparseBlockStructuredMatrix* const m, const char *filename)
{

}
void SBM_free_pred(SparseBlockStructuredMatrixPred *blmatpred)
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


unsigned int * SBM_diagonal_block_indices(SparseBlockStructuredMatrix* const M)
{
  assert(M);
  if (M->diagonal_blocks) return M->diagonal_blocks;

  unsigned int * diagonal_blocks = (unsigned int*) malloc(M->blocksize0[M->blocknumber0-1]* sizeof(unsigned int));
  for (size_t currentRowNumber = 0 ; currentRowNumber < M->filled1 - 1; ++currentRowNumber)
  {
    for (size_t blockNum = M->index1_data[currentRowNumber];
         blockNum < M->index1_data[currentRowNumber + 1]; ++blockNum)
    {
      if (M->index2_data[blockNum] == currentRowNumber)
      {
        diagonal_blocks[currentRowNumber] = blockNum;
        break;
      }
      /* Warning: no check if the diagonal block does not exists */
    }
  }
  M->diagonal_blocks = diagonal_blocks;
  return diagonal_blocks;
}


unsigned int SBM_diagonal_block_index(SparseBlockStructuredMatrix* const M, unsigned int row)
{

  unsigned int * diagonal_blocks = SBM_diagonal_block_indices(M);
  return diagonal_blocks[row];
  
  /* /\* Look for the first block of row number num *\/ */
  /* unsigned int pos; */
  /* size_t firstBlockOfRow = M->index1_data[row]; */

  /* /\* Look at the diagonal block *\/ */
  /* for (pos = (unsigned int)firstBlockOfRow; M->index2_data[pos] != row; ++pos, assert(pos < M->filled2)); */

  /* return pos; */
}

int SBM_zentry(const SparseBlockStructuredMatrix* const M, unsigned int row, unsigned int col, double val)
{
  DEBUG_BEGIN("SBM_zentry(...)\n");
  DEBUG_PRINTF("row= %i, col =% i, val =%e\n", row, col, val);


  /* find the correct row number of blocks */
  unsigned int rowNumber =0;
  while (M->blocksize0[rowNumber] <= row)
  {
    rowNumber++;
    if (rowNumber >= M->blocknumber0)
    {
      numerics_warning("SBM_zentry", "The row number exceeds the size of the matrix");
      DEBUG_END("SBM_zentry(...)\n");
      return 0;
    }
  }
  DEBUG_PRINTF("rowNumber= %i\n", rowNumber);


  unsigned int colNumber =0;

  while (M->blocksize1[colNumber] <= col)
  {
    colNumber++;
    if (colNumber >= M->blocknumber1)
    {
      numerics_warning("SBM_zentry", "The col number exceeds the size of the matrix");
      DEBUG_END("SBM_zentry(...)\n");
      return 0;
    }
  }
  DEBUG_PRINTF("colNumber= %i\n", colNumber);

  size_t nbRows, nbColumns;
  size_t row_pos=0, col_pos=0;
  for (size_t blockNum = M->index1_data[rowNumber];
       blockNum < M->index1_data[rowNumber + 1]; ++blockNum)
  {
    assert(blockNum < M->filled2);
    if (colNumber == M->index2_data[blockNum])
    {
      DEBUG_PRINTF("blockNum = %zu\n", blockNum);
      /* Get dim. of the current block */
      nbRows = M->blocksize0[rowNumber];
      row_pos = row;
      if (rowNumber != 0)
      {
        nbRows -= M->blocksize0[rowNumber - 1];
        row_pos -= M->blocksize0[rowNumber - 1];
      }
      assert(nbRows);
      assert(row_pos < nbRows);

      nbColumns = M->blocksize1[colNumber];
      col_pos=col;
      if (colNumber != 0)
      {
        nbColumns -= M->blocksize1[colNumber - 1];
        col_pos -= M->blocksize1[colNumber - 1];
      }
      assert(nbColumns);
      assert(col_pos < nbColumns);

      DEBUG_PRINTF("row_pos = %zu, col_pos = %zu \n", row_pos, col_pos);

      M->block[blockNum][row_pos + col_pos * nbRows] = val;
      DEBUG_END("SBM_zentry(...)\n");
      return 1;
    }
  }
  numerics_warning("SBM_zentry", "no existing block for inserting entry\n");
  DEBUG_END("SBM_zentry(...)\n");
  return 0;
}

double SBM_get_value(const SparseBlockStructuredMatrix* const M, unsigned int row, unsigned int col)
{
  /*      Search the row of blocks and the column of blocks */
  unsigned int rowblocknumber = M->blocknumber0;
  unsigned int colblocknumber = M->blocknumber1;
  int rowpos = -1;
  int colpos = -1;
  int blockNumber = -1;
  size_t colNumber;
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
      for (size_t blockNum = M->index1_data[rowblocknumber];
           blockNum < M->index1_data[rowblocknumber + 1]; ++blockNum)
      {
        assert(blockNum < M->filled2);

        colNumber = M->index2_data[blockNum];
        if (colblocknumber == colNumber)
        {
          blockNumber = (int)blockNum;
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

int SBM_copy(const SparseBlockStructuredMatrix* const A, SparseBlockStructuredMatrix*  B, unsigned int copyBlock)
{
  assert(A);
  assert(B);

  int need_blocks = 0;

  if (B->nbblocks < A->nbblocks)
  {
    need_blocks = 1;
    for (unsigned i=0; i<B->nbblocks; ++i)
    {
      free(B->block [i]);
      B->block [i] = NULL;
    }
    B->block = (double **) realloc(B->block, A->nbblocks * sizeof(double *));
  }
  B->nbblocks = A->nbblocks;

  if (B->blocknumber0 < A->blocknumber0)
  {
    B->blocksize0 = (unsigned int*) realloc(B->blocksize0, A->blocknumber0 * sizeof(unsigned int));
  }
  B->blocknumber0 = A->blocknumber0;

  if (B->blocknumber1 < A->blocknumber1)
  {
    B->blocksize1 = (unsigned int*) realloc(B->blocksize1, A->blocknumber1 * sizeof(unsigned int));
  }
  B->blocknumber1 = A->blocknumber1;

  if (B->filled1 < A->filled1)
  {
    B->index1_data = (size_t*) realloc(B->index1_data, A->filled1 * sizeof(size_t));
  }
  B->filled1 = A->filled1;

  if (B->filled2 < A->filled2)
  {
    B->index2_data = (size_t*) realloc(B->index2_data, A->filled2 * sizeof(size_t));
  }
  B->filled2 = A->filled2;


  memcpy(B->blocksize0, A->blocksize0, A->blocknumber0 * sizeof(unsigned int));
  memcpy(B->blocksize1, A->blocksize1, A->blocknumber1 * sizeof(unsigned int));
  memcpy(B->index1_data, A->index1_data, A->filled1 * sizeof(size_t));
  memcpy(B->index2_data, A->index2_data, A->filled2 * sizeof(size_t));

  if (copyBlock)
  {
    unsigned int currentRowNumber ;
    size_t colNumber;
    unsigned int nbRows, nbColumns;
    for (currentRowNumber = 0 ; currentRowNumber < A->filled1 - 1; ++currentRowNumber)
    {
      for (size_t blockNum = A->index1_data[currentRowNumber];
           blockNum < A->index1_data[currentRowNumber + 1]; ++blockNum)
      {
        assert(blockNum < A->filled2);
        colNumber = A->index2_data[blockNum];
        /* Get dim. of the current block */
        nbRows = A->blocksize0[currentRowNumber];
        if (currentRowNumber != 0)
          nbRows -= A->blocksize0[currentRowNumber - 1];
        nbColumns = A->blocksize1[colNumber];

        if (colNumber != 0)
          nbColumns -= A->blocksize1[colNumber - 1];

        if (need_blocks)
        {
          B->block[blockNum] = (double*)malloc(nbRows * nbColumns * sizeof(double));
        }

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

  if (B->diagonal_blocks)
  {
    free(B->diagonal_blocks);
    B->diagonal_blocks=NULL;
  }
 
  if (A->diagonal_blocks)
  {
     B->diagonal_blocks = (unsigned int*) malloc(A->blocksize0[A->blocknumber0-1]* sizeof(unsigned int));
     memcpy(B->diagonal_blocks,A->diagonal_blocks, A->blocksize0[A->blocknumber0-1]* sizeof(unsigned int));
  }



  
  return 0;
}











int SBM_transpose(const SparseBlockStructuredMatrix* const A, SparseBlockStructuredMatrix*  B)
{
  assert(A);
  assert(B);
  B->nbblocks = A->nbblocks;
  B->blocknumber0 = A->blocknumber1;
  B->blocknumber1 = A->blocknumber0;
  B->blocksize0 = (unsigned int*)malloc(B->blocknumber0 * sizeof(unsigned int));
  for (unsigned int i = 0; i < B->blocknumber0; i++) B->blocksize0[i] = A->blocksize1[i];
  B->blocksize1 = (unsigned int*)malloc(B->blocknumber1 * sizeof(unsigned int));
  for (unsigned int i = 0; i < B->blocknumber1; i++) B->blocksize1[i] = A->blocksize0[i];

  B->filled1 = A->blocknumber1 + 1;;
  B->filled2 = A->filled2;
  B->index1_data = (size_t *)malloc(B->filled1 * sizeof(size_t));
  B->index2_data = (size_t *)malloc(B->filled2 * sizeof(size_t));


  unsigned int currentRowNumberofA;
  unsigned int currentColNumberofA ;
  size_t colNumberofA;

  size_t * blockMap  = (size_t *)malloc(B->filled2 * sizeof(size_t));



  int blockNumB = -1;
  B->index1_data[0] = 0;
  for (currentColNumberofA = 0 ; currentColNumberofA < A->blocknumber1; ++currentColNumberofA)
  {
    assert(currentColNumberofA + 1 < B->filled1);
    B->index1_data[currentColNumberofA + 1] = B->index1_data[currentColNumberofA];
    for (currentRowNumberofA = 0 ; currentRowNumberofA < A->blocknumber0; ++currentRowNumberofA)
    {
      for (size_t blockNum = A->index1_data[currentRowNumberofA];
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
  size_t colNumber;
  unsigned int nbRows, nbColumns;
  for (currentRowNumber = 0 ; currentRowNumber < B->filled1 - 1; ++currentRowNumber)
  {
    for (size_t blockNum = B->index1_data[currentRowNumber];
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
int SBM_inverse_diagonal_block_matrix_in_place(const SparseBlockStructuredMatrix*  M,  int* ipiv)
{
  for (unsigned int i = 0; i < M->filled1 - 1; i++)
  {
    size_t numberofblockperrow = M->index1_data[i + 1] - M->index1_data[i];
    if (numberofblockperrow != 1)
    {
      fprintf(stderr, "SparseBlockMatrix : SBM_inverse_diagonal_block_matrix: Not a diagonal block matrix\n");
      exit(EXIT_FAILURE);
    }
  }
  for (unsigned int i = 0; i < M->filled2; i++)
  {
    if (M->index2_data[i] != i)
    {
      fprintf(stderr, "SparseBlockMatrix : SBM_inverse_diagonal_block_matrix: Not a diagonal block matrix\n");
      exit(EXIT_FAILURE);
    }
  }

  unsigned int currentRowNumber ;
  size_t colNumber;
  unsigned int nbRows, nbColumns;
  lapack_int infoDGETRF = 0;
  lapack_int infoDGETRI = 0;
  int info = 0;

  lapack_int* lapack_ipiv = (lapack_int *) ipiv;


  for (currentRowNumber = 0 ; currentRowNumber < M->filled1 - 1; ++currentRowNumber)
  {
    for (size_t blockNum = M->index1_data[currentRowNumber];
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

      DGETRF(nbRows, nbColumns, M->block[blockNum], nbRows, lapack_ipiv, &infoDGETRF);
      assert(!infoDGETRF);

      DGETRI(nbRows, M->block[blockNum], nbRows, lapack_ipiv, &infoDGETRI);
      assert(!infoDGETRI);

    }
  }
  if ((!infoDGETRF) || (!infoDGETRI)) info = 0;

  return info;


}

void SBM_to_dense(const SparseBlockStructuredMatrix* const A, double *denseMat)
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
      denseMat[i + j * n] = SBM_get_value(A, i, j);
    }
  }
}

int SBM_to_sparse_init_memory(const SparseBlockStructuredMatrix* const A, CSparseMatrix *sparseMat)
{
  assert(A);
  assert(A->blocksize0);
  assert(A->blocksize1);
  int n = A->blocksize0[A->blocknumber0 - 1];
  int m = A->blocksize1[A->blocknumber1 - 1];

  sparseMat->m = n;
  sparseMat->n = m;

  sparseMat->nz = -2; /* csr */
  sparseMat->nzmax = 0;
  sparseMat->p = (CS_INT*)malloc((sparseMat->m + 1) * sizeof(CS_INT));

  /* Row (block) position of the current block */
  unsigned int currentRowNumber ;
  /* Column (block) position of the current block*/
  size_t colNumber;
  /* Number of rows/columns of the current block */
  int nbRows, nbColumns;

  for (currentRowNumber = 0 ; currentRowNumber < A->filled1 - 1; ++currentRowNumber)
  {
    for (size_t blockNum = A->index1_data[currentRowNumber];
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
  sparseMat->i = (CS_INT*)malloc((sparseMat->nzmax) * sizeof(CS_INT));
  sparseMat->x = (double*)malloc((sparseMat->nzmax) * sizeof(double));

  return 0;
}

sparse_matrix_iterator sparseMatrixBegin(const CSparseMatrix* const sparseMat)
{
  if (sparseMat->nz >= 0)
  {
    sparse_matrix_iterator i = { 0, -1, 0, 0, NAN, sparseMat};
    return i;
  }
  else if (sparseMat->nz == -1)
  {
    sparse_matrix_iterator i = { 0, sparseMat->p[0], 0, 0, NAN, sparseMat };
    return i;
  }
  else if (sparseMat->nz == -2)
  {
    sparse_matrix_iterator i = { 0, sparseMat->p[0], 0, 0, NAN, sparseMat };
    return i;
  }

  assert(0);
  sparse_matrix_iterator z = { 0, 0, 0, 0, NAN, sparseMat };
  return z;
}

int sparseMatrixNext(sparse_matrix_iterator* it)
{
  if (it->mat->nz >= 0) /* triplet */
  {
    assert(it->counter2 == -1);

    if (it->counter1 < it->mat->nz)
    {
      assert(it->counter1 >= 0);

      it->first = it->mat->i[it->counter1];
      it->second = it->mat->p[it->counter1];
      it->third = it->mat->x[it->counter1];

      it->counter1++;

      return 1;
    }
    else
    {
      return 0; /* stop */
    }
  }
  else if (it->mat->nz == -1) /* csc */
  {
    if (it->counter1 < it->mat->n)
    {
      if (it->counter2 < it->mat->p[it->counter1 + 1])
      {
        /* next line */
        assert(it->counter2 >= 0);
        assert(it->counter2 < it->mat->nzmax);

        it->first = it->mat->i[it->counter2];
        it->second = it->counter1;
        it->third = it->mat->x[it->counter2];
        it->counter2++;
        return 1;
      }
      else
      {
        /* next column */
        it->counter1++;

        assert(it->counter1 >= 0);
        assert(it->counter1 < (it->mat->n + 1));
        it->counter2 = it->mat->p[it->counter1];

        DEBUG_PRINTF("it->counter1 = %ld, it->counter2 = %ld\n", it->counter1, it->counter2);

        if (it->counter2 < it->mat->nzmax)
        {

          assert(it->counter2 >= 0);
          assert(it->counter2 < it->mat->nzmax);

          it->first = it->mat->i[it->counter2];
          it->second = it->counter1;
          it->third = it->mat->x[it->counter2];

          it->counter2++;
          return 1;
        }
        else
        {
          return 0; /* stop */
        }
      }
    }
    return 0; /* stop */
  }
  else if (it->mat->nz == -2) /* csr */
  {
    if (it->counter1 < it->mat->m)
    {
      if (it->counter2 < it->mat->p[it->counter1 + 1])
      {
        /* next column */

        assert(it->counter2 >= 0);
        assert(it->counter2 < it->mat->nzmax);

        it->first = it->counter1;
        it->second = it->mat->i[it->counter2];
        it->third = it->mat->x[it->counter2];
        it->counter2++;
        return 1;
      }
      else
      {
        /* next line */
        it->counter1++;

        assert(it->counter1 >= 0);
        assert(it->counter1 < (it->mat->m + 1));

        it->counter2 = it->mat->p[it->counter1];

        DEBUG_PRINTF("it->counter1 = %ld, it->counter2 = %ld\n", it->counter1, it->counter2);

        if (it->counter2 < it->mat->nzmax)
        {
          assert(it->counter2 >= 0);
          assert(it->counter2 < it->mat->nzmax);

          it->first = it->counter1;
          it->second = it->mat->i[it->counter2];
          it->third = it->mat->x[it->counter2];

          it->counter2++;
          return 1;
        }
        else
        {
          return 0; /* stop */
        }
      }
    }
    else
    {
      return 0; /* stop */
    }
  }

  assert(0);
  return 0;
}


void SBCM_null(SparseBlockCoordinateMatrix* MC)
{
  MC->blocknumber0 = 0;
  MC->blocknumber1 = 0; 
  MC->nbblocks = 0;
  
  MC->row = NULL;
  MC->column = NULL;
  MC->blocksize0 = NULL;
  MC->blocksize1 = NULL;
  MC->block=NULL ;
  
}

SparseBlockCoordinateMatrix* SBCM_new(void)
{
  SparseBlockCoordinateMatrix* MC = (SparseBlockCoordinateMatrix*) malloc(sizeof(SparseBlockCoordinateMatrix));
  SBCM_null(MC);
  return MC;
}



SparseBlockCoordinateMatrix*  SBCM_new_3x3(unsigned int m, unsigned int n,
    unsigned int nbblocks,
    unsigned int *row,
    unsigned int *column,
    double *block)
{
  
  SparseBlockCoordinateMatrix* MC = SBCM_new();
  MC->blocknumber0 = m; /* block row */
  MC->blocknumber1 = n; /* block column */
  MC->nbblocks = nbblocks;
  MC->row = (unsigned int *)  malloc(sizeof(unsigned int) * nbblocks);
  MC->column = (unsigned int *)  malloc(sizeof(unsigned int) * nbblocks);
  for (unsigned int i = 0; i < nbblocks; ++i)
  {
    MC->row[i] = row[i] - 1;
  }
  for (unsigned int i = 0; i < nbblocks; ++i)
  {
    MC->column[i] = column[i] - 1;
  }
  MC->block = (double **) malloc(sizeof(double*)*nbblocks);
  MC->blocksize0 = (unsigned int *) malloc(sizeof(unsigned int) * nbblocks);
  MC->blocksize1 = (unsigned int *) malloc(sizeof(unsigned int) * nbblocks);
  for (unsigned int i = 0; i < nbblocks; ++i)
  {
    MC->blocksize0[i] = 3 * (i + 1);
    MC->blocksize1[i] = 3 * (i + 1);
    MC->block[i] = &block[i * 9];
  }

  return MC;
}

void  SBCM_free_3x3(SparseBlockCoordinateMatrix *MC)
{
  free(MC->block);
  free(MC->blocksize0);
  free(MC->blocksize1);
  free(MC->row);
  free(MC->column);
}

/* quite obvious alg but in case of incomprehension see Nathan Bell coo_tocsr */
/* i.e coo.h file under scipy sparsetools */
SparseBlockStructuredMatrix* SBCM_to_SBM(SparseBlockCoordinateMatrix* MC)
{
  SparseBlockStructuredMatrix* M = SBM_new();

  M->nbblocks = MC->nbblocks;
  M->filled2 = MC->nbblocks;
  M->blocknumber0 = MC->blocknumber0;
  M->blocknumber1 = MC->blocknumber1;
  M->blocksize0 = MC->blocksize0;
  M->blocksize1 = MC->blocksize1;

  M->index1_data = (size_t *) calloc(M->blocknumber0 + 1, sizeof(size_t));
  M->index2_data = (size_t *) malloc(sizeof(size_t) * M->nbblocks);
  M->filled1 = 2;
  for (unsigned int i = 0; i < M->nbblocks; ++i)
  {
    M->index1_data[MC->row[i] + 1]++;
    M->filled1 = ((MC->row[i] + 2) > M->filled1) ? MC->row[i] + 2 : M->filled1;
  }

  for (unsigned int i = 1; i < M->blocknumber0 + 1; ++i)
  {
    M->index1_data[i] += M->index1_data[i - 1];
  }

  M->block = (double **) malloc(sizeof(double*)*M->nbblocks);

  for (unsigned int i = 0; i < M->nbblocks; ++i)
  {
    unsigned int row  = MC->row[i];
    size_t dest  = M->index1_data[row];

    assert(dest < M->nbblocks);
    M->index2_data[dest] = MC->column[i];

    assert(dest < M->nbblocks);
    M->block[dest] = MC->block[i];

    assert(row < (M->blocknumber0 + 1));
    M->index1_data[row]++;
  }

  size_t last = 0;
  for (unsigned int i = 0; i <= M->blocknumber0; i++)
  {
    size_t temp = M->index1_data[i];
    M->index1_data[i]  = last;
    last = temp;
  }

  return M;
}

void SBM_free_from_SBCM(SparseBlockStructuredMatrix* M)
{
  free(M->index1_data);
  free(M->index2_data);
  free(M->block);
  free(M);
  M = NULL;
}

int SBM_from_csparse(int blocksize, const CSparseMatrix* const sparseMat, SparseBlockStructuredMatrix* A)
{
  DEBUG_PRINT("SBM_from_csparse start\n")
  assert(sparseMat);
  assert(sparseMat->p);
  assert(sparseMat->i);
  assert(sparseMat->x);
  /* assert(sparseMat->nz == -2); */

  assert(sparseMat->m % blocksize == 0);
  assert(sparseMat->n % blocksize == 0);

  CS_INT bnrow = sparseMat->m / blocksize;
  CS_INT bncol = sparseMat->n / blocksize;
  DEBUG_PRINTF("SBM_from_csparse. bnrow =%li\n", bnrow);
  DEBUG_PRINTF("SBM_from_csparse. bncol =%li\n", bncol);
  A->blocknumber0 = (int) bnrow;
  A->blocknumber1 = (int) bncol;

  //  assert(A->blocksize0 == NULL);
  A->blocksize0 = (unsigned int*) malloc(A->blocknumber0 * sizeof(unsigned int));

  //  assert(A->blocksize1 == NULL);
  A->blocksize1 = (unsigned int*) malloc(A->blocknumber1 * sizeof(unsigned int));

  for (unsigned int i = 0; i < A->blocknumber0; i++)
  {
    A->blocksize0[i] = (i + 1) * blocksize;
  }

  for (unsigned int i = 0; i < A->blocknumber1; i++)
  {
    A->blocksize1[i] = (i + 1) * blocksize;
  }

  /* we have to find non empty blocks */

  CS_INT blockindexmax = -1;
  CS_INT blocklinemax = -1;
  int* blockline;
  int* blocknum;

  /* 1: find blockindexmax (<= bnrow + bncol * bnrow) & blocklinemax */
  for (sparse_matrix_iterator it = sparseMatrixBegin(sparseMat);
       sparseMatrixNext(&it);)
  {

    CS_INT row = it.first;
    CS_INT col = it.second;

    DEBUG_PRINTF("it.first = %li, it.second = %li \n", row, col );
    DEBUG_PRINTF("it.third = %g,  \n", it.third );

    CS_INT brow = row / blocksize;
    CS_INT bcol = col / blocksize;
    CS_INT blockindex = brow * bncol + bcol;

    if ((fabs(it.third) > 0.0) && (blockindex > blockindexmax - 1))
    {
      blockindexmax = blockindex + 1;
    };

    if ((fabs(it.third) > 0.0) && brow > (blocklinemax - 1))
    {
      blocklinemax = brow + 1;
    }
  }
  DEBUG_PRINTF("SBM_from_csparse. blockindexmax =%li\n",blockindexmax);
  DEBUG_PRINTF("SBM_from_csparse. blocklinemax=%li\n", blocklinemax);

  // assert(blockindexmax <= bnrow + bncol * bnrow + 1);
  // assert(blocklinemax <= bnrow + 1);

  /* 2: allocate temporary memory for blocknumbers & blocklines */
  blocknum = (int *) malloc(blockindexmax * sizeof(int));
  blockline = (int *) malloc(blocklinemax * sizeof(int));
  for (int i = 0; i < blockindexmax; i++)
  {
    blocknum[i] = -1;
  }

  for (int i = 0; i < blocklinemax; i++)
  {
    blockline[i] = -1;
  }

  /* 3: flag non empty blocks & lines */
  for (sparse_matrix_iterator it = sparseMatrixBegin(sparseMat);
       sparseMatrixNext(&it);)
  {
    CS_INT row = it.first;
    CS_INT col = it.second;

    CS_INT brow = row / blocksize;
    CS_INT bcol = col / blocksize;

    CS_INT blockindex = brow * bncol + bcol;

    if (fabs(it.third) > 0.)
    {
      assert(blockindex < blockindexmax);
      blocknum[blockindex] = -2;
      assert(brow < blocklinemax);
      blockline[brow] = -2;
    }
  }

  /* 4: count non empty blocks */
  A->nbblocks = 0;
  for (int i = 0; i < blockindexmax; i++)
  {
    assert(blocknum[i] == -1 || blocknum[i] == -2);

    if (blocknum[i] == -2)
    {
      DEBUG_PRINTF("blocknum[%d] = %d\n", i, A->nbblocks);
      blocknum[i] = A->nbblocks++;
    }
  }

  /* 5: allocate memory for contiguous blocks */
  assert(A->nbblocks>0);

  if (!A->block)
    A->block = (double **) malloc(A->nbblocks * sizeof(double *));
  for (unsigned int i = 0; i < A->nbblocks; i++)
  {
    A->block[i] = (double *) malloc(blocksize * blocksize * sizeof(double));

    /* fill block with 0 */
    for (int birow = 0; birow < blocksize; ++birow)
    {
      for (int bicol = 0; bicol < blocksize; ++bicol)
      {
        A->block[i][birow + bicol * blocksize] = 0.;
      }
    }
  }

  A->filled2 = (size_t) A->nbblocks; /* one of them should be deprecated! */

  /* 6: count non empty lines */
  A->filled1 = 1; /* A->filled1 = number of non empty lines + 1 */

  for (int i = 0; i < blocklinemax; i++)
  {
    assert(blockline[i] == -1 || blockline[i] == -2);

    if (blockline[i] == -2)
    {
      A->filled1++;
    }
  }
  DEBUG_PRINTF("A->filled1 =%zu\n",A->filled1);
  DEBUG_PRINTF("A->filled2 =%zu\n",A->filled2);
  /* 7: allocate memory for index data */

  if (!A->index1_data)
    A->index1_data = (size_t*) malloc(A->filled1 * sizeof(size_t));

  if (!A->index2_data)
    A->index2_data = (size_t*) malloc(A->filled2 * sizeof(size_t));


  for (size_t i = 0; i < A->filled1; i++)
  {
    A->index1_data[i] = blockindexmax + 1;
  }

  for (size_t i = 0; i < A->filled2; i++)
  {
    A->index2_data[i] = 0;
  }

  /* 8: fill index1_data & index2_data & copy values in contiguous
   * blocks */
  for (sparse_matrix_iterator it = sparseMatrixBegin(sparseMat);
       sparseMatrixNext(&it);)
  {

    CS_INT row = it.first;
    CS_INT col = it.second;

    assert(row < sparseMat->m);
    assert(col < sparseMat->n);

    CS_INT brow = row / blocksize;
    CS_INT bcol = col / blocksize;

    assert(brow < bnrow);
    assert(bcol < bncol);

    CS_INT blockindex = brow * bncol + bcol;

    CS_INT birow = row % blocksize; /* block inside row */
    CS_INT bicol = col % blocksize; /* block inside column */

    A->index1_data[A->filled1 - 1] = A->filled2;

    if ((blockindex < blockindexmax) && (blocknum[blockindex] >= 0))
    {

      /* this is an non empty block */

      /* index1_data[rowNumber]<= blockNumber */

      assert(brow < (CS_INT)A->filled1);
      if (A->index1_data[brow] > (size_t)blocknum[blockindex])
      {
        A->index1_data[brow] = blocknum[blockindex];
      }

      A->index2_data[blocknum[blockindex]] = bcol;

      assert(birow + bicol * blocksize <= blocksize * blocksize);

      assert(blockindex < blockindexmax);
      assert(blocknum[blockindex] < (CS_INT)A->nbblocks);
      assert(blocknum[blockindex] < (CS_INT)A->filled2);

      DEBUG_PRINTF("A->block[blocknum[blockindex=%ld]=%d][birow=%ld + bicol=%ld * blocksize=%d] = it.third=%g\n",
                   blockindex, blocknum[blockindex], birow, bicol, blocksize, it.third);
      A->block[blocknum[blockindex]][birow + bicol * blocksize] = it.third;
    }
  }

  /* 9: free temp memory */
  free(blocknum);
  free(blockline);
  DEBUG_PRINT("SBM_from_csparse end\n")

  return 0;

}

int  SBM_to_sparse(const SparseBlockStructuredMatrix* const A, CSparseMatrix *outSparseMat)
{
  assert(A);
  assert(A->blocksize0);
  assert(A->blocksize1);

  assert(outSparseMat);
  assert(outSparseMat->p);
  assert(outSparseMat->i);
  assert(outSparseMat->x);

  /* Row (block) position of the current block */
  unsigned int currentRowNumber ;
  /* Column (block) position of the current block*/
  size_t colNumber;
  /* Number of rows/columns of the current block */
  int nbRows, nbColumns;

  int nnz = 0;
  int isparserowend, isparserowstart, isparserow;
  int isparsecolumnend, isparsecolumnstart, isparsecolumn;
  outSparseMat->p[0] = 0; /* We assume that the first row is non empty */
  for (currentRowNumber = 0 ; currentRowNumber < (unsigned int) A->filled1 - 1; ++currentRowNumber)
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


      outSparseMat->p[isparserow + 1] =  outSparseMat->p[isparserow];

      for (size_t blockNum = A->index1_data[currentRowNumber];
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

        outSparseMat->p[isparserow + 1] += nbColumns  ;

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
          outSparseMat->i[nnz] = isparsecolumn;

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

          outSparseMat->x[nnz] = A->block[blockNum][colintheblock * nbRows + rowintheblock];
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
  {
    printf("val1 : %d", NUMERICS_SBM_FREE_SBM);
    printf("val2 : %d", level);
    free(A);
  }
}

//#define SBM_DEBUG_SBM_row_to_dense
void SBM_row_to_dense(const SparseBlockStructuredMatrix* const A, int row, double *denseMat, int rowPos, int rowNb)
{
  assert(A);
  int BlockRowNb = 0;
  int ColNb = A->blocksize1[A->blocknumber1 - 1];
  if (row)
    BlockRowNb = A->blocksize0[row] - A->blocksize0[row - 1];
  else
    BlockRowNb = A->blocksize0[row];
#ifdef SBM_DEBUG_SBM_row_to_dense
  printf("SBM_row_to_dense : copi block row %i, containing %i row and %i col.\n", row, BlockRowNb, ColNb);
#endif

  //zero memory
  for (int numRow = rowPos; numRow < rowPos + BlockRowNb; numRow++)
    for (int numCol = 0; numCol < ColNb; numCol++)
      denseMat[numRow + numCol * rowNb] = 0.0;

  //index1_data[rowNumber]<= blockNumber <index1_data[rowNumber+1]
  for (size_t numBlock = A->index1_data[row]; numBlock < A->index1_data[row + 1]; numBlock++)
  {
    double * beginBlock = A->block[numBlock];
    size_t colNumber = A->index2_data[numBlock];
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
#ifdef SBM_DEBUG_SBM_row_to_dense
  printf("SBM_row_to_dense : res in file SBM_row_to_dense.txt.");
  FILE * titi  = fopen("SBM_row_to_dense.txt", "w");
  SBM_write_in_fileForScilab(A, titi);
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
void SBM_row_permutation(unsigned int *rowIndex, SparseBlockStructuredMatrix* A, SparseBlockStructuredMatrix*  C)
{
#ifdef SBM_DEBUG_SBM_ROW_PERM
  FILE * titi  = fopen("SBM_row_permutation_input.txt", "w");
  SBM_write_in_fileForScilab(A, titi);
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
    printf("SBM_row_permutation rowA=%i, rowC=%i\n", rowA, rowC);
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
  titi  = fopen("SBM_row_permutation_output.txt", "w");
  SBM_write_in_fileForScilab(C, titi);
  fclose(titi);
#endif
}
//#define SBM_DEBUG_SBM_COL_PERM
void SBM_column_permutation(unsigned int *colIndex, SparseBlockStructuredMatrix* A, SparseBlockStructuredMatrix*  C)
{
#ifdef SBM_DEBUG_SBM_COL_PERM
  FILE * titi  = fopen("SBM_column_permutation_input.txt", "w");
  SBM_write_in_fileForScilab(A, titi);
  fclose(titi);
#endif
  SBM_copy(A, C, 0);
  for (unsigned int n = 0; n < C->nbblocks; n++)
  {
    C->index2_data[n] = colIndex[C->index2_data[n]];
  }
  //int curColnb=0;
  int nbBlockCol = A->blocknumber1;
  for (int numCol = 0; numCol < nbBlockCol; numCol++)
  {
#ifdef SBM_DEBUG_SBM_COL_PERM
    printf("SBM_column_permutation colA=%i, colC=%i\n", numCol, colIndex[numCol]);
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
  titi  = fopen("SBM_column_permutation_output.txt", "w");
  SBM_write_in_fileForScilab(C, titi);
  fclose(titi);
#endif
}
