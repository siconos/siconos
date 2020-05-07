
/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

/*
  Tests functions for NumericsMatrix structure

 */

#include "numericsMatrixTestFunction.h"
#include <math.h>                  // for fabs
#include <stdio.h>                 // for printf
#include <stdlib.h>                // for malloc, size_t, free, rand, calloc
#include "NumericsMatrix.h"        // for NumericsMatrix, NM_new, NM_clearDense
#include "NumericsSparseMatrix.h"  // for NSM_CSC, NumericsSparseMatrix
#include "SparseBlockMatrix.h"     // for SparseBlockStructuredMatrix, SBM_new

NumericsMatrix * test_matrix_1()
{
  NumericsMatrix * M1 = NM_new();
  int n = 8;
  /* Double * storage (column-major) */
  double m0[64] = {1, 2, 0, 5, 0, 0, 0, 0,
                   2, 1, 0, 0, 0, 0, 0, 0,
                   0, 0, 1, -1, 0, 0, 2, 2,
                   4, 0, -1, 6, 0, 0, 1, 2,
                   3, 4, 0, 0, 1, 0, 0, 0,
                   -1, 1, 0, 6, 0, 2, 0, 0,
                   0, 0, 0, 0, 0, 0, 2, -1,
                   0, 0, 0, 0, 5, 2, 2, 2
                  };
  M1->storageType = 0;
  M1->size0 = n;
  M1->size1 = n;
  M1->matrix0 = (double *)malloc(n * n * sizeof(double));


  for(int i = 0; i < n * n; i++)
    M1->matrix0[i] = m0[i];
  return M1;
}

NumericsMatrix * test_matrix_3()
{
  NumericsMatrix * M3 = NM_new();
  int n = 8;

  int nn = 4;
  /* Double * storage (column-major) */
  double m00[] = {1, 2, 0, 5, 0, 0, 0, 0,
                  2, 1, 0, 0, 0, 0, 0, 0,
                  0, 0, 1, -1, 0, 0, 2, 2,
                  4, 0, -1, 6, 0, 0, 1, 2
                 };
  M3->storageType = 0;
  M3->size0 = n;
  M3->size1 = nn;
  M3->matrix0 = (double *)malloc(nn * n * sizeof(double));
  for(int i = 0; i < nn * n; i++)
    M3->matrix0[i] = m00[i];

  return M3;
}
NumericsMatrix * test_matrix_2()
{
  NumericsMatrix * M2 = NM_new();
  int n = 8;

  /* Build a NumericsMatrix with sparse-block storage */
  M2->storageType = 1;
  M2->size0 = n;
  M2->size1 = n;

  SparseBlockStructuredMatrix * SBM = SBM_new();
  M2->matrix1 = SBM;
  SBM->nbblocks = 6;
  SBM->blocknumber0 = 3;
  SBM->blocknumber1 = SBM->blocknumber0 ;

  SBM->blocksize0 = (unsigned int*)malloc(3 * sizeof(unsigned int));
  SBM->blocksize0[0] = 4;
  SBM->blocksize0[1] = 6;
  SBM->blocksize0[2] = 8;
  SBM->blocksize1 = SBM->blocksize0;


  SBM->filled1 = 4;
  SBM->filled2 = SBM->nbblocks;

  SBM->index1_data = (size_t*)malloc((SBM->filled1) * sizeof(size_t));
  SBM->index1_data[0] = 0;
  SBM->index1_data[1] = 2;
  SBM->index1_data[2] = 4;
  SBM->index1_data[3] = 6;

  SBM->index2_data = (size_t*)malloc((SBM->filled2) * sizeof(size_t));
  SBM->index2_data[0] =  0;
  SBM->index2_data[1] =  1;
  SBM->index2_data[2] =  1;
  SBM->index2_data[3] =  2;
  SBM->index2_data[4] =  0;
  SBM->index2_data[5] =  2;

  SBM->block = (double **)malloc(SBM->nbblocks * sizeof(*(SBM->block)));
  double block0[16] = {1, 2, 0, 5, 2, 1, 0, 0, 0, 0, 1, -1, 4, 0, -1, 6};
  double block1[8] = {3, 4, 0, 0, -1, 1, 0, 6};
  double block2[4] = {1, 0, 0, 2};
  double block3[4] = {0, 0, 5, 2};
  double block4[8] = {0, 0, 0, 0, 2, 2, 1, 2};
  double block5[4] = {2, -1, 2, 2};
  SBM->block[0] = (double *)malloc(16 * sizeof(double));
  SBM->block[1] = (double *)malloc(8 * sizeof(double));
  SBM->block[2] = (double *)malloc(4 * sizeof(double));
  SBM->block[3] = (double *)malloc(4 * sizeof(double));
  SBM->block[4] = (double *)malloc(8 * sizeof(double));
  SBM->block[5] = (double *)malloc(4 * sizeof(double));
  for(int i = 0; i < 16; i++)
    SBM->block[0][i] = block0[i];
  for(int i = 0; i < 8; i++)
    SBM->block[1][i] = block1[i];
  for(int i = 0; i < 4; i++)
    SBM->block[2][i] = block2[i];
  for(int i = 0; i < 4; i++)
    SBM->block[3][i] = block3[i];
  for(int i = 0; i < 8; i++)
    SBM->block[4][i] = block4[i];
  for(int i = 0; i < 4; i++)
    SBM->block[5][i] = block5[i];
  /* SBM_print(SBM); */

  return M2;

}



NumericsMatrix * test_matrix_4()
{
  NumericsMatrix * M4 = NM_new();
  int n = 8;

  /* Build a NumericsMatrix with sparse-block storage */
  M4->storageType = 1;
  M4->size0 = n;
  M4->size1 = 4;

  M4->matrix1 = SBM_new();
  SparseBlockStructuredMatrix * SBM4 = M4->matrix1;

  SBM4->nbblocks = 2;
  SBM4->blocknumber0 = 3;

  SBM4->blocksize0 = (unsigned int*)malloc(SBM4->blocknumber0 * sizeof(unsigned int));
  SBM4->blocksize0[0] = 4;
  SBM4->blocksize0[1] = 6;
  SBM4->blocksize0[2] = 8;

  SBM4->blocknumber1 = 1;
  SBM4->blocksize1 = (unsigned int*)malloc(SBM4->blocknumber1 * sizeof(unsigned int));
  SBM4->blocksize1[0] = 4;

  SBM4->filled1 = 4;
  SBM4->filled2 = SBM4->nbblocks;

  SBM4->index1_data = (size_t*)malloc((SBM4->filled1) * sizeof(size_t));
  SBM4->index1_data[0] = 0;
  SBM4->index1_data[1] = 1;
  SBM4->index1_data[2] = 1;
  SBM4->index1_data[3] = 2;

  SBM4->index2_data = (size_t*)malloc((SBM4->filled2) * sizeof(size_t));
  SBM4->index2_data[0] =  0;
  SBM4->index2_data[1] =  0;

  SBM4->block = (double **)malloc(SBM4->nbblocks * sizeof(*(SBM4->block)));
  double block00[16] = {1, 2, 0, 5, 2, 1, 0, 0, 0, 0, 1, -1, 4, 0, -1, 6};
  double block40[8] = {0, 0, 0, 0, 2, 2, 1, 2};
  SBM4->block[0] = (double *)malloc(16 * sizeof(double));
  SBM4->block[1] = (double *)malloc(8 * sizeof(double));
  for(int i = 0; i < 16; i++)
    SBM4->block[0][i] = block00[i];
  for(int i = 0; i < 8; i++)
    SBM4->block[1][i] = block40[i];

  return M4;
}

NumericsMatrix * test_matrix_5()
{
  NumericsMatrix * M2 = test_matrix_2();

  NM_csc(M2);

  NM_clearDense(M2);
  NM_clearSparseBlock(M2);
  M2->storageType=NM_SPARSE;
  numericsSparseMatrix(M2)->origin = NSM_CSC;
  return M2;
}
NumericsMatrix * test_matrix_6()
{
  NumericsMatrix * M4 = test_matrix_4();

  NM_csc(M4);

  NM_clearDense(M4);
  NM_clearSparseBlock(M4);
  M4->storageType=NM_SPARSE;
  numericsSparseMatrix(M4)->origin = NSM_CSC;

  return M4;
}

NumericsMatrix * test_matrix_9()
{
  NumericsMatrix * M = NM_new();
  int n = 8;
  /* Double * storage (column-major) */
  double m0[64] = {1, 2, 0, 5,   0, 0,   0, 0,
                   2, 1, 0, 0,   0, 0,   0, 0,
                   0, 0, 1, -1,  0, 0,   2, 2,
                   4, 0, -1, 6,  0, 0,   1, 2,

                   0, 0, 0, 0,   0, 0,   0, 0,
                   0, 0, 0, 0,   0, 0,   0, 0,

                   0, 0, 0, 0,   0, 0,   2, -1,
                   0, 0, 0, 0,   5, 2,   2, 2
                  };
  M->storageType = 0;
  M->size0 = n;
  M->size1 = n;
  M->matrix0 = (double *)malloc(n * n * sizeof(double));


  for(int i = 0; i < n * n; i++)
    M->matrix0[i] = m0[i];
  return M;
}

NumericsMatrix * test_matrix_10()
{
  NumericsMatrix * M = NM_new();
  int n = 8;

  /* Build a NumericsMatrix with sparse-block storage */
  M->storageType = 1;
  M->size0 = n;
  M->size1 = n;

  SparseBlockStructuredMatrix * SBM = SBM_new();
  M->matrix1 = SBM;
  SBM->nbblocks = 4;
  SBM->blocknumber0 = 3;
  SBM->blocknumber1 = SBM->blocknumber0 ;

  SBM->blocksize0 = (unsigned int*)malloc(3 * sizeof(unsigned int));
  SBM->blocksize0[0] = 4;
  SBM->blocksize0[1] = 6;
  SBM->blocksize0[2] = 8;
  SBM->blocksize1 = SBM->blocksize0;


  SBM->filled1 = 4;
  SBM->filled2 = SBM->nbblocks;

  SBM->index1_data = (size_t*)malloc((SBM->filled1) * sizeof(size_t));
  SBM->index1_data[0] = 0;
  SBM->index1_data[1] = 1;
  SBM->index1_data[2] = 2;
  SBM->index1_data[3] = 4;

  SBM->index2_data = (size_t*)malloc((SBM->filled2) * sizeof(size_t));
  SBM->index2_data[0] =  0;
  SBM->index2_data[1] =  2;
  SBM->index2_data[2] =  0;
  SBM->index2_data[3] =  2;

  SBM->block = (double **)malloc(SBM->nbblocks * sizeof(*(SBM->block)));
  double block0[16] = {1, 2, 0, 5, 2, 1, 0, 0, 0, 0, 1, -1, 4, 0, -1, 6};
  /* double block1[8] = {3, 4, 0, 0, -1, 1, 0, 6}; */
  /* double block2[4] = {1, 0, 0, 2}; */
  double block3[4] = {0, 0, 5, 2};
  double block4[8] = {0, 0, 0, 0, 2, 2, 1, 2};
  double block5[4] = {2, -1, 2, 2};

  SBM->block[0] = (double *)malloc(16 * sizeof(double));
  /* SBM->block[1] = (double *)malloc(8 * sizeof(double)); */
  /* SBM->block[2] = (double *)malloc(4 * sizeof(double)); */
  SBM->block[1] = (double *)malloc(4 * sizeof(double));
  SBM->block[2] = (double *)malloc(8 * sizeof(double));
  SBM->block[3] = (double *)malloc(4 * sizeof(double));
  for(int i = 0; i < 16; i++)
    SBM->block[0][i] = block0[i];
  /* for (int i = 0; i < 8; i++) */
  /*   SBM->block[1][i] = block1[i]; */
  /* for (int i = 0; i < 4; i++) */
  /*   SBM->block[2][i] = block2[i]; */
  for(int i = 0; i < 4; i++)
    SBM->block[1][i] = block3[i];
  for(int i = 0; i < 8; i++)
    SBM->block[2][i] = block4[i];
  for(int i = 0; i < 4; i++)
    SBM->block[3][i] = block5[i];
  /* SBM_print(SBM); */

  return M;
}
NumericsMatrix * test_matrix_20()
{
  NumericsMatrix * M = NM_new();
  int n = 8;

  /* Build a NumericsMatrix with sparse-block storage */
  M->storageType = 1;
  M->size0 = n;
  M->size1 = n;

  SparseBlockStructuredMatrix * SBM = SBM_new();
  M->matrix1 = SBM;
  SBM->nbblocks = 9;
  SBM->blocknumber0 = 3;
  SBM->blocknumber1 = SBM->blocknumber0 ;

  SBM->blocksize0 = (unsigned int*)malloc(3 * sizeof(unsigned int));
  SBM->blocksize0[0] = 4;
  SBM->blocksize0[1] = 6;
  SBM->blocksize0[2] = 8;
  SBM->blocksize1 = SBM->blocksize0;

  SBM->filled1 = 4;
  SBM->filled2 = SBM->nbblocks;

  SBM->index1_data = (size_t*)malloc((SBM->filled1) * sizeof(size_t));
  SBM->index1_data[0] = 0;
  SBM->index1_data[1] = 3;
  SBM->index1_data[2] = 6;
  SBM->index1_data[3] = 9;

  SBM->index2_data = (size_t*)malloc((SBM->filled2) * sizeof(size_t));
  SBM->index2_data[0] =  0;
  SBM->index2_data[1] =  1;
  SBM->index2_data[2] =  2;
  SBM->index2_data[3] =  0;
  SBM->index2_data[4] =  1;
  SBM->index2_data[5] =  2;
  SBM->index2_data[6] =  0;
  SBM->index2_data[7] =  1;
  SBM->index2_data[8] =  2;

  SBM->block = (double **)malloc(SBM->nbblocks * sizeof(*(SBM->block)));


  for(size_t currentRowNumber = 0 ; currentRowNumber < SBM->filled1 - 1; ++currentRowNumber)
  {
    for(size_t blockNum = SBM->index1_data[currentRowNumber];
        blockNum < SBM->index1_data[currentRowNumber + 1]; ++blockNum)
    {
      unsigned int blocksize0 = SBM->blocksize0[currentRowNumber];
      if(currentRowNumber != 0)
        blocksize0  -= SBM->blocksize0[currentRowNumber - 1];

      unsigned int colNumber = SBM->index2_data[blockNum];
      unsigned int blocksize1 = SBM->blocksize1[colNumber];
      if(colNumber != 0)
        blocksize1 -= SBM->blocksize1[colNumber - 1];

      SBM->block[blockNum] = (double *)calloc(blocksize0*blocksize1,sizeof(double));
    }
  }


  return M;
}




int SBM_dense_equal(SparseBlockStructuredMatrix * M, double * mat, double tol)
{
  int info =0;

  int n = M->blocksize0[M->blocknumber0-1];
  int m = M->blocksize1[M->blocknumber1-1];

  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < m; j++)
    {
      if(fabs(mat[i + j * n] - SBM_get_value(M, i, j)) > tol)
      {
        printf("Matrices are not equal on element %i %i ", i, j);
        printf(" with values %e, %e\n", mat[i + j * n], SBM_get_value(M, i, j));
        info = 1;
      }
      if(info == 1) break;
    }
    if(info == 1) break ;
  }
  return info;
}
int NM_dense_equal(NumericsMatrix * M, double * mat, double tol)
{
  int info =0;

  int n = M->size0;
  int m = M->size1;

  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < m; j++)
    {
      if(fabs(mat[i + j * n] - NM_get_value(M, i, j)) > tol)
      {
        printf("Matrices are not equal on element %i %i ", i, j);
        printf(" with values %e, %e\n", mat[i + j * n], NM_get_value(M, i, j));
        info = 1;
      }
      if(info == 1) break;
    }
    if(info == 1) break ;
  }
  return info;
}
int test_build_first_4_NM(NumericsMatrix** MM)
{

  MM[0] = test_matrix_1();
  MM[1] = test_matrix_2();
  MM[2] = test_matrix_3();
  MM[3] = test_matrix_4();

  NumericsMatrix* M1 =  MM[0];
  NumericsMatrix* M2 =  MM[1];
  NumericsMatrix* M3 =  MM[2];
  NumericsMatrix* M4 =  MM[3];

  int info = 0;
  /* Build two equal Numerics Matrices, one with double* storage (MM1), the other with sparse storage (MM2)*/

  /* M1 and M2 must have the same values.*/
  double tol = 1e-12;
  info = SBM_dense_equal(M2->matrix1, M1->matrix0, tol);

  if(info)
    return info;

  /*   SBM_print(SBM2); */
  /* M3 and M4 must have the same values.*/
  info = SBM_dense_equal(M4->matrix1, M3->matrix0, tol);

  return info;
}
/* ============================================================================================================================== */


int test_SBM_row_to_dense(SparseBlockStructuredMatrix *M)
{
  double * denseRes = (double*) malloc(M->blocksize0[M->blocknumber0 - 1] * M->blocksize1[M->blocknumber1 - 1] * sizeof(double));
  unsigned int curRow = 0;
  unsigned int nbCol = M->blocksize1[M->blocknumber1 - 1];
  for(unsigned int i = 0; i < M->blocknumber0; i++)
  {
    unsigned int lLin = 0;
    unsigned int nbBlockRow = M->blocksize0[i] - curRow;
    SBM_row_to_dense(M, i, denseRes, 0, nbBlockRow);
    for(unsigned int lin = curRow; lin < M->blocksize0[i]; lin++)
    {
      unsigned int lCol = 0;
      for(unsigned int col = 0; col < nbCol; col++)
      {
        if(fabs(SBM_get_value(M, lin, col) - denseRes[lLin + lCol * (nbBlockRow)]) > 10e-12)
        {
          free(denseRes);
          return 1;
        }
        lCol++;
      }
      lLin++;
    }
    curRow = M->blocksize0[i];
  }
  curRow = 0;
  for(unsigned int i = 0; i < M->blocknumber0; i++)
  {

    //    int lLin=0;
    //    int nbBlockRow=M->blocksize0[i]-curRow;
    SBM_row_to_dense(M, i, denseRes, curRow, M->blocksize0[M->blocknumber0 - 1]);
    curRow = M->blocksize0[i];
  }

  double * denseRes2 = (double*) malloc(M->blocksize0[M->blocknumber0 - 1] * M->blocksize1[M->blocknumber1 - 1] * sizeof(double));

  SBM_to_dense(M, denseRes2);
  for(unsigned int n = 0; n < M->blocksize0[M->blocknumber0 - 1]*M->blocksize1[M->blocknumber1 - 1]; n++)
    if(fabs(denseRes2[n] - denseRes[n]) > 10e-12)
    {
      free(denseRes);
      free(denseRes2);
      return 1;
    }

  free(denseRes);
  free(denseRes2);
  return 0;
}
int test_SBM_row_permutation(SparseBlockStructuredMatrix *M)
{
  SparseBlockStructuredMatrix MRes;
  unsigned int nbRow = M->blocknumber0;
  unsigned int * rowBlockIndex = (unsigned int*) malloc(nbRow * sizeof(unsigned int));
  unsigned int * mark = (unsigned int*) malloc(nbRow * sizeof(unsigned int));
  for(unsigned int i = 0; i < nbRow; i++)
  {
    mark[i] = 0;
  }
  for(unsigned int i = 0; i < nbRow; i++)
  {
    int candidate = rand() % nbRow;
    while(mark[candidate])
      candidate = rand() % nbRow;
    rowBlockIndex[i] = candidate;
    mark[candidate] = 1;
  }
  SBM_row_permutation(rowBlockIndex, M, &MRes);
  double * denseMRes = (double*) malloc(M->blocksize0[M->blocknumber0 - 1] * M->blocksize1[M->blocknumber1 - 1] * sizeof(double));
  SBM_to_dense(&MRes, denseMRes);
  double * denseM = (double*) malloc(M->blocksize0[M->blocknumber0 - 1] * M->blocksize1[M->blocknumber1 - 1] * sizeof(double));
  unsigned int curRow = 0;
  unsigned int nbRowInM = M->blocksize0[M->blocknumber0 - 1];
  for(unsigned int i = 0; i < nbRow; i++)
  {
    unsigned int rowInM = rowBlockIndex[i];
    unsigned int nbRow = 0;
    if(rowInM)
      nbRow = M->blocksize0[rowInM] - M->blocksize0[rowInM - 1];
    else
      nbRow = M->blocksize0[rowInM];
    SBM_row_to_dense(M, rowInM, denseM, curRow, nbRowInM);
    curRow += nbRow;
  }
  for(unsigned int n = 0; n < M->blocksize0[M->blocknumber0 - 1]*M->blocksize1[M->blocknumber1 - 1]; n++)
    if(fabs(denseMRes[n] - denseM[n]) > 10e-12)
    {
      free(denseM);
      free(denseMRes);
      free(rowBlockIndex);
      free(mark);
      SBMfree(&MRes, 0);
      return 1;
    }
  free(denseM);
  free(denseMRes);
  free(rowBlockIndex);
  free(mark);
  SBMfree(&MRes, 0);
  return 0;
}
int test_SBM_column_permutation(SparseBlockStructuredMatrix *M)
{
  SparseBlockStructuredMatrix MRes;
  SBM_null(&MRes);
  int nbCol = M->blocknumber1;
  unsigned int * colBlockIndex = (unsigned int*) malloc(nbCol * sizeof(unsigned int));
  int * mark = (int*) malloc(nbCol * sizeof(int));
  for(int i = 0; i < nbCol; i++)
    mark[i] = 0;
  for(int i = 0; i < nbCol; i++)
  {
    int candidate = rand() % nbCol;
    while(mark[candidate])
      candidate = rand() % nbCol;
    colBlockIndex[i] = candidate;
    mark[candidate] = 1;
  }
  SBM_column_permutation(colBlockIndex, M, &MRes);
  free(colBlockIndex);
  free(mark);
  SBMfree(&MRes, 0);
  return 0;
}
