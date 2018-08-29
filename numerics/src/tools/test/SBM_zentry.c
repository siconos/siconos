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

/*
  Tests functions for NumericsMatrix structure

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NumericsMatrix.h"
#include <math.h>
#include "numericsMatrixTestFunction.h"
#include "SparseBlockMatrix.h"
#include "SiconosLapack.h"
#include "numerics_verbose.h"
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"

static void add_initial_value_square_1(NumericsMatrix * M)
{
  int i=0, j=0;
  for (i=0; i < 4 ; i++)
  {
    for (j=0; j < 4 ; j++)
      NM_zentry(M,i,j,1.0);
  }
  for (i=0; i < 4 ; i++)
  {
    for (j=4; j < 6 ; j++)
      NM_zentry(M,i,j,2.0);
  }
  for (i=4; i < 6 ; i++)
  {
    for (j=4; j < 6 ; j++)
      NM_zentry(M,i,j,3.0);
  }
  for (i=4; i < 6 ; i++)
  {
    for (j=6; j < 8 ; j++)
      NM_zentry(M,i,j,4.0);
  }
  for (i=6; i < 8 ; i++)
  {
    for (j=0; j < 4 ; j++)
      NM_zentry(M,i,j,5.0);
  }
  for (i=6; i < 8 ; i++)
  {
    for (j=6; j < 8 ; j++)
      NM_zentry(M,i,j,6.0);
  }
}
static void add_initial_value_square_SBM_1(SparseBlockStructuredMatrix * M)
{
  int i=0, j=0;
  for (i=0; i < 4 ; i++)
  {
    for (j=0; j < 4 ; j++)
      CHECK_RETURN(SBM_zentry(M,i,j,1.0));
  }
  for (i=0; i < 4 ; i++)
  {
    for (j=4; j < 6 ; j++)
      CHECK_RETURN(SBM_zentry(M,i,j,2.0));
  }
  for (i=4; i < 6 ; i++)
  {
    for (j=4; j < 6 ; j++)
      SBM_zentry(M,i,j,3.0);
  }
  for (i=4; i < 6 ; i++)
  {
    for (j=6; j < 8 ; j++)
      SBM_zentry(M,i,j,4.0);
  }
  for (i=6; i < 8 ; i++)
  {
    for (j=0; j < 4 ; j++)
      SBM_zentry(M,i,j,5.0);
  }
  for (i=6; i < 8 ; i++)
  {
    for (j=6; j < 8 ; j++)
      SBM_zentry(M,i,j,6.0);
  }
}


static int SBM_zentry_test1(double tol)
{

  int info=0;

  NumericsMatrix * M2 = test_matrix_2();

  NumericsMatrix * C= NM_create(NM_DENSE, M2->size0, M2->size1);
  add_initial_value_square_1(C);


  NumericsMatrix * C2= NM_create(NM_SPARSE_BLOCK, M2->size0, M2->size1);
  C2->matrix1 = SBM_zero_matrix_for_multiply(M2->matrix1, M2->matrix1);
  add_initial_value_square_SBM_1(C2->matrix1);
  
  DEBUG_EXPR(NM_display(C));
  DEBUG_EXPR(NM_display(C2));
  info = NM_dense_equal(C2,C->matrix0,tol);
  NM_free(M2);
  NM_free(C2);
  NM_free(C);
  return info;
}

static int SBM_zentry_test2(double tol)
{

  int info=0;

  NumericsMatrix * M2 = test_matrix_2();
  
  CHECK_RETURN(SBM_zentry(M2->matrix1,0,8,1.0));
  CHECK_RETURN(SBM_zentry(M2->matrix1,8,0,1.0));
  
  CHECK_RETURN(SBM_zentry(M2->matrix1,0,7,1.0));

  NM_free(M2);

  return info;
}


int main(void)
{

  int info =0;
  double tol = 1e-14;

  printf("========= Starts SBM tests SBM_zentry  ========= \n");

  info = SBM_zentry_test1(tol);
  if (info == 1)
  {
    printf("========= Ends SBM tests SBM_zentry  :  Unsuccessfull ========= \n");
    return info;
  }
  info = SBM_zentry_test2(tol);
  if (info == 1)
  {
    printf("========= Ends SBM tests SBM_zentry  :  Unsuccessfull ========= \n");
    return info;
  }
  printf("========= Ends SBM tests SBM_zentry  :  successfull ========= \n");

  return info;

}
