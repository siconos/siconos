#include <stdio.h>                 // for printf, fclose, fopen, FILE, NULL
/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include <stdlib.h>                // for malloc
#include "CSparseMatrix_internal.h"         // for cs_dl_entry, CS_INT, cs_dl_print
#include "NumericsFwd.h"           // for NumericsMatrix, NumericsSparseMatrix
#include "NumericsMatrix.h"        // for NM_entry, NM_display, NM_create
#include "NumericsSparseMatrix.h"  // for NumericsSparseMatrix, NSM_TRIPLET
#include "NumericsVector.h"        // for NV_display

#ifdef SICONOS_HAS_MPI
#include <mpi.h>
#endif
int add_square_triplet(void);
int add_square_csc(void);
int add_square_triplet_into_csc(void);
int add_rectangle_triplet(void);

int add_square_triplet()
{


  int size0 =3;
  int size1 =3;

  // product of triplet matrices into triplet matrix
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;
  NM_entry(A, 0, 0, 1);
  NM_entry(A, 0, 1, 2);
  NM_entry(A, 0, 2, 3);
  NM_entry(A, 1, 1, 2);
  NM_entry(A, 1, 2, 3);
  NM_display(A);


  NumericsMatrix * B  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(B,0);
  B->matrix2->origin= NSM_TRIPLET;
  NM_entry(B, 0, 0, 1);
  NM_entry(B, 1, 1, 2);
  NM_entry(B, 2, 2, 3);
  NM_display(B);

  double alpha = 2.0;
  double beta =2.0;
  NumericsMatrix * C = NM_add(alpha, A, beta, B);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(Cref,0);
  Cref->matrix2->origin= NSM_TRIPLET;
  NM_entry(Cref, 0, 0, 4);
  NM_entry(Cref, 0, 1, 4);
  NM_entry(Cref, 0, 2, 6);
  NM_entry(Cref, 1, 1, 8);
  NM_entry(Cref, 1, 2, 6);
  NM_entry(Cref, 2, 2, 6);
  NM_display(Cref);

  printf("add_square_triplet: NM_equal(C,Cref) =%i \n", NM_equal(C,Cref));
  int info = (int)!NM_equal(C,Cref);
  
  NM_clear(A);
  NM_clear(B);
  NM_clear(C);
  NM_clear(Cref);

  return info;

}

int add_square_csc()
{


  int size0 =3;
  int size1 =3;

  // product of csc matrices into csc matrix
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(A,0);
  A->matrix2->origin= NSM_CSC;
  NM_entry(A, 0, 0, 1);
  NM_entry(A, 0, 1, 2);
  NM_entry(A, 0, 2, 3);
  NM_entry(A, 1, 1, 2);
  NM_entry(A, 1, 2, 3);
  /* NM_display(A); */


  NumericsMatrix * B  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(B,0);
  B->matrix2->origin= NSM_CSC;
  NM_entry(B, 0, 0, 1);
  NM_entry(B, 1, 1, 2);
  NM_entry(B, 2, 2, 3);
  /* NM_display(B); */

  double alpha = 2.0;
  double beta =  2.0;
  NumericsMatrix * C = NM_add(alpha, A, beta, B);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(Cref,0);
  Cref->matrix2->origin= NSM_CSC;
  NM_entry(Cref, 0, 0, 4);
  NM_entry(Cref, 0, 1, 4);
  NM_entry(Cref, 0, 2, 6);
  NM_entry(Cref, 1, 1, 8);
  NM_entry(Cref, 1, 2, 6);
  NM_entry(Cref, 2, 2, 6);
  NM_display(Cref);
  printf("add_square_csc: NM_equal(C,Cref) =%i \n", NM_equal(C,Cref));
  int info = (int)!NM_equal(C,Cref);
  
  NM_clear(A);
  NM_clear(B);
  NM_clear(C);
  NM_clear(Cref);

  return info;

}

int add_rectangle_triplet()
{


  int size0 =3;
  int size1 =9;

  // product of triplet matrices into triplet matrix
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;
  NM_entry(A, 0, 0, 1);
  NM_entry(A, 0, 1, 2);
  NM_entry(A, 0, 2, 3);
  NM_entry(A, 1, 1, 2);
  NM_entry(A, 1, 2, 3);
  NM_entry(A, 2, 6, 2);
  NM_entry(A, 2, 5, 22);
  NM_display(A);


  NumericsMatrix * B  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(B,0);
  B->matrix2->origin= NSM_TRIPLET;
  NM_entry(B, 0, 0, 1);
  NM_entry(B, 1, 1, 2);
  NM_entry(B, 2, 2, 3);
  NM_entry(B, 0, 3, 1);
  NM_entry(B, 1, 4, 2);
  NM_entry(B, 2, 5, 3);
  NM_display(B);

  double alpha = 1.0;
  double beta  = 2.0;
  NumericsMatrix * C = NM_add(alpha, A, beta, B);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(Cref,0);
  Cref->matrix2->origin= NSM_TRIPLET;
  NM_entry(Cref, 0, 0, 3);
  NM_entry(Cref, 0, 1, 2);
  NM_entry(Cref, 0, 2, 3);
  NM_entry(Cref, 0, 3, 2);

  NM_entry(Cref, 1, 1, 6);
  NM_entry(Cref, 1, 2, 3);
  NM_entry(Cref, 1, 4, 4);

  NM_entry(Cref, 2, 2, 6);


  NM_entry(Cref, 2, 5, 28);
  NM_entry(Cref, 2, 6, 2);
  NM_display(Cref);

  printf("add_rectangle_triplet : NM_equal(C,Cref) =%i \n", NM_equal(C,Cref));
  int info = (int)!NM_equal(C,Cref);
  
  NM_clear(A);
  NM_clear(B);
  NM_clear(C);
  NM_clear(Cref);

  return info;



}

static int add_test(void)
{

  int info = add_square_triplet();
  info += add_square_csc();
  info +=  add_rectangle_triplet();

  return info;
}


/* create an empty triplet matrix, insert 2 elements, print and free */
static int test_CSparseMatrix_alloc(void)
{
  CSparseMatrix *m = cs_spalloc(0,0,0,0,1); /* coo format */

  CS_INT info1 = 1-cs_entry(m, 3, 4, 1.0);
  CS_INT info2 = 1-cs_entry(m, 1, 2, 2.0);

  CS_INT info3 = 1-cs_print(m, 0);

  m=cs_spfree(m);

  CS_INT info4 = 1-(m==NULL);

  return (int)(info1+info2+info3+info4);
}

static int test_CSparseMatrix_spsolve_unit(CSparseMatrix *M )
{
  //cs_print(M, 0);

  CSparseMatrix *b_triplet = cs_spalloc(M->m, M->n,M->n, 1, 1); /* coo format */
  for (int i =0; i < M->n; i++)
    cs_entry(b_triplet, i, i, 1.0);

  CSparseMatrix *B = cs_compress(b_triplet);

  CSparseMatrix *X = cs_spalloc(M->m, M->n, M->nzmax, 1,0); /* csr format */

  CSparseMatrix_factors* cs_lu_M = (CSparseMatrix_factors*) malloc(sizeof(CSparseMatrix_factors));

  int info = 1-CSparseMatrix_lu_factorization(1, M, 1e-14, cs_lu_M);


  if (info)
  {
    printf("problem in Lu factor\n");
    return info;
  }
  /* printf(" L:\n"); */
  /* cs_print(cs_lu_M->N->L, 0); */
  /* printf(" U:\n"); */
  /* cs_print(cs_lu_M->N->U, 0); */


  info = 1-CSparseMatrix_spsolve(cs_lu_M, X, B);
  if (info)
  {
    printf("problem in spsolve\n");
    return info;
  }
  CSparseMatrix* I = cs_multiply(M, B);
  //printf(" M * M^-1:\n");
  //cs_print(I, 0);

  CSparseMatrix *Id = cs_compress(b_triplet);

  CSparseMatrix* check = cs_add(I, Id, 1.0, -1.0);
  //cs_print(check, 0);

  double error = cs_norm(check);
  printf("residual =%e\n", error);

  cs_spfree(b_triplet);
  cs_spfree(B);
  cs_spfree(X);
  CSparseMatrix_free_lu_factors(cs_lu_M);
  cs_spfree(I);
  cs_spfree(Id);
  cs_spfree(check);

  
  if (error > 1e-12)
  {
    return 1;

  }

  return  info;
}


static int test_CSparseMatrix_spsolve(void)
{
  printf("start - test_CSparseMatrix_spsolve\n");
  CSparseMatrix *m_triplet = cs_spalloc(3,3,3,1,1); /* coo format */
  cs_entry(m_triplet, 0, 0, 1.0);
  cs_entry(m_triplet, 1, 1, 2.0);
  cs_entry(m_triplet, 2, 2, 4.0);
//  CS_INT info4 = 1-cs_print(m_triplet, 0);
  CSparseMatrix *M = cs_compress(m_triplet);

  int info =  test_CSparseMatrix_spsolve_unit(M);

  cs_entry(m_triplet, 0, 1, 3.0);
  cs_entry(m_triplet, 0, 2, 6.0);
  cs_entry(m_triplet, 1, 2, 5.0);

  CSparseMatrix *M1 = cs_compress(m_triplet);
  info +=  test_CSparseMatrix_spsolve_unit(M1);

  cs_entry(m_triplet, 1, 0, 7.0);
  cs_entry(m_triplet, 2, 0, 8.0);
  cs_entry(m_triplet, 2, 1, 9.0);

  CSparseMatrix *M2 = cs_compress(m_triplet);
  info +=  test_CSparseMatrix_spsolve_unit(M2);
  cs_spfree(M);
  cs_spfree(M1);
  cs_spfree(M2);
  cs_spfree(m_triplet);
  

  int size0 =10;
  int size1 =10;
  CSparseMatrix *a_triplet = cs_spalloc(size0,size1,size0,1,1); /* coo format */
  for(int i =0; i < size0; i++)
  {
    for(int j =i; j < size1; j++)
    {
      cs_entry(a_triplet, i, j, i+j+1);
    }
  }
  CSparseMatrix *A = cs_compress(a_triplet);
  info +=  test_CSparseMatrix_spsolve_unit(A);
  printf("end - test_CSparseMatrix_spsolve\n");
  cs_spfree(A);
  cs_spfree(a_triplet);
  return  info;
}
static int test_CSparseMatrix_chol_spsolve_unit(CSparseMatrix *M )
{
  //cs_print(M, 0);


  CSparseMatrix *b_triplet = cs_spalloc(M->m, M->n,M->n, 1, 1); /* coo format */
  for (int i =0; i < M->n; i++)
    cs_entry(b_triplet, i, i, 1.0);


  CSparseMatrix *B = cs_compress(b_triplet);

  CSparseMatrix *X = cs_spalloc(M->m, M->n, M->nzmax, 1,0); /* csc format */

  CSparseMatrix_factors* cs_chol_M = (CSparseMatrix_factors*) malloc(sizeof(CSparseMatrix_factors));

  int info = 1-CSparseMatrix_chol_factorization(1, M, cs_chol_M);

  if (info)
  {
    printf("problem in Cholesky factor\n");
    return info;
  }
  /* printf(" L:\n"); */
  /* cs_print(cs_chol_M->N->L, 0); */

  info = 1-CSparseMatrix_chol_spsolve(cs_chol_M, X, B);
  if (info)
  {
    printf("problem in chol_spsolve\n");
    return info;
  }
  CSparseMatrix* I = cs_multiply(M, B);
  //printf(" M * M^-1:\n");
  //cs_print(I, 0);

  CSparseMatrix *Id = cs_compress(b_triplet);

  CSparseMatrix* check = cs_add(I, Id, 1.0, -1.0);
  //cs_print(check, 0);

  double error = cs_norm(check);
  printf("residual =%e\n", error);
  cs_spfree(b_triplet);
  cs_spfree(B);
  cs_spfree(X);
  CSparseMatrix_free_lu_factors(cs_chol_M);
  cs_spfree(I);
  cs_spfree(Id);
  cs_spfree(check);

  if (error > 1e-12)
  {
    return 1;

  }
  //printf("info =%i\n", info);
  return  info;
}


static int test_CSparseMatrix_chol_spsolve(void)
{
  printf("start - test_CSparseMatrix_chol_spsolve \n");
  CSparseMatrix *m_triplet = cs_spalloc(3,3,3,1,1); /* coo format */
  cs_entry(m_triplet, 0, 0, 1.0);
  cs_entry(m_triplet, 1, 1, 2.0);
  cs_entry(m_triplet, 2, 2, 4.0);
//  CS_INT info4 = 1-cs_print(m_triplet, 0);
  CSparseMatrix *M = cs_compress(m_triplet);

  int info =  test_CSparseMatrix_chol_spsolve_unit(M);
  if (info) return 1;

  cs_entry(m_triplet, 0, 1, 3.0);
  cs_entry(m_triplet, 0, 2, 6.0);
  cs_entry(m_triplet, 1, 2, 5.0);
  cs_entry(m_triplet, 1, 0, 3.0);
  cs_entry(m_triplet, 2, 0, 6.0);
  cs_entry(m_triplet, 2, 1, 5.0);
  CSparseMatrix *M2 = cs_compress(m_triplet);
  CSparseMatrix *M2T = cs_transpose(M2,1);
  CSparseMatrix *M2M2T = cs_multiply(M2,M2T);
  info =  test_CSparseMatrix_chol_spsolve_unit(M2M2T);
  cs_spfree(M);
  cs_spfree(M2);
  cs_spfree(M2T);
  cs_spfree(M2M2T);
  cs_spfree(m_triplet);
  if (info) return 1;

  int size0 =10;
  int size1 =10;
  CSparseMatrix *a_triplet = cs_spalloc(size0,size1,size0,1,1); /* coo format */
  for(int i =0; i < size0; i++)
  {
    for(int j =i; j < size1; j++)
    {
      cs_entry(a_triplet, i, j, i+j+1);
    }
  }
  CSparseMatrix *A = cs_compress(a_triplet);
  CSparseMatrix *AT = cs_transpose(A,1);
  CSparseMatrix *AAT = cs_multiply(A,AT);

  info =  test_CSparseMatrix_chol_spsolve_unit(AAT);
  cs_spfree(a_triplet);
  cs_spfree(A);
  cs_spfree(AT);
  cs_spfree(AAT);
  printf("end - test_CSparseMatrix_chol_spsolve \n");
  return  info;
}
static int test_CSparseMatrix_ldlt_solve_unit(CSparseMatrix *M )
{
  //cs_print(M, 0);

  double * b = (double * )calloc(M->n, sizeof(double));
  for (int k = 0 ; k < M->n; k++)
  {
    b[k] = 3*k+3;
  }
  double * b_backup = (double * )calloc(M->n, sizeof(double));
  for (int k = 0 ; k < M->n; k++)
  {
    b_backup[k] = b[k];
  }
  double * x = (double * )calloc(M->n, sizeof(double));

  CSparseMatrix_factors* cs_ldlt_M = (CSparseMatrix_factors*) malloc(sizeof(CSparseMatrix_factors));

  int info = 1-CSparseMatrix_ldlt_factorization(1, M, cs_ldlt_M);

  if (info)
  {
    printf("problem in LDLT factor\n");
    return info;
  }
  /* printf(" L:\n"); */
  /* cs_print(cs_ldlt_M->N->L, 0); */

  info = 1-CSparseMatrix_ldlt_solve(cs_ldlt_M, x, b);
  if (info)
  {
    printf("problem in ldlt_spsolve\n");
    return info;
  }


  /* compute residual */
  CSparseMatrix_aaxpby(1.0, M, b, -1.0, b_backup);
  
  double error =0.0;
  for (int k =0; k < M->n ; k++)
  {
    error += b_backup[k]*b_backup[k];
  }
  error = sqrt(error);
  
  free(b);
  free(x);
  free(b_backup);

  printf("residual =%e\n", error);
  if (error > 1e-12)
  {
    return 1;

  }
  printf("info =%i\n", info);
  return  info;
}


static int test_CSparseMatrix_ldlt_solve(void)
{
  printf("\nstart - test_CSparseMatrix_ldlt_solve \n");
  CSparseMatrix *m_triplet = cs_spalloc(3,3,3,1,1); /* coo format */
  cs_entry(m_triplet, 0, 0, 1.0);
  cs_entry(m_triplet, 1, 1, 2.0);
  cs_entry(m_triplet, 2, 2, 4.0);
//  CS_INT info4 = 1-cs_print(m_triplet, 0);
  CSparseMatrix *M = cs_compress(m_triplet);
  printf("test 1 ....");
  int info =  test_CSparseMatrix_ldlt_solve_unit(M);
  if (info) return 1;
  printf("ok\n");

  cs_entry(m_triplet, 0, 1, 3.0);
  cs_entry(m_triplet, 0, 2, 6.0);
  cs_entry(m_triplet, 1, 2, 5.0);
  cs_entry(m_triplet, 1, 0, 3.0);
  cs_entry(m_triplet, 2, 0, 6.0);
  cs_entry(m_triplet, 2, 1, 5.0);
  CSparseMatrix *M2 = cs_compress(m_triplet);
  CSparseMatrix *M2T = cs_transpose(M2,1);
  CSparseMatrix *M2M2T = cs_multiply(M2,M2T);
  printf("test 2 ....");
  info =  test_CSparseMatrix_ldlt_solve_unit(M2M2T);
  if (info) return 1;
  printf("ok\n");
  int size0 =10;
  int size1 =10;
  CSparseMatrix *a_triplet = cs_spalloc(size0,size1,size0,1,1); /* coo format */
  for(int i =0; i < size0; i++)
  {
    for(int j =i; j < size1; j++)
    {
      cs_entry(a_triplet, i, j, i+j+1);
    }
  }
  CSparseMatrix *A = cs_compress(a_triplet);
  CSparseMatrix *AT = cs_transpose(A,1);
  CSparseMatrix *AAT = cs_multiply(A,AT);
  printf("test 3 ....");
  info =  test_CSparseMatrix_ldlt_solve_unit(AAT);
  printf("ok\n");
  printf("end - test_CSparseMatrix_ldlt_solve \n");
  return  info;
}


int main()
{
#ifdef SICONOS_HAS_MPI
  MPI_Init(NULL, NULL);
#endif
  int info = add_test();

  info += test_CSparseMatrix_alloc();
  info += test_CSparseMatrix_spsolve();
  printf("info : %i\n", info);
  info += test_CSparseMatrix_chol_spsolve();
  info += test_CSparseMatrix_ldlt_solve();
  printf("info : %i\n", info);

#ifdef SICONOS_HAS_MPI
  MPI_Finalize();
#endif
  return info;
}
