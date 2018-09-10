#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NumericsMatrix.h"
#include <math.h>
#include "numericsMatrixTestFunction.h"
#include "SiconosLapack.h"
#include "sanitizer.h"
#include "SparseBlockMatrix.h"
#include "NumericsSparseMatrix.h"
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"



static int NM_add_to_diag3_test(NumericsMatrix* M, double alpha)
{
  printf("\n == Numerics tests: NM_add_to_diag3(...) == \n");
  printf("Starts NM_add_to_diag3_test for alpha = %e\n",alpha);
  
  int info =-1;

  /***********************************************************/
  /* C = C + alpha +I NM_SPARSE_BLOCK storage, square matrix  */
  /***********************************************************/
  int n = M->size0;
  int m = M->size1;
  NumericsMatrix * C1= NM_new();
  NM_copy(M, C1);

  NumericsMatrix * Cref = NM_create(NM_DENSE, n, m);
  SBM_to_dense(M->matrix1,Cref->matrix0);

  
  NM_add_to_diag3(C1, alpha);
  DEBUG_EXPR(NM_display(C1););

  
  double * Id = (double * ) calloc (n*m, sizeof(double));
  for (int i = 0; i < n; i++)
  {
    Id[i + i  *n  ] =1.0;
  }

  cblas_daxpy(n*m, alpha, Id, 1, Cref->matrix0, 1);
  DEBUG_EXPR(NM_display(Cref););
  info = NM_dense_equal(C1, Cref->matrix0, 1e-14 );

  if (info == 0)
    printf("Step 0 ( C = C + alpha*I,  NM_SPARSE_BLOCK storage, square matrix) ok ...\n");
  else
  {
    printf("Step 0 (C = C + alpha*I, NM_SPARSE_BLOCK storage, square matrix) failed ...\n");
    goto exit_0;
  }

  /***********************************************************/
  /* C = C + alpha +I NM_SPARSE storage, square matrix  */
  /***********************************************************/
  NumericsMatrix * C2= NM_create(NM_SPARSE, n, m);
  NM_copy_to_sparse(M,C2);
  NM_add_to_diag3(C2, alpha);

  info = NM_dense_equal(C2, Cref->matrix0, 1e-14 );

  if (info == 0)
    printf("Step 1 ( C = C + alpha*I,  NM_SPARSE storage, square matrix) ok ...\n");
  else
  {
    printf("Step 1 (C = C + alpha*I, NM_SPARSE storage, square matrix) failed ...\n");
    goto exit_1;
  }

exit_1:
  NM_free(C2);
exit_0:
  free(Id);
  NM_free(Cref);
  
  return info;
}


int main(void)
{

  printf("========= Starts Numerics tests for NM_add_to_diag3 ========= \n");

  SparseBlockStructuredMatrix * SBM = SBM_new();
  FILE *file = fopen("data/SBM1.dat", "r");
  SBM_new_from_file(SBM, file);
  fclose(file);

  NumericsMatrix * M = NM_create(NM_SPARSE_BLOCK, SBM->blocksize0[SBM->blocknumber0-1],SBM->blocksize1[SBM->blocknumber1-1] );
  M->matrix1=SBM;

  
  int info = NM_add_to_diag3_test(M,1.0);
  if (info != 0)
  {
    printf("End of  : Numerics tests for NM_add_to_diag3unsucessfull\n");
    return info;
  }

  printf("End of Numerics tests for NM_add_to_diag3 ...\n");
  if (info != 0) return info;
  /* free memory */


  NM_free(M);

  printf("========= End Numerics tests for NM_add_to_diag3 ========= \n");
  return info;
}
