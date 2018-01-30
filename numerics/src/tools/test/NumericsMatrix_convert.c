#include "NumericsMatrix.h"
#include "NumericsVector.h"
#include "NumericsSparseMatrix.h"
#include "CSparseMatrix.h"
int triplet_to_dense(void);
int csc_to_dense(void);

int triplet_to_dense(void)
{
  int info =1;
  NumericsMatrix *A;
  char * filename =  "./data/NSM_triplet_162x162.dat";
  A = NM_create_from_filename(filename);
  /* NM_display(A); */

  NumericsMatrix *B = NM_create(NM_DENSE,A->size0,A->size1);

  info =  NM_to_dense(A, B);

  /* NM_display(B); */

  return (int) NM_equal(A,B) -1 ;
}

int csc_to_dense(void)
{
  int info =1;
  NumericsMatrix *A;
  /* char * filename =  "./data/NSM_triplet_162x162.dat"; */
  /* A = NM_create_from_filename(filename);   */
  /* /\* NM_display(A); *\/ */
  /* CSparseMatrix* A_csc = NM_csc(A); */
  /* NumericsSparseMatrix * C_NSM =  numericsSparseMatrix_new(); */
  /* C_NSM->csc = A_csc ; */
  /* C_NSM->origin= NS_CSC; */
  /* NumericsMatrix*  A_CSC  = NM_create_from_data(NM_SPARSE, A->size0, A->size1, (void*)C_NSM); */
  /* NM_display(A_CSC); */
  /* NM_write_in_filename(A_CSC, "./data/NSM_csc_162x162.dat"); */


  char * filename =  "./data/NSM_csc_162x162.dat";
  A = NM_new_from_filename(filename);
  /* NM_display(A); */


  NumericsMatrix *B = NM_create(NM_DENSE,A->size0,A->size1);
  info =  NM_to_dense(A, B);

  /* NM_display(B);   */

  return (int) NM_equal(A,B) -1 ; 

}

int main()
{

  int info = triplet_to_dense();
  printf("triplet_to_dense() :  info = %i\n", info);
  info += csc_to_dense();
  printf("csc_to_dense() :  info = %i\n", info);

  return info;
}
