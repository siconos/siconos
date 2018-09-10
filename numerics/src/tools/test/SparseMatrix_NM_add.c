#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "CSparseMatrix.h"
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
  NM_zentry(A, 0, 0, 1);
  NM_zentry(A, 0, 1, 2);
  NM_zentry(A, 0, 2, 3);
  NM_zentry(A, 1, 1, 2);
  NM_zentry(A, 1, 2, 3);
  NM_display(A);


  NumericsMatrix * B  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(B,0);
  B->matrix2->origin= NSM_TRIPLET;
  NM_zentry(B, 0, 0, 1);
  NM_zentry(B, 1, 1, 2);
  NM_zentry(B, 2, 2, 3);
  NM_display(B);

  double alpha = 2.0;
  double beta =2.0;
  NumericsMatrix * C = NM_add(alpha, A, beta, B);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(Cref,0);
  Cref->matrix2->origin= NSM_TRIPLET;
  NM_zentry(Cref, 0, 0, 4);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 6);
  NM_zentry(Cref, 1, 1, 8);
  NM_zentry(Cref, 1, 2, 6);
  NM_zentry(Cref, 2, 2, 6);
  NM_display(Cref);

  printf("add_square_triplet: NM_equal(C,Cref) =%i \n", NM_equal(C,Cref));
  return (int)!NM_equal(C,Cref);


}

int add_square_csc()
{


  int size0 =3;
  int size1 =3;

  // product of csc matrices into csc matrix
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(A,0);
  A->matrix2->origin= NSM_CSC;
  NM_zentry(A, 0, 0, 1);
  NM_zentry(A, 0, 1, 2);
  NM_zentry(A, 0, 2, 3);
  NM_zentry(A, 1, 1, 2);
  NM_zentry(A, 1, 2, 3);
   /* NM_display(A); */


  NumericsMatrix * B  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(B,0);
  B->matrix2->origin= NSM_CSC;
  NM_zentry(B, 0, 0, 1);
  NM_zentry(B, 1, 1, 2);
  NM_zentry(B, 2, 2, 3);
  /* NM_display(B); */

  double alpha = 2.0;
  double beta =  2.0;
  NumericsMatrix * C = NM_add(alpha, A, beta, B);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(Cref,0);
  Cref->matrix2->origin= NSM_CSC;
  NM_zentry(Cref, 0, 0, 4);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 6);
  NM_zentry(Cref, 1, 1, 8);
  NM_zentry(Cref, 1, 2, 6);
  NM_zentry(Cref, 2, 2, 6);
  NM_display(Cref);
  printf("add_square_csc: NM_equal(C,Cref) =%i \n", NM_equal(C,Cref));
  return  (int)!NM_equal(C,Cref);;

}

int add_rectangle_triplet()
{


  int size0 =3;
  int size1 =9;

  // product of triplet matrices into triplet matrix
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;
  NM_zentry(A, 0, 0, 1);
  NM_zentry(A, 0, 1, 2);
  NM_zentry(A, 0, 2, 3);
  NM_zentry(A, 1, 1, 2);
  NM_zentry(A, 1, 2, 3);
  NM_zentry(A, 2, 6, 2);
  NM_zentry(A, 2, 5, 22);
  NM_display(A);


  NumericsMatrix * B  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(B,0);
  B->matrix2->origin= NSM_TRIPLET;
  NM_zentry(B, 0, 0, 1);
  NM_zentry(B, 1, 1, 2);
  NM_zentry(B, 2, 2, 3);
  NM_zentry(B, 0, 3, 1);
  NM_zentry(B, 1, 4, 2);
  NM_zentry(B, 2, 5, 3);
  NM_display(B);

  double alpha = 1.0;
  double beta  = 2.0;
  NumericsMatrix * C = NM_add(alpha, A, beta, B);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(Cref,0);
  Cref->matrix2->origin= NSM_TRIPLET;
  NM_zentry(Cref, 0, 0, 3);
  NM_zentry(Cref, 0, 1, 2);
  NM_zentry(Cref, 0, 2, 3);
  NM_zentry(Cref, 0, 3, 2);

  NM_zentry(Cref, 1, 1, 6);
  NM_zentry(Cref, 1, 2, 3);
  NM_zentry(Cref, 1, 4, 4);

  NM_zentry(Cref, 2, 2, 6);


  NM_zentry(Cref, 2, 5, 28);
  NM_zentry(Cref, 2, 6, 2);
  NM_display(Cref);

  printf("add_rectangle_triplet : NM_equal(C,Cref) =%i \n", NM_equal(C,Cref));
  return (int)!NM_equal(C,Cref);


}




int main()
{

  int info = add_square_triplet();
  info += add_square_csc();
  info +=  add_rectangle_triplet();

  return info;
}
