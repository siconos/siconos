#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "CSparseMatrix.h"

#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "CSparseMatrix.h"

#include "CSparseMatrix_internal.h"

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




static int add_test(void)
{

  int info = add_square_triplet();
  info += add_square_csc();
  info +=  add_rectangle_triplet();

  return info;
}

int gemm_square_triplet(void);
int gemm_square_csc(void);
int gemm_square_triplet_into_csc(void);

int gemm_rectangle_triplet(void);

int gemm_square_triplet()
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

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(C,0);
  C->matrix2->origin= NSM_TRIPLET;

  double alpha = 1.0;
  double beta =0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(Cref,0);
  Cref->matrix2->origin= NSM_TRIPLET;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */


  printf("NM_equal(C,Cref) =%i \n", NM_equal(C,Cref));
  return (int)!NM_equal(C,Cref);


}

int gemm_square_csc()
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

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size1);
  int nzmax = 10;
  NM_csc_empty_alloc(C,nzmax);
  C->matrix2->origin= NSM_CSC;

  NM_display(C);
  double alpha = 1.0;
  double beta = 0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(Cref,0);
  Cref->matrix2->origin= NSM_CSC;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */
  return  (int)!NM_equal(C,Cref);;

}
int gemm_square_triplet_into_csc()
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

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size1);
  int nzmax = 10;
  NM_csc_empty_alloc(C,nzmax);
  C->matrix2->origin= NSM_CSC;

  double alpha = 1.0;
  double beta =0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);
  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(Cref,0);
  Cref->matrix2->origin= NSM_CSC;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */
  return  (int)!NM_equal(C,Cref);;
}

int gemm_rectangle_triplet()
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




  NumericsMatrix * B  = NM_create(NM_SPARSE, size1, size0);
  NM_triplet_alloc(B,0);
  B->matrix2->origin= NSM_TRIPLET;
  NM_zentry(B, 0, 0, 1);
  NM_zentry(B, 1, 1, 2);
  NM_zentry(B, 2, 2, 3);
  NM_zentry(B, 3, 0, 1);
  NM_zentry(B, 4, 1, 2);
  NM_zentry(B, 5, 2, 3);
  NM_display(B);

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size0);
  NM_triplet_alloc(C,0);
  C->matrix2->origin= NSM_TRIPLET;

  double alpha = 1.0;
  double beta =0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size0);
  NM_triplet_alloc(Cref,0);
  Cref->matrix2->origin= NSM_TRIPLET;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_zentry(Cref, 2, 2, 66);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */


  printf("NM_equal(C,Cref) =%i \n", NM_equal(C,Cref));
  return (int)!NM_equal(C,Cref);


}




static int gemm_test(void)
{

  int info = gemm_square_triplet();
  info += gemm_square_csc();
  info +=  gemm_square_triplet_into_csc();
  info +=  gemm_rectangle_triplet();

  return info;
}

int gemm_square_triplet_1(void);
int gemm_square_csc_1(void);
int gemm_square_triplet_into_csc_1(void);
int gemm_rectangle_triplet_1(void);

int gemm_square_triplet_1()
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

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(C,0);
  C->matrix2->origin= NSM_TRIPLET;

  double alpha = 1.0;
  double beta =0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(Cref,0);
  Cref->matrix2->origin= NSM_TRIPLET;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */


  printf("NM_equal(C,Cref) =%i \n", NM_equal(C,Cref));
  return (int)!NM_equal(C,Cref);


}

int gemm_square_csc_1()
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

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size1);
  int nzmax = 10;
  NM_csc_empty_alloc(C,nzmax);
  C->matrix2->origin= NSM_CSC;

  NM_display(C);
  double alpha = 1.0;
  double beta = 0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(Cref,0);
  Cref->matrix2->origin= NSM_CSC;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */
  return  (int)!NM_equal(C,Cref);;

}
int gemm_square_triplet_into_csc_1()
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

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size1);
  int nzmax = 10;
  NM_csc_empty_alloc(C,nzmax);
  C->matrix2->origin= NSM_CSC;

  double alpha = 1.0;
  double beta =0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);
  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(Cref,0);
  Cref->matrix2->origin= NSM_CSC;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */
  return  (int)!NM_equal(C,Cref);;
}

int gemm_rectangle_triplet_1()
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




  NumericsMatrix * B  = NM_create(NM_SPARSE, size1, size0);
  NM_triplet_alloc(B,0);
  B->matrix2->origin= NSM_TRIPLET;
  NM_zentry(B, 0, 0, 1);
  NM_zentry(B, 1, 1, 2);
  NM_zentry(B, 2, 2, 3);
  NM_zentry(B, 3, 0, 1);
  NM_zentry(B, 4, 1, 2);
  NM_zentry(B, 5, 2, 3);
  NM_display(B);

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size0);
  NM_triplet_alloc(C,0);
  C->matrix2->origin= NSM_TRIPLET;

  double alpha = 1.0;
  double beta =0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size0);
  NM_triplet_alloc(Cref,0);
  Cref->matrix2->origin= NSM_TRIPLET;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_zentry(Cref, 2, 2, 66);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */


  printf("NM_equal(C,Cref) =%i \n", NM_equal(C,Cref));
  return (int)!NM_equal(C,Cref);


}




static int gemm_test2(void)
{

  int info = gemm_square_triplet();
  info += gemm_square_csc();
  info +=  gemm_square_triplet_into_csc();
  info +=  gemm_rectangle_triplet();

  return info;
}



/* create an empty triplet matrix, insert 2 elements, print and free */
static int alloc_test(void)
{
  CSparseMatrix *m = cs_spalloc(0,0,0,0,1); /* coo format */

  CS_INT info1 = 1-cs_entry(m, 3, 4, 1.0);
  CS_INT info2 = 1-cs_entry(m, 1, 2, 2.0);

  CS_INT info3 = 1-cs_print(m, 0);

  m=cs_spfree(m);

  CS_INT info4 = 1-(m==NULL);

  return (int)(info1+info2+info3+info4);
}

/* create an empty triplet matrix, insert 2 elements, print and free */
static int inv_test(void)
{
  int size0 =10;
  int size1 =10;
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;
  for (int i =0; i < size0; i++)
  {
    for (int j =i; j < size1; j++)
    {
      NM_zentry(A, i, j, i+j+1);
    }
  }

  //NM_zentry(A, size0-1, size0-1, 10);

  NM_display(A);
  FILE * fileout = fopen("dataA.py", "w");
  NM_write_in_file_python(A, fileout);
  fclose(fileout);
 
  NumericsMatrix * Ainv  = NM_new();
  Ainv->size0 = size0;
  Ainv->size1 = size1;
  Ainv->storageType = NM_SPARSE;

  NM_inv(A, Ainv);

  NumericsMatrix* AAinv = NM_multiply(A,Ainv);
  //NM_display(AAinv);

  NumericsMatrix * Id  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(Id,0);
  Id->matrix2->origin= NSM_TRIPLET;
  for (int i =0; i < size0; i++)
  {
    NM_zentry(Id, i, i, 1.0);
  }
  //NM_display(Id);

  
  return  !NM_equal(AAinv, Id);
}



int main(int argc, char *argv[])
{

#ifdef SICONOS_HAS_MPI
  MPI_Init(&argc, &argv);
#endif


  int info = add_test();

  info += gemm_test();

  info += gemm_test2();

  info += alloc_test();
  
  info += inv_test();

#ifdef SICONOS_HAS_MPI
    MPI_Finalize();
#endif

  
  return info;
}
