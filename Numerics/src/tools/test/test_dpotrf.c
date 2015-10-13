#include <stdlib.h>
#include <stdio.h>
#include "SiconosLapack.h"
#include <assert.h>
#include "test_utils.h"

/* Parameters */
#define N 3
#define LDA N

/* Main program */
int main() {
  /* Locals */
  int n = N,lda = LDA;
  int info = 0;
  double a[LDA*N] = {
    2, -1, 0,
    -1, 2, -1,
    0, -1,  2
  };
        
  /* Executable statements */
  printf( "Dpotrf  Example Program Results\n" );
  DPOTRF(LA_UP, n, a, LDA, &info );
  /* Check for the exact singularity */
  if( info > 0 ) {
    printf( "The diagonal element of the triangular factor of A,\n" );
    printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
    printf( "the solution could not be computed.\n" );
    exit( 1 );
  }
        
  /* Print details of LU factorization */
  print_matrix( "Details of LU factorization", n, n, a, lda );
        
  exit( 0 );
} 
