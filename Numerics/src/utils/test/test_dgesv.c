#include <stdlib.h>
#include <stdio.h>
#include "SiconosLapack.h"
#include <assert.h>
#include "test_utils.h"
#include "math.h"

/* Parameters */
#define N 5
#define NRHS 3
#define LDA N
#define LDB N

/* Main program */
int main() {
        /* Locals */
        int n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
        /* Local arrays */
        int ipiv[N];
        double a[LDA*N] = {
            6.80, -2.11,  5.66,  5.97,  8.23,
           -6.05, -3.30,  5.36, -4.44,  1.08,
           -0.45,  2.58, -2.70,  0.27,  9.04,
            8.32,  2.71,  4.35, -7.17,  2.14,
           -9.67, -5.14, -7.26,  6.08, -6.87
        };
        double b[LDB*NRHS] = {
            4.02,  6.19, -8.22, -7.57, -3.03,
           -1.56,  4.00, -8.67,  1.75,  2.86,
            9.81, -4.09, -4.57, -8.61,  8.99
        };
        /* Executable statements */
        printf( "LAPACKE_dgesv (column-major, high-level) Example Program Results\n" );
        /* Solve the equations A*X = B */
        DGESV(n, nrhs, a, lda, ipiv, b, ldb, &info );
        /* Check for the exact singularity */
        if( info > 0 ) {
                printf( "The diagonal element of the triangular factor of A,\n" );
                printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
                printf( "the solution could not be computed.\n" );
                exit( 1 );
        }

        double sol[LDB*NRHS] = {
          -0.80,  -0.7, 0.59, 1.32, 0.57,
          -0.39,  -0.55, 0.84, -0.1, 0.11,
          0.96, 0.22, 1.90, 5.36, 4.04
        };

       unsigned int err = 0;
       for (int i = 0; i<LDB*NRHS;++i)
          if( fabs(b[i]- sol[i]) > 1e-8)
          {
              printf("err at index %i: %6.10f\n", i, b[i]);
              err = 1;
          }
        //return 1;

        /* Print solution */
        print_matrix( "Solution", n, nrhs, b, ldb );
        /* Print details of LU factorization */
        print_matrix( "Details of LU factorization", n, n, a, lda );
        /* Print pivot indices */
        print_int_vector( "Pivot indices", n, ipiv );
        exit( err );
} /* End of LAPACKE_dgesv Example */

