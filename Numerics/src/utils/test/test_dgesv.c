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

        double sol[LDB*NRHS] = {-0.80071403, -0.69524338, 0.59391499, 1.32172561, 0.5657562, 
                                -0.38962139, -0.55442713,  0.84222739, -0.10380185,  0.10571095,
                                0.95546491, 0.22065963, 1.90063673, 5.35766149, 4.04060266
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

