#include <stdlib.h>
#include <stdio.h>
#include "SiconosLapack.h"
#include "test_utils.h"

/* Parameters */
#define M 6
#define N 4
#define NRHS 2
#define LDA M
#define LDB M

/* Main program */
int main() {
        /* Locals */
        int m = M, n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
        /* Local arrays */
        double a[LDA*N] = {
            1.44, -9.96, -7.55,  8.34,  7.08, -5.45,
           -7.84, -0.28,  3.24,  8.09,  2.52, -5.70,
           -4.39, -3.24,  6.27,  5.28,  0.74, -1.19,
            4.53,  3.83, -6.64,  2.06, -2.47,  4.70
        };
        double b[LDB*NRHS] = {
            8.58,  8.26,  8.48, -5.28,  5.72,  8.93,
            9.35, -4.43, -0.70, -0.26, -7.36, -2.52
        };
        /* Executable statements */
        printf( "LAPACKE_dgels (column-major, high-level) Example Program Results\n" );
        /* Solve the equations A*X = B */

        DGELS('N', m, n, nrhs, a, lda, b, ldb , &info);
        /* Check for the full rank */
        if( info > 0 ) {
                printf( "The diagonal element %i of the triangular factor ", info );
                printf( "of A is zero, so that A does not have full rank;\n" );
                printf( "the least squares solution could not be computed.\n" );
                exit( 1 );
        }
        /* Print least squares solution */
        print_matrix( "Least squares solution", n, nrhs, b, ldb );

        /* Print residual sum of squares for the solution */
        print_vector_norm( "Residual sum of squares for the solution", m-n, nrhs,
                        &b[n], ldb );
        /* Print details of QR factorization */
        print_matrix( "Details of QR factorization", m, n, a, lda );
        exit( 0 );
} /* End of LAPACKE_dgels Example */
