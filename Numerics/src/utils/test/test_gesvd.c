#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "SiconosLapack.h"
#include "test_utils.h"
/* Parameters */
#define M 6
#define N 5
#define MACRO_MIN(a,b) ((a)>(b) ? (b):(a))
#define LDA M
#define LDU M
#define LDVT N

/* Main program */
int main() {
        /* Locals */
        int m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT, info = 0;
        double superb[MACRO_MIN(M,N)-1];
        double tol = 1e-2;
        /* Local arrays */
        double s[N], u[LDU*M], vt[LDVT*N];
        double a[LDA*N] = {
            8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
            9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
            9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
            5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
            3.16,  7.98,  3.01,  5.80,  4.27, -5.31
        };
        /* Executable statements */
       
        /* Compute SVD */
        //info = LAPACKE_dgesvd(CblasColMajor,'A', 'A', m, n, a, lda,s, u, ldu, vt, ldvt, superb);
        char JOBU = 'S', JOBVT = 'S';
        DGESVD(JOBU,JOBVT, m, n, a, lda,s, u, ldu, vt, ldvt, superb,&info);
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm computing SVD failed to converge.\n" );
                exit( 1 );
        }
        // printf( "LAPACKE_dgesvd (column-major, high-level) Example Program Results\n" );
        double singularValues[N] = {27.47, 22.64, 8.56, 5.99, 2.01};
        for (int i = 0; i<N;++i)
          if(fabs(s[i]- singularValues[i])>tol) exit(1);
 
        double LeftSingularVectors[LDU*N] =
          { -0.59 , -0.40, -0.03, -0.43, -0.47, 0.29,
            0.26, 0.24, -0.6, 0.24, -0.35, 0.58,
            0.36, -0.22, -0.45, -0.69, 0.39, -0.02,
            0.31, -0.75, 0.23, 0.33, 0.16, 0.38,
            0.23, -0.36, -0.31, 0.16, -0.52, -0.65
          };
 
        for (int i = 0; i<LDU*N;++i)
          if(fabs(u[i]- LeftSingularVectors[i])>tol) exit(1);

        double RightSingularVectors[LDVT*N] =
          { -0.25, 0.81, -0.26, 0.4, -0.22,
            -0.4, 0.36, 0.7, -0.45, 0.14,
            -0.69, -0.25, -0.22, 0.25, 0.59,
            -0.37, -0.37, 0.39, 0.43, -0.63,
            -0.41, -0.1, -0.49, -0.62, -0.44
          };
 
        for (int i = 0; i<LDVT*N;++i)
          if(fabs(vt[i]- RightSingularVectors[i])>tol) exit(1);
         
        /* Print singular values */
        print_matrix( "Singular values", 1, n, s, 1 );
        /* Print left singular vectors */
        print_matrix( "Left singular vectors (stored columnwise)", m, n, u, ldu );
        /* Print right singular vectors */
        print_matrix( "Right singular vectors (stored rowwise)", n, n, vt, ldvt );
        exit( 0 );
} /* End of LAPACKE_dgesvd Example */

