#include <stdlib.h>
#include <stdio.h>
#include "SiconosLapack.h"
#include <assert.h>
#include "test_utils.h"

/* Parameters */
#define N 5
#define LDA N

/* Main program */
int main() {
        /* Locals */
        int n = N, lda = LDA,  info;
        /* Local arrays */
        int ipiv[N];
        double a[LDA*N] = {
            6.80, -2.11,  5.66,  5.97,  8.23,
           -6.05, -3.30,  5.36, -4.44,  1.08,
           -0.45,  2.58, -2.70,  0.27,  9.04,
            8.32,  2.71,  4.35, -7.17,  2.14,
           -9.67, -5.14, -7.26,  6.08, -6.87
        };
        double b[LDA*N];
        cblas_dcopy (LDA*N, a, 1, b, 1);
        /* Executable statements */
        printf( "LAPACKE_dgetr (column-major, high-level) Example Program Results\n" );
        /* Solve the equations A*X = B */
        DGETRF(n, n, a, lda, ipiv, &info );
        if( info > 0 ) {
          printf( "DGETRF failed.\n" );
          exit( 1 );
        }

        DGETRI(n, a, lda, ipiv, &info );
        if( info > 0 ) {
          printf( "DGETRI failed.\n" );
          exit( 1 );
        }
        double tol = 1e-9;
        double c[LDA*N];
        cblas_dgemm(CblasColMajor, CblasNoTrans,CblasNoTrans, N, N, N, 1.0, a, N, b, N, 0.0, c, N);
        for(int i=0;i<N;++i)
          for(int j = 0; j<N; ++j)
            if(i==j)
            {if (fabs(c[i+N*j]-1.0)>tol) exit(1);}
            else 
            {if (fabs(c[i+N*j]) > tol) exit(1);}

        print_matrix("Id?", N,N,b,N);
        
        exit( 0 );
} 

