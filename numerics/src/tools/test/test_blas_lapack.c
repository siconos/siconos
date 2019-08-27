#include <stdlib.h>
#include <stdio.h>
#include "math.h"
#include<SiconosBlas.h>

#include "SiconosLapack.h"
#include "test_utils.h"


#include <assert.h>

/* Parameters */
#define M 2
#define N 2
#define K 2

/* Main program */
static int test_dgemm(void)
{
  double alpha = 1.0;
  double A[M*K];
  double B[K*N];
  int lda = M;
  int ldb = N;
  double beta = 1.0;
  double C[M*N];
  int ldc = M;
  
  for(int i=0;i<M;i++)
    for(int j=0;j<N;j++)
    {
      A[i+j*M] = i + 10*j;
      B[i+j*M] = 2.0;
      C[i+j*M] = 0.0;
    }

  for(int i=0;i<M;i++)
  {
    for(int j=0;j<N;j++)
    {
      printf("%f\t", A[i+j*M]);
    }
    printf("\n");
  }
  for(int i=0;i<M;i++)
  {
    for(int j=0;j<N;j++)
    {
      printf("%f\t", B[i+j*M]);
    }
    printf("\n");
  }

  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);

  for(int i=0;i<M;i++)
  {
    for(int j=0;j<N;j++)
    {
      printf("%f\t", C[i+j*M]);
    }
    printf("\n");
  }
  
  return 0;
} 


/* Parameters */
#undef M
#define M 6
#undef N
#define N 4
#define NRHS 2
#undef LDA
#define LDA M
#undef LDB
#define LDB M

/* Main program */
#ifdef HAS_LAPACK_dgels
static int test_dgels(void)
{
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
#endif

/* Parameters */
#undef N
#define N 5
#undef NRHS
#define NRHS 3
#undef LDA
#define LDA N
#undef LDB
#define LDB N

/* Main program */
static int test_dgesv(void)
{
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

/* Parameters */
#undef N
#define N 5
#undef LDA
#define LDA N

/* Main program */
static int test_dgetrf(void)
{
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

/* Parameters */
#undef N
#define N 3
#undef LDA
#define LDA N

/* Main program */
static int test_dpotrf(void)
{
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

/* Parameters */
#undef M
#define M 6
#undef N
#define N 5
#define MACRO_MIN(a,b) ((a)>(b) ? (b):(a))
#undef LDA
#define LDA M
#undef LDU
#define LDU M
#undef LDVT
#define LDVT N

/* Main program */
#ifdef HAS_LAPACK_dgesvd
static int test_dgesvd(void)
{
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
#endif



int main(void)
{
  int info = test_dgemm();

  info += test_dgesv();
  info += test_dgetrf();
  info += test_dpotrf();
#ifdef HAS_LAPACK_dgesvd
  info += test_dgesvd();
#endif
#ifdef  HAS_LAPACK_dgels
  info += test_dgels();
#endif
  return info;
}
