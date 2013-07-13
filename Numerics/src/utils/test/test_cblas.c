#include <stdlib.h>
#include <stdio.h>

#include<SiconosBlas.h>

/* Parameters */
#define M 2
#define N 2
#define K 2

/* Main program */
int main() {
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

