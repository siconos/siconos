///////////////////////////////////////////////////////////////////////////
//  This main file allows the resolution of LCP (Linear Complementary Problem):
//  try (z,w) such that:
//
//                        M z - w = q        (1)
//                        0<= z , 0<= w      (2)
//                        z.w= O             (3)
//
//  here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional
//  vector and w an n-dimensional vector.
//
//  This system of equations and inequalities is solved thanks to subroutine_lcp:
//        lemke_lcp (M,q,n,itermax,z,w,it_end,res,info)
//        gsnl_lcp(M,q,n,itermax,tol,z,w,it_end,res,info)
//        gcp_lcp(M,q,n,itermax,tol,z,w,it_end,res,info)
//        latin_lcp (M,q,n,k_latin,itermax,tol,z,w,it_end,res,info)
//        qp_lcp(M,q,n,res,z,w,fail) (not available for the moment)
//
//  where _ itermax is the maximum iterations required, it's an integer
//        _ res is the residue, it's a float
//        _ it_end is the number of iterations carried out, it's an integer
//        _ tol is the tolerance value desired, it's a float
//        _ k_latin is the parameter research of the latin, it's a float
//        _ fail shows the termination reason, it's an integer
//        _ info shows the trmination reason (0 successful otherwise 1), it's an integer.
//
//  The subroutine's call is due to the function solve_lcp:
//
//  int solve_lcp (float (*M)[maxcols],float * q, int n, methode *pt, float *z, float * w)
//
//  where M is an n by n matrix, q an n-dimensional vector, n is the row dimension of M, and pt a pointer other a structure (methode), z and w are n-dimensional vector, the solutions of the lcp.
//  methode is a variable with a structure type; this structure gives to the function solve_lcp, the name and the parameters (itermax, tol, k_latin) of the method we want to use.
//  This function return an interger:  0 successful return otherwise 1.
//
//
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "solverpack.h"


void test_lcp_series(int n, double ** M, double *q)
{

  int i, j;
  int info;
  double *z, *w;
  double *vec;



  z = malloc(n * sizeof(double));
  w = malloc(n * sizeof(double));

  vec = (double*)malloc(n * n * sizeof(double));

  //// Fortran compatibility
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      vec[j * n + i] = M[i][j];

  /////////////////////////////////
  // Premier test //////////////////
  /////////////////////////////////
  printf("\n Test de GSNL\n");

  static methode_lcp meth_lcp  = {"Gsnl", 1001, 0.000001, 0.6};

  printf("\n we go in the function solve\n");
  // for (i=0;i<n*n;i++)  printf("vec[%i] = %g\n",i,vec[i]);


  info = solve_lcp(vec, q, &n, &meth_lcp, z, w);

  for (i = 0; i < n; i++)
    printf("z %10.4e\t\t w %10.4e\n", z[i], w[i]);


  printf("\n we go out the function and info is %d\n", info);

  /////////////////////////////////
  // second test //////////////////
  /////////////////////////////////
  printf("\n Test de Gcp\n");

  static methode_lcp meth_lcp2  = {"Gcp", 1000, -0.0000001, 0.6};

  printf("\n we go in the function solve\n");
  //for (i=0;i<n*n;i++) printf("vec[%i] = %g\n",i,vec[i]);

  info = solve_lcp(vec, q, &n, &meth_lcp2, z, w);

  for (i = 0; i < n; i++)
    printf("z %10.4e\t\t w %10.4e\n", z[i], w[i]);


  printf("\n we go out the function and info is %d\n", info);

  /////////////////////////////////
  // third test //////////////////
  /////////////////////////////////
  printf("\n Test de Latin\n");
  static methode_lcp meth_lcp3  = {"Latin", 1000, -0.0000001, 0.7};
  printf("\n we go in the function\n");
  // for (i=0;i<n*n;i++)  printf("vec[%i] = %g\n",i,vec[i]);

  info = solve_lcp(vec, q, &n, &meth_lcp3, z, w);

  for (i = 0; i < n; i++)
    printf("z %10.4e\t\t w %10.4e\n", z[i], w[i]);



  printf("\n we go out the function and info is %d\n", info);

  /////////////////////////////////
  // Lemke test //////////////////
  /////////////////////////////////
  printf("\n Test de Lemke\n");
  static methode_lcp meth_lcp4 = {"Lemke", 1000, -0.0000001, 0.7};
  //for (i=0;i<n*n;i++) printf("vec[%i] = %g\n",i,vec[i]);
  info = solve_lcp(vec, q, &n, &meth_lcp4, z, w);
  for (i = 0; i < n; i++)
    printf("z %10.4e\t\t w %10.4e\n", z[i], w[i]);
  printf("\n Info is %d\n", info);


  /////////////////////////////////
  // Qp test //////////////////
  /////////////////////////////////


  printf("\n Test du Qp\n");
  static methode_lcp meth_lcp5 = {"Qp", 1000, 0.001, 0.7};
  //for (i=0;i<n*n;i++) printf("vec[%i] = %g\n",i,vec[i]);
  info = solve_lcp(vec, q, &n, &meth_lcp5, z, w);
  for (i = 0; i < n; i++)
    printf("z %10.4e\t\t w %10.4e\n", z[i], w[i]);



  printf("\n Info is %d\n", info);

  /////////////////////////////////
  // Qpnonsym test //////////////////
  /////////////////////////////////


  printf("\n Test du Qpnonsym\n");
  static methode_lcp meth_lcp6 = {"Qpnonsym", 1000, 0.001, 0.7};


  printf("\n we go in the function solve_lcp\n");
  //for (i=0;i<n*n;i++) printf("vec[%i] = %g\n",i,vec[i]);

  info = solve_lcp(vec, q, &n, &meth_lcp6, z, w);

  for (i = 0; i < n; i++)
    printf("z %10.4e\t\t w %10.4e\n", z[i], w[i]);



  printf("\n we go out the function and info is %d\n", info);
  /////////////////////////////////
  // LexicoLemke test //////////////////
  /////////////////////////////////


  printf("\n Test du LexicoLemke\n");

  static methode_lcp meth_lcp7 = {"LexicoLemke", 1000, 1e-8, 0.7};

  //for (i=0;i<n*n;i++) printf("vec[%i] = %g\n",i,vec[i]);
  info = solve_lcp(vec, q, &n, &meth_lcp7, z, w);
  for (i = 0; i < n; i++)
    printf("z %10.4e\t\t w %10.4e\n", z[i], w[i]);
  printf("\n Info is %d\n", info);

  printf("\n we go out the function and info is %d\n", info);

  free(z);
  free(w);


}



void testmmc(void)
{

  printf("/////////////////////////////////////////////\n");
  printf("////////////  TEST MMC             //////////\n");
  printf("/////////////////////////////////////////////\n");

  FILE *f1, *f2, *f3, *f4, *f5, *f6;
  int i, j, nl, nc, nll, it, n, dimM;
  double *q, *vec;
  int info;
  double qi, Mij;
  char val[14], vall[14];
  int symmetric;
  double **M;
  // computation of the size of the problem
  n = 0;
  if ((f1 = fopen("MM_mmc_mu2.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }
  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    n = nl;
  }

  printf("size of the problem n= %i\n", n);
  fclose(f1);

  dimM = n;

  // Memory allocation

  //double (*M)[n];
  //M =malloc(dimM*dimM*sizeof(double));


  M = (double **)malloc(dimM * sizeof(double*));
  for (i = 0; i < n; i++)
    M[i] = (double*)malloc(dimM * sizeof(double));




  // Data loading of M and q

  if ((f1 = fopen("MM_mmc_mu2.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }

  if ((f2 = fopen("qq_mmc_mu2.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }







  for (i = 0; i < dimM; i++)
    for (j = 0; j < dimM; j++)
      M[i][j] = 0.;

  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);

    /////////////       on met la transpos       ////////////////
    /*fprintf(f3,"%d %d %.14e\n",nc,nl,Mij);
      fscanf(f3,"%d %d %.14e\n", &nc, &nl, &Mij);*/
    *(*(M + nc - 1) + nl - 1) = Mij;
    //////////////         fin transpos         ////////////////////

  }



  ////////////////////////////////////////////////////////////////////////

  q = malloc(dimM * sizeof(double));

  while (!feof(f2))
  {
    fscanf(f2, "%d", &nll);
    fscanf(f2, "%s", vall);
    qi = atof(vall);

    *(q + nll - 1) = -qi;
  }
  fclose(f2);
  fclose(f1);

  // Test on M and q
  // Is  M symmetric ?
  symmetric = 0;
  for (i = 0; i < dimM; i++)
  {
    for (j = 0; j < i; j++)
    {
      if (abs(M[i][j] - M[j][i]) > 1e-10)
      {
        symmetric = 1;
        break;
      }
    }
  }



  if (symmetric)
  {
    printf("M is a not symmetrix matrix\n");
  }
  else
  {
    printf("M is a symmetrix matrix\n");
  }



  test_lcp_series(n, M, q);

  free(M);
  free(q);
}

void testmathieuold(void)
{

  printf("/////////////////////////////////////////////\n");
  printf("////////////  TEST MATHIEU OLD     //////////\n");
  printf("/////////////////////////////////////////////\n");

  FILE *f1, *f2, *f3, *f4, *f5, *f6;
  int i, j, k, nl, nc, nll, it, n, dimM;
  double *q, *vec, *data;
  int info;
  double qi, Mij;
  char val[14], vall[14], totochar[100];
  int symmetric;
  double **M;
  // computation of the size of the problem
  n = 6;
  dimM = n;

  // Memory allocation


  M = (double **)malloc(dimM * sizeof(double*));
  for (i = 0; i < n; i++)
    M[i] = (double*)malloc(dimM * sizeof(double));

  q = malloc(dimM * sizeof(double));

  // Data loading of M and q

  if ((f1 = fopen("matrix.test.old", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }

  for (i = 0; i < dimM; i++)
    for (j = 0; j < dimM; j++)
      M[i][j] = 0.;

  fscanf(f1, "%s", val);
  k = 0;

  //fscanf(f1,"%f %f %f %f %f %f",M[i][0],M[i][1],M[i][2],M[i][3],M[i][4],M[i][5]);
  data = malloc(dimM * dimM * dimM * sizeof(double));
  while (!feof(f1))
  {
    fscanf(f1, "%s", val);
    data[k] = atof(val);
    k++;
  }
  for (i = 0; i < k; i++)
  {
    printf("data[%i] = %g\n", i, data[i]);
  }

  fclose(f1);


  for (i = 0; i < dimM; i++)
    for (j = 0; j < dimM; j++)
      M[i][j] = data[i * dimM + j];

  for (i = 0; i < dimM; i++)
    q[i] = -data[dimM * dimM + i];




  for (i = 0; i < dimM; i++)
  {
    printf("q[%i] = %g\n", i, q[i]);
    for (j = 0; j < dimM; j++) printf("M[%i,%i ] = %g\n", i, j, M[i][j]);
  }


  // Test on M and q
  // Is  M symmetric ?
  symmetric = 0;
  for (i = 0; i < dimM; i++)
  {
    for (j = 0; j < i; j++)
    {
      if (abs(M[i][j] - M[j][i]) > 1e-10)
      {
        symmetric = 1;
        break;
      }
    }
  }



  if (symmetric)
  {
    printf("M is a not symmetrix matrix\n");
  }
  else
  {
    printf("M is a symmetrix matrix\n");
  }



  test_lcp_series(n, M, q);

  free(M);
  free(q);
}

void testmathieu1(void)
{

  printf("/////////////////////////////////////////////\n");
  printf("////////////  TEST MATHIEU1        //////////\n");
  printf("/////////////////////////////////////////////\n");

  FILE *f1, *f2, *f3, *f4, *f5, *f6;
  int i, j, k, nl, nc, nll, it, n, dimM;
  double *q, *vec, *data;
  int info;
  double qi, Mij;
  char val[14], vall[14], totochar[100];
  int symmetric;
  double **M;
  // computation of the size of the problem
  n = 6;
  dimM = n;

  // Memory allocation


  M = (double **)malloc(dimM * sizeof(double*));
  for (i = 0; i < n; i++)
    M[i] = (double*)malloc(dimM * sizeof(double));

  q = malloc(dimM * sizeof(double));
  for (i = 0; i < dimM; i++)
  {
    q[i] = 0.0;
    for (j = 0; j < dimM; j++) M[i][j] = 0.;
  }

  // Data loading of M and q
  data = malloc(dimM * dimM * dimM * sizeof(double));
  if ((f1 = fopen("mathieu.m1.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }
  k = 0;

  while (!feof(f1))
  {
    fscanf(f1, "%s", val);
    data[k] = atof(val);
    k++;
  }
  for (i = 0; i < k; i++)
  {
    printf("data[%i] = %g\n", i, data[i]);
  }

  fclose(f1);


  for (i = 0; i < dimM; i++)
    for (j = 0; j < dimM; j++)
      M[i][j] = data[i * dimM + j];

  if ((f1 = fopen("mathieu.q1.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }
  k = 0;

  while (!feof(f1))
  {
    fscanf(f1, "%s", val);
    data[k] = atof(val);
    k++;
  }
  for (i = 0; i < k; i++)
  {
    printf("data[%i] = %g\n", i, data[i]);
  }

  fclose(f1);




  for (i = 0; i < dimM; i++)
    q[i] = -data[i];




  for (i = 0; i < dimM; i++)
  {
    printf("q[%i] = %g\n", i, q[i]);
    for (j = 0; j < dimM; j++) printf("M[%i,%i ] = %g\n", i, j, M[i][j]);
  }


  // Test on M and q
  // Is  M symmetric ?
  symmetric = 0;
  for (i = 0; i < dimM; i++)
  {
    for (j = 0; j < i; j++)
    {
      if (abs(M[i][j] - M[j][i]) > 1e-10)
      {
        symmetric = 1;
        break;
      }
    }
  }



  if (symmetric)
  {
    printf("M is a not symmetrix matrix\n");
  }
  else
  {
    printf("M is a symmetrix matrix\n");
  }



  test_lcp_series(n, M, q);

  free(M);
  free(q);
}

void test1(void)
{
  FILE *f1, *f2, *f3, *f4, *f5, *f6;
  int i, j, nl, nc, nll, it, n, dimM;
  double *q, *vec;
  int info;
  double qi, Mij;
  char val[14], vall[14];
  int symmetric;
  double **M;

  printf("/////////////////////////////////////////////\n");
  printf("////////////  TEST 1               //////////\n");
  printf("/////////////////////////////////////////////\n");


  dimM = 2;
  n = 2;

  // Memory allocation

  M = (double **)malloc(dimM * sizeof(double*));
  for (i = 0; i < n; i++)
    M[i] = (double*)malloc(dimM * sizeof(double));

  q = (double*)malloc(dimM * sizeof(double));

  M[0][0] = 2;
  M[0][1] = 1;
  M[1][0] = 1;
  M[1][1] = 2;

  q[0] = 5.0;
  q[1] = 6.0;

  test_lcp_series(n, M, q);

  q[0] = -5.0;
  q[1] = -6.0;

  test_lcp_series(n, M, q);


  free(M);
  free(q);
}


void test2(void)
{
  FILE *f1, *f2, *f3, *f4, *f5, *f6;
  int i, j, nl, nc, nll, it, n, dimM;
  double *q, *vec;
  int info;
  double qi, Mij;
  char val[14], vall[14];
  int symmetric;
  double **M;
  printf("/////////////////////////////////////////////\n");
  printf("////////////  TEST 2               //////////\n");
  printf("/////////////////////////////////////////////\n");
  //M=[[1 -1 -1 -1];[-1 1 -1 -1 ];[1 1 2 0];[1 1 0 2]]
  //q=[-3 ;-5 ; 9 ; 5]

  dimM = 4;
  n = 4;

  // Memory allocation

  M = (double **)malloc(dimM * sizeof(double*));
  for (i = 0; i < n; i++)
    M[i] = (double*)malloc(dimM * sizeof(double));

  q = (double*)malloc(dimM * sizeof(double));

  M[0][0] = 1;
  M[0][1] = -1;
  M[0][2] = -1;
  M[0][3] = -1;

  M[1][0] = -1;
  M[1][1] = 1;
  M[1][2] = -1;
  M[1][3] = -1;

  M[2][0] = 1;
  M[2][1] = 1;
  M[2][2] = 2;
  M[2][3] = 0;

  M[3][0] = 1;
  M[3][1] = 1;
  M[3][2] = 0;
  M[3][3] = 2;





  q[0] = -3.0;
  q[1] = -5.0;
  q[2] = 9.0;
  q[3] = 5.0;

  test_lcp_series(n, M, q);


  free(M);
  free(q);
}



void test3(void)
{
  FILE *f1, *f2, *f3, *f4, *f5, *f6;
  int i, j, nl, nc, nll, it, n, dimM;
  double *q, *vec;
  int info;
  double qi, Mij;
  char val[14], vall[14];
  int symmetric;
  double **M;
  printf("/////////////////////////////////////////////\n");
  printf("////////////  TEST 3               //////////\n");
  printf("/////////////////////////////////////////////\n");
  //M=[[1 0 2 -1];[0 1 1 4 ];[-2 -1 0 0 ];[1 -4 0 0 ]]
  //q=[1/2 ; -1/2 ; 6 ;6 ]

  dimM = 4;
  n = 4;

  // Memory allocation

  M = (double **)malloc(dimM * sizeof(double*));
  for (i = 0; i < n; i++)
    M[i] = (double*)malloc(dimM * sizeof(double));

  q = (double*)malloc(dimM * sizeof(double));

  M[0][0] = 1;
  M[0][1] = 0;
  M[0][2] = 2;
  M[0][3] = -1;

  M[1][0] = 0;
  M[1][1] = 1;
  M[1][2] = 1;
  M[1][3] = 4;

  M[2][0] = -2;
  M[2][1] = -1;
  M[2][2] = 0;
  M[2][3] = 0;

  M[3][0] = 1;
  M[3][1] = -4;
  M[3][2] = 0;
  M[3][3] = 0;





  q[0] = -1. / 2;
  q[1] = 1. / 2;
  q[2] = -6.0;
  q[3] = -6.0;

  test_lcp_series(n, M, q);


  free(M);
  free(q);
}





void test4diodes(void)
{
  FILE *f1, *f2, *f3, *f4, *f5, *f6;
  int i, j, nl, nc, nll, it, n, dimM;
  double *q, *vec;
  int info;
  double qi, Mij;
  char val[14], vall[14];
  int symmetric;
  double **M;
  printf("/////////////////////////////////////////////\n");
  printf("////////////  TEST 4 diodes       ///////////\n");
  printf("/////////////////////////////////////////////\n");


  //M = [0.001  0.001  -1  0; 0.001  0.001  0  -1; 1  0  9.90099  -9.90099; 0  1  -9.90099  9.90099]
  //q = [-0 ;-0; 9.90099; -9.90099]

  dimM = 4;
  n = 4;

  // Memory allocation

  M = (double **)malloc(dimM * sizeof(double*));
  for (i = 0; i < n; i++)
    M[i] = (double*)malloc(dimM * sizeof(double));

  q = (double*)malloc(dimM * sizeof(double));


  M[0][0] = 0.001;
  M[0][1] = 0.001;
  M[0][2] = -1;
  M[0][3] = 0;

  M[1][0] = 0.001;
  M[1][1] = 0.001;
  M[1][2] = 0;
  M[1][3] = -1;

  M[2][0] = 1;
  M[2][1] = 0;
  M[2][2] = 9.90099;
  M[2][3] = -9.90099;

  M[3][0] = 0;
  M[3][1] = 1;
  M[3][2] = -9.90099;
  M[3][3] = 9.90099;





  q[0] = 0;
  q[1] = 0;
  q[2] = -9.90099;
  q[3] = 9.90099;

  test_lcp_series(n, M, q);


  free(M);
  free(q);
}






int main(void)

{


  test1();
  test2();
  test3();
  testmmc();
  testmathieu1();
  //testmathieuold();
  test4diodes();




  printf("\n The End \n");
}

