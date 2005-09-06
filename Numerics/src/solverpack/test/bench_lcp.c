/*
 *  Benchmark matrices for LCP (Linear Complementary Problem)
 *  This file generate the reading of data file contains in the folder MATRIX
 *  and try to found (z,w) such that:
 *
 *                        M z + q = w        (1)
 *                        0 <= z , 0 <= w    (2)
 *                        z.w= O             (3)
 *
 *  where M is an (n x n)-matrix, q an n-vector, z an n-vector and w an n-vector.
 *
 *  This system of equations and inequalities is solved thanks to subroutine_lcp:
 *
 *        lexicolemke_lcp (M,q,n,itermax,z,w,it_end,res,info)
 *        gsnl_lcp(M,q,n,itermax,tol,z,w,it_end,res,info)
 *        gcp_lcp(M,q,n,itermax,tol,z,w,it_end,res,info)
 *        qp_lcp(M,q,n,res,z,w,fail)
 *
 *  where _ itermax is the maximum iterations required, it's an integer
 *        _ res is the residue, it's a float
 *        _ it_end is the number of iterations carried out, it's an integer
 *        _ tol is the tolerance value desired, it's a float
 *        _ k_latin is the parameter research of the latin, it's a float
 *        _ fail shows the termination reason, it's an integer
 *        _ info shows the trmination reason (0 successful otherwise 1), it's an integer.
 *
 *  The subroutine's call is due to the function solve_lcp:
 *
 *  int solve_lcp (float (*M)[maxcols],float * q, int n, methode *pt, float *z, float * w)
 *
 *  where M is an n by n matrix, q an n-dimensional vector, n is the row dimension of M,
 *  and pt a pointer other a structure (methode), z and w are n-dimensional vector, the solutions of the lcp.
 *  methode is a variable with a structure type; this structure gives to the function solve_lcp,
 *  the name and the parameters (itermax, tol, k_latin) of the method we want to use.
 *  This function return an interger:  0 successful return otherwise 1.
 *
 * author : Mathieu Renouf
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "solverpack.h"
#include "blaslapack.h"

int main(void)
{

  FILE *LCPfile;

  int i, j;
  int incx = 1 , incy = 1;
  int info ;
  int dim , dim2;

  double a1;

  double RES;
  double *q , *z , *w , *sol;
  double *vecM;

  char val[20];

  /*
   * methode declaration :
   *
   * methode 1 for Gsnl
   * methode 2 for Cpg
   * methode 3 for lexicolemke
   *
   */

  static methode_lcp meth_lcp1 = { "Gsnl"        , 1000 , 1e-8 , 0.6 , 1.0 , 1 , "N2" };
  static methode_lcp meth_lcp2 = { "Gcp"         , 1000 , 1e-8 , 0.6 , 1.0 , 1 , "N2" };
  static methode_lcp meth_lcp3 = { "LexicoLemke" , 1000 , 1e-8 , 0.6 , 1.0 , 1 , "N2" };

  /*************************************************************
   *
   *                       TEST NUMBER ZERO
   *
   *************************************************************/

  if ((LCPfile = fopen("MATRIX/deudeu.dat", "r")) == NULL)
  {
    perror("fopen LCPfile: deudeu.dat");
    exit(1);
  }
  printf("\n SCAN 2x2 MATRIX FOR LCP \n");

  fscanf(LCPfile , "%d" , &dim);

  //dim = atoi(val);
  dim2 = dim * dim;

  vecM = (double*)malloc(dim2 * sizeof(double));
  q    = (double*)malloc(dim * sizeof(double));
  sol  = (double*)malloc(dim * sizeof(double));

  for (i = 0 ; i < dim ; ++i)  q[i]    = 0.0;
  for (i = 0 ; i < dim ; ++i)  sol[i]  = 0.0;
  for (i = 0 ; i < dim2 ; ++i) vecM[i] = 0.0;

  for (i = 0 ; i < dim ; ++i)
  {
    for (j = 0 ; j < dim ; ++j)
    {
      fscanf(LCPfile, "%s", val);
      vecM[ dim * j + i ] = atof(val);
    }
  }

  for (i = 0 ; i < dim ; ++i)
  {
    fscanf(LCPfile , "%s" , val);
    q[i] = atof(val);
  }

  /* Print INFO */

  for (i = 0 ; i < dim ; ++i)
  {
    for (j = 0 ; j < dim ; ++j) printf(" %10.4g " , vecM[ dim * j + i ]);
    printf("\n");
  }
  printf("\n");
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , q[i]);

  /* Print INFO */

  fscanf(LCPfile , "%s" , val);

  if (!feof(LCPfile))
  {

    sol[0] = atof(val);

    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }

    printf("\n SOLUTION OF LCP :\n");
    for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
    printf("\n");
  }
  else
  {
    printf("\n NO EXACT SOLUTION FOR LCP :\n");
    for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
  }

  fclose(LCPfile);

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp1 , z , w);

  printf("***********************\n");
  printf(" With NLGS:            \n");
  printf(" Solution :            \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  /* CPG */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp2 , z , w);

  printf("**********************\n");
  printf(" With CPG:            \n");
  printf(" Solution :           \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  /* Lemke */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp3 , z , w);

  printf("**********************\n");
  printf(" With Lemke:          \n");
  printf(" Solution :           \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  free(q);
  free(sol);
  free(vecM);

  /*************************************************************
   *
   *                       TEST NUMBER ONE
   *
   *************************************************************/

  if ((LCPfile = fopen("MATRIX/ortiz.dat", "r")) == NULL)
  {
    perror("fopen LCPfile:ortiz.dat");
    exit(1);
  }
  printf("\n SCAN ORTIZ MATRIX FOR LCP \n");

  fscanf(LCPfile , "%d" , &dim);

  //dim = atoi(val);
  dim2 = dim * dim;

  vecM = (double*)malloc(dim2 * sizeof(double));
  q    = (double*)malloc(dim * sizeof(double));
  sol  = (double*)malloc(dim * sizeof(double));

  for (i = 0 ; i < dim2 ; ++i) vecM[i] = 0.0;
  for (i = 0 ; i < dim ; ++i)  q[i]    = 0.0;
  for (i = 0 ; i < dim ; ++i)  sol[i]  = 0.0;

  for (i = 0 ; i < dim ; ++i)
  {
    for (j = 0 ; j < dim ; ++j)
    {
      fscanf(LCPfile, "%s", val);
      vecM[ dim * j + i ] = atof(val);
    }
  }

  for (i = 0 ; i < dim ; ++i)
  {
    fscanf(LCPfile , "%s" , val);
    q[i] = atof(val);
  }

  /* Print INFO */

  for (i = 0 ; i < dim ; ++i)
  {
    for (j = 0 ; j < dim ; ++j) printf(" %10.4g " , vecM[ dim * j + i ]);
    printf("\n");
  }
  printf("\n");
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , q[i]);

  /* Print INFO */

  fscanf(LCPfile , "%s" , val);

  if (!feof(LCPfile))
  {

    sol[0] = atof(val);

    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }

    printf("\n SOLUTION OF LCP :\n");
    for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
    printf("\n");
  }
  else
  {
    printf("\n NO EXACT SOLUTION FOR LCP :\n");
    for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
  }

  fclose(LCPfile);

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp1 , z , w);

  printf("***********************\n");
  printf(" With NLGS:            \n");
  printf(" Solution :            \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  /* CPG */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp2 , z , w);

  printf("**********************\n");
  printf(" With CPG:            \n");
  printf(" Solution :           \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  /* Lemke */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp3 , z , w);

  printf("**********************\n");
  printf(" With Lemke:          \n");
  printf(" Solution :           \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  free(q);
  free(sol);
  free(vecM);

  /*************************************************************
   *
   *                       TEST NUMBER TWO
   *
   *************************************************************/

  if ((LCPfile = fopen("MATRIX/diodes.dat", "r")) == NULL)
  {
    perror("fopen LCPfile: diodes.dat");
    exit(1);
  }
  printf("\n SCAN DIODE MATRIX FOR LCP \n");

  fscanf(LCPfile , "%d" , &dim);

  //dim = atoi(val);
  dim2 = dim * dim;

  vecM = (double*)malloc(dim2 * sizeof(double));
  q    = (double*)malloc(dim * sizeof(double));
  sol  = (double*)malloc(dim * sizeof(double));

  for (i = 0 ; i < dim ; ++i)  q[i]    = 0.0;
  for (i = 0 ; i < dim ; ++i)  sol[i]  = 0.0;
  for (i = 0 ; i < dim2 ; ++i) vecM[i] = 0.0;

  for (i = 0 ; i < dim ; ++i)
  {
    for (j = 0 ; j < dim ; ++j)
    {
      fscanf(LCPfile, "%s", val);
      vecM[ dim * j + i ] = atof(val);
    }
  }

  for (i = 0 ; i < dim ; ++i)
  {
    fscanf(LCPfile , "%s" , val);
    q[i] = atof(val);
  }

  /* Print INFO */

  for (i = 0 ; i < dim ; ++i)
  {
    for (j = 0 ; j < dim ; ++j) printf(" %10.4g " , vecM[ dim * j + i ]);
    printf("\n");
  }
  printf("\n");
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , q[i]);

  /* Print INFO */

  fscanf(LCPfile , "%s" , val);

  if (!feof(LCPfile))
  {

    sol[0] = atof(val);

    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }

    printf("\n SOLUTION OF LCP :\n");
    for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
    printf("\n");
  }
  else
  {
    printf("\n NO EXACT SOLUTION FOR LCP :\n");
    for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
  }

  fclose(LCPfile);

  printf("\n NOTE : DET(M) = 1 & COND(M) = 432 \n");
  printf("          M NON SYMETRIC             \n");
  printf("          Ker(M) = nul               \n");

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp1 , z , w);

  printf("***********************\n");
  printf(" With NLGS:            \n");
  printf(" Solution :            \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  /* CPG */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp2 , z , w);

  printf("**********************\n");
  printf(" With CPG:            \n");
  printf(" Solution :           \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  /* Lemke */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp3 , z , w);

  printf("**********************\n");
  printf(" With Lemke:          \n");
  printf(" Solution :           \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  free(q);
  free(sol);
  free(vecM);

  /*************************************************************
   *
   *                       TEST NUMBER THREE
   *
   *************************************************************/

  if ((LCPfile = fopen("MATRIX/pang.dat", "r")) == NULL)
  {
    perror("fopen LCPfile: pang.dat");
    exit(1);
  }
  printf("\n SCAN PANG MATRIX FOR LCP \n");

  fscanf(LCPfile , "%d" , &dim);

  //dim = atoi(val);
  dim2 = dim * dim;

  vecM = (double*)malloc(dim2 * sizeof(double));
  q    = (double*)malloc(dim * sizeof(double));
  sol  = (double*)malloc(dim * sizeof(double));

  for (i = 0 ; i < dim ; ++i)  q[i]    = 0.0;
  for (i = 0 ; i < dim ; ++i)  sol[i]  = 0.0;
  for (i = 0 ; i < dim2 ; ++i) vecM[i] = 0.0;

  for (i = 0 ; i < dim ; ++i)
  {
    for (j = 0 ; j < dim ; ++j)
    {
      fscanf(LCPfile, "%s", val);
      vecM[ dim * j + i ] = atof(val);
    }
  }

  for (i = 0 ; i < dim ; ++i)
  {
    fscanf(LCPfile , "%s" , val);
    q[i] = atof(val);
  }

  /* Print INFO */

  for (i = 0 ; i < dim ; ++i)
  {
    for (j = 0 ; j < dim ; ++j) printf(" %10.4g " , vecM[ dim * j + i ]);
    printf("\n");
  }
  printf("\n");
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , q[i]);

  /* Print INFO */

  fscanf(LCPfile , "%s" , val);

  if (!feof(LCPfile))
  {

    sol[0] = atof(val);

    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }

    printf("\n SOLUTION OF LCP :\n");
    for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
    printf("\n");
  }
  else
  {
    printf("\n NO EXACT SOLUTION FOR LCP \n");
    for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
  }

  fclose(LCPfile);

  printf("\n NOTE : DET(M) = 0        \n");
  printf("          M NON SYMETRIC    \n");
  printf("          Ker(M)={(-1,1,0)} \n");

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp1 , z , w);

  printf("***********************\n");
  printf(" With NLGS:            \n");
  printf(" Solution :            \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  /* CPG */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp2 , z , w);

  printf("**********************\n");
  printf(" With CPG:            \n");
  printf(" Solution :           \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  /* Lemke */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp3 , z , w);

  printf("**********************\n");
  printf(" With Lemke:          \n");
  printf(" Solution :           \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  free(q);
  free(sol);
  free(vecM);

  /*************************************************************
   *
   *                       TEST NUMBER FOUR
   *
   *************************************************************/

  if ((LCPfile = fopen("MATRIX/murty1.dat", "r")) == NULL)
  {
    perror("fopen LCPfile: murty1.dat");
    exit(1);
  }
  printf("\n SCAN MURTY MATRIX N°1 FOR LCP \n");

  fscanf(LCPfile , "%d" , &dim);

  //dim = atoi(val);
  dim2 = dim * dim;

  vecM = (double*)malloc(dim2 * sizeof(double));
  q    = (double*)malloc(dim * sizeof(double));
  sol  = (double*)malloc(dim * sizeof(double));

  for (i = 0 ; i < dim ; ++i)  q[i]    = 0.0;
  for (i = 0 ; i < dim ; ++i)  sol[i]  = 0.0;
  for (i = 0 ; i < dim2 ; ++i) vecM[i] = 0.0;

  for (i = 0 ; i < dim ; ++i)
  {
    for (j = 0 ; j < dim ; ++j)
    {
      fscanf(LCPfile, "%s", val);
      vecM[ dim * j + i ] = atof(val);
    }
  }

  for (i = 0 ; i < dim ; ++i)
  {
    fscanf(LCPfile , "%s" , val);
    q[i] = atof(val);
  }

  /* Print INFO */

  for (i = 0 ; i < dim ; ++i)
  {
    for (j = 0 ; j < dim ; ++j) printf(" %10.4g " , vecM[ dim * j + i ]);
    printf("\n");
  }
  printf("\n");
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , q[i]);

  /* Print INFO */

  fscanf(LCPfile , "%s" , val);

  if (!feof(LCPfile))
  {

    sol[0] = atof(val);

    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }

    printf("\n SOLUTION OF LCP :\n");
    for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
    printf("\n");
  }
  else
  {
    printf("\n NO EXACT SOLUTION FOR LCP :\n");
    for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
  }

  fclose(LCPfile);

  printf("\n NOTE : DET(M) = 16 & COND(M) = 4 \n");
  printf("          M NON SYMETRIC            \n");
  printf("          Ker(M) = nul              \n");

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp1 , z , w);

  printf("***********************\n");
  printf(" With NLGS:            \n");
  printf(" Solution :            \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  /* CPG */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp2 , z , w);

  printf("**********************\n");
  printf(" With CPG:            \n");
  printf(" Solution :           \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  /* Lemke */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp3 , z , w);

  printf("**********************\n");
  printf(" With Lemke:          \n");
  printf(" Solution :           \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  free(q);
  free(sol);
  free(vecM);

  /*************************************************************
    *
    *                       TEST TRIVIAL
    *
    *************************************************************/

  if ((LCPfile = fopen("MATRIX/trivial.dat", "r")) == NULL)
  {
    perror("fopen LCPfile: trivial.dat");
    exit(1);
  }
  printf("\n SCAN TRIVIAL MATRIX FOR LCP \n");

  fscanf(LCPfile , "%d" , &dim);

  //dim = atoi(val);
  dim2 = dim * dim;

  vecM = (double*)malloc(dim2 * sizeof(double));
  q    = (double*)malloc(dim * sizeof(double));
  sol  = (double*)malloc(dim * sizeof(double));

  for (i = 0 ; i < dim ; ++i)  q[i]    = 0.0;
  for (i = 0 ; i < dim ; ++i)  sol[i]  = 0.0;
  for (i = 0 ; i < dim2 ; ++i) vecM[i] = 0.0;

  for (i = 0 ; i < dim ; ++i)
  {
    for (j = 0 ; j < dim ; ++j)
    {
      fscanf(LCPfile, "%s", val);
      vecM[ dim * j + i ] = atof(val);
    }
  }

  for (i = 0 ; i < dim ; ++i)
  {
    fscanf(LCPfile , "%s" , val);
    q[i] = atof(val);
  }

  /* Print INFO */

  for (i = 0 ; i < dim ; ++i)
  {
    for (j = 0 ; j < dim ; ++j) printf(" %10.4g " , vecM[ dim * j + i ]);
    printf("\n");
  }
  printf("\n");
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , q[i]);

  /* Print INFO */

  fscanf(LCPfile , "%s" , val);

  if (!feof(LCPfile))
  {

    sol[0] = atof(val);

    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }

    printf("\n SOLUTION OF LCP :\n");
    for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
    printf("\n");
  }
  else
  {
    printf("\n NO EXACT SOLUTION FOR LCP :\n");
    for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
  }

  fclose(LCPfile);

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp1 , z , w);

  printf("***********************\n");
  printf(" With NLGS:            \n");
  printf(" Solution :            \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  /* CPG */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp2 , z , w);

  printf("**********************\n");
  printf(" With CPG:            \n");
  printf(" Solution :           \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  /* Lemke */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = solve_lcp(vecM , q , &dim , &meth_lcp3 , z , w);

  printf("**********************\n");
  printf(" With Lemke:          \n");
  printf(" Solution :           \n");

  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);

  a1 = -1.;

  daxpy_(&dim , &a1 , sol , &incx , z , &incy);
  RES = dnrm2_(&dim , z , &incx);

  printf("\n ||z-sol|| = %g\n" , RES);

  printf("\n We go out the function and info is %d\n", info);

  free(q);
  free(sol);
  free(vecM);

  return 0;
}




