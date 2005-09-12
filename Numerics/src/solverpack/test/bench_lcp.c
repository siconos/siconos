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
 *        lcp_lexicolemke(M,q,n,itermax,z,w,it_end,res,info)
 *        lcp_gsnl(M,q,n,itermax,tol,z,w,it_end,res,info)
 *        lcp_gcp(M,q,n,itermax,tol,z,w,it_end,res,info)
 *        lcp_qp(M,q,n,res,z,w,fail)
 *
 *  where _ itermax is the maximum iterations required, it's an integer
 *        _ res is the residue, it's a float
 *        _ it_end is the number of iterations carried out, it's an integer
 *        _ tol is the tolerance value desired, it's a float
 *        _ k_latin is the parameter research of the latin, it's a float
 *        _ fail shows the termination reason, it's an integer
 *        _ info shows the trmination reason (0 successful otherwise 1), it's an integer.
 *
 *  The subroutine's call is due to the function lcp_solver:
 *
 *  int lcp_solver (float (*M)[maxcols],float * q, int n, method *pt, float *z, float * w)
 *
 *  where M is an n by n matrix, q an n-dimensional vector, n is the row dimension of M,
 *  and pt a pointer other a structure (method), z and w are n-dimensional vector, the solutions of the lcp.
 *  method is a variable with a structure type; this structure gives to the function lcp_solver,
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

#define BAVARD

int main(void)
{

  FILE *LCPfile;

  int i, j;
  int info ;
  int dim , dim2;
  int dn, db, db2, nblock, itmp;

  int *inb, *iid;
  double *q , *z , *w , *sol;
  double *vecM;

  char val[20];

  int iter, titer;
  double criteria;

  /*
   * method declaration :
   *
   * method 1 for Gsnl
   * method 2 for Cpg
   * method 3 for lexicolemke
   *
   */

  static method_lcp meth_lcp1 = { "NLGS"        , 1000 , 1e-8 , 0.6 , 1.0 , 0 , "N2" };
  static method_lcp meth_lcp2 = { "CPG"         , 1000 , 1e-8 , 0.6 , 1.0 , 0 , "N2" };
  static method_lcp meth_lcp3 = { "LexicoLemke" , 1000 , 1e-8 , 0.6 , 1.0 , 0 , "N2" };
  static method_lcp meth_lcp4 = { "QP"          , 1000 , 1e-8 , 0.6 , 1.0 , 0 , "N2" };
  static method_lcp meth_lcp5 = { "NSQP"        , 1000 , 1e-8 , 0.6 , 1.0 , 0 , "N2" };

  /****************************************************************/

  iter  = 0;
  titer = 0;

  /****************************************************************/
#ifdef BAVARD
  printf("\n ********** BENCHMARK FOR LCP_SOLVER ********** \n\n");
#endif
  /****************************************************************/

  if ((LCPfile = fopen("MATRIX/deudeu.dat", "r")) == NULL)
  {
    perror("fopen LCPfile: deudeu.dat");
    exit(1);
  }

  fscanf(LCPfile , "%d" , &dim);

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

  fscanf(LCPfile , "%s" , val);

  if (!feof(LCPfile))
  {

    sol[0] = atof(val);

    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }
  }
  else
  {
    for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
  }

  fclose(LCPfile);

#ifdef BAVARD
  printf("\n solution deudeu  : ");
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
  printf("\n");
#endif

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp1 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           NLGS   :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* CPG */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp2 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           CPG    :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* Lemke */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp3 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           Lemke  :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* QP */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp4 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           QP     :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  free(q);
  free(w);
  free(z);
  free(sol);
  free(vecM);

  /****************************************************************/
  /****************************************************************/

  if ((LCPfile = fopen("MATRIX/trivial.dat", "r")) == NULL)
  {
    perror("fopen LCPfile: trivial.dat");
    exit(1);
  }

  fscanf(LCPfile , "%d" , &dim);

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

  fscanf(LCPfile , "%s" , val);

  if (!feof(LCPfile))
  {

    sol[0] = atof(val);

    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }
  }
  else
  {
    for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
  }

  fclose(LCPfile);

#ifdef BAVARD
  printf("\n");
  printf("\n solution trivial : ");
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
  printf("\n");
#endif

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp1 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           NLGS   :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* CPG */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp2 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           CPG    :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* Lemke */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp3 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           Lemke  :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* QP */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp4 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           QP     :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  free(q);
  free(w);
  free(z);
  free(sol);
  free(vecM);

  /****************************************************************/
  /****************************************************************/

  if ((LCPfile = fopen("MATRIX/ortiz.dat", "r")) == NULL)
  {
    perror("fopen LCPfile: ortiz.dat");
    exit(1);
  }

  fscanf(LCPfile , "%d" , &dim);

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

  fscanf(LCPfile , "%s" , val);

  if (!feof(LCPfile))
  {

    sol[0] = atof(val);

    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }
  }
  else
  {
    for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
  }

  fclose(LCPfile);

#ifdef BAVARD
  printf("\n");
  printf("\n solution ortiz   : ");
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
  printf("\n");
#endif

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp1 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           NLGS   :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* CPG */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp2 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           CPG    :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* Lemke */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp3 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           Lemke  :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* NSQP */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp5 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           NSQP   :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  free(q);
  free(sol);
  free(w);
  free(z);
  free(vecM);

  /****************************************************************/
  /****************************************************************/

  if ((LCPfile = fopen("MATRIX/pang.dat", "r")) == NULL)
  {
    perror("fopen LCPfile: pang.dat");
    exit(1);
  }

  fscanf(LCPfile , "%d" , &dim);

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

  fscanf(LCPfile , "%s" , val);

  if (!feof(LCPfile))
  {

    sol[0] = atof(val);

    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }
  }
  else
  {
    for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
  }

  fclose(LCPfile);

#ifdef BAVARD
  printf("\n");
  printf("\n solution pang    : ");
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
  printf("\n");
#endif

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp1 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           NLGS   :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* CPG */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp2 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           CPG    :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* Lemke */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp3 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           Lemke  :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* NSQP */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp5 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           NSQP   :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  free(q);
  free(sol);
  free(w);
  free(z);
  free(vecM);

  /****************************************************************/
  /****************************************************************/

  if ((LCPfile = fopen("MATRIX/diodes.dat", "r")) == NULL)
  {
    perror("fopen LCPfile: diodes.dat");
    exit(1);
  }

  fscanf(LCPfile , "%d" , &dim);

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

  fscanf(LCPfile , "%s" , val);

  if (!feof(LCPfile))
  {

    sol[0] = atof(val);

    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }
  }
  else
  {
    for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
  }

  fclose(LCPfile);

#ifdef BAVARD
  printf("\n");
  printf("\n solution diodes  : ");
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
  printf("\n");
#endif

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp1 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           NLGS   :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* CPG */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp2 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           CPG    :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* Lemke */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp3 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           Lemke  :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* NSQP */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp5 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           NSQP   :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  free(q);
  free(sol);
  free(w);
  free(z);
  free(vecM);

  /****************************************************************/
  /****************************************************************/

  if ((LCPfile = fopen("MATRIX/murty1.dat", "r")) == NULL)
  {
    perror("fopen LCPfile: murty1.dat");
    exit(1);
  }

  fscanf(LCPfile , "%d" , &dim);

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

  fscanf(LCPfile , "%s" , val);

  if (!feof(LCPfile))
  {

    sol[0] = atof(val);

    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }
  }
  else
  {
    for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
  }

  fclose(LCPfile);

#ifdef BAVARD
  printf("\n");
  printf("\n solution murty1  : ");
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
  printf("\n");
#endif

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp1 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           NLGS   :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* CPG */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp2 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           CPG    :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* Lemke */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp3 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           Lemke  :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* NSQP */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp5 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           NSQP   :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  free(q);
  free(w);
  free(z);
  free(sol);
  free(vecM);

  /****************************************************************/
  /****************************************************************/

  if ((LCPfile = fopen("MATRIX/confeti.dat", "r")) == NULL)
  {
    perror("fopen LCPfile: confeti.dat");
    exit(1);
  }

  fscanf(LCPfile , "%d" , &dim);

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

  fscanf(LCPfile , "%s" , val);

  if (!feof(LCPfile))
  {

    sol[0] = atof(val);

    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }
  }
  else
  {
    for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
  }

  fclose(LCPfile);

#ifdef BAVARD
  printf("\n");
  printf("\n solution confeti : ");
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
  printf("\n");
#endif

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp1 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           NLGS   :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* CPG */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp2 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           CPG    :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* Lemke */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp3 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           Lemke  :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* NSQP */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp5 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           NSQP   :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  free(q);
  free(sol);
  free(w);
  free(z);
  free(vecM);

  /****************************************************************/
  /****************************************************************/

  if ((LCPfile = fopen("MATRIX/mathieu2.dat", "r")) == NULL)
  {
    perror("fopen LCPfile: mathieu2.dat");
    exit(1);
  }

  fscanf(LCPfile , "%d" , &dim);

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

  fscanf(LCPfile , "%s" , val);

  if (!feof(LCPfile))
  {

    sol[0] = atof(val);

    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }
  }
  else
  {
    for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
  }

  fclose(LCPfile);

#ifdef BAVARD
  printf("\n");
  printf("\n solution mathieu2 : ");
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
  printf("\n");
#endif

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp1 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           NLGS   :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* CPG */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp2 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           CPG    :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* Lemke */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp3 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           Lemke  :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  /* NSQP */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver(vecM , q , &dim , &meth_lcp5 , z , w , &iter , &criteria);

#ifdef BAVARD
  printf("\n           NSQP   :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d ) " , iter);
#endif

  free(q);
  free(sol);
  free(w);
  free(z);
  free(vecM);

  /****************************************************************/
#ifdef BAVARD
  printf("\n\n ******** BENCHMARK FOR LCP_SOLVER_BLOCK ******** \n\n");
#endif
  /****************************************************************/

  if ((LCPfile = fopen("BLOCKMATRIX/mathieu1.dat", "r")) == NULL)
  {
    perror("fopen LCPfile: mathieu1.dat");
    exit(1);
  }

  fscanf(LCPfile , "%d" , &dn);
  fscanf(LCPfile , "%d" , &db);

  nblock = 0;

  inb  = (int*)malloc(dn * sizeof(int));

  for (i = 0 ; i < dn ; ++i)
  {
    fscanf(LCPfile , "%d" , &itmp);
    inb[i] = itmp;
    nblock += itmp;
  }

  iid  = (int*)malloc(nblock * sizeof(int));

  for (i = 0 ; i < nblock ; ++i)
  {
    fscanf(LCPfile , "%d" , &itmp);
    iid[i] = itmp;
  }

  dim  = db * dn;
  db2 = db * db;
  dim2 = nblock * db2;

  q    = (double*)malloc(dim * sizeof(double));
  sol  = (double*)malloc(dim * sizeof(double));
  vecM = (double*)malloc(dim2 * sizeof(double));

  for (i = 0 ; i < dim ; ++i)  q[i]    = 0.0;
  for (i = 0 ; i < dim ; ++i)  sol[i]  = 0.0;
  for (i = 0 ; i < dim2 ; ++i) vecM[i] = 0.0;


  for (i = 0 ; i < nblock ; ++i)
  {
    for (j = 0 ; j < db2 ; ++j)
    {
      fscanf(LCPfile, "%s", val);
      vecM[i * db2 + j] = atof(val);
    }
  }

  for (i = 0 ; i < dim ; ++i)
  {
    fscanf(LCPfile , "%s" , val);
    q[i] = atof(val);
  }

  fscanf(LCPfile , "%s" , val);

  if (!feof(LCPfile))
  {

    sol[0] = atof(val);

    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }
  }
  else
  {
    for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
  }

  fclose(LCPfile);

#ifdef BAVARD
  printf("\n");
  printf("\n solution mathieu1 : ");
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
  printf("\n");
#endif

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver_block(inb , iid , vecM , q , &dn , &db , &meth_lcp1 , z , w , &iter , &titer , &criteria);

#ifdef BAVARD
  printf("\n           NLGS   :%1d ", info);
  for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , z[i]);
  printf(" ( %6d -- %6d ) " , iter , titer);
#endif

  free(q);
  free(w);
  free(z);
  free(sol);
  free(vecM);


  printf(" \n");

  return 0;
}




