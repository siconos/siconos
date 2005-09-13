/*!
 ******************************************************************************
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).\n
 *
 * Try \f$(z,w)\f$ such that:\n
 * \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *    M z + q= w\\
 *    0 \le z \perp w \ge 0\\
 *   \end{array}
 *  \right.
 * \f$
 *
 * where M is an (n x n)-matrix, q , w and z n-vectors.\n
 *
 *  This system of equations and inequalities is solved thanks to lcp solvers\n
 *
 *        lemke_lcp (M,q,n,itermax,z,w,it_end,res,info)

 *        lcp_nlgs( n , M , q , z , w , info , iparam , dparam )
 *        lcp_cpg ( n , M , q , z , w , info , iparam , dparam )
 *        lcp_lexicolemke( n , M , q , z , w , info , iparam , dparam )
 *        lcp_qp( n , M , q , z , w , info , iparam , dparam )
 *        lcp_nsqp( n , M , q , z , w , info , iparam , dparam )
 *
 *  where info shows the termination result (0 for success) and iparam and dparam are respectivelly
 *  pointer over integer and pointer over double which contain specific parameters of each solver.
 *
 *  The solver's call is performed via the function lcp_solver:
 *
 *  int lcp_solver( double *vec , double *q , int *nn , method *pt , double *z , double *w , int *it_end , double *res )
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "solverpack.h"
#include "blaslapack.h"

#define BAVARD

void test_lcp_series(int n , double *vec , double *q)
{

  int i, j, info, iter;
  int nonsymmetric;
  int incx = 1, incy = 1;

  double criteria, comp;
  double *z, *w;

  z = malloc(n * sizeof(double));
  w = malloc(n * sizeof(double));

  /* Method definition */

  static method_lcp method_lcp1 = { "NLGS"       , 1001 , 1e-8 , 0.6 , 1.0 , 0 , "N2"};
  static method_lcp method_lcp2 = { "CPG"        , 1000 , 1e-8 , 0.6 , 1.0 , 0 , "N2"};
  static method_lcp method_lcp3 = { "Latin"      , 1000 , 1e-8 , 0.7 , 1.0 , 0 , "N2"};
  static method_lcp method_lcp4 = { "Lemke"      , 1000 , 1e-8 , 0.7 , 1.0 , 0 , "N2"};
  static method_lcp method_lcp5 = { "QP"         , 1000 , 1e-8 , 0.7 , 1.0 , 0 , "N2"};
  static method_lcp method_lcp6 = { "NSQP"       , 1000 , 1e-8 , 0.7 , 1.0 , 0 , "N2"};
  static method_lcp method_lcp7 = { "LexicoLemke", 1000 , 1e-8 , 0.7 , 1.0 , 0 , "N2"};

  info     = -1;
  iter     = 0;
  criteria = 0.0;
  comp     = 0.0;

  nonsymmetric = 0;

  /* Is M symmetric ? */

  for (i = 0 ; i < n ; ++i)
  {
    for (j = 0 ; j < i ; ++j)
    {
      if (abs(vec[i * n + j] - vec[j * n + i]) > 1e-16)
      {
        nonsymmetric = 1;
        break;
      }
    }
  }

  if (nonsymmetric) printf("\n !! WARNING !!\n M is a non symmetric matrix \n");
  else printf(" M is a symmetric matrix \n");

  /* #1 NLGS TEST */

  printf("\n ***** NLGS TEST ************* \n");
  for (i = 0 ; i < n ; ++i) z[i] = 0.0;

  info = lcp_solver(vec , q , &n , &method_lcp1 , z , w , &iter , &criteria);

  printf("\n NLGS LOG         : %d ( ITER/PIVOT= %d - RES= %g )\n", info, iter, criteria);
  printf(" SOLUTION: ");
  for (i = 0 ; i < n ; ++i) printf(" %10.4g " , z[i]);
  comp = ddot_(&n , z , &incx , w , &incy);
  printf("\n COMPLEMENTARITY : %g \n", comp);

  /* #2 CPG TEST */

  printf("\n ***** CPG TEST ************** \n");
  for (i = 0 ; i < n ; ++i) z[i] = 0.0;

  info = lcp_solver(vec , q , &n , &method_lcp2 , z , w , &iter , &criteria);

  printf("\n CPG LOG          : %d ( ITER/PIVOT= %d - RES= %g )\n", info, iter, criteria);
  printf(" SOLUTION: ");
  for (i = 0 ; i < n ; ++i) printf(" %10.4g " , z[i]);
  comp = ddot_(&n , z , &incx , w , &incy);
  printf("\n COMPLEMENTARITY : %g \n", comp);

  /* #3 QP TEST */

  printf("\n ***** QP TEST *************** \n");
  for (i = 0 ; i < n ; ++i) z[i] = 0.0;

  info = lcp_solver(vec , q , &n , &method_lcp5 , z , w , &iter , &criteria);

  printf("\n QP LOG           : %d ( ITER/PIVOT= %d - RES= %g )\n", info, iter, criteria);
  printf(" SOLUTION: ");
  for (i = 0 ; i < n ; ++i) printf(" %10.4g " , z[i]);
  comp = ddot_(&n , z , &incx , w , &incy);
  printf("\n COMPLEMENTARITY : %g \n", comp);

  /* #4 NSQP TEST */
  printf("\n ***** NSQP TEST ************* \n");
  for (i = 0 ; i < n ; ++i) z[i] = 0.0;

  info = lcp_solver(vec , q , &n , &method_lcp6 , z , w , &iter , &criteria);

  printf("\n NSQP LOG        : %d ( ITER/PIVOT= %d - RES= %g )\n", info, iter, criteria);
  printf(" SOLUTION: ");
  for (i = 0 ; i < n ; ++i) printf(" %10.4g " , z[i]);
  comp = ddot_(&n , z , &incx , w , &incy);
  printf("\n COMPLEMENTARITY : %g \n", comp);

  /* #5 LEXICO LEMKE TEST */
  printf("\n ***** LEXICO LEMKE TEST ***** \n");
  for (i = 0 ; i < n ; ++i) z[i] = 0.0;

  info = lcp_solver(vec , q , &n , &method_lcp7 , z , w , &iter , &criteria);

  printf("\n LEXICO LEMKE LOG: %d ( ITER/PIVOT= %d - RES= %g )\n", info, iter, criteria);
  printf(" SOLUTION: ");
  for (i = 0 ; i < n ; ++i) printf(" %10.4g " , z[i]);
  comp = ddot_(&n , z , &incx , w , &incy);
  printf("\n COMPLEMENTARITY : %g \n", comp);

  /* #6 LATIN TEST */
  printf("\n ***** LATIN TEST ************ \n");
  for (i = 0 ; i < n ; ++i) z[i] = 0.0;

  info = lcp_solver(vec , q , &n , &method_lcp3 , z , w , &iter , &criteria);

  printf("\n LATIN LOG        : %d ( ITER/PIVOT= %d - RES= %g )\n", info, iter, criteria);
  printf(" SOLUTION: ");
  for (i = 0 ; i < n ; ++i) printf(" %10.4g " , z[i]);
  comp = ddot_(&n , z , &incx , w , &incy);
  printf("\n COMPLEMENTARITY : %g \n", comp);

  /* #7 LEMKE TEST */
  printf("\n ***** LEMKE TEST ************ \n");
  for (i = 0 ; i < n ; ++i) z[i] = 0.0;

  info = lcp_solver(vec , q , &n , &method_lcp4 , z , w , &iter , &criteria);

  printf("\n LEMKE LOG        : %d ( ITER/PIVOT= %d - RES= %g )\n", info, iter, criteria);
  printf(" SOLUTION: ");
  for (i = 0 ; i < n ; ++i) printf(" %10.4g " , z[i]);
  comp = ddot_(&n , z , &incx , w , &incy);
  printf("\n COMPLEMENTARITY : %g \n", comp);

  free(z);
  free(w);

}

void test_mmc(void)
{

  FILE *f1, *f2;

  int i, nl, nc, n;

  double *q, *vecM;
  double qi, Mij;

  char val[14], vall[14];

  printf("* *** ******************** *** * \n");
  printf("* ***        TEST MMC      *** * \n");
  printf("* *** ******************** *** * \n");

  // computation of the size of the problem

  n = 0;

  if ((f1 = fopen("DATA/MM_mmc_mu2.dat", "r")) == NULL)
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

  printf("\n SIZE OF THE PROBLEM : %d \n", n);
  fclose(f1);

  vecM = (double *)malloc(n * n * sizeof(double));

  /* Data loading of M and q */

  if ((f1 = fopen("DATA/MM_mmc_mu2.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }

  if ((f2 = fopen("DATA/qq_mmc_mu2.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }

  for (i = 0 ; i < n * n ; ++i) vecM[i] = 0.0;

  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);

    Mij = atof(val);
    vecM[(nc - 1)*n + nl - 1 ] = Mij;
  }


  q = (double *)malloc(n * sizeof(double));

  for (i = 0 ; i < n ; ++i) q[i] = 0.0;

  while (!feof(f2))
  {
    fscanf(f2, "%d", &nl);
    fscanf(f2, "%s", vall);
    qi = atof(vall);

    q[nl - 1] = qi;
  }

  fclose(f2);
  fclose(f1);

  test_lcp_series(n , vecM , q);

  free(vecM);
  free(q);

}

void test_matrix(void)
{

  FILE *LCPfile;

  int i, j, itest, NBTEST;
  int isol;
  int dim , dim2;

  double *q, *sol;
  double *vecM;

  char val[20];

  int iter;
  double criteria;

  iter  = 0;
  criteria = 0.0;

  NBTEST = 8;

  /****************************************************************/
#ifdef BAVARD
  printf("\n ********** BENCHMARK FOR LCP_SOLVER ********** \n\n");
#endif
  /****************************************************************/

  for (itest = 0 ; itest < NBTEST ; ++itest)
  {

    switch (itest)
    {
    case 0:
      printf("\n 2x2 LCP ");
      if ((LCPfile = fopen("MATRIX/deudeu.dat", "r")) == NULL)
      {
        perror("fopen LCPfile: deudeu.dat");
        exit(1);
      }
      break;
    case 1:
      printf("\n DIAGONAL LCP ");
      if ((LCPfile = fopen("MATRIX/trivial.dat", "r")) == NULL)
      {
        perror("fopen LCPfile: trivial.dat");
        exit(1);
      }
      break;
    case 2:
      printf("\n ORTIZ LCP ");
      if ((LCPfile = fopen("MATRIX/ortiz.dat", "r")) == NULL)
      {
        perror("fopen LCPfile: ortiz.dat");
        exit(1);
      }
      break;
    case 3:
      printf("\n PANG LCP ");
      if ((LCPfile = fopen("MATRIX/pang.dat", "r")) == NULL)
      {
        perror("fopen LCPfile: pang.dat");
        exit(1);
      }
      break;
    case 4:
      printf("\n DIODES RLC LCP ");
      if ((LCPfile = fopen("MATRIX/diodes.dat", "r")) == NULL)
      {
        perror("fopen LCPfile: diodes.dat");
        exit(1);
      }
      break;
    case 5:
      printf("\n MURTY1 LCP ");
      if ((LCPfile = fopen("MATRIX/murty1.dat", "r")) == NULL)
      {
        perror("fopen LCPfile: murty1.dat");
        exit(1);
      }
      break;
    case 6:
      printf("\n APPROXIMATED FRICTION CONE LCP ");
      if ((LCPfile = fopen("MATRIX/confeti.dat", "r")) == NULL)
      {
        perror("fopen LCPfile: confeti.dat");
        exit(1);
      }
      break;
    case 7:
      printf("\n MATHIEU2 LCP ");
      if ((LCPfile = fopen("MATRIX/mathieu2.dat", "r")) == NULL)
      {
        perror("fopen LCPfile: mathieu2.dat");
        exit(1);
      }
      break;
    }

    fscanf(LCPfile , "%d" , &dim);

    dim2 = dim * dim;
    isol = 1;

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
      isol = 0;
      for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
    }

    fclose(LCPfile);

#ifdef BAVARD
    printf("\n exact solution : ");
    if (isol) for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
    else printf(" unknown ");
    printf("\n");
#endif

    test_lcp_series(dim , vecM , q);

    free(vecM);
    free(q);
  }
}

void test_blockmatrix(void)
{

  FILE *LCPfile;

  int i, j;
  int info, isol;
  int dim , dim2;
  int dn, db, db2, nblock, itmp;

  int *inb, *iid;
  double *q , *z , *w , *sol;
  double *vecM;

  char val[20];

  int iter, titer;
  double criteria;

  iter  = 0;
  titer = 0;
  criteria = 0.0;

  /****************************************************************/
#ifdef BAVARD
  printf("\n\n ******** BENCHMARK FOR LCP_SOLVER_BLOCK ******** \n\n");
#endif
  /****************************************************************/

  /* Method definition */

  static method_lcp method_lcp1 = { "NLGS"       , 1001 , 1e-8 , 0.6 , 1.0 , 0 , "N2"};
  /*   static method_lcp method_lcp2 = { "CPG"        , 1000 , 1e-8 , 0.6 , 1.0 , 0 , "N2"}; */
  /*   static method_lcp method_lcp3 = { "Latin"      , 1000 , 1e-8 , 0.7 , 1.0 , 0 , "N2"}; */
  /*   static method_lcp method_lcp4 = { "Lemke"      , 1000 , 1e-8 , 0.7 , 1.0 , 0 , "N2"}; */
  /*   static method_lcp method_lcp5 = { "QP"         , 1000 , 1e-8 , 0.7 , 1.0 , 0 , "N2"}; */
  /*   static method_lcp method_lcp6 = { "NSQP"       , 1000 , 1e-8 , 0.7 , 1.0 , 0 , "N2"}; */
  /*   static method_lcp method_lcp7 = { "LexicoLemke", 1000 , 1e-8 , 0.7 , 1.0 , 0 , "N2"}; */


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
  printf("\n exact solution : ");
  if (isol) for (i = 0 ; i < dim ; ++i) printf(" %10.4g " , sol[i]);
  else printf(" unknown ");
  printf("\n");
#endif

  z = malloc(dim * sizeof(double));
  w = malloc(dim * sizeof(double));

  /* NLGS */

  for (i = 0 ; i < dim ; ++i) z[i] = 0.0;

  info = lcp_solver_block(inb , iid , vecM , q , &dn , &db , &method_lcp1 , z , w , &iter , &titer , &criteria);

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

}

int main(void)
{

  //test_mmc();
  test_matrix();
  test_blockmatrix();

  return 0;
}

