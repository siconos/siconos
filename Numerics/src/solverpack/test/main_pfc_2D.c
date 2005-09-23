/*!
 *
 *  This main file allows the primal resolution of contact problems with friction:
 *  try (z,w) such that:
 *
 *
 *                    M z  + q = w                 (1)
 *                    0<= zn , 0<= wn , zn.wn= O   (2) (Signorini contact law)
 *                    -wt in diff(PSI)   (zt)      (3) (Coulomb friction law)
 *                                 [-mu*zn, mu*zn]
 *
 *
 *  here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional
 *  vector and w an n-dimensional vector.
 *
 *  This system of equations and inequalities is solved thanks to pfc_2D_subroutine:
 *        pfc_2D_nlgs (M,q,n,mu,itermax,tol,z,w,it_end,res,info)
 *        pfc_2D_cpg (M,q,n,mu,itermax,tol,z,w,it_end,res,info)
 *        pfc_2D_latin (M,q,n,k_latin,mu,itermax,tol,z,w,it_end,res,info)
 *
 *  where _ itermax is the maximum iterations required, it's an integer
 *        _ mu is the friction coefficient, it's a float
 *        _ res is the residue, it's a float
 *        _ it_end is the number of iterations carried out, it's an integer
 *        _ tol is the tolerance value desired, it's a float
 *        _ k_latin is the parameter research of the latin, it's a float
 *        _ z and w are the solutions of the problem
 *        _ info shows the termination reason,0 is successful otherwise 1, it's an integer.
 *
 *
 * \author Mathieu Renouf
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "solverpack.h"
#include "blaslapack.h"

#define BAVARD

void pfc_2D_series(int n , double *vec , double *q)
{

  int i, j;
  int nonsymmetric;
  int info[3];
  int n2;
  int incx = 1, incy = 1;

  double comp;

  double *z1, *w1, *z2, *w2, *z3, *w3;
  printf(" START SERIES \n");
  for (i = 0 ; i < 3 ; ++i) info[i] = -1;

  /* Methods*/

  static method_pfc meth_pfc1  = { "NLGS"  , 1000 , 1e-08 , 0.3 , 0.7 , 1 , "N2" , 0 , 0.0 };
  static method_pfc meth_pfc2  = { "CPG"   , 1000 , 1e-08 , 0.3 , 0.7 , 1 , "N2" , 0 , 0.0 };
  static method_pfc meth_pfc3  = { "Latin" , 1000 , 1e-08 , 0.3 , 35. , 1 , "N2" , 0 , 0.0 };

  z1 = malloc(n * sizeof(double));
  w1 = malloc(n * sizeof(double));
  z2 = malloc(n * sizeof(double));
  w2 = malloc(n * sizeof(double));
  z3 = malloc(n * sizeof(double));
  w3 = malloc(n * sizeof(double));

  nonsymmetric = 0;
  n2 = n / 2;

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

#ifdef BAVARD
  if (nonsymmetric) printf("\n !! WARNING !!\n M is a non symmetric matrix \n");
  else printf(" M is a symmetric matrix \n");
#endif

  /* #1 NLGS TEST */
#ifdef BAVARD
  printf("**** NLGS TEST ****\n");
#endif
  for (i = 0 ; i < n ; ++i) z1[i] = 0.0;

  info[0] = pfc_2D_solver(vec , q , &n2 , &meth_pfc1 , z1 , w1);

  /* #1 NLGS TEST */
#ifdef BAVARD
  printf("**** CPG TEST *****\n");
#endif
  for (i = 0 ; i < n ; ++i) z2[i] = 0.0;

  info[1] = pfc_2D_solver(vec , q , &n2 , &meth_pfc2 , z2 , w2);
  /* #1 NLGS TEST */
#ifdef BAVARD
  printf("**** Latin TEST ***\n");
#endif
  for (i = 0 ; i < n ; ++i) z3[i] = 0.0;

  //  info[2] = pfc_2D_solver( vec , q , &n2 , &meth_pfc3 , z3 , w3 );


#ifdef BAVARD
  printf(" *** ************************************** ***\n");
  printf("\n   NLGS RESULT : ");
  for (i = 0 ; i < n ; ++i) printf(" %10.4g " , z1[i]);
  printf("\n    CPG RESULT : ");
  for (i = 0 ; i < n ; ++i) printf(" %10.4g " , z2[i]);
  printf("\n  LATIN RESULT : ");
  for (i = 0 ; i < n ; ++i) printf(" %10.4g " , z3[i]);
  printf("\n\n");

  printf(" -- SOLVEUR ------ ITER/PIVOT ----- ERR ----- COMP");

  comp = ddot_(&n , z1 , &incx , w1 , &incy);

  printf("\n    NLGS   (LOG:%1d)|      %5d | %10.4g | %10.4g |" , info[0] , meth_pfc1.iter , meth_pfc1.err , comp);

  comp = ddot_(&n , z2 , &incx , w2 , &incy);

  printf("\n    CPG    (LOG:%1d)|      %5d | %10.4g | %10.4g |" , info[1] , meth_pfc2.iter , meth_pfc2.err , comp);

  comp = ddot_(&n , z3 , &incx , w3 , &incy);

  printf("\n    LATIN  (LOG:%1d)|      %5d | %10.4g | %10.4g |" , info[2] , meth_pfc3.iter , meth_pfc3.err , comp);

#endif

  free(z1);
  free(w1);
  free(z2);
  free(w2);
  free(z3);
  free(w3);

}

/* ****************************************************** */
/*          READ MATRICES OF MATRIX DIRECTORY             */
/* ****************************************************** */

void test_matrix(void)
{

  FILE *LCPfile;

  int i, j, itest, NBTEST;
  int isol;
  int dim , dim2;

  double *q, *sol;
  double *vecM;

  char val[20];

  NBTEST = 2;

  for (itest = 0 ; itest < NBTEST ; ++itest)
  {

    switch (itest)
    {
    case 0:
      printf("\n\n 2x2 LCP \n");
      if ((LCPfile = fopen("MATRIX/deudeu.dat", "r")) == NULL)
      {
        perror("fopen LCPfile: deudeu.dat");
        exit(1);
      }
      break;
    case 1:
      printf("\n\n DIAGONAL LCP \n");
      if ((LCPfile = fopen("MATRIX/trivial2.dat", "r")) == NULL)
      {
        perror("fopen LCPfile: trivial2.dat");
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

    pfc_2D_series(dim , vecM , q);

    free(sol);
    free(vecM);
    free(q);
  }
}

/* ****************************************************** *
 *          READ MATRICES OF DATA DIRECTORY               *
 * ****************************************************** */

void test_data(void)
{

  FILE *f1, *f2;

  int nl, nc, n;

  double qi, Mij;
  double *q, *vec;

  char val[14];

  /* WARNING STATIC SIZE */

  n = 152;

  printf("\n\n GRANUL TEST \n");
  if ((f1 = fopen("DATA/MM_gran_mu12.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }

  vec = (double*)malloc(n * n * sizeof(double));
  q   = (double*)malloc(n * sizeof(double));

  while (!feof(f1))
  {

    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);
    vec[(nc - 1)*n + (nl - 1) ] = Mij;

  }

  if ((f2 = fopen("DATA/qq_gran_mu12.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }

  while (!feof(f2))
  {
    fscanf(f2, "%d", &nl);
    fscanf(f2, "%s", val);
    qi = atof(val);
    q[nl - 1] = qi;
  }

  fclose(f2);
  fclose(f1);

  pfc_2D_series(n , vec , q);

  free(vec);
  free(q);

}

int main(void)
{

  /****************************************************************/
#ifdef BAVARD
  printf("\n ********** BENCHMARK FOR PFC 2D SOLVERS ********** \n\n");
#endif
  /****************************************************************/

  test_data();
  test_matrix();

}
