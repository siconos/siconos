/*!
 *
 *  This main file allows the primal resolution of contact problems with friction:
 *  try (z,w) such that:
 *
 *
 *                    M z  - w = q                 (1)
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
 *    For more information about the methods see the chapter 4 of the Siconos manual theory.
 *
 *
 *
 *  The subroutine's call is due to the function pfc_2D_solver:
 *
 *  int pfc_2D_solver(double (*M)[maxcols],double *q,int n,method *pt, double *z,double *w)
 *
 *  where M is an n by n matrix, q an n-dimensional vector, n is the row dimension
 *  of M, and pt a pointer other a structure ( method). z and w are n-dimensional
 *  vectors solution.
 *  method is a variable with a structure type; this structure gives to the function
 *  pfc_2D_solver, the name and the parameters (itermax, tol, k_latin,...) of the method we want to use.
 *  This function return an interger:  0 successful return otherwise 1.
 *
 * \author Mathieu Renouf
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "solverpack.h"

#define BAVARD

void main(void)
{

  FILE *f1, *f2;

  int i, nl, nc, n, dimM;

  int info[3];

  double qi, Mij;
  double *q, *vec;

  double *z1, *w1, *z2, *w2, *z3, *w3;

  char val[14];

  for (i = 0 ; i < 3 ; ++i) info[i] = -1;

  /****************************************************************/
#ifdef BAVARD
  printf("\n ********** BENCHMARK FOR PFC 2D SOLVERS ********** \n\n");
#endif
  /****************************************************************/

  /* WARNING STATIC SIZE */

  n = 152;
  dimM = n;

  /* Methods*/

  static method_pfc_2D meth_pfc1  = { "NLGS"  , 1000 , 1e-08 , 0.3 , 0.7 };
  static method_pfc_2D meth_pfc2  = { "Latin" , 1000 , 1e-08 , 0.3 , 35. };
  static method_pfc_2D meth_pfc3  = { "CPG"   , 1000 , 1e-08 , 0.3 , 0.7 };

  if ((f1 = fopen("DATA/MM_gran_mu12.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }

  vec = (double*)malloc(dimM * dimM * sizeof(double));
  q   = (double*)malloc(dimM * sizeof(double));

  z1 = malloc(dimM * sizeof(double));
  w1 = malloc(dimM * sizeof(double));
  z2 = malloc(dimM * sizeof(double));
  w2 = malloc(dimM * sizeof(double));
  z3 = malloc(dimM * sizeof(double));
  w3 = malloc(dimM * sizeof(double));

  while (!feof(f1))
  {

    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);
    vec[(nc - 1)*dimM + (nl - 1) ] = Mij;

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



  for (i = 0 ; i < dimM ; ++i) z1[i] = 0.0;

  info[0] = pfc_2D_solver(vec , q , &n , &meth_pfc1 , z1 , w1);

  for (i = 0 ; i < dimM ; ++i) z2[i] = 0.0;

  info[1] = pfc_2D_solver(vec , q , &n , &meth_pfc2 , z2 , w2);

  for (i = 0 ; i < dimM ; ++i) z3[i] = 0.0;

  info[2] = pfc_2D_solver(vec , q , &n , &meth_pfc3 , z3 , w3);


#ifdef BAVARD
  printf(" *** ************************************** ***\n");
  printf("\n   NLGS RESULT (LOG:%1d)|: " , info[0]);
  for (i = 0 ; i < dimM ; ++i) printf(" %10.4g " , z1[i]);
  printf("\n    CPG RESULT (LOG:%1d)|: " , info[1]);
  for (i = 0 ; i < dimM ; ++i) printf(" %10.4g " , z3[i]);
  printf("\n  LATIN RESULT (LOG:%1d)|: " , info[2]);
  for (i = 0 ; i < dimM ; ++i) printf(" %10.4g " , z2[i]);
  printf("\n\n");
#endif

  free(q);
  free(z1);
  free(w1);
  free(z2);
  free(w2);
  free(z3);
  free(w3);
  free(vec);

}


