/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
/*!
 ******************************************************************************
 *
 * This subroutine allows the resolution of MLCP (Mixed Linear Complementary Problem).\n
 *
 *
 *    \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *   A u + Cv +a =0\\
 *   D u + Bv +b = w
 *   0 \le v \perp  w \ge 0\\
 *   \end{array}
 *  \right.
 * \f$
 *
 *  This system of equations and inequalities is solved thanks to mlcp solvers\n
 *
 *        mlcp_pgs( n ,m, A , B , C , D , a  , b, u, v, w , &info , iparamMLCP , dparamMLCP )
 *
 *  where info shows the termination result (0 for success) and iparam and dparam are respectivelly
 *  pointer over integer and pointer over double which contain specific parameters of each solver.
 *
 *  The solver's call is performed via the function mlcp_solver:
 *
 *  int mlcp_solver(double *A , double *B , double *C , double *D , double *a , double *b, int *n , int* m, method *pt ,  double *u, double *v, double *w  )
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NSSpack.h"
#include "LA.h"

#define BAVARD

/*
 ******************************************************************************
 */

void test_mlcp_series(int n , int m, double *A , double *B , double *C , double *D , double *a , double *b)
{

  int i, j;
  int nonsymmetric;
  int incx = 1, incy = 1;
  int info1 = -1;

  double comp, diff, alpha, beta;

  double *u1;
  double *v1;
  double *w1;


  char NT = 'N';

  u1 = malloc(n * sizeof(double));
  v1 = malloc(m * sizeof(double));
  w1 = malloc(m * sizeof(double));

  /* Method definition */

  static method_mlcp method_mlcp1 = { "PGS"       , 1001 , 1e-8 , 0.6 , 1.0 , 1 , 0 , 0.0 };

  /* Method definition */

  static method_mlcp method_mlcp2 = { "RPGS"       , 1001 , 1e-8 , 0.6 , 1.0 , 1 , 0 , 0.0 };

  /*   nonsymmetric = 0; */

  /*   /\* Is M symmetric ? *\/ */

  /*   for( i = 0 ; i < n ; ++i ){ */
  /*     for( j = 0 ; j < i ; ++j ){ */
  /*       if( abs( vec[i*n+j] - vec[j*n+i] ) > 1e-16 ){ */
  /*  nonsymmetric = 1; */
  /*  break; */
  /*       } */
  /*     } */
  /*   } */

#ifdef BAVARD
  if (nonsymmetric) printf("\n !! WARNING !!\n M is a non symmetric matrix \n");
  else printf(" M is a symmetric matrix \n");
#endif

  /* #1 PGS TEST */
#ifdef BAVARD
  printf("**** PGS TEST ****\n");
#endif
  for (i = 0 ; i < m ; ++i)
  {
    v1[i] = 0.0;
    w1[i] = 0.0;
  }
  for (i = 0 ; i < n ; ++i)
  {
    u1[i] = 0.0;
  }

  info1 = mlcp_solver(A, B, C, D, a, b, &n , &m, &method_mlcp1 , u1 , v1, w1);


#ifdef BAVARD
  /*  printf(" *** ************************************** ***\n"); */

  printf(" ****** z = ********************************\n");
  for (i = 0 ; i < n ; i++)
    printf("PGS: %14.7e  \n", u1[i]);

  for (i = 0 ; i < m ; i++)
    printf("PGS: %14.7e  \n", v1[i]);

  printf(" ****** w = ********************************\n");
  for (i = 0 ; i < m ; i++)
    printf("PGS: %14.7e  \n", w1[i]);


  /*   nonsymmetric = 0; */

  /*   /\* Is M symmetric ? *\/ */

  /*   for( i = 0 ; i < n ; ++i ){ */
  /*     for( j = 0 ; j < i ; ++j ){ */
  /*       if( abs( vec[i*n+j] - vec[j*n+i] ) > 1e-16 ){ */
  /*  nonsymmetric = 1; */
  /*  break; */
  /*       } */
  /*     } */
  /*   } */

#ifdef BAVARD
  if (nonsymmetric) printf("\n !! WARNING !!\n M is a non symmetric matrix \n");
  else printf(" M is a symmetric matrix \n");
#endif


  /* #1 RPGS TEST */
#ifdef BAVARD
  printf("**** RPGS TEST ****\n");
#endif
  for (i = 0 ; i < m ; ++i)
  {
    v1[i] = 0.0;
    w1[i] = 0.0;
  }
  for (i = 0 ; i < n ; ++i)
  {
    u1[i] = 0.0;
  }

  info1 = mlcp_solver(A, B, C, D, a, b, &n , &m, &method_mlcp2 , u1 , v1, w1);


#ifdef BAVARD
  /*  printf(" *** ************************************** ***\n"); */

  printf(" ****** z = ********************************\n");
  for (i = 0 ; i < n ; i++)
    printf("RPGS: %14.7e  \n", u1[i]);

  for (i = 0 ; i < m ; i++)
    printf("RPGS: %14.7e  \n", v1[i]);

  printf(" ****** w = ********************************\n");
  for (i = 0 ; i < m ; i++)
    printf("RPGS: %14.7e  \n", w1[i]);

#endif
#endif

  free(u1);
  free(v1);
  free(w1);



}

void test_matrix(void)
{

  FILE *MLCPfile;

  int i, j, itest, NBTEST;
  int isol;
  int dim , dim2;
  int n , n2;
  int m, m2;

  double *a, *b, *sol;
  double *vecA, *vecB, *vecC, *vecD;

  char val[20];

  int iter;
  double criteria;

  printf("* *** ******************** *** * \n");
  printf("* *** STARTING TEST MATRIX *** * \n");
  printf("* *** ******************** *** * \n");


  iter  = 0;
  criteria = 0.0;

  NBTEST = 1;

  /****************************************************************/
#ifdef BAVARD
  printf("\n ********** BENCHMARK FOR MLCP_SOLVER ********** \n\n");
#endif
  /****************************************************************/

  for (itest = 0 ; itest < NBTEST ; ++itest)
  {

    switch (itest)
    {
    case 0:
      printf("\n\n 2x2 MLCP ");
      if ((MLCPfile = fopen("MATRIX/deudeu_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: deudeu.dat");
        exit(1);
      }
      break;
    case 1:
      printf("\n\n DIAGONAL MLCP ");
      if ((MLCPfile = fopen("MATRIX/trivial.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: trivial.dat");
        exit(1);
      }
      break;
    }

    fscanf(MLCPfile , "%d" , &n);
    fscanf(MLCPfile , "%d" , &m);

    n2 = n * n;
    m2 = m * m;
    isol = 1;

    vecA = (double*)malloc(n2 * sizeof(double));
    vecB = (double*)malloc(m2 * sizeof(double));
    vecC = (double*)malloc(n * m * sizeof(double));
    vecD = (double*)malloc(m * m * sizeof(double));
    a    = (double*)malloc(n * sizeof(double));
    b    = (double*)malloc(m * sizeof(double));
    sol  = (double*)malloc((n + m) * sizeof(double));

    for (i = 0 ; i < n ; ++i)  a[i]    = 0.0;
    for (i = 0 ; i < n ; ++i)  b[i]    = 0.0;

    for (i = 0 ; i < n2 ; ++i) vecA[i] = 0.0;
    for (i = 0 ; i < m2 ; ++i) vecB[i] = 0.0;
    for (i = 0 ; i < n * m ; ++i) vecC[i] = 0.0;
    for (i = 0 ; i < m * n ; ++i) vecD[i] = 0.0;



    for (i = 0 ; i < m ; ++i)  sol[i]  = 0.0;

    for (i = 0 ; i < n ; ++i)
    {
      for (j = 0 ; j < n ; ++j)
      {
        fscanf(MLCPfile, "%s", val);
        vecA[ n * j + i ] = atof(val);
      }
    }
    for (i = 0 ; i < m ; ++i)
    {
      for (j = 0 ; j < m ; ++j)
      {
        fscanf(MLCPfile, "%s", val);
        vecB[ n * j + i ] = atof(val);
      }
    }
    for (i = 0 ; i < n ; ++i)
    {
      for (j = 0 ; j < m ; ++j)
      {
        fscanf(MLCPfile, "%s", val);
        vecC[ n * j + i ] = atof(val);
      }
    }
    for (i = 0 ; i < m ; ++i)
    {
      for (j = 0 ; j < n ; ++j)
      {
        fscanf(MLCPfile, "%s", val);
        vecD[ n * j + i ] = atof(val);
      }
    }

    for (i = 0 ; i < n ; ++i)
    {
      fscanf(MLCPfile , "%s" , val);
      a[i] = atof(val);
    }
    for (i = 0 ; i < m ; ++i)
    {
      fscanf(MLCPfile , "%s" , val);
      b[i] = atof(val);
    }

    fscanf(MLCPfile , "%s" , val);

    if (!feof(MLCPfile))
    {

      sol[0] = atof(val);

      for (i = 1 ; i < dim ; ++i)
      {
        fscanf(MLCPfile , "%s" , val);
        sol[i] = atof(val);
      }
    }
    else
    {
      isol = 0;
      for (i = 0 ; i < (n + m) ; ++i) sol[i] = 0.0;
    }

    fclose(MLCPfile);

#ifdef BAVARD
    printf("\n exact solution : ");
    if (isol) for (i = 0 ; i < (n + m) ; ++i) printf(" %10.4g " , sol[i]);
    else printf(" unknown ");
    printf("\n");
#endif

    test_mlcp_series(n, m , vecA, vecB, vecC , vecD, a, b);

    free(sol);
    free(vecA);
    free(vecB);
    free(vecC);
    free(vecD);
    free(a);
    free(b);
  }
  printf("* *** ******************** *** * \n");
  printf("* *** END OF TEST MATRIX   *** * \n");
  printf("* *** ******************** *** * \n");


}

int main(void)
{


  test_matrix();

  return 0;
}

