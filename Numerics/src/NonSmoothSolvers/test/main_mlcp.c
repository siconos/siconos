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
 *  The solver's call is performed via the function mlcp_driver:
 *
 *  int mlcp_driver(double *A , double *B , double *C , double *D , double *a , double *b, int *n , int* m, method *pt ,  double *u, double *v, double *w  )
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NonSmoothDrivers.h"
#include "LA.h"

#define BAVARD
#define NBTEST 10
#define NBMETHODS 5

#define PATH_DRIVER

#ifdef PATH_DRIVER
const unsigned short int *__ctype_b;
const __int32_t *__ctype_tolower ;
#endif /*PATH_DRIVER*/


typedef struct
{
  char file[124];
  char cv[5][NBMETHODS];
  int nbSteps[NBMETHODS];
} dataSummary;

static dataSummary summary[NBTEST];
static int itest;

/*
 ******************************************************************************
 */
void printSolution(char *name, int n, int m, double *u1, double *v1, double *w1)
{
  int i;
#ifdef BAVARD
  printf(" *** ************************************** ***\n");
  printf(" ****** z = ********************************\n");
  for (i = 0 ; i < n ; i++)
    printf("%s: %14.7e  \n", name, u1[i]);
  for (i = 0 ; i < m ; i++)
    printf("%s: %14.7e  \n", name, v1[i]);
  printf(" ****** w = ********************************\n");
  for (i = 0 ; i < m ; i++)
    printf("%s: %14.7e  \n", name, w1[i]);
#endif
}
void test_mlcp_series(int n , int m, double *A , double *B , double *C , double *D , double *a , double *b, double *sol)
{

  int i;
  int info1 = -1;

  double *u1;
  double *v1;
  double *w1;


  u1 = malloc(n * sizeof(double));
  v1 = malloc(m * sizeof(double));
  w1 = malloc(m * sizeof(double));

  /* Method definition */
  static method_mlcp method_mlcp1 = { "PGS"       , 101 , 1e-8 , 0.6 , 1.0 , 1 , 0 , 0.0 };
  static method_mlcp method_mlcp2 = { "RPGS"       , 101 , 1e-8 , 0.6 , 1.0 , 1 , 0 , 0.0 };
  static method_mlcp method_mlcp3 = { "PSOR"       , 101 , 1e-8 , 2.0 , 1.0 , 1 , 0 , 0.0 };
  static method_mlcp method_mlcp4 = { "RPSOR"       , 101 , 1e-8 , 2.0 , 1.0 , 1 , 0 , 0.0 };
  static method_mlcp method_mlcp5 = { "PATH"       , 0 , 1e-8 , 0.0 , 0.0 , 0 , 0 , 0.0 };

  method myMethod;

  /* #1 PGS TEST */
#ifdef BAVARD
  printf("**** PGS TEST ****\n");
#endif
  for (i = 0 ; i < m ; ++i)
  {
    v1[i] = sol[i + n];
    w1[i] = sol[n + m + i];
  }
  for (i = 0 ; i < n ; ++i)
  {
    u1[i] = sol[i];
  }
  myMethod.mlcp =  method_mlcp1;
  info1 = mlcp_driver(A, B, C, D, a, b, &n , &m, &myMethod , u1 , v1, w1);

  strcpy(summary[itest].cv[0], "CV");
  if (info1 > 0)
    strcpy(summary[itest].cv[0], "NO");
  printSolution("PGS", n, m, u1, v1, w1);



  /* #1 RPGS TEST */
#ifdef BAVARD
  printf("**** RPGS TEST ****\n");
#endif
  for (i = 0 ; i < m ; ++i)
  {
    v1[i] = sol[i + n];
    w1[i] = sol[n + m + i];
  }
  for (i = 0 ; i < n ; ++i)
  {
    u1[i] = sol[i];
  }

  myMethod.mlcp =  method_mlcp2;
  info1 = mlcp_driver(A, B, C, D, a, b, &n , &m, &myMethod , u1 , v1, w1);
  strcpy(summary[itest].cv[1], "CV");
  if (info1 > 0)
    strcpy(summary[itest].cv[1], "NO");
  printSolution("RPGS", n, m, u1, v1, w1);

  /* #1 PSOR TEST */
#ifdef BAVARD
  printf("**** PSOR TEST ****\n");
#endif
  for (i = 0 ; i < m ; ++i)
  {
    v1[i] = sol[i + n];
    w1[i] = sol[n + m + i];
  }
  for (i = 0 ; i < n ; ++i)
  {
    u1[i] = sol[i];
  }

  myMethod.mlcp =  method_mlcp3;
  info1 = mlcp_driver(A, B, C, D, a, b, &n , &m, &myMethod , u1 , v1, w1);
  strcpy(summary[itest].cv[2], "CV");
  if (info1 > 0)
    strcpy(summary[itest].cv[2], "NO");
  printSolution("PSOR", n, m, u1, v1, w1);

  /* #1 RPSOR TEST */
#ifdef BAVARD
  printf("**** RPSOR TEST ****\n");
#endif
  for (i = 0 ; i < m ; ++i)
  {
    v1[i] = sol[i + n];
    w1[i] = sol[n + m + i];
  }
  for (i = 0 ; i < n ; ++i)
  {
    u1[i] = sol[i];
  }

  myMethod.mlcp =  method_mlcp4;
  info1 = mlcp_driver(A, B, C, D, a, b, &n , &m, &myMethod , u1 , v1, w1);
  strcpy(summary[itest].cv[3], "CV");
  if (info1 > 0)
    strcpy(summary[itest].cv[3], "NO");
  printSolution("RPSOR", n, m, u1, v1, w1);
  /* PATH TEST */
#ifdef BAVARD
  printf("**** PATH TEST ****\n");
#endif
  for (i = 0 ; i < m ; ++i)
  {
    v1[i] = sol[i + n];
    w1[i] = sol[n + m + i];
  }
  for (i = 0 ; i < n ; ++i)
  {
    u1[i] = sol[i];
  }

  myMethod.mlcp =  method_mlcp5;
  info1 = mlcp_driver(A, B, C, D, a, b, &n , &m, &myMethod , u1 , v1, w1);
  strcpy(summary[itest].cv[4], "CV");
  if (info1 > 0)
    strcpy(summary[itest].cv[4], "NO");
  printSolution("PATH", n, m, u1, v1, w1);

  free(u1);
  free(v1);
  free(w1);



}

void test_matrix(void)
{

  FILE *MLCPfile;

  int i, j;
  int isol;
  int n , n2;
  int m, m2;
  int withSol = 0;

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
      strcpy(summary[itest].file, "2x2 MLCP");
      if ((MLCPfile = fopen("MATRIX/deudeu_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: deudeu_mlcp.dat");
        exit(1);
      }
      break;
    case 1:
      printf("\n\n PD ");
      strcpy(summary[itest].file, "PD");
      if ((MLCPfile = fopen("MATRIX/PD_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: PD_mlcp.dat");
        exit(1);
      }
      break;
    case 2:
      printf("\n\n m2n1_mlcp.dat ");
      strcpy(summary[itest].file, "m2n1_mlcp");
      if ((MLCPfile = fopen("MATRIX/m2n1_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: m2n1SOL_mlcp.dat");
        exit(1);
      }
      break;
    case 3:
      printf("\n\n m3n2_mlcp.dat ");
      strcpy(summary[itest].file, "m3n2_mlcp");
      if ((MLCPfile = fopen("MATRIX/m3n2_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: m3n2SOL_mlcp.dat");
        exit(1);
      }
      break;
    case 4:
      printf("\n\n PDSym_mlcp.dat ");
      strcpy(summary[itest].file, "PDSym_mlcp");
      if ((MLCPfile = fopen("MATRIX/PDSym_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: PDSym_mlcp.dat");
        exit(1);
      }
      break;
    case 5:
      printf("\n\n diodeBridge MLCP ");
      strcpy(summary[itest].file, "diodeBridge MLCP");
      if ((MLCPfile = fopen("MATRIX/diodeBridge_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: diodeBridge_mlcp.dat");
        exit(1);
      }
      break;
    case 7:
      printf("\n\n diodeBridge 20 MLCP ");
      strcpy(summary[itest].file, "diodeBridge 20 MLCP");
      if ((MLCPfile = fopen("MATRIX/diodeBridge20_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: diodeBridge20_mlcp.dat");
        exit(1);
      }
      break;
    case 9:
      printf("\n\n diodeBridge 40 MLCP ");
      strcpy(summary[itest].file, "diodeBridge 40 MLCP");
      if ((MLCPfile = fopen("MATRIX/diodeBridge40_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: diodeBridge40_mlcp.dat");
        exit(1);
      }
      break;
    case 8:
      printf("\n\n Buck converter MLCP ");
      strcpy(summary[itest].file, "Buck converter MLCP");
      if ((MLCPfile = fopen("MATRIX/buckconverterregul2.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: buckconverterregul2.dat");
        exit(1);
      }
      break;

    case 6:
      printf("\n\n RCD ");
      strcpy(summary[itest].file, "RCD");
      if ((MLCPfile = fopen("MATRIX/RCD_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: RCD_mlcp.dat");
        exit(1);
      }
      break;
    default :
      exit(1);
    }

    fscanf(MLCPfile , "%d" , &n);
    fscanf(MLCPfile , "%d" , &m);

    n2 = n * n;
    m2 = m * m;
    isol = 1;

    vecA = (double*)calloc(n2, sizeof(double));
    vecB = (double*)calloc(m2, sizeof(double));
    vecC = (double*)calloc(n * m, sizeof(double));
    vecD = (double*)calloc(m * n, sizeof(double));
    a    = (double*)calloc(n, sizeof(double));
    b    = (double*)calloc(m, sizeof(double));
    sol  = (double*)calloc((n + m + m), sizeof(double));




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
        vecB[ m * j + i ] = atof(val);


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
        vecD[ m * j + i ] = atof(val);
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
    withSol = 0;
    if (!feof(MLCPfile))
    {
      withSol = 1;

      sol[0] = atof(val);

      for (i = 1 ; i < n + m + m ; ++i)
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
    printf("\n With exact solution : ");
    printf("\n ------------------- : \n");
    if (isol) for (i = 0 ; i < (n + m + m) ; ++i) printf(" %10.4g " , sol[i]);
    else printf(" unknown ");
    printf("\n");
#endif
    if (withSol)
    {
      test_mlcp_series(n, m , vecA, vecB, vecC , vecD, a, b, sol);
      printf("\n Without exact solution : ");
      printf("\n ---------------------- : \n");
    }
    for (i = 0; i < n + m + m; i++)
      sol[i] = 0;
    test_mlcp_series(n, m , vecA, vecB, vecC , vecD, a, b, sol);

    free(sol);
    free(vecA);
    free(vecB);
    free(vecC);
    free(vecD);
    free(a);
    free(b);
  }
  printf("* *** ******************** *** * \n");
  printf("* *** SUMMARY             *** * \n");
  printf("* *** ******************** *** * \n");

  for (itest = 0; itest < NBTEST ; ++itest)
  {
    printf(" test %s  :\n", summary[itest].file);
    printf(" PGS %s  \t", summary[itest].cv[0]);
    printf("| RPGS %s  \t", summary[itest].cv[1]);
    printf("| PSOR %s  \t", summary[itest].cv[2]);
    printf("| RPSOR %s  \t", summary[itest].cv[3]);
    printf("| PATH %s  \n", summary[itest].cv[4]);
  }
  printf("* *** ******************** *** * \n");
  printf("* *** END OF TEST MATRIX   *** * \n");
  printf("* *** ******************** *** * \n");


}

int main(void)
{

  /*   test_matrix(); */

  return 1;
}

