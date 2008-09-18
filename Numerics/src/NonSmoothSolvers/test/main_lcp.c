/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NonSmoothDrivers.h"
#include "LA.h"
#include "blaslapack.h" /* for dlamch */

//#define PATH_DRIVER

#ifdef PATH_DRIVER
const unsigned short int *__ctype_b;
const __int32_t *__ctype_tolower ;
#endif /*PATH_DRIVER*/
/*
******************************************************************************
*/

int test_lcp_series(LinearComplementarity_Problem * problem, int* solversList)
{
  /*
     solversList[i] = 1 => the corresponding solver will be applied to the input problem
     solversList[i] = 0 => solver ignored.

     0: PGS
     1: RPGS
     2: CPG
     3: Lemke
     4: Latin
     5: Latin_w
     6: path
     7: QP
     8: NSQP
     9: Newton_Min
     10: Newton FB
     11: Enumeratif

     MIND TO CHANGE totalNBSolver if you add a new solver in the list

  */

  int totalNBSolver = 12;
  int nbSolvers = 0; /* Real number of solvers called in the tests */
  int i, j;
  for (i = 0; i < totalNBSolver; ++i)
  {
    if (solversList[i] == 1)
      nbSolvers++;
  }

  int nonsymmetric;
  int incx = 1, incy = 1;
  double comp, diff;

  int n = problem->size;
  double **z = malloc(nbSolvers * sizeof(*z));
  double **w = malloc(nbSolvers * sizeof(*w));

  for (i = 0; i < nbSolvers; i++)
  {
    z[i] = malloc(n * sizeof(double));
    w[i] = malloc(n * sizeof(double));

    for (j = 0; j < n; ++j)
    {
      z[i][j] = 0.0;
      w[i][j] = 0.0;
    }
  }
  char nameList[300] = "";

  /* Buffer for w, to check if w = Mz+q  */
  double *wBuffer = malloc(n * sizeof(double));

  Numerics_Options global_options;
  global_options.verboseMode = 0;

  nonsymmetric = 0;
  int info = 0;
  double alpha = -1, beta  = 1;
  int maxIter = 1001;
  double tolerance = 1e-8;

  Solver_Options * options ;
  Solver_Options * local_options = NULL;
  int numberOfSolvers ;
  int isSparse = 0;
  /* Dense storage */
  if (problem->M->storageType == 0)
  {
    printf("\n\n  ");
    printf("The matrix of the LCP is dense (ie double* storage) ");
    printf("\n\n  ");
    numberOfSolvers = 1;
    options = malloc(numberOfSolvers * sizeof(*options));
    local_options = options;
    /* Is M symmetric ? */
    for (i = 0 ; i < n ; ++i)
    {
      for (j = 0 ; j < i ; ++j)
      {
        if (abs(problem->M->matrix0[i * n + j] - problem->M->matrix0[j * n + i]) > 1e-16)
        {
          nonsymmetric = 1;
          break;
        }
      }
    }

    if (nonsymmetric) printf("\n !! WARNING !!\n M is a non symmetric matrix \n");
    else printf(" M is a symmetric matrix \n");

  }
  /* Sparse Block storage */
  else
  {
    printf("\n\n  ");
    printf("The matrix of the LCP is a SparseBlockStructuredMatrix.");
    printf("\n\n  ");
    isSparse = 1;
    numberOfSolvers = 2;
    options = malloc(numberOfSolvers * sizeof(*options));
    strcpy(options[0].solverName, "GaussSeidel_SBM");
    int iparam[3] = {maxIter, 0, 0};
    double dparam[3] = {tolerance, 0.0, 0.0};
    options[0].iSize = 3;
    options[0].dSize = 3;
    options[0].iparam = iparam;
    options[0].dparam = dparam;
    options[0].isSet = 1;
    options[0].filterOn = 0;
    local_options = &options[1];
  }

  printf("\n\n  ");
  printf("The following solvers are called:");
  printf("\n\n  ");
  printf("      SOLVER     | ITER/PIVOT |   ERROR    |     w.z     | ||w-Mz-q|| |");
  printf("\n\n  ");

  /* Current solver number */
  int k = 0;

  /* PGS */
  if (solversList[0] == 1)
  {
    strcat(nameList, "    PGS     |");
    strcpy(local_options->solverName, "PGS");
    int iparam[2] = {maxIter, 0};
    double dparam[2] = {tolerance, 0.0};
    local_options->iSize = 2;
    local_options->dSize = 2;
    local_options->iparam = iparam;
    local_options->dparam = dparam;
    local_options->isSet = 1;
    local_options->filterOn = 0;
    int info1 = lcp_driver(problem, z[k] , w[k], options, numberOfSolvers, &global_options);
    comp = DDOT(n , z[k] , incx , w[k] , incy);
    DCOPY(n , w[k], incx, wBuffer , incy);
    DAXPY(n , alpha , problem->q , incx , wBuffer , incy);
    prod(n, n, beta, problem->M, z[k], alpha, wBuffer);
    diff = DNRM2(n , wBuffer  , incx);

    if (isSparse == 0)
      printf("  PGS     (LOG:%1d)|      %5d | %10.4g | %10.4g | %10.4g |\n", info1, local_options->iparam[1], local_options->dparam[1], comp, diff);
    else
      printf("Gauss-Seidel/PGS     (LOG:%1d)|      %5d | %10.4g | %10.4g | %10.4g |\n", info1, options[0].iparam[1], options[0].dparam[1], comp, diff);

    if (info1 != 0)
      info = info1;
    k++;
  }
  /* RPGS */
  if (solversList[1] == 1)
  {
    strcat(nameList, "    RPGS     |");
    strcpy(local_options->solverName, "RPGS");
    int iparam[2] = {maxIter, 0};
    double dparam[3] = {tolerance, 0.0, 1.0};
    local_options->iSize = 2;
    local_options->dSize = 3;
    local_options->iparam = iparam;
    local_options->dparam = dparam;
    local_options->isSet = 1;
    int info1 = lcp_driver(problem, z[k] , w[k], options, numberOfSolvers, &global_options);
    comp = DDOT(n , z[k] , incx , w[k] , incy);
    DCOPY(n , w[k], incx, wBuffer , incy);
    DAXPY(n , alpha , problem->q , incx , wBuffer  , incy);
    prod(n, n, beta, problem->M, z[k], alpha, wBuffer);
    diff = DNRM2(n , wBuffer , incx);

    if (isSparse == 0)
      printf("    RPGS    (LOG:%1d)|      %5d | %10.4g | %10.4g | %10.4g |\n", info1, local_options->iparam[1], local_options->dparam[1], comp, diff);
    else
      printf("  Gauss-Seidel/RPGS    (LOG:%1d)|      %5d | %10.4g | %10.4g | %10.4g |\n", info1, options[0].iparam[1], options[0].dparam[1], comp, diff);

    if (info1 != 0)
      info = info1;
    k++;
  }
  /* CPG */
  if (solversList[2] == 1)
  {
    strcat(nameList, "     CPG     |");
    strcpy(local_options->solverName, "CPG");
    int iparam[2] = {maxIter, 0};
    double dparam[2] = {tolerance, 0.0};
    local_options->iSize = 2;
    local_options->dSize = 2;
    local_options->iparam = iparam;
    local_options->dparam = dparam;
    local_options->isSet = 1;
    int info1 = lcp_driver(problem, z[k] , w[k], options, numberOfSolvers, &global_options);
    comp = DDOT(n , z[k] , incx , w[k] , incy);
    DCOPY(n , w[k], incx, wBuffer , incy);
    DAXPY(n , alpha , problem->q , incx , wBuffer  , incy);
    prod(n, n, beta, problem->M, z[k], alpha, wBuffer);
    diff = DNRM2(n , wBuffer , incx);

    printf("    CPG     (LOG:%1d)|      %5d | %10.4g | %10.4g | %10.4g |\n", info1, local_options->iparam[1], local_options->dparam[1], comp, diff);

    if (info1 != 0)
      info = info1;
    k++;
  }
  /* Lemke */
  if (solversList[3] == 1)
  {
    strcat(nameList, "    Lemke    |");
    //double tol = F77NAME(dlamch)("e");
    double tol = tolerance;
    strcpy(local_options->solverName, "Lemke");
    int iparam[2] = {maxIter, 0};
    double dparam[2] = {tol, 0.0};
    local_options->iSize = 2;
    local_options->dSize = 2;
    local_options->iparam = iparam;
    local_options->dparam = dparam ;
    local_options->isSet = 1;
    local_options->filterOn = 1;
    int info1 = lcp_driver(problem, z[k] , w[k], options, numberOfSolvers, &global_options);
    comp = DDOT(n , z[k] , incx , w[k] , incy);
    DCOPY(n , w[k], incx, wBuffer , incy);
    DAXPY(n , alpha , problem->q , incx , wBuffer , incy);
    prod(n, n, beta, problem->M, z[k], alpha, wBuffer);
    diff = DNRM2(n , wBuffer , incx);

    if (isSparse == 0)
      printf("    Lemke   (LOG:%1d)|      %5d |  not set   | %10.4g | %10.4g |\n", info1, local_options->iparam[1], comp, diff);
    else
      printf("  Gauss-Seidel/Lemke   (LOG:%1d)|      %5d |  not set   | %10.4g | %10.4g |\n", info1, options[0].iparam[1], comp, diff);


    if (info1 != 0)
      info = info1;
    k++;
  }

  /* Latin */
  if (solversList[4] == 1)
  {
    strcat(nameList, "   LATIN    |");
    strcpy(local_options->solverName, "Latin");
    int iparam[2] = {maxIter, 0};
    double dparam[3] = {tolerance, 0.0, 0.3};
    local_options->iSize = 2;
    local_options->dSize = 3;
    local_options->iparam = iparam;
    local_options->dparam = dparam;
    local_options->isSet = 1;
    int info1 = lcp_driver(problem, z[k] , w[k], options, numberOfSolvers, &global_options);
    comp = DDOT(n , z[k] , incx , w[k] , incy);
    DCOPY(n , w[k], incx, wBuffer , incy);
    DAXPY(n , alpha , problem->q , incx , wBuffer , incy);
    prod(n, n, beta, problem->M, z[k], alpha, wBuffer);
    diff = DNRM2(n , wBuffer , incx);

    if (isSparse == 0)
      printf("    LATIN   (LOG:%1d)|      %5d | %10.4g | %10.4g | %10.4g |\n", info1, local_options->iparam[1], local_options->dparam[1], comp, diff);
    else
      printf("    Gauss-Seidel/LATIN   (LOG:%1d)|      %5d | %10.4g | %10.4g | %10.4g |\n", info1, options[0].iparam[1], options[0].dparam[1], comp, diff);

    if (info1 != 0)
      info = info1;
    k++;
  }
  /* Latin_w */
  if (solversList[5] == 1)
  {
    strcat(nameList, "  LATIN_w    |");
    strcpy(local_options->solverName, "Latin_w");
    int iparam[2] = {maxIter, 0};
    double dparam[4] = {tolerance, 0.0, 0.3, 1.0};
    local_options->iSize = 2;
    local_options->dSize = 4;
    local_options->iparam = iparam;
    local_options->dparam = dparam;
    local_options->isSet = 1;
    int info1 = lcp_driver(problem, z[k] , w[k], options, numberOfSolvers, &global_options);
    comp = DDOT(n , z[k] , incx , w[k] , incy);
    DCOPY(n , w[k], incx, wBuffer , incy);
    DAXPY(n , alpha , problem->q , incx , wBuffer , incy);
    prod(n, n, beta, problem->M, z[k], alpha, wBuffer);
    diff = DNRM2(n , wBuffer , incx);

    if (isSparse == 0)
      printf("    LATIN_W (LOG:%1d)|      %5d | %10.4g | %10.4g | %10.4g |\n", info1, local_options->iparam[1], local_options->dparam[1], comp, diff);
    else
      printf("    Gauss-Seidel/LATIN_W (LOG:%1d)|      %5d | %10.4g | %10.4g | %10.4g |\n", info1, options[0].iparam[1], options[0].dparam[1], comp, diff);

    if (info1 != 0)
      info = info1;
    k++;
  }

  /* Path */
#ifdef HAVE_PATHFERRIS
  if (solversList[6] == 1)
  {
    strcat(nameList, "    PATH     |");
    strcpy(local_options->solverName, "Path");
    double dparam[2] = {tolerance, 0.0};
    local_options->iSize = 0;
    local_options->dSize = 2;
    local_options->iparam = NULL;
    local_options->dparam = dparam;
    local_options->isSet = 1;
    int info1 = lcp_driver(problem, z[k] , w[k], options, numberOfSolvers, &global_options);
    comp = DDOT(n , z[k] , incx , w[k] , incy);
    DCOPY(n , w[k], incx, wBuffer , incy);
    DAXPY(n , alpha , problem->q , incx , wBuffer , incy);
    prod(n, n, beta, problem->M, z[k], alpha, wBuffer);
    diff = DNRM2(n , wBuffer , incx);

    if (isSparse == 0)
      printf("    Path    (LOG:%1d)|  not set   | %10.4g | %10.4g | %10.4g |\n", info1, local_options->dparam[1], comp, diff);
    else
      printf("    Gauss-Seidel/Path    (LOG:%1d)|  not set   | %10.4g | %10.4g | %10.4g |\n", info1, options[0].dparam[1], comp, diff);

    if (info1 != 0)
      info = info1;
    k++;
  }
#endif
  /* QP */
  if (solversList[7] == 1)
  {
    strcat(nameList, "     QP      |");
    strcpy(local_options->solverName, "QP");
    double dparam[2] = {tolerance, 0.0};
    local_options->iSize = 0;
    local_options->dSize = 2;
    local_options->iparam = NULL;
    local_options->dparam = dparam;
    local_options->isSet = 1;
    int info1 = lcp_driver(problem, z[k] , w[k], options, numberOfSolvers, &global_options);
    comp = DDOT(n , z[k] , incx , w[k] , incy);
    DCOPY(n , w[k], incx, wBuffer , incy);
    DAXPY(n , alpha , problem->q , incx , wBuffer , incy);
    prod(n, n, beta, problem->M, z[k], alpha, wBuffer);
    diff = DNRM2(n , wBuffer , incx);

    printf("    QP      (LOG:%1d)|  not set   |  not set   | %10.4g | %10.4g |\n", info1, comp, diff);

    if (info1 != 0)
      info = info1;
    k++;
  }
  /* NSQP */
  if (solversList[8] == 1)
  {
    strcat(nameList, "    NSQP     |");
    strcpy(local_options->solverName, "NSQP");
    double dparam[2] = {tolerance, 0.0};
    local_options->iSize = 0;
    local_options->dSize = 2;
    local_options->iparam = NULL;
    local_options->dparam = dparam;
    local_options->isSet = 1;
    int info1 = lcp_driver(problem, z[k] , w[k], options, numberOfSolvers, &global_options);
    comp = DDOT(n , z[k] , incx , w[k] , incy);
    DCOPY(n , w[k], incx, wBuffer , incy);
    DAXPY(n , alpha , problem->q , incx , wBuffer , incy);
    prod(n, n, beta, problem->M, z[k], alpha, wBuffer);
    diff = DNRM2(n , wBuffer , incx);

    printf("    NSQP    (LOG:%1d)|  not set   |  not set   | %10.4g | %10.4g |\n", info1, comp, diff);

    if (info1 != 0)
      info = info1;
    k++;
  }
  /* Newton Min */
  if (solversList[9] == 1)
  {
    strcat(nameList, "   NewtonMin |");
    strcpy(local_options->solverName, "NewtonMin");
    int iparam[2] = {maxIter, 0};
    double dparam[2] = {tolerance, 0.0};
    local_options->iSize = 2;
    local_options->dSize = 2;
    local_options->iparam = iparam;
    local_options->dparam = dparam;
    local_options->isSet = 1;
    int info1 = lcp_driver(problem, z[k] , w[k], options, numberOfSolvers, &global_options);
    comp = DDOT(n , z[k] , incx , w[k] , incy);
    DCOPY(n , w[k], incx, wBuffer , incy);
    DAXPY(n , alpha , problem->q , incx , wBuffer , incy);
    prod(n, n, beta, problem->M, z[k], alpha, wBuffer);
    diff = DNRM2(n , wBuffer , incx);

    if (isSparse == 0)
      printf(" Newton Min (LOG:%1d)|      %5d | %10.4g | %10.4g | %10.4g |\n", info1, local_options->iparam[1], local_options->dparam[1], comp, diff);
    else
      printf(" Gauss-Seidel/Newton Min (LOG:%1d)|      %5d | %10.4g | %10.4g | %10.4g |\n", info1, options[0].iparam[1], options[0].dparam[1], comp, diff);

    if (info1 != 0)
      info = info1;
    k++;
  }

  /*Enumeratif solver*/
  if (solversList[10] == 1)
  {
    strcat(nameList, "   ENUM |");
    strcpy(local_options->solverName, "ENUM");
    int iparam[2] = {0, 0};
    double dparam[1] = {tolerance};
    local_options->iSize = 2;
    local_options->dSize = 1;
    local_options->iparam = iparam;
    local_options->dparam = dparam;
    int info1 = 1;
    if (problem->size < 20)
    {
      local_options->dWork = (double*) malloc((3 * problem->size + problem->size * problem->size) * sizeof(double));
      local_options->iWork = (int*) malloc(2 * problem->size * sizeof(int));
      info1 = lcp_driver(problem, z[k] , w[k], options, numberOfSolvers, &global_options);
      free(local_options->dWork);
      free(local_options->iWork);
      comp = DDOT(n , z[k] , incx , w[k] , incy);
      DCOPY(n , w[k], incx, wBuffer , incy);
      DAXPY(n , alpha , problem->q , incx , wBuffer , incy);
      prod(n, n, beta, problem->M, z[k], alpha, wBuffer);
      diff = DNRM2(n , wBuffer , incx);
    }
    else
    {
      comp = 1;
      diff = 1;
    }

    printf("    ENUM    (LOG:%1d)|      %5d | %10.4g | %10.4g | %10.4g |\n", info1, local_options->iparam[1], local_options->dparam[1], comp, diff);

    if (info1 != 0)
      info = info1;
    k++;
  }


  /* Newton Fischer-Burmeister */
  if (solversList[11] == 1)
  {
    strcat(nameList, "   Newton FB |");
    strcpy(local_options->solverName, "NewtonFB");
    int iparam[2] = {10, 0};
    double dparam[2] = {tolerance, 0.0};
    local_options->iSize = 2;
    local_options->dSize = 2;
    local_options->iparam = iparam;
    local_options->dparam = dparam;
    local_options->isSet = 1;
    int info1 = lcp_driver(problem, z[k] , w[k], options, numberOfSolvers, &global_options);
    comp = DDOT(n , z[k] , incx , w[k] , incy);
    DCOPY(n , w[k], incx, wBuffer , incy);
    DAXPY(n , alpha , problem->q , incx , wBuffer , incy);
    prod(n, n, beta, problem->M, z[k], alpha, wBuffer);
    diff = DNRM2(n , wBuffer , incx);

    if (isSparse == 0)
      printf("\n    Newton FB   (LOG:%1d)|      %5d | %10.4g | %10.4g | %10.4g |", info1, local_options->iparam[1], local_options->dparam[1], comp, diff);
    else
      printf("\n    Gauss-Seidel/Newton FB   (LOG:%1d)|      %5d | %10.4g | %10.4g | %10.4g |", info1, options[0].iparam[1], options[0].dparam[1], comp, diff);

    if (info1 != 0)
      info = info1;
    k++;
  }


  printf("\n\n");

  /* =========================== Ouput: comparison between the different methods =========================== */

  strcat(nameList, "\n");
  printf(" *****   z = \n\n");
  printf(nameList);
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < nbSolvers; j++)
      printf("%11.5g | ", z[j][i]);
    printf("\n");
  }
  printf("\n\n");
  printf(" *****   w = \n\n");
  printf(nameList);
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < nbSolvers; j++)
      printf("%11.5g | ", w[j][i]);
    printf("\n");
  }
  printf("\n\n");

  for (i = 0; i < nbSolvers; i++)
  {
    free(z[i]);
    free(w[i]);
  }
  free(z);
  free(w);

  return info;

  free(wBuffer);
  free(options);

}

int test_mmc(void)
{
  printf("========================================================================================================== \n");
  printf("                                     LCP Solvers tests (function: test_mmc)  \n");
  printf("==========================================================================================================\n");
  FILE *f1, *f2;
  int i, nl, nc;
  double qi, Mij;
  char val[20], vall[20];

  /* Building of the LCP */
  LinearComplementarity_Problem * problem = malloc(sizeof(*problem));

  // computation of the size of the problem

  int n = 0;
  // Open files ...
  if ((f1 = fopen("DATA/M_rectangle1.dat", "r")) == NULL)
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

  problem->size = n;

  // M
  double * vecM = (double*)malloc(n * n * sizeof(double));

  /* Data loading for M and q */

  if ((f1 = fopen("DATA/M_rectangle1.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }

  if ((f2 = fopen("DATA/q_rectangle1.dat", "r")) == NULL)
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
  NumericsMatrix * MM = malloc(sizeof(*MM));
  MM->matrix0 = vecM;
  MM->size0 = n;
  MM->size1 = n;
  MM->storageType = 0;

  double * q = (double*)malloc(n * sizeof(double));
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

  // Set M and q of the problem
  problem->q = q;
  problem->M = MM;

  // Call tests

  printf(" ----------------------------------------------------------\n");
  printf("Run working tests ...\n");
  /* Stable: */
  //  int solversList[11] ={1,1,1,1,0,0,0,1,1,0,0};
  int solversList[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0};
  int info = test_lcp_series(problem, solversList);
  printf(" ----------------------------------------------------------\n");

  /* Fail or unstable: */
  printf("---------------------------------------------------------- \n");
  printf("\n Run unstable tests (results may be wrong or log !=0)...\n");
  int solversList2[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1};
  int infoFail = test_lcp_series(problem, solversList2);
  printf("--------- End of unstable tests --------------------------- \n");

  // Release memory
  problem->M = NULL;
  problem->q = NULL;
  free(MM);
  free(problem);
  free(vecM);
  free(q);

  printf("========================================================================================================== \n");
  printf("                                           END OF TEST MMC     \n");
  printf("==========================================================================================================\n");
  return info;

}

/* To read in a file a LinearComplementarity_Problem with a "double*" storage for M */
void getProblem(char* name, LinearComplementarity_Problem *  problem)
{

  FILE * LCPfile =  fopen(name, "r");
  if (LCPfile == NULL)
  {
    fprintf(stderr, "fopen LCPfile: %s\n", name);
    exit(1);
  }
  printf("\n\n******************************************************\n");
  printf("Read Linear Complementarity Problem in file %s\n", name);
  printf("******************************************************\n");

  /* Dim of the LCP */
  int dim;
  fscanf(LCPfile , "%d" , &dim);
  int dim2 = dim * dim;

  problem->M->matrix0 = malloc(dim2 * sizeof(double));
  problem->q = (double*)malloc(dim * sizeof(double));

  double * vecM = problem->M->matrix0;
  int i, j;
  char val[20];
  /* fill M */
  for (i = 0 ; i < dim ; ++i)
  {
    for (j = 0 ; j < dim ; ++j)
    {
      fscanf(LCPfile, "%s", val);
      vecM[ dim * j + i ] = atof(val);
    }
  }

  /* fill q */
  for (i = 0 ; i < dim ; ++i)
  {
    fscanf(LCPfile , "%s" , val);
    problem->q[i] = atof(val);
  }

  /* fill sol */
  double* sol = NULL;
  fscanf(LCPfile , "%s" , val);
  if (!feof(LCPfile))
  {
    sol  = (double*)malloc(dim * sizeof(double));
    sol[0] = atof(val);
    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }
  }
  printf("\n exact solution : ");
  if (sol != NULL) for (i = 0 ; i < problem->size ; ++i) printf(" %10.4g " , sol[i]);
  else printf(" unknown ");
  printf("\n");

  problem->size = dim;
  problem->M->size0 = dim;
  problem->M->size1 = dim;
  fclose(LCPfile);
  if (sol != NULL)
    free(sol);
}

/* To read in a file a LinearComplementarity_Problem with a "SparseBlockStructuredMatrix*" storage for M */
void getProblemSBM(char* name, LinearComplementarity_Problem *  problem)
{

  FILE * LCPfile =  fopen(name, "r");
  if (LCPfile == NULL)
  {
    fprintf(stderr, "fopen LCPfile: %s\n", name);
    exit(1);
  }
  printf("\n\n******************************************************\n");
  printf("Read Linear Complementarity Problem in file %s\n", name);
  printf("******************************************************\n");

  printf("\n The matrix M of the LCP is a SparseBlockStructuredMatrix.\n");

  SparseBlockStructuredMatrix * blmat =  problem->M->matrix1;

  int i, j;
  char val[20];

  /***** M *****/
  fscanf(LCPfile , "%d" , &blmat->nbblocks);
  fscanf(LCPfile , "%d" , &blmat->size);
  blmat->blocksize = (int*)malloc(blmat->size * sizeof(int));
  for (i = 0 ; i < blmat->size ; i++) fscanf(LCPfile , "%d" , &blmat->blocksize[i]);
  blmat->RowIndex = (int*)malloc(blmat->nbblocks * sizeof(int));
  blmat->ColumnIndex = (int*)malloc(blmat->nbblocks * sizeof(int));
  blmat->index1_data = (size_t*) malloc((blmat->size + 1) * sizeof(size_t));
  for (i = 0 ; i < blmat->nbblocks ; i++)
  {
    fscanf(LCPfile , "%d" , &blmat->RowIndex[i]);
    fscanf(LCPfile , "%d" , &blmat->ColumnIndex[i]);
  }


  /* Boost sparse compressed : index1_data construction */
  /* see SparseBlockMatrix.h */
  int current_block;
  int current_row;
  int k;

  current_block = 0;
  current_row = 0;
  while (current_block < blmat->nbblocks)
  {
    for (k = current_block; blmat->RowIndex[current_block] == blmat->RowIndex[k] && k < blmat->nbblocks; ++k) ;
    for (i = current_row; i < blmat->RowIndex[current_block] + 1; ++i)
      blmat->index1_data[i] = current_block;

    blmat->index1_data[i] = k;
    current_row = blmat->RowIndex[current_block] + 1;
    current_block = k;
  };

  int k1, k2;

  for (i = 0, k1 = 0; k1 < blmat->nbblocks; k1 = k2, i++)
  {
    for (k2 = k1; blmat->RowIndex[k1] == blmat->RowIndex[k2] && k2 < blmat->nbblocks; k2++);
    blmat->index1_data[i] = k1;
    blmat->index1_data[i + 1] = k2;
  };

  blmat->block = (double**)malloc(blmat->nbblocks * sizeof(double*));
  int pos, sizebl, numberOfRows, numberOfColumns;
  for (i = 0 ; i < blmat->nbblocks ; i++)
  {
    pos = blmat->RowIndex[i];
    numberOfRows = blmat->blocksize[pos];
    if (pos > 0)
      numberOfRows -= blmat->blocksize[pos - 1];
    pos = blmat->ColumnIndex[i];
    numberOfColumns = blmat->blocksize[pos];
    if (pos > 0)
      numberOfColumns -= blmat->blocksize[pos - 1];
    sizebl = numberOfRows * numberOfColumns;
    blmat->block[i] = (double*)malloc(sizebl * sizeof(double));
    for (j = 0 ; j < sizebl ; j++)
    {
      fscanf(LCPfile, "%s", val);
      blmat->block[i][j] = atof(val);
    }
  }

  int dim = blmat->blocksize[blmat->size - 1];
  /**** q ****/
  problem->q = (double*)malloc(dim * sizeof(double));
  for (i = 0 ; i < dim ; i++)
  {
    fscanf(LCPfile , "%s" , val);
    problem->q[i] = atof(val);
  }

  fscanf(LCPfile , "%s" , val);

  double* sol = NULL;
  if (!feof(LCPfile))
  {
    sol  = (double*)malloc(dim * sizeof(double));
    sol[0] = atof(val);
    for (i = 1 ; i < dim ; i++)
    {
      fscanf(LCPfile , "%s" , val);
      sol[i] = atof(val);
    }
  }
  printf("\n exact solution : ");
  if (sol != NULL) for (i = 0 ; i < problem->size ; ++i) printf(" %10.4g " , sol[i]);
  else printf(" unknown ");
  printf("\n");

  problem->size = dim;
  problem->M->size0 = dim;
  problem->M->size1 = dim;
  fclose(LCPfile);
  if (sol != NULL)
    free(sol);
}

int test_matrix(void)
{
  /*
    Two problems: one with dense (double*) storage for M, "problem", and the other
    with SparseBlockStructuredMatrix storage, "problemSBM".

    For each case, one or both problems are read in a dat file. \n
    The according to the value of the lists solversList and solversList2 (for problem) ,\n
    solversListSBM (for problemSBM).
    the problems are solved with different solvers. \n
    list[i] = 1 => solver is called, list[i] = 0, solver is ignored.
    Check on top of test_lcp_series() function for the list of available solvers and their corresponding indices.

    "stable" tests, that must succeed, are those defined in solversList, while "unstable", that may fail or
    return an unexpected termination value are those defined in solversList2.

   */

  printf("========================================================================================================== \n");
  printf("                         LCP Solvers tests (function: test_matrix)  \n");
  printf("==========================================================================================================\n");

  FILE *LCPfile = NULL, *LCPfileBlock = NULL;
  int i, j, itest;
  int iter = 0;
  double criteria = 0.0;
  double *sol = NULL;

  int NBTEST = 15;

  /* === Building of the LCPs === */

  /* LCP with dense storage */
  LinearComplementarity_Problem * problem = malloc(sizeof(*problem));
  problem->M = malloc(sizeof(*(problem->M)));
  problem->M->storageType = 0;
  problem->M->matrix1 = NULL;

  /* LCP with sparse-block storage */
  LinearComplementarity_Problem * problemSBM = malloc(sizeof(*problemSBM));
  problemSBM->M = malloc(sizeof(*(problemSBM->M)));
  problemSBM->M->storageType = 1;
  problemSBM->M->matrix0 = NULL;
  problemSBM->M->matrix1 = malloc(sizeof(*(problemSBM->M->matrix1)));

  /* List of working solvers */
  int * solversList = NULL; /* for dense */
  int * solversListSBM = NULL; /* for sparse */

  /* List of unstable solvers (failed or wrong termination value) */
  int * solversList2 = NULL;
  int hasSBM = 0; /* =1 if the read matrix exists also in sparse */
  int hasDense = 0;
  int hasUnstable = 0;
  int info = -1;
  for (itest = 0 ; itest < NBTEST ; ++itest)
  {
    hasSBM = 0;
    hasUnstable = 0;
    hasDense = 0;
    /* Read the LCP */
    switch (itest)
    {
    case 0:
      getProblem("MATRIX/deudeu.dat", problem);
      getProblemSBM("MATRIX/deudeu_block.dat", problemSBM);
      hasSBM = 1;
      hasDense = 1;
      hasUnstable = 1;
      {
        int l1[12] = {1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0};
        int l2[12] = {0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1};
        //int l3[11] = {1,1,0,1,0,0,0,1,1,0,0};
        int * l3 = l1;
        solversList = l1;
        solversList2 = l2;
        solversListSBM = l3;
      }
      break;
    case 1:
      getProblem("MATRIX/trivial.dat", problem);
      getProblemSBM("MATRIX/trivial_block.dat", problemSBM);
      hasSBM = 1;
      hasDense = 1;
      hasUnstable = 1;
      {
        int l1[12] = {1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0};
        int l2[12] = {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
        int l3[12] = {1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0};
        solversList = l1;
        solversList2 = l2;
        solversListSBM = l3;
      }
      break;
    case 2:
      getProblem("MATRIX/ortiz.dat", problem);
      getProblemSBM("MATRIX/ortiz_block.dat", problemSBM);
      hasSBM = 1;
      hasDense = 1;
      hasUnstable = 1;
      {
        int l1[12] = {1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0};
        int l2[12] = {0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1};
        int l3[12] = {1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0};
        solversList = l1;
        solversList2 = l2;
        solversListSBM = l3;
      }
      break;
    case 3: /* Ortiz again ... */
      getProblemSBM("MATRIX/ortiz_monoblock.dat", problemSBM);
      hasSBM = 1;
      hasDense = 0;
      hasUnstable = 0;
      {
        int l3[12] = {1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0};
        solversListSBM = l3;
      }
      break;
    case 4: /* Ortiz again ... */
      getProblemSBM("MATRIX/ortiz_block31.dat", problemSBM);
      hasSBM = 1;
      hasDense = 0;
      hasUnstable = 0;
      {
        int l3[12] = {1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0};
        solversListSBM = l3;
      }
      break;
    case 5:
      getProblem("MATRIX/pang.dat", problem);
      hasSBM = 0;
      hasDense = 1;
      hasUnstable = 1;
      {
        int l1[12] = {0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0};
        int l2[12] = {1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1};
        solversList = l1;
        solversList2 = l2;
      }
      break;
    case 6:
      getProblem("MATRIX/diodes.dat", problem);
      hasSBM = 0;
      hasDense = 1;
      hasUnstable = 1;
      {
        int l1[12] = {0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0};
        int l2[12] = {1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1};
        solversList = l1;
        solversList2 = l2;
      }
      break;
    case 7:
      getProblem("MATRIX/murty1.dat", problem);
      hasSBM = 0;
      hasDense = 1;
      hasUnstable = 1;
      {
        int l1[12] = {0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0};
        int l2[12] = {1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1};
        solversList = l1;
        solversList2 = l2;
      }
      break;
    case 8:
      printf("\n\n APPROXIMATED FRICTION CONE LCP ");
      getProblem("MATRIX/confeti.dat", problem);
      getProblemSBM("MATRIX/confeti_block.dat", problemSBM);
      hasSBM = 1;
      hasDense = 1;
      hasUnstable = 1;
      {
        int l1[12] = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0};
        int l2[12] = {1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1};
        int l3[12] = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
        solversList = l1;
        solversList2 = l2;
        solversListSBM = l3;
      }
      break;
    case 9:
      getProblem("MATRIX/mathieu1.dat", problem);
      getProblemSBM("MATRIX/bloc22_mathieu1.dat", problemSBM);
      hasSBM = 1;
      hasDense = 1;
      hasUnstable = 1;
      {
        int l1[12] = {1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0};
        int l2[12] = {0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1};
        int l3[12] = {1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0};
        solversList = l1;
        solversList2 = l2;
        solversListSBM = l3;
      }
      break;
    case 10: /* Same as above with sparse storage (one block) */
      getProblemSBM("MATRIX/monobloc_mathieu1.dat", problemSBM);
      hasSBM = 1;
      hasDense = 0;
      hasUnstable = 0;
      {
        int l1[12] = {1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0};
        solversListSBM = l1;
      }
      break;

    case 11:
      getProblem("MATRIX/mathieu2.dat", problem);
      getProblemSBM("MATRIX/bloc3333_mathieu2.dat", problemSBM);
      hasSBM = 1;
      hasDense = 1;
      hasUnstable = 1;
      {
        int l1[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        int l2[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        int l3[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        solversList = l1;
        solversList2 = l2;
        solversListSBM = l3;
      }
      break;

    case 12:
      getProblem("MATRIX/trivial3.dat", problem);
      getProblemSBM("MATRIX/trivial3_block.dat", problemSBM);
      hasSBM = 1;
      hasDense = 1;
      hasUnstable = 1;
      {
        int l1[12] = {1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1};
        int l2[12] = {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0};
        int l3[12] = {1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1};
        solversList = l1;
        solversList2 = l2;
        solversListSBM = l3;
      }
      break;
    case 13:
      getProblem("MATRIX/buckconverterregul2.dat", problem);
      getProblemSBM("MATRIX/buckconverterregul2_block.dat", problemSBM);
      hasSBM = 1;
      hasDense = 1;
      hasUnstable = 1;
      {
        int l1[12] = {0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0};
        int l2[12] = {1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1};
        int l3[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        solversList = l1;
        solversList2 = l2;
        solversListSBM = l3;
      }
    case 14:
      getProblem("MATRIX/relay.dat", problem);
      hasSBM = 0;
      hasDense = 1;
      hasUnstable = 1;
      {
        int l1[12] = {0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0};
        int l2[12] = {1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1};
        solversList = l1;
        solversList2 = l2;
      }
      break;
    }

    printf(" ----------------------------------------------------------\n");
    printf("Run working tests ...\n");
    /* Stable: */
    int infoTmp = -1;
    if (hasDense == 1)
      infoTmp = test_lcp_series(problem, solversList);

    if (hasSBM == 1)
    {
      printf("Run working tests for sparse storage ...\n");
      infoTmp = test_lcp_series(problemSBM, solversListSBM);
    }

    printf(" ----------------------------------------------------------\n\n");

    if (hasUnstable == 1)
    {
      /* Fail or unstable: */
      printf("---------------------------------------------------------- \n");
      printf("\n Run unstable tests (results may be wrong or log !=0)...\n");
      int infoFail = test_lcp_series(problem, solversList2);
      printf("--------- End of unstable tests --------------------------- \n\n");
    }

    /* Free Memory */
    if (problem->M->matrix0 != NULL)
      free(problem->M->matrix0);
    problem->M->matrix0 = NULL;
    if (problem->q != NULL)
      free(problem->q);
    problem->q = NULL;
    if (problemSBM->q != NULL)
      free(problemSBM->q);
    problemSBM->q = NULL;
    info = infoTmp;
    if (hasSBM == 1)
      freeSBM(problemSBM->M->matrix1);

    if (infoTmp != 0)
      break;
  }
  free(problemSBM->M->matrix1);
  free(problemSBM);
  free(problem->M);
  free(problem);
  printf("========================================================================================================== \n");
  printf("                                  END OF TEST MATRIX   \n");
  printf("========================================================================================================== \n");

  return info;

}


int main(void)
{

  /* In each test function, two series of tests are called:
     - working tests => must return info = 0 , else error, test fail.
     - unstable tests => just print results on screen

  */

  int info1 = test_mmc();

  if (info1 != 0)
  {
    // printf("Warning: test_mmc log output different from 0 for some solvers.\n");
    return info1;
  }

  int info2 = test_matrix();
  if (info2 != 0)
  {
    // printf("Warning: test_mmc log output different from 0 for some solvers.\n");
    return info2;
  }

  return 0;
}

