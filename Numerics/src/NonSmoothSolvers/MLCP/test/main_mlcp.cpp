/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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
#include "sys/time.h"
#include "LA.h"
/*only for debug functions*/
#include "NonSmoothNewtonNeighbour.h"
#include "MixedLinearComplementarityProblem.h"

#define NAX_NBTESTS 20

//#define BAVARD
//#define NBTEST 19
#define NBTEST 10

#define ENUM_ID 0
#define PGS_EX_ID 1
#define PGS_IM_ID 2
#define RPGS_ID 3
#define PSOR_05_ID 4
#define PSOR_1_ID 5
#define PSOR_15_ID 6
#define PSOR_2_ID 7
#define RPSOR_ID 8
#define PATH_ID 9
#define SIMPLEX_ID 10
#define DIRECT_ENUM_ID 11
#define FB_ID 12
#define DIRECT_FB_ID 13
#define NBMETHODS 14

#define PATH_DRIVER
#define MAX_DIM_ENUM 19

/*#ifdef PATH_DRIVER
const unsigned short int *__ctype_b;
const __int32_t *__ctype_tolower ;
#endif
*/

static int sNbNOCV = 0;
static int sIdWithSol = 0;
typedef struct
{
  char file[124];
  char cv[NBMETHODS][2][5];
  int nbSteps[NBMETHODS][2];
  long times[NBMETHODS][2];
} dataSummary;

typedef struct
{
  char Name[124];
  struct timeval mStart;
  long mCumul;
} dataTime;
static dataTime sDt;
void startTimer()
{
  gettimeofday(&(sDt.mStart), NULL);
  sDt.mCumul = 0;
}
void stopTimer()
{
  struct timeval aux;
  gettimeofday(&aux, NULL);
  sDt.mCumul += (aux.tv_sec - sDt.mStart.tv_sec) * 1000000 + (aux.tv_usec - sDt.mStart.tv_usec) ;
}
static dataSummary summary[NBTEST];
static int sRunMethod[NBMETHODS];
static int itest;


/*
 ******************************************************************************
 */
void printSolution(const char *name, int n, int m, int NbLines, double *z, double *w)
{
  int i;
#ifdef BAVARD
  printf(" *** ************************************** ***\n");
  printf(" ****** z = ********************************\n");
  for (i = 0 ; i < n + m ; i++)
    printf("%s: %14.7e  \n", name, z[i]);
  printf(" ****** w = ********************************\n");
  for (i = 0 ; i < m ; i++)
    printf("%s: %14.7e  \n", name, w[i + (NbLines - m)]);
#endif
}
/*
 *sol = (z,w)
 *
 */
void solTozw(int n, int m, double *z, double *w, double *sol)
{
  int i;
  for (i = 0; i < n + m; i++)
    z[i] = sol[i];
  for (i = 0; i < m; i++)
    w[i] = sol[n + m + i];

}

void test_mlcp_series(MixedLinearComplementarityProblem* problem, double *z, double *w, double *sol)
{
  int info;
  SolverOptions mlcpOptions;
  NumericsOptions global_options;
  double tol1 = 1e-15;
  double tol2 = 1e-6;
  double error = 0;
  int n = problem->n;
  int m = problem->m;
  double ohmega = 0;
  int _ID = PSOR_05_ID - 1;
  int ndW = 0;
  int niW = 0;
  int aux = 0;
  long laux;
  int NbLines = problem->M->size0;

  //   mlcpOptions.isSet=0;
  //   mlcpOptions.iSize=0;
  //   mlcpOptions.iparam=0;
  //   mlcpOptions.dSize=0;
  //   mlcpOptions.dparam=0;
  //   mlcpOptions.filterOn=0;
  //   mlcpOptions.dWork=0;
  //   mlcpOptions.iWork=0;
  //   mlcpOptions.iparam=(int*)malloc(10*sizeof(int));
  //   mlcpOptions.dparam=(double*)malloc(10*sizeof(double));

  //  global_options.verboseMode=1;
  setNumericsOptions(&global_options);

  //   mlcpOptions.iparam[5]=3;/*Number of registered configurations*/
  //   mlcpOptions.iparam[8]=0;/*Prb nedd a update*/
  //   mlcpOptions.dparam[5]=1e-12;
  //   mlcpOptions.dparam[6]=1e-12;

  //   /*iwork and dwork memory*/
  //   strcpy(mlcpOptions.solverName,"DIRECT_ENUM");
  //   niW=mlcp_driver_get_iwork(problem,&mlcpOptions);
  //   ndW=mlcp_driver_get_dwork(problem,&mlcpOptions);

  //   strcpy(mlcpOptions.solverName,"DIRECT_FB");
  //   aux= mlcp_driver_get_iwork(problem,&mlcpOptions);
  //   if (aux > niW)
  //     niW = aux;
  //   aux=mlcp_driver_get_dwork(problem,&mlcpOptions);
  //   if (aux > ndW)
  //     ndW = aux;


  //   strcpy(mlcpOptions.solverName,"DIRECT_PATH");
  //   aux= mlcp_driver_get_iwork(problem,&mlcpOptions);
  //   if (aux > niW)
  //     niW = aux;
  //   aux=mlcp_driver_get_dwork(problem,&mlcpOptions);
  //   if (aux > ndW)
  //     ndW = aux;



  //   mlcpOptions.dWork=(double*)malloc( ndW*sizeof(double));
  //   mlcpOptions.iWork=(int*)malloc(niW*sizeof(int));
  //   if (mlcpOptions.dWork == 0 || mlcpOptions.iWork ==0){
  //     printf("test_mlcp_series allocate memory failed %d %d %d %d.\n",niW, ndW,problem->n,problem->m);
  //     exit(1);
  //   }



  /*SOLVER ENUM*/
  if (sRunMethod[ENUM_ID])
  {
    solTozw(n, m, z, w, sol);
    printf("TRY SOLVER %s\n", "ENUM");
    strcpy(mlcpOptions.solverName, "ENUM");
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    //    mlcpOptions.iSize=1;
    //    mlcpOptions.dSize=1;
    mlcpOptions.dparam[0] = tol1;

    mlcp_driver_init(problem, &mlcpOptions);
    startTimer();
    info = 1;
    if (n + m < MAX_DIM_ENUM)
      info = mlcp_driver(problem, z, w, &mlcpOptions, &global_options);
    stopTimer();
    summary[itest].times[ENUM_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[ENUM_ID][sIdWithSol], "CV");
    if (info > 0)
    {
      printf("Can't find a solution\n");
      strcpy(summary[itest].cv[ENUM_ID][sIdWithSol], "NO");
      sNbNOCV++;
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      printf("find a solution with error %lf \n", error);
      printSolution("ENUM", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }

  /*SOLVER PGS*/
  if (sRunMethod[PGS_IM_ID])
  {
    printf("TRY SOLVER %s\n", "PGS");
    solTozw(n, m, z, w, sol);
    strcpy(mlcpOptions.solverName, "PGS");
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    //"PGS"       , 101 , 1e-8 , 0.6 , 1.0 , 1 , 0 , 0.0 }; */
    mlcpOptions.dparam[0] = tol2;
    mlcp_driver_init(problem, &mlcpOptions);

    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions, &global_options);
    stopTimer();
    summary[itest].times[PGS_IM_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[PGS_IM_ID][sIdWithSol], "CV");
    if (info > 0)
    {
      printf("Can't find a solution\n");
      strcpy(summary[itest].cv[PGS_IM_ID][sIdWithSol], "NO");
      sNbNOCV++;
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      printf("find a solution with error %lf \n", error);
      printSolution("PGS", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }

  /*SOLVER PGS*/
  if (sRunMethod[PGS_EX_ID])
  {

    printf("TRY SOLVER %s\n", "EX PGS");
    solTozw(n, m, z, w, sol);
    strcpy(mlcpOptions.solverName, "PGS");
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    //"PGS"       , 101 , 1e-8 , 0.6 , 1.0 , 1 , 0 , 0.0 }; */
    mlcpOptions.dparam[0] = tol2;

    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions, &global_options);
    stopTimer();
    summary[itest].times[PGS_EX_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[PGS_EX_ID][sIdWithSol], "CV");
    if (info > 0)
    {
      printf("Can't find a solution\n");
      strcpy(summary[itest].cv[PGS_EX_ID][sIdWithSol], "NO");
      sNbNOCV++;
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      printf("find a solution with error %lf \n", error);
      printSolution("PGS", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER RPGS*/
  if (sRunMethod[RPGS_ID])
  {
    printf("TRY SOLVER %s\n", "RPGS");
    solTozw(n, m, z, w, sol);
    strcpy(mlcpOptions.solverName, "RPGS");
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    mlcpOptions.dparam[0] = tol2;
    mlcpOptions.dparam[2] = 0.5; /*rho*/
    mlcp_driver_init(problem, &mlcpOptions);

    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions, &global_options);
    stopTimer();
    summary[itest].times[RPGS_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[RPGS_ID][sIdWithSol], "CV");
    if (info > 0)
    {
      printf("Can't find a solution\n");
      strcpy(summary[itest].cv[RPGS_ID][sIdWithSol], "NO");
      sNbNOCV++;
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      printf("find a solution with error %lf \n", error);
      printSolution("RPGS", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER PSOR*/
  if (sRunMethod[_ID])
  {
    printf("TRY SOLVER %s\n", "PSOR");
    ohmega = 0;
    _ID = PSOR_05_ID - 1;
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    for (int cmp = 0; cmp < 4; cmp++)
    {
      _ID++;
      ohmega += 0.5;

      solTozw(n, m, z, w, sol);
      strcpy(mlcpOptions.solverName, "PSOR");
      mlcpOptions.dparam[0] = tol2;
      mlcpOptions.dparam[2] = ohmega;
      mlcp_driver_init(problem, &mlcpOptions);

      startTimer();
      info = mlcp_driver(problem, z, w, &mlcpOptions, &global_options);
      stopTimer();
      summary[itest].times[_ID][sIdWithSol] = sDt.mCumul;
      strcpy(summary[itest].cv[_ID][sIdWithSol], "CV");
      if (info > 0)
      {
        printf("Can't find a solution\n");
        strcpy(summary[itest].cv[_ID][sIdWithSol], "NO");
        sNbNOCV++;
      }
      else
      {
        mlcp_compute_error(problem, z, w, tol1,  &error);
        printf("find a solution with error %lf \n", error);
        printSolution("PSOR", n, m, NbLines, z, w);
      }
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER RPSOR*/
  if (sRunMethod[RPSOR_ID])
  {
    printf("TRY SOLVER %s\n", "RPSOR");
    solTozw(n, m, z, w, sol);
    strcpy(mlcpOptions.solverName, "RPSOR");
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    mlcpOptions.dparam[0] = tol2;
    mlcpOptions.dparam[2] = 0.5; /*rho*/
    mlcpOptions.dparam[3] = 2; /*ohmega*/
    mlcp_driver_init(problem, &mlcpOptions);

    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions, &global_options);
    stopTimer();
    summary[itest].times[RPSOR_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[RPSOR_ID][sIdWithSol], "CV");
    if (info > 0)
    {
      printf("Can't find a solution\n");
      strcpy(summary[itest].cv[RPSOR_ID][sIdWithSol], "NO");
      sNbNOCV++;
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      printf("find a solution with error %lf \n", error);
      printSolution("RPSOR", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER PATH*/
  if (sRunMethod[PATH_ID])
  {
    printf("TRY SOLVER %s\n", "PATH");
    solTozw(n, m, z, w, sol);
    strcpy(mlcpOptions.solverName, "PATH");
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_init(problem, &mlcpOptions);

    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions, &global_options);
    stopTimer();
    summary[itest].times[PATH_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[PATH_ID][sIdWithSol], "CV");
    if (info > 0)
    {
      printf("Can't find a solution\n");
      strcpy(summary[itest].cv[PATH_ID][sIdWithSol], "NO");
      sNbNOCV++;
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      if (error > tol2)
      {
        strcpy(summary[itest].cv[PATH_ID][sIdWithSol], "NO");
        sNbNOCV++;
      }
      printf("find a solution with error %lf \n", error);
      printSolution("PATH", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER SIMPLEX*/
  if (sRunMethod[SIMPLEX_ID])
  {
    printf("TRY SOLVER %s\n", "SIMPLEX");


    solTozw(n, m, z, w, sol);
    strcpy(mlcpOptions.solverName, "SIMPLEX");
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    mlcpOptions.iparam[0] = 1000000;
    mlcpOptions.iparam[1] = 1;
    mlcpOptions.dparam[0] = 1e-12;
    mlcpOptions.dparam[1] = 1e-12;
    mlcpOptions.dparam[2] = 1e-9;
    mlcp_driver_init(problem, &mlcpOptions);

    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions, &global_options);
    stopTimer();
    summary[itest].times[SIMPLEX_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[SIMPLEX_ID][sIdWithSol], "CV");
    if (info > 0)
    {
      printf("Can't find a solution\n");
      strcpy(summary[itest].cv[SIMPLEX_ID][sIdWithSol], "NO");
      sNbNOCV++;
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      if (error > 1e-9)
      {
        strcpy(summary[itest].cv[SIMPLEX_ID][sIdWithSol], "NO");
        sNbNOCV++;
      }
      printf("find a solution with error %lf \n", error);
      printSolution("SIMPLEX", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER DIRECT ENUM*/
  if (sRunMethod[DIRECT_ENUM_ID])
  {
    printf("TRY SOLVER %s\n", "DIRECT_ENUM");
    solTozw(n, m, z, w, sol);
    strcpy(mlcpOptions.solverName, "DIRECT_ENUM");
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    mlcpOptions.dparam[0] = 1e-12;
    mlcp_driver_init(problem, &mlcpOptions);

    printf("Solver direct, precompute solution.\n");
    if (n + m < MAX_DIM_ENUM)
      info = mlcp_driver(problem, z, w, &mlcpOptions, &global_options);
    if (info == 0)
    {
      printf("Solver direct, solution computed. ");
      printf("Now, run direct solver.\n");
      startTimer();
      if (n + m < MAX_DIM_ENUM)
        info = mlcp_driver(problem, z, w, &mlcpOptions, &global_options);
      stopTimer();
      summary[itest].times[DIRECT_ENUM_ID][sIdWithSol] = sDt.mCumul;
      strcpy(summary[itest].cv[DIRECT_ENUM_ID][sIdWithSol], "CV");
    }
    if (info > 0)
    {
      printf("Can't find a solution\n");
      strcpy(summary[itest].cv[DIRECT_ENUM_ID][sIdWithSol], "NO");
      sNbNOCV++;
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      if (error > 1e-9)
      {
        strcpy(summary[itest].cv[DIRECT_ENUM_ID][sIdWithSol], "NO");
        sNbNOCV++;
      }
      printf("find a solution with error %lf \n", error);
      printSolution("DIRECT_ENUM_ID", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER FB*/
  if (sRunMethod[FB_ID])
  {

    printf("TRY SOLVER %s\n", "FB");
    solTozw(n, m, z, w, sol);
    strcpy(mlcpOptions.solverName, "FB");
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    mlcpOptions.dparam[0] = 1e-10;
    mlcpOptions.dparam[1] = 0;


    mlcp_driver_init(problem, &mlcpOptions);
    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions, &global_options);
    stopTimer();
    summary[itest].times[FB_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[FB_ID][sIdWithSol], "CV");
    if (info > 0)
    {
      printf("Can't find a solution\n");
      strcpy(summary[itest].cv[FB_ID][sIdWithSol], "NO");
      sNbNOCV++;
      mlcp_compute_error(problem, z, w, tol1,  &error);
      printf("current point with error %lf \n", error);
      printSolution("FB", n, m, NbLines, z, w);
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      //    if (error > 1e-9)
      strcpy(summary[itest].cv[FB_ID][sIdWithSol], "CV");
      printf("find a solution with error %lf \n", error);
      printSolution("FB", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER DIRECT_FB*/
  if (sRunMethod[DIRECT_FB_ID])
  {

    printf("TRY SOLVER %s\n", "DIRECT_FB");
    solTozw(n, m, z, w, sol);
    strcpy(mlcpOptions.solverName, "DIRECT_FB");
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    mlcpOptions.iparam[0] = 500;
    mlcpOptions.iparam[1] = 0;

    mlcpOptions.dparam[0] = 1e-11;
    mlcpOptions.dparam[1] = 0;
    info = 1;

    mlcp_driver_init(problem, &mlcpOptions);
    info = mlcp_driver(problem, z, w, &mlcpOptions, &global_options);
    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions, &global_options);
    stopTimer();
    laux = sDt.mCumul;
    summary[itest].times[DIRECT_FB_ID][sIdWithSol] = sDt.mCumul;

    if (info > 0)
    {
      printf("Can't find a solution\n");
      strcpy(summary[itest].cv[DIRECT_FB_ID][sIdWithSol], "NO");
      sNbNOCV++;
      mlcp_compute_error(problem, z, w, tol1,  &error);
      printf("current point with error %lf \n", error);
      printSolution("DIRECT_FB", n, m, NbLines, z, w);
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      strcpy(summary[itest].cv[DIRECT_FB_ID][sIdWithSol], "CV");
      printf("find a solution with error %lf \n", error);
      printSolution("DIRECT_FB", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);

  }


}

void test_matrix(void)
{

  FILE *MLCPfile;

  int i, j;
  int isol;
  int n , n2;
  int m, m2;
  int NbLines;
  int withSol = 0;

  double *a, *b, *sol, *z, *w;
  double *vecA, *vecB, *vecC, *vecD, *vecM, *vecQ;
  NumericsMatrix M;

  char val[128];

  int iter;
  double criteria;
  MixedLinearComplementarityProblem problem;
  problem.n = 0;
  problem.m = 0;
  problem.q = 0;
  problem.A = 0;
  problem.B = 0;
  problem.C = 0;
  problem.D = 0;
  problem.a = 0;
  problem.b = 0;
  problem.problemType = 0;


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

  for (itest = 0 ; itest < NBTEST ; itest++)
  {

    switch (itest)
    {
    case 0:
      printf("\n\n 2x2 MLCP **************************************************************************");
      strcpy(summary[itest].file, "2x2 MLCP");
      if ((MLCPfile = fopen("data/deudeu_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: deudeu_mlcp.dat");
        exit(1);
      }
      break;
    case 1:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n PD **************************************************************************");
      strcpy(summary[itest].file, "PD");
      if ((MLCPfile = fopen("data/PD_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: PD_mlcp.dat");
        exit(1);
      }
      break;
    case 2:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n m2n1_mlcp.dat **************************************************************************");
      strcpy(summary[itest].file, "m2n1_mlcp");
      if ((MLCPfile = fopen("data/m2n1_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: m2n1SOL_mlcp.dat");
        exit(1);
      }
      break;
    case 3:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n relay2_mlcp.dat **************************************************************************");
      strcpy(summary[itest].file, "relay2_mlcp");
      if ((MLCPfile = fopen("data/relay2_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: relay2_mlcp.dat");
        exit(1);
      }
      break;
      //     case 4:
      //       printf("BEGIN A NEWTEST  **************************************************************************");
      //       printf("\n\n simple_mlcp.dat **************************************************************************");
      //       strcpy(summary[itest].file,"simple_mlcp");
      //       if( ( MLCPfile = fopen( "data/simple_mlcp.dat","r" ) ) == NULL ){
      //  perror("fopen MLCPfile: simple_mlcp.dat");
      //  exit(1);
      //       }
      //       break;
    case 11:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n deltasigma2_mlcp.dat **************************************************************************");
      strcpy(summary[itest].file, "deltasigma2_mlcp");
      if ((MLCPfile = fopen("data/deltasigma2_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: deltasigma2_mlcp.dat");
        exit(1);
      }
      break;
    case 12:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n relay3_mlcp.dat **************************************************************************");
      strcpy(summary[itest].file, "relay3_mlcp");
      if ((MLCPfile = fopen("data/relay3_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: relay3_mlcp.dat");
        exit(1);
      }
      break;
    case 13:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n deltasigma_mlcp.dat **************************************************************************");
      strcpy(summary[itest].file, "deltasigma_mlcp");
      if ((MLCPfile = fopen("data/deltasigma_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: deltasigma_mlcp.dat");
        exit(1);
      }
      break;
    case 14:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n m3n2_mlcp.dat **************************************************************************");
      strcpy(summary[itest].file, "m3n2_mlcp");
      if ((MLCPfile = fopen("data/m3n2_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: m3n2SOL_mlcp.dat");
        exit(1);
      }
      break;
    case 9:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n PDSym_mlcp.dat **************************************************************************");
      strcpy(summary[itest].file, "PDSym_mlcp");
      if ((MLCPfile = fopen("data/PDSym_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: PDSym_mlcp.dat");
        exit(1);
      }
      break;
    case 4:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n RLCD_mlcp MLCP **************************************************************************");
      strcpy(summary[itest].file, "RLCD_mlcp MLCP");
      if ((MLCPfile = fopen("data/RLCD_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: RLCD_mlcp.dat");
        exit(1);
      }
      break;
    case 16:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n Rectangular_1_n0_m4_mlcp MLCP **************************************************************************");
      strcpy(summary[itest].file, "Rectangular_1_n0_m4_mlcp MLCP");
      if ((MLCPfile = fopen("data/Rectangular_1_n0_m4_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: Rectangular_1_n0_m4_mlcp.dat");
        exit(1);
      }
      break;
    case 17:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n Rectangular_mlcp MLCP **************************************************************************");
      strcpy(summary[itest].file, "Rectangular_mlcp MLCP");
      if ((MLCPfile = fopen("data/Rectangular_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: Rectangular_mlcp.dat");
        exit(1);
      }
      break;
    case 6:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n diodeBridge MLCP **************************************************************************");
      strcpy(summary[itest].file, "diodeBridge MLCP");
      if ((MLCPfile = fopen("data/diodeBridge_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: diodeBridge_mlcp.dat");
        exit(1);
      }
      break;
    case 7:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n RCD **************************************************************************");
      strcpy(summary[itest].file, "RCD");
      if ((MLCPfile = fopen("data/RCD_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: RCD_mlcp.dat");
        exit(1);
      }
      break;
    case 5:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n diodeBridge 20 MLCP **************************************************************************");
      strcpy(summary[itest].file, "diodeBridge 20 MLCP");
      if ((MLCPfile = fopen("data/diodeBridge20_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: diodeBridge20_mlcp.dat");
        exit(1);
      }
      break;
    case 18:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n diodeBridge 40 MLCP **************************************************************************");
      strcpy(summary[itest].file, "diodeBridge 40 MLCP");
      if ((MLCPfile = fopen("data/diodeBridge40_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: diodeBridge40_mlcp.dat");
        exit(1);
      }
      break;
    case 8:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n Buck converter **************************************************************************");
      strcpy(summary[itest].file, "BuckConverter_mlcp");
      if ((MLCPfile = fopen("data/BuckConverter_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: BuckConverter_mlcp.dat");
        exit(1);
      }
      break;
    case 19:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n Buck converter FB failed **************************************************************************");
      strcpy(summary[itest].file, "Buck2");
      if ((MLCPfile = fopen("data/Buck2_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: Buck2_mlcp.dat");
        exit(1);
      }
      break;
    case 20:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n Buck converter First step **************************************************************************");
      strcpy(summary[itest].file, "BuckFirstStep");
      if ((MLCPfile = fopen("data/BuckFirstStep_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: BuckFirstStep_mlcp.dat");
        exit(1);
      }
      break;
    case 15:
      printf("BEGIN A NEWTEST  **************************************************************************");
      printf("\n\n Relay **************************************************************************");
      strcpy(summary[itest].file, "relay_mlcp");
      if ((MLCPfile = fopen("data/relay_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: relay_mlcp.dat");
        exit(1);
      }
      break;

    default :
      exit(1);
    }

    fscanf(MLCPfile , "%d" , &n);
    fscanf(MLCPfile , "%d" , &m);
    fscanf(MLCPfile , "%d" , &NbLines);

    n2 = n * n;
    m2 = m * m;
    isol = 1;

    vecM = (double*)malloc((n + m) * (NbLines) * sizeof(double));
    vecQ = (double*)malloc((NbLines) * sizeof(double));
    z = (double*)calloc((n + m), sizeof(double));
    w = (double*)calloc((NbLines), sizeof(double));
    vecA = (double*)malloc(n * (NbLines - m) * sizeof(double));
    vecB = (double*)malloc(m2 * sizeof(double));
    vecC = (double*)malloc((NbLines - m) * m * sizeof(double));
    vecD = (double*)malloc(m * n * sizeof(double));
    a    = (double*)malloc((NbLines - m) * sizeof(double));
    b    = (double*)malloc(m * sizeof(double));
    sol  = (double*)malloc((n + m + m) * sizeof(double));

    M.storageType = 0;
    M.matrix0 = vecM;
    problem.M = &M;
    problem.problemType = 0;
    problem.q = vecQ;
    problem.A = vecA;
    problem.B = vecB;
    problem.C = vecC;
    problem.D = vecD;
    problem.a = a;
    problem.b = b;

    problem.n = n;
    problem.m = m;
    M.size0 = NbLines;
    M.size1 = n + m;



    for (i = 0 ; i < NbLines - m ; ++i)
    {
      for (j = 0 ; j < n ; ++j)
      {
        fscanf(MLCPfile, "%s", val);
        vecA[(NbLines - m)*j + i ] = atof(val);
        vecM[(NbLines)*j + i ] = atof(val);
      }
    }
    for (i = 0 ; i < m ; ++i)
    {
      for (j = 0 ; j < m ; ++j)
      {
        fscanf(MLCPfile, "%s", val);
        vecB[ m * j + i ] = atof(val);
        /*  vecM[ n*(m+n)+(n+m)*j+n+i ] = atof(val);*/
        vecM[ n * (NbLines) + (NbLines)*j + (NbLines - m) + i ] = atof(val);

      }
    }
    for (i = 0 ; i < NbLines - m ; ++i)
    {
      for (j = 0 ; j < m ; ++j)
      {
        fscanf(MLCPfile, "%s", val);
        vecC[(NbLines - m)*j + i ] = atof(val);
        vecM[(NbLines) * (n + j) + i ] = atof(val);
      }
    }
    for (i = 0 ; i < m ; ++i)
    {
      for (j = 0 ; j < n ; ++j)
      {
        fscanf(MLCPfile, "%s", val);
        vecD[ m * j + i ] = atof(val);
        vecM[(NbLines)*j + i + (NbLines - m) ] = atof(val);
      }
    }

    for (i = 0 ; i < NbLines - m ; ++i)
    {
      fscanf(MLCPfile , "%s" , val);
      a[i] = atof(val);
      vecQ[i] = atof(val);
    }
    for (i = 0 ; i < m ; ++i)
    {
      fscanf(MLCPfile , "%s" , val);
      b[i] = atof(val);
      vecQ[i + NbLines - m] = atof(val);
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
    printf("\n");
    for (i = 0; i < NbLines; i++)
    {
      for (j = 0; j < n + m; j++)
        printf("%f ", vecM[NbLines * j + i]);
      printf("\n");
    }
    fclose(MLCPfile);

#ifdef BAVARD
    printf("\n With exact solution : ");
    printf("\n ------------------- : \n");
    if (isol) for (i = 0 ; i < (n + m + m) ; ++i) printf(" %10.4g " , sol[i]);
    else printf(" unknown ");
    printf("\n");
#endif
    sIdWithSol = 1;
    /*ONLY FOR DEBUG    NSNN_thisIsTheSolution(n+m,sol);*/
    if (withSol)
    {
      test_mlcp_series(&problem, z, w, sol);
      printf("\n Without exact solution : ");
      printf("\n ---------------------- : \n");
    }

    /*ONLY FOR DEBUG        NSNN_thisIsTheSolution(n+m,sol);*/
    for (i = 0; i < n + m + m; i++)
      sol[i] = 0;
    sIdWithSol = 0;

    test_mlcp_series(&problem, z, w, sol);

    free(sol);
    free(vecQ);
    free(vecM);
    free(vecA);
    free(vecB);
    free(vecC);
    free(vecD);
    free(a);
    free(b);
    free(z);
    free(w);
  }
  for (j = 0; j < 2; j++)
  {
    sIdWithSol = j;
    printf("* *** ************************ *** * \n");
    printf("* *** SUMMARY with solution=%d  *** * \n", j);
    printf("* *** ************************ *** * \n");

    for (itest = 0; itest < NBTEST ; itest++)
    {
      printf(" test %s  :\n", summary[itest].file);
      printf(" ===================================================  \n");
      if (sRunMethod[ENUM_ID])
        printf(" ENUM %s %ld \t", summary[itest].cv[ENUM_ID][sIdWithSol], summary[itest].times[ENUM_ID][sIdWithSol]);
      if (sRunMethod[PGS_IM_ID])
        printf(" PGS IM %s %ld \t", summary[itest].cv[PGS_IM_ID][sIdWithSol], summary[itest].times[PGS_IM_ID][sIdWithSol]);
      if (sRunMethod[PGS_EX_ID])
        printf(" PGS EX %s %ld \t", summary[itest].cv[PGS_EX_ID][sIdWithSol], summary[itest].times[PGS_EX_ID][sIdWithSol]);
      if (sRunMethod[RPGS_ID])
        printf(" RPGS %s %ld \t", summary[itest].cv[RPGS_ID][sIdWithSol], summary[itest].times[RPGS_ID][sIdWithSol]);
      if (sRunMethod[PSOR_05_ID])
        printf(" PSOR 05 %s %ld \t", summary[itest].cv[PSOR_05_ID][sIdWithSol], summary[itest].times[PSOR_05_ID][sIdWithSol]);
      if (sRunMethod[PSOR_1_ID])
        printf(" PSOR 1 %s %ld \t", summary[itest].cv[PSOR_1_ID][sIdWithSol], summary[itest].times[PSOR_1_ID][sIdWithSol]);
      if (sRunMethod[PSOR_15_ID])
        printf(" PSOR 15 %s %ld \t", summary[itest].cv[PSOR_15_ID][sIdWithSol], summary[itest].times[PSOR_15_ID][sIdWithSol]);
      if (sRunMethod[PSOR_2_ID])
        printf(" PSOR 2 %s %ld \t \n", summary[itest].cv[PSOR_2_ID][sIdWithSol], summary[itest].times[PSOR_2_ID][sIdWithSol]);
      if (sRunMethod[RPSOR_ID])
        printf(" RPSOR %s %ld \t", summary[itest].cv[RPSOR_ID][sIdWithSol], summary[itest].times[RPSOR_ID][sIdWithSol]);
      if (sRunMethod[PATH_ID])
        printf(" PATH %s %ld \t", summary[itest].cv[PATH_ID][sIdWithSol], summary[itest].times[PATH_ID][sIdWithSol]);
      if (sRunMethod[SIMPLEX_ID])
        printf(" SIMPLEX %s %ld \t", summary[itest].cv[SIMPLEX_ID][sIdWithSol], summary[itest].times[SIMPLEX_ID][sIdWithSol]);
      if (sRunMethod[DIRECT_ENUM_ID])
        printf(" DIR_ENUM %s %ld ", summary[itest].cv[DIRECT_ENUM_ID][sIdWithSol], summary[itest].times[DIRECT_ENUM_ID][sIdWithSol]);
      if (sRunMethod[FB_ID])
        printf(" FB %s %ld ", summary[itest].cv[FB_ID][sIdWithSol], summary[itest].times[FB_ID][sIdWithSol]);
      if (sRunMethod[DIRECT_FB_ID])
        printf(" DIR_FB %s %ld ", summary[itest].cv[DIRECT_FB_ID][sIdWithSol], summary[itest].times[DIRECT_FB_ID][sIdWithSol]);
      printf("\n");
    }
  }
  printf("* *** ******************** *** * \n");
  printf("* *** END OF TEST MATRIX   *** * \n");
  printf("* *** ******************** *** * \n");


}

int main(void)
{
  verbose = 0;
  int i;
  for (i = 0; i < NBMETHODS; i++)
    sRunMethod[i] = 1;
  sRunMethod[PATH_ID] = 0;
  //  sRunMethod[ENUM_ID]=1;
  sRunMethod[SIMPLEX_ID] = 0;

  test_matrix();
  printf("nb no cv %d\n", sNbNOCV);
  if (sNbNOCV > 90)
    return 0;
  else
    return 1;
}

