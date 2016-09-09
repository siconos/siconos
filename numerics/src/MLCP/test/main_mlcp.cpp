/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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

#pragma GCC diagnostic ignored "-Wmissing-declarations"

#include "SiconosConfig.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NonSmoothDrivers.h"
#include <sys/time.h>
#include <cassert>
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
  int cvState[NBMETHODS][2];
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

void initDataSummary()
{
  for (int i = 0; i < NBTEST; i++)
    for (int j = 0; j < NBMETHODS; j++)
    {
      summary[i].cvState[j][0] = 0;
      summary[i].cvState[j][1] = 0;
    }
}
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
  int info = -1;
  SolverOptions mlcpOptions;
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
    mlcpOptions.solverId = SICONOS_MLCP_ENUM;
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    //    mlcpOptions.iSize=1;
    //    mlcpOptions.dSize=1;
    mlcpOptions.dparam[0] = tol1;

    mlcp_driver_init(problem, &mlcpOptions);
    startTimer();
    info = 1;
    if (n + m < MAX_DIM_ENUM)
      info = mlcp_driver(problem, z, w, &mlcpOptions);
    stopTimer();
    summary[itest].times[ENUM_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[ENUM_ID][sIdWithSol], "CV");
    summary[itest].cvState[ENUM_ID][sIdWithSol] = 1;
    if (info > 0)
    {
      strcpy(summary[itest].cv[ENUM_ID][sIdWithSol], "NO");
      sNbNOCV++;
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      printSolution("ENUM", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }

  /*SOLVER PGS*/
  if (sRunMethod[PGS_IM_ID])
  {
    solTozw(n, m, z, w, sol);
    mlcpOptions.solverId = SICONOS_MLCP_PGS;
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    //"PGS"       , 101 , 1e-8 , 0.6 , 1.0 , 1 , 0 , 0.0 }; */
    mlcpOptions.dparam[0] = tol2;
    mlcp_driver_init(problem, &mlcpOptions);

    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions);
    stopTimer();
    summary[itest].times[PGS_IM_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[PGS_IM_ID][sIdWithSol], "CV");
    summary[itest].cvState[PGS_IM_ID][sIdWithSol] = 1;
    if (info > 0)
    {
      strcpy(summary[itest].cv[PGS_IM_ID][sIdWithSol], "NO");
      sNbNOCV++;
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      printSolution("PGS", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }

  /*SOLVER PGS*/
  if (sRunMethod[PGS_EX_ID])
  {

    solTozw(n, m, z, w, sol);
    mlcpOptions.solverId = SICONOS_MLCP_PGS;
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    //"PGS"       , 101 , 1e-8 , 0.6 , 1.0 , 1 , 0 , 0.0 }; */
    mlcpOptions.dparam[0] = tol2;

    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions);
    stopTimer();
    summary[itest].times[PGS_EX_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[PGS_EX_ID][sIdWithSol], "CV");
    summary[itest].cvState[PGS_EX_ID][sIdWithSol] = 1;
    if (info > 0)
    {
      strcpy(summary[itest].cv[PGS_EX_ID][sIdWithSol], "NO");
      sNbNOCV++;
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      printSolution("PGS", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER RPGS*/
  if (sRunMethod[RPGS_ID])
  {
    solTozw(n, m, z, w, sol);
    mlcpOptions.solverId = SICONOS_MLCP_RPGS;
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    mlcpOptions.dparam[0] = tol2;
    mlcpOptions.dparam[2] = 0.5; /*rho*/
    mlcp_driver_init(problem, &mlcpOptions);

    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions);
    stopTimer();
    summary[itest].times[RPGS_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[RPGS_ID][sIdWithSol], "CV");
    summary[itest].cvState[RPGS_ID][sIdWithSol] = 1;
    if (info > 0)
    {
      strcpy(summary[itest].cv[RPGS_ID][sIdWithSol], "NO");
      sNbNOCV++;
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      printSolution("RPGS", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER PSOR*/
  if (sRunMethod[_ID])
  {
    ohmega = 0;
    _ID = PSOR_05_ID - 1;
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    for (int cmp = 0; cmp < 4; cmp++)
    {
      _ID++;
      ohmega += 0.5;

      solTozw(n, m, z, w, sol);
      mlcpOptions.solverId = SICONOS_MLCP_PSOR;
      mlcpOptions.dparam[0] = tol2;
      mlcpOptions.dparam[2] = ohmega;
      mlcp_driver_init(problem, &mlcpOptions);

      startTimer();
      info = mlcp_driver(problem, z, w, &mlcpOptions);
      stopTimer();
      summary[itest].times[_ID][sIdWithSol] = sDt.mCumul;
      strcpy(summary[itest].cv[_ID][sIdWithSol], "CV");
      summary[itest].cvState[_ID][sIdWithSol] = 1;
      if (info > 0)
      {
        strcpy(summary[itest].cv[_ID][sIdWithSol], "NO");
        sNbNOCV++;
      }
      else
      {
        mlcp_compute_error(problem, z, w, tol1,  &error);
        printSolution("PSOR", n, m, NbLines, z, w);
      }
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER RPSOR*/
  if (sRunMethod[RPSOR_ID])
  {
    solTozw(n, m, z, w, sol);
    mlcpOptions.solverId = SICONOS_MLCP_RPSOR;
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    mlcpOptions.dparam[0] = tol2;
    mlcpOptions.dparam[2] = 0.5; /*rho*/
    mlcpOptions.dparam[3] = 2; /*ohmega*/
    mlcp_driver_init(problem, &mlcpOptions);

    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions);
    stopTimer();
    summary[itest].times[RPSOR_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[RPSOR_ID][sIdWithSol], "CV");
    summary[itest].cvState[RPSOR_ID][sIdWithSol] = 1;
    if (info > 0)
    {
      strcpy(summary[itest].cv[RPSOR_ID][sIdWithSol], "NO");
      sNbNOCV++;
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      printSolution("RPSOR", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER PATH*/
  if (sRunMethod[PATH_ID])
  {
    solTozw(n, m, z, w, sol);
    mlcpOptions.solverId = SICONOS_MLCP_PATH;
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_init(problem, &mlcpOptions);

    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions);
    stopTimer();
    summary[itest].times[PATH_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[PATH_ID][sIdWithSol], "CV");
    summary[itest].cvState[PATH_ID][sIdWithSol] = 1;
    if (info > 0)
    {
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
      printSolution("PATH", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER SIMPLEX*/
  if (sRunMethod[SIMPLEX_ID])
  {


    solTozw(n, m, z, w, sol);
    mlcpOptions.solverId = SICONOS_MLCP_SIMPLEX;
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    mlcpOptions.iparam[0] = 1000000;
    mlcpOptions.iparam[1] = 1;
    mlcpOptions.dparam[0] = 1e-12;
    mlcpOptions.dparam[1] = 1e-12;
    mlcpOptions.dparam[2] = 1e-9;
    mlcp_driver_init(problem, &mlcpOptions);

    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions);
    stopTimer();
    summary[itest].times[SIMPLEX_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[SIMPLEX_ID][sIdWithSol], "CV");
    summary[itest].cvState[SIMPLEX_ID][sIdWithSol] = 1;
    if (info > 0)
    {
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
      printSolution("SIMPLEX", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER DIRECT ENUM*/
  if (sRunMethod[DIRECT_ENUM_ID])
  {
    solTozw(n, m, z, w, sol);
    mlcpOptions.solverId = SICONOS_MLCP_DIRECT_ENUM;
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    mlcpOptions.dparam[0] = 1e-12;
    mlcp_driver_init(problem, &mlcpOptions);

    if (n + m < MAX_DIM_ENUM)
      info = mlcp_driver(problem, z, w, &mlcpOptions);
    if (info == 0)
    {
      startTimer();
      if (n + m < MAX_DIM_ENUM)
        info = mlcp_driver(problem, z, w, &mlcpOptions);
      stopTimer();
      summary[itest].times[DIRECT_ENUM_ID][sIdWithSol] = sDt.mCumul;
      strcpy(summary[itest].cv[DIRECT_ENUM_ID][sIdWithSol], "CV");
      summary[itest].cvState[DIRECT_ENUM_ID][sIdWithSol] = 1;
    }
    if (info > 0)
    {
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
      printSolution("DIRECT_ENUM_ID", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER FB*/
  if (sRunMethod[FB_ID])
  {

    solTozw(n, m, z, w, sol);
    mlcpOptions.solverId = SICONOS_MLCP_FB;
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    mlcpOptions.dparam[0] = 1e-10;
    mlcpOptions.dparam[1] = 0;


    mlcp_driver_init(problem, &mlcpOptions);
    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions);
    stopTimer();
    summary[itest].times[FB_ID][sIdWithSol] = sDt.mCumul;
    strcpy(summary[itest].cv[FB_ID][sIdWithSol], "CV");
    summary[itest].cvState[FB_ID][sIdWithSol] = 1;
    if (info > 0)
    {
      strcpy(summary[itest].cv[FB_ID][sIdWithSol], "NO");
      sNbNOCV++;
      mlcp_compute_error(problem, z, w, tol1,  &error);
      printSolution("FB", n, m, NbLines, z, w);
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      //    if (error > 1e-9)
      strcpy(summary[itest].cv[FB_ID][sIdWithSol], "CV");
      summary[itest].cvState[FB_ID][sIdWithSol] = 1;
      printSolution("FB", n, m, NbLines, z, w);
    }
    mixedLinearComplementarity_deleteDefaultSolverOptions(problem, &mlcpOptions);
    mlcp_driver_reset(problem, &mlcpOptions);
  }
  /*SOLVER DIRECT_FB*/
  if (sRunMethod[DIRECT_FB_ID])
  {

    solTozw(n, m, z, w, sol);
    mlcpOptions.solverId = SICONOS_MLCP_DIRECT_FB;
    mixedLinearComplementarity_setDefaultSolverOptions(problem, &mlcpOptions);
    mlcpOptions.iparam[0] = 500;
    mlcpOptions.iparam[1] = 0;

    mlcpOptions.dparam[0] = 1e-11;
    mlcpOptions.dparam[1] = 0;
    info = 1;

    mlcp_driver_init(problem, &mlcpOptions);
    info = mlcp_driver(problem, z, w, &mlcpOptions);
    startTimer();
    info = mlcp_driver(problem, z, w, &mlcpOptions);
    stopTimer();
    laux = sDt.mCumul;
    summary[itest].times[DIRECT_FB_ID][sIdWithSol] = sDt.mCumul;

    if (info > 0)
    {
      strcpy(summary[itest].cv[DIRECT_FB_ID][sIdWithSol], "NO");
      sNbNOCV++;
      mlcp_compute_error(problem, z, w, tol1,  &error);
      printSolution("DIRECT_FB", n, m, NbLines, z, w);
    }
    else
    {
      mlcp_compute_error(problem, z, w, tol1,  &error);
      strcpy(summary[itest].cv[DIRECT_FB_ID][sIdWithSol], "CV");
      summary[itest].cvState[DIRECT_FB_ID][sIdWithSol] = 1;
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
  // int n , n2;
  // int m, m2;
  // int NbLines;

  int n, m ;
  int withSol = 0;

  // double *a, *b,*sol,*z,*w;
  // double *vecA, *vecB, *vecC, *vecD, *vecM, *vecQ;
  // NumericsMatrix M;

  double *sol, *z, *w;


  char val[128];

  int iter;
  double criteria;



#ifdef BAVARD
  printf("* *** ******************** *** * \n");
  printf("* *** STARTING TEST MATRIX *** * \n");
  printf("* *** ******************** *** * \n");
#endif

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
    case 5:
      strcpy(summary[itest].file, "2x2 MLCP");
      if ((MLCPfile = fopen("data/deudeu_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: deudeu_mlcp.dat");
        exit(1);
      }
      break;
    case 1:
      strcpy(summary[itest].file, "PD");
      if ((MLCPfile = fopen("data/PD_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: PD_mlcp.dat");
        exit(1);
      }
      break;
    case 2:
      strcpy(summary[itest].file, "m2n1_mlcp");
      if ((MLCPfile = fopen("data/m2n1_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: m2n1SOL_mlcp.dat");
        exit(1);
      }
      break;
    case 3:
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
      strcpy(summary[itest].file, "deltasigma2_mlcp");
      if ((MLCPfile = fopen("data/deltasigma2_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: deltasigma2_mlcp.dat");
        exit(1);
      }
      break;
    case 12:
      strcpy(summary[itest].file, "relay3_mlcp");
      if ((MLCPfile = fopen("data/relay3_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: relay3_mlcp.dat");
        exit(1);
      }
      break;
    case 13:
      strcpy(summary[itest].file, "deltasigma_mlcp");
      if ((MLCPfile = fopen("data/deltasigma_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: deltasigma_mlcp.dat");
        exit(1);
      }
      break;
    case 14:
      strcpy(summary[itest].file, "m3n2_mlcp");
      if ((MLCPfile = fopen("data/m3n2_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: m3n2SOL_mlcp.dat");
        exit(1);
      }
      break;
    case 9:
      strcpy(summary[itest].file, "PDSym_mlcp");
      if ((MLCPfile = fopen("data/PDSym_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: PDSym_mlcp.dat");
        exit(1);
      }
      break;
    case 4:
      strcpy(summary[itest].file, "RLCD_mlcp MLCP");
      if ((MLCPfile = fopen("data/RLCD_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: RLCD_mlcp.dat");
        exit(1);
      }
      break;
    case 6:
      strcpy(summary[itest].file, "diodeBridge MLCP");
      if ((MLCPfile = fopen("data/diodeBridge_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: diodeBridge_mlcp.dat");
        exit(1);
      }
      break;
    case 7:
      strcpy(summary[itest].file, "RCD");
      if ((MLCPfile = fopen("data/RCD_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: RCD_mlcp.dat");
        exit(1);
      }
      break;
    case 0:
      strcpy(summary[itest].file, "diodeBridge 20 MLCP");
      if ((MLCPfile = fopen("data/deudeu_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: diodeBridge20_mlcp.dat");
        exit(1);
      }
      break;
    case 8:
      strcpy(summary[itest].file, "BuckConverter_mlcp");
      if ((MLCPfile = fopen("data/BuckConverter_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: BuckConverter_mlcp.dat");
        exit(1);
      }
      break;
    case 15:
      strcpy(summary[itest].file, "Rectangular_1_n0_m4_mlcp MLCP");
      if ((MLCPfile = fopen("data/Rectangular_1_n0_m4_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: Rectangular_1_n0_m4_mlcp.dat");
        exit(1);
      }
      break;
    case 16:
      strcpy(summary[itest].file, "Rectangular_mlcp MLCP");
      if ((MLCPfile = fopen("data/Rectangular_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: Rectangular_mlcp.dat");
        exit(1);
      }
      break;

    case 17:
      strcpy(summary[itest].file, "diodeBridge 40 MLCP");
      if ((MLCPfile = fopen("data/diodeBridge40_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: diodeBridge40_mlcp.dat");
        exit(1);
      }
      break;

    case 18:
      strcpy(summary[itest].file, "Buck2");
      if ((MLCPfile = fopen("data/Buck2_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: Buck2_mlcp.dat");
        exit(1);
      }
      break;
    case 19:
      strcpy(summary[itest].file, "BuckFirstStep");
      if ((MLCPfile = fopen("data/BuckFirstStep_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: BuckFirstStep_mlcp.dat");
        exit(1);
      }
      break;
      // case 20:
      //   strcpy(summary[itest].file,"relay_mlcp");
      //   if((MLCPfile = fopen("data/relay_mlcp.dat","r")) == NULL)
      //   {
      //     perror("fopen MLCPfile: relay_mlcp.dat");
      //     exit(1);
      //   }
      //   break;
    default :
      exit(1);
    }

    MixedLinearComplementarityProblem * problem = (MixedLinearComplementarityProblem *) malloc(sizeof(MixedLinearComplementarityProblem));

    mixedLinearComplementarity_newFromFile(problem, MLCPfile);
    //mixedLinearComplementarity_newFromFileOld(problem, MLCPfile);
    printf("\n");
    printf("====================== \n")     ;
    printf("test on  = %s\n", summary[itest].file);

    // n2 = n*n;
    // m2 = m*m;
    // isol = 1;

    n = problem->n;
    m = problem->m;
    assert(n>0);
    assert(m>0);

    z = (double*)calloc((n + m), sizeof(double));
    w = (double*)calloc((n + m), sizeof(double));

    sol  = (double*)malloc((n + m + m) * sizeof(double));




    // for(i = 0 ; i < NbLines-m ; ++i)
    // {
    //   for(j = 0 ; j < n ; ++j)
    //   {
    //     fscanf(MLCPfile,"%s",val);
    //     vecA[(NbLines-m)*j+i ] = atof(val);
    //     vecM[(NbLines)*j+i ] = atof(val);
    //   }
    // }
    // for(i = 0 ; i < m ; ++i)
    // {
    //   for(j = 0 ; j < m ; ++j)
    //   {
    //     fscanf(MLCPfile,"%s",val);
    //     vecB[ m*j+i ] = atof(val);
    //     /* vecM[ n*(m+n)+(n+m)*j+n+i ] = atof(val);*/
    //     vecM[ n*(NbLines)+(NbLines)*j+(NbLines-m)+i ] = atof(val);

    //   }
    // }
    // for(i = 0 ; i < NbLines-m ; ++i)
    // {
    //   for(j = 0 ; j < m ; ++j)
    //   {
    //     fscanf(MLCPfile,"%s",val);
    //     vecC[(NbLines-m)*j+i ] = atof(val);
    //     vecM[(NbLines)*(n+j)+i ] = atof(val);
    //   }
    // }
    // for(i = 0 ; i < m ; ++i)
    // {
    //   for(j = 0 ; j < n ; ++j)
    //   {
    //     fscanf(MLCPfile,"%s",val);
    //     vecD[ m*j+i ] = atof(val);
    //     vecM[(NbLines)*j+i+(NbLines-m) ] = atof(val);
    //   }
    // }

    // for(i = 0 ; i < NbLines-m ; ++i)
    // {
    //   fscanf(MLCPfile , "%s" , val);
    //   a[i] = atof(val);
    //   vecQ[i] = atof(val);
    // }
    // for(i = 0 ; i < m ; ++i)
    // {
    //   fscanf(MLCPfile , "%s" , val);
    //   b[i] = atof(val);
    //   vecQ[i+NbLines-m] = atof(val);
    // }



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
#ifdef BAVARD
    printf("\n");
    for (i = 0; i < n + m; i++)
    {
      for (j = 0; j < n + m; j++)
        printf("%f ", problem->M->matrix0[(n + m)*j + i]);
      printf("\n");
    }
#endif
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
      test_mlcp_series(problem, z, w, sol);
#ifdef BAVARD
      printf("\n Without exact solution : ");
      printf("\n ---------------------- : \n");
#endif
    }

    /*ONLY FOR DEBUG        NSNN_thisIsTheSolution(n+m,sol);*/
    for (i = 0; i < n + m + m; i++)
      sol[i] = 0;
    sIdWithSol = 0;

    test_mlcp_series(problem, z, w, sol);

    free(sol);
    free(z);
    free(w);
    freeMixedLinearComplementarityProblem(problem);
  }
#ifdef BAVARD
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
  for (itest = 0; itest < NBTEST ; itest++)
  {
    for (int imethod = 0; imethod < NBMETHODS; imethod++)
      printf("%i %i ", summary[itest].cvState[imethod][0], summary[itest].cvState[imethod][1]);
    printf("\n");
  }
#endif

}

int main(void)
{
  verbose = 0;
  int i;
  initDataSummary();
  for (i = 0; i < NBMETHODS; i++)
    sRunMethod[i] = 1;
  sRunMethod[PATH_ID] = 0;
  //  sRunMethod[ENUM_ID]=1;
  sRunMethod[SIMPLEX_ID] = 0;

  test_matrix();
#ifdef BAVARD
  printf("nb no cv %d\n", sNbNOCV);
#endif
  int T[280] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
                1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
                1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
                1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0
               };
  int cmp = 0;
  for (itest = 0; itest < NBTEST ; itest++)
  {
    for (int imethod = 0; imethod < NBMETHODS; imethod++)
    {
      if (T[cmp] < summary[itest].cvState[imethod][0])
      {
#ifdef BAVARD
        printf("Warning, test=%i method=%i failed.", itest, imethod);
#endif
        return 1;//failed
      }
      cmp++;
      if (T[cmp] < summary[itest].cvState[imethod][1])
      {
#ifdef BAVARD
        printf("Warning, test=%i method=%i failed.", itest, imethod);
#endif
        return 1;//failed
      }
      cmp++;
    }
  }
#ifdef BAVARD
  printf("Test succed\n");
#endif
  return 0;
}

