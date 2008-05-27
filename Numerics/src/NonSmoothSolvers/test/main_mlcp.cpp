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

#define BAVARD
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
#define NBMETHODS 12

#define PATH_DRIVER
#define MAX_DIM_ENUM 15
/*#ifdef PATH_DRIVER
const unsigned short int *__ctype_b;
const __int32_t *__ctype_tolower ;
#endif
*/

typedef struct
{
  char file[124];
  char cv[5][NBMETHODS];
  int nbSteps[NBMETHODS];
  long times[NBMETHODS];
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
static int itest;

/*
 ******************************************************************************
 */
void printSolution(char *name, int n, int m, double *z, double *w)
{
  int i;
#ifdef BAVARD
  printf(" *** ************************************** ***\n");
  printf(" ****** z = ********************************\n");
  for (i = 0 ; i < n + m ; i++)
    printf("%s: %14.7e  \n", name, z[i]);
  printf(" ****** w = ********************************\n");
  for (i = 0 ; i < m ; i++)
    printf("%s: %14.7e  \n", name, w[i + n]);
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

void test_mlcp_series(MixedLinearComplementarity_Problem* problem, double *z, double *w, double *sol)
{
  int info;
  Solver_Options mlcpOptions;
  Numerics_Options global_options;
  double tol1 = 1e-16;
  double tol2 = 1e-6;
  double error = 0;
  int n = problem->n;
  int m = problem->m;
  double ohmega = 0;
  int _ID = PSOR_05_ID - 1;


  mlcpOptions.iparam = (int*)malloc(10 * sizeof(int));
  mlcpOptions.dparam = (double*)malloc(10 * sizeof(double));

  /*SOLVER ENUM*/
  //   solTozw(n,m,z,w,sol);
  //   strcpy(mlcpOptions.solverName,"ENUM");
  //   mlcpOptions.iSize=1;
  //   mlcpOptions.iparam[0]=0;/*verbose 1*/
  //   mlcpOptions.dSize=1;
  //   mlcpOptions.dparam[0]=tol1;

  //   mlcpOptions.dWork=(double*)malloc( mlcp_enum_getNbDWork(problem,0)*sizeof(double));
  //   mlcpOptions.iWork=(int*)malloc(mlcp_enum_getNbIWork(problem,0)*sizeof(int));
  //   startTimer();
  //   info=1;
  //   if ( n+m < MAX_DIM_ENUM)
  //     info = mlcp_driver( problem,z,w, &mlcpOptions,&global_options);
  //   stopTimer();
  //   summary[itest].times[ENUM_ID]=sDt.mCumul;
  //   strcpy(summary[itest].cv[ENUM_ID],"CV");
  //   if (info>0) {
  //     printf("Can't find a solution\n");
  //     strcpy(summary[itest].cv[ENUM_ID],"NO");
  //   }else{
  //     mlcp_compute_error(problem, z, w, tol1,  &error);
  //     printf("find a solution with error %lf \n",error);
  //     printSolution("ENUM",n,m,z,w);
  //   }
  //   free(mlcpOptions.dWork);
  //   free(mlcpOptions.iWork);
  // /*SOLVER PGS*/
  //   solTozw(n,m,z,w,sol);
  //   strcpy(mlcpOptions.solverName,"PGS");
  //   //"PGS"       , 101 , 1e-8 , 0.6 , 1.0 , 1 , 0 , 0.0 }; */
  //   mlcpOptions.iSize=2;
  //   mlcpOptions.iparam[0]=101;
  //   mlcpOptions.iparam[2]=0;//implicit
  //   mlcpOptions.dSize=2;
  //   mlcpOptions.dparam[0]=tol2;

  //   startTimer();
  //   info = mlcp_driver( problem,z,w, &mlcpOptions,&global_options);
  //   stopTimer();
  //   summary[itest].times[PGS_IM_ID]=sDt.mCumul;
  //   strcpy(summary[itest].cv[PGS_IM_ID],"CV");
  //   if (info>0) {
  //     printf("Can't find a solution\n");
  //     strcpy(summary[itest].cv[PGS_IM_ID],"NO");
  //   }else{
  //     mlcp_compute_error(problem, z, w, tol1,  &error);
  //     printf("find a solution with error %lf \n",error);
  //     printSolution("PGS",n,m,z,w);
  //   }
  // /*SOLVER PGS*/
  //   solTozw(n,m,z,w,sol);
  //   strcpy(mlcpOptions.solverName,"PGS");
  //   //"PGS"       , 101 , 1e-8 , 0.6 , 1.0 , 1 , 0 , 0.0 }; */
  //   mlcpOptions.iSize=2;
  //   mlcpOptions.iparam[0]=101;
  //   mlcpOptions.iparam[2]=1;//explicit
  //   mlcpOptions.dSize=2;
  //   mlcpOptions.dparam[0]=tol2;

  //   startTimer();
  //   info = mlcp_driver( problem,z,w, &mlcpOptions,&global_options);
  //   stopTimer();
  //   summary[itest].times[PGS_EX_ID]=sDt.mCumul;
  //   strcpy(summary[itest].cv[PGS_EX_ID],"CV");
  //   if (info>0) {
  //     printf("Can't find a solution\n");
  //     strcpy(summary[itest].cv[PGS_EX_ID],"NO");
  //   }else{
  //     mlcp_compute_error(problem, z, w, tol1,  &error);
  //     printf("find a solution with error %lf \n",error);
  //     printSolution("PGS",n,m,z,w);
  //   }
  // /*SOLVER RPGS*/
  //   solTozw(n,m,z,w,sol);
  //   strcpy(mlcpOptions.solverName,"RPGS");
  //   mlcpOptions.iSize=2;
  //   mlcpOptions.iparam[0]=101;
  //   mlcpOptions.dSize=3;
  //   mlcpOptions.dparam[0]=tol2;
  //   mlcpOptions.dparam[2]=0.5;/*rho*/

  //   startTimer();
  //   info = mlcp_driver( problem,z,w, &mlcpOptions,&global_options);
  //   stopTimer();
  //   summary[itest].times[RPGS_ID]=sDt.mCumul;
  //   strcpy(summary[itest].cv[RPGS_ID],"CV");
  //   if (info>0) {
  //     printf("Can't find a solution\n");
  //     strcpy(summary[itest].cv[RPGS_ID],"NO");
  //   }else{
  //     mlcp_compute_error(problem, z, w, tol1,  &error);
  //     printf("find a solution with error %lf \n",error);
  //     printSolution("RPGS",n,m,z,w);
  //   }

  // /*SOLVER PSOR*/
  //   ohmega=0;
  //   _ID=PSOR_05_ID-1;
  //   for(int cmp=0;cmp<4;cmp++){
  //     _ID++;
  //     ohmega+=0.5;

  //     solTozw(n,m,z,w,sol);
  //     strcpy(mlcpOptions.solverName,"PSOR");
  //     mlcpOptions.iSize=2;
  //     mlcpOptions.iparam[0]=101;
  //     mlcpOptions.dSize=3;
  //     mlcpOptions.dparam[0]=tol2;
  //     mlcpOptions.dparam[2]=ohmega;

  //     startTimer();
  //     info = mlcp_driver( problem,z,w, &mlcpOptions,&global_options);
  //     stopTimer();
  //     summary[itest].times[_ID]=sDt.mCumul;
  //     strcpy(summary[itest].cv[_ID],"CV");
  //     if (info>0) {
  //       printf("Can't find a solution\n");
  //       strcpy(summary[itest].cv[_ID],"NO");
  //     }else{
  //       mlcp_compute_error(problem, z, w, tol1,  &error);
  //       printf("find a solution with error %lf \n",error);
  //       printSolution("PSOR",n,m,z,w);
  //     }
  //   }
  // /*SOLVER RPSOR*/
  //   solTozw(n,m,z,w,sol);
  //   strcpy(mlcpOptions.solverName,"RPSOR");
  //   mlcpOptions.iSize=2;
  //   mlcpOptions.iparam[0]=101;
  //   mlcpOptions.dSize=4;
  //   mlcpOptions.dparam[0]=tol2;
  //   mlcpOptions.dparam[2]=0.5;/*rho*/
  //   mlcpOptions.dparam[3]=2;/*ohmega*/

  //   startTimer();
  //   info = mlcp_driver( problem,z,w, &mlcpOptions,&global_options);
  //   stopTimer();
  //   summary[itest].times[RPSOR_ID]=sDt.mCumul;
  //   strcpy(summary[itest].cv[RPSOR_ID],"CV");
  //   if (info>0) {
  //     printf("Can't find a solution\n");
  //     strcpy(summary[itest].cv[RPSOR_ID],"NO");
  //   }else{
  //     mlcp_compute_error(problem, z, w, tol1,  &error);
  //     printf("find a solution with error %lf \n",error);
  //     printSolution("RPSOR",n,m,z,w);
  //   }
  /*SOLVER PATH*/
  solTozw(n, m, z, w, sol);
  strcpy(mlcpOptions.solverName, "PATH");
  mlcpOptions.iSize = 1;
  mlcpOptions.iparam[0] = 101;
  mlcpOptions.dSize = 1;
  mlcpOptions.dparam[0] = tol2;

  startTimer();
  info = mlcp_driver(problem, z, w, &mlcpOptions, &global_options);
  stopTimer();
  summary[itest].times[PATH_ID] = sDt.mCumul;
  strcpy(summary[itest].cv[PATH_ID], "CV");
  if (info > 0)
  {
    printf("Can't find a solution\n");
    strcpy(summary[itest].cv[PATH_ID], "NO");
  }
  else
  {
    mlcp_compute_error(problem, z, w, tol1,  &error);
    if (error > tol2)
      strcpy(summary[itest].cv[PATH_ID], "NO");
    printf("find a solution with error %lf \n", error);
    printSolution("PATH", n, m, z, w);
  }
  /*SOLVER SIMPLEX*/
  //   solTozw(n,m,z,w,sol);
  //   strcpy(mlcpOptions.solverName,"SIMPLEX");
  //   mlcpOptions.iSize=2;
  //   mlcpOptions.iparam[0]=1000000;
  //   mlcpOptions.iparam[1]=1;/*VERBOSE*/
  //   mlcpOptions.dSize=3;
  //   mlcpOptions.dparam[0]=1e-12;
  //   mlcpOptions.dparam[1]=1e-12;
  //   mlcpOptions.dparam[2]=1e-9;

  //   startTimer();
  //   info = mlcp_driver( problem,z,w, &mlcpOptions,&global_options);
  //   stopTimer();
  //   summary[itest].times[SIMPLEX_ID]=sDt.mCumul;
  //   strcpy(summary[itest].cv[SIMPLEX_ID],"CV");
  //   if (info>0) {
  //     printf("Can't find a solution\n");
  //     strcpy(summary[itest].cv[SIMPLEX_ID],"NO");
  //   }else{
  //     mlcp_compute_error(problem, z, w, tol1,  &error);
  //     if (error > 1e-9)
  //       strcpy(summary[itest].cv[SIMPLEX_ID],"NO");
  //     printf("find a solution with error %lf \n",error);
  //     printSolution("SIMPLEX",n,m,z,w);
  //   }
  // /*SOLVER DIRECT ENUM*/
  //   solTozw(n,m,z,w,sol);
  //   strcpy(mlcpOptions.solverName,"DIRECT_ENUM");
  //   mlcpOptions.iSize=6;
  //   mlcpOptions.iparam[5]=3;/*Number of registered configurations*/
  //   mlcpOptions.iparam[0]=0;/*VERBOSE*/
  //   mlcpOptions.dSize=6;
  //   mlcpOptions.dparam[0]=1e-12;
  //   mlcpOptions.dparam[5]=1e-12;
  //   mlcpOptions.dparam[6]=1e-12;
  //  int aux = mlcp_direct_enum_getNbDWork(problem,&mlcpOptions);
  //  double * daux = (double*)malloc( aux*sizeof(double));
  //   mlcpOptions.dWork=daux;
  //   aux = mlcp_direct_enum_getNbIWork(problem,&mlcpOptions);
  //   int * iaux = (int*)malloc(aux*sizeof(int));
  //   mlcpOptions.iWork=iaux;

  //   mlcp_direct_enum_init(problem,&mlcpOptions);
  //   printf("Solver direct, precompute solution.\n");
  //   info = mlcp_driver( problem,z,w, &mlcpOptions,&global_options);
  //   if (info==0) {
  //     printf("Solver direct, solution computed. ");
  //     printf("Now, run direct solver.\n");
  //     startTimer();
  //     info = mlcp_driver( problem,z,w, &mlcpOptions,&global_options);
  //     stopTimer();
  //     summary[itest].times[DIRECT_ENUM_ID]=sDt.mCumul;
  //     strcpy(summary[itest].cv[DIRECT_ENUM_ID],"CV");
  //   }
  //   if (info>0) {
  //     printf("Can't find a solution\n");
  //     strcpy(summary[itest].cv[DIRECT_ENUM_ID],"NO");
  //   }else{
  //     mlcp_compute_error(problem, z, w, tol1,  &error);
  //     if (error > 1e-9)
  //       strcpy(summary[itest].cv[DIRECT_ENUM_ID],"NO");
  //     printf("find a solution with error %lf \n",error);
  //     printSolution("DIRECT_ENUM_ID",n,m,z,w);
  //   }

  //  free( daux);
  //  free( iaux);
  mlcp_direct_enum_reset();
  deleteSolverOptions(&mlcpOptions);




}

void test_matrix(void)
{

  FILE *MLCPfile;

  int i, j;
  int isol;
  int n , n2;
  int m, m2;
  int withSol = 0;

  double *a, *b, *sol, *z, *w;
  double *vecA, *vecB, *vecC, *vecD, *vecM, *vecQ;
  NumericsMatrix M;

  char val[20];

  int iter;
  double criteria;
  MixedLinearComplementarity_Problem problem;

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
    case 9:
      printf("\n\n Buck converter ");
      strcpy(summary[itest].file, "BuckConverter_mlcp");
      if ((MLCPfile = fopen("MATRIX/BuckConverter_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: BuckConverter_mlcp.dat");
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
    case 8:
      printf("\n\n diodeBridge 40 MLCP ");
      strcpy(summary[itest].file, "diodeBridge 40 MLCP");
      if ((MLCPfile = fopen("MATRIX/diodeBridge40_mlcp.dat", "r")) == NULL)
      {
        perror("fopen MLCPfile: diodeBridge40_mlcp.dat");
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

    vecM = (double*)malloc((n + m) * (n + m) * sizeof(double));
    vecQ = (double*)malloc((n + m) * sizeof(double));
    z = (double*)malloc((n + m) * sizeof(double));
    w = (double*)malloc((n + m) * sizeof(double));
    vecA = (double*)malloc(n2 * sizeof(double));
    vecB = (double*)malloc(m2 * sizeof(double));
    vecC = (double*)malloc(n * m * sizeof(double));
    vecD = (double*)malloc(m * n * sizeof(double));
    a    = (double*)malloc(n * sizeof(double));
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
    M.size0 = n + m;
    M.size1 = n + m;



    for (i = 0 ; i < n ; ++i)
    {
      for (j = 0 ; j < n ; ++j)
      {
        fscanf(MLCPfile, "%s", val);
        vecA[ n * j + i ] = atof(val);
        vecM[(n + m)*j + i ] = atof(val);
      }
    }
    for (i = 0 ; i < m ; ++i)
    {
      for (j = 0 ; j < m ; ++j)
      {
        fscanf(MLCPfile, "%s", val);
        vecB[ m * j + i ] = atof(val);
        vecM[ n * (m + n) + (n + m)*j + n + i ] = atof(val);


      }
    }
    for (i = 0 ; i < n ; ++i)
    {
      for (j = 0 ; j < m ; ++j)
      {
        fscanf(MLCPfile, "%s", val);
        vecC[ n * j + i ] = atof(val);
        vecM[(n + m) * (n + j) + i ] = atof(val);
      }
    }
    for (i = 0 ; i < m ; ++i)
    {
      for (j = 0 ; j < n ; ++j)
      {
        fscanf(MLCPfile, "%s", val);
        vecD[ m * j + i ] = atof(val);
        vecM[(n + m)*j + i + n ] = atof(val);
      }
    }

    for (i = 0 ; i < n ; ++i)
    {
      fscanf(MLCPfile , "%s" , val);
      a[i] = atof(val);
      vecQ[i] = atof(val);
    }
    for (i = 0 ; i < m ; ++i)
    {
      fscanf(MLCPfile , "%s" , val);
      b[i] = atof(val);
      vecQ[i + n] = atof(val);
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
      //      test_mlcp_series( n,m , vecA,vecB,vecC ,vecD,a,b,sol);
      printf("\n Without exact solution : ");
      printf("\n ---------------------- : \n");
    }
    for (i = 0; i < n + m + m; i++)
      sol[i] = 0;
    test_mlcp_series(&problem, z, w, sol);

    free(sol);
    free(vecA);
    free(vecB);
    free(vecC);
    free(vecD);
    free(a);
    free(b);
    free(z);
    free(w);
  }
  printf("* *** ******************** *** * \n");
  printf("* *** SUMMARY             *** * \n");
  printf("* *** ******************** *** * \n");

  for (itest = 0; itest < NBTEST ; ++itest)
  {
    printf(" test %s  :\n", summary[itest].file);
    printf(" ===================================================  \n");
    printf(" ENUM %s %d \t", summary[itest].cv[ENUM_ID], summary[itest].times[ENUM_ID]);
    printf(" PGS IM %s %d \t", summary[itest].cv[PGS_IM_ID], summary[itest].times[PGS_IM_ID]);
    printf(" PGS EX %s %d \t", summary[itest].cv[PGS_EX_ID], summary[itest].times[PGS_EX_ID]);
    printf(" RPGS %s %d \t", summary[itest].cv[RPGS_ID], summary[itest].times[RPGS_ID]);
    /*printf(" PSOR 05 %s %d \t",summary[itest].cv[PSOR_05_ID],summary[itest].times[PSOR_05_ID]);
    printf(" PSOR 1 %s %d \t",summary[itest].cv[PSOR_1_ID],summary[itest].times[PSOR_1_ID]);
    printf(" PSOR 15 %s %d \t",summary[itest].cv[PSOR_15_ID],summary[itest].times[PSOR_15_ID]);
    printf(" PSOR 2 %s %d \t",summary[itest].cv[PSOR_2_ID],summary[itest].times[PSOR_2_ID]);
    printf(" RPSOR %s %d \t",summary[itest].cv[RPSOR_ID],summary[itest].times[RPSOR_ID]);*/
    printf(" PATH %s %d \t", summary[itest].cv[PATH_ID], summary[itest].times[PATH_ID]);
    printf(" SIMPLEX %s %d \t", summary[itest].cv[SIMPLEX_ID], summary[itest].times[SIMPLEX_ID]);
    printf(" DIR_ENUM %s %d \n", summary[itest].cv[DIRECT_ENUM_ID], summary[itest].times[DIRECT_ENUM_ID]);
  }
  printf("* *** ******************** *** * \n");
  printf("* *** END OF TEST MATRIX   *** * \n");
  printf("* *** ******************** *** * \n");


}

int main(void)
{

  test_matrix();

  return 1;
}

