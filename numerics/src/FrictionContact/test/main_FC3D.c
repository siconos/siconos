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

/*
  Tests functions for NumericsMatrix structure

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NonSmoothDrivers.h"

#include <math.h>


int test_Series_nsgs(FrictionContactProblem* problem,
                     double tolerance, double localtolerance,
                     int nLocalSolver, int * localSolverList)
{
  int info = -1, i;


  SolverOptions numerics_solver_options;
  numerics_solver_options.filterOn = 0;
  numerics_solver_options.isSet = 1;

  (numerics_solver_options.solverId = SICONOS_FRICTION_3D_NSGS;

   numerics_solver_options.iSize = 5;
   numerics_solver_options.iparam = (int*)malloc(numerics_solver_options.iSize * sizeof(int));
   numerics_solver_options.dSize = 5;
   numerics_solver_options.dparam = (double*)malloc(numerics_solver_options.dSize * sizeof(double));


   int NC = problem->numberOfContacts;
   double *reaction = (double*)malloc(3 * NC * sizeof(double));
   double *velocity = (double*)malloc(3 * NC * sizeof(double));

   int nmax ; // Max number of iteration

   for (i = 0; i < nLocalSolver; i++)
{
  nmax = 100;
  // 0: projection on Cone, 1: Newton/AlartCurnier,  2: projection on Cone with local iteration, 3: projection on Disk  with diagonalization,
  for (int k = 0 ; k < 3 * NC; k++)
    {
      velocity[k] = 0.0;
      reaction[k] = 0.0;
    }
    numerics_solver_options.iparam[0] = nmax ;
    numerics_solver_options.iparam[1] = 0 ;
    numerics_solver_options.iparam[2] = 0 ;
    numerics_solver_options.iparam[3] = 0 ;
    numerics_solver_options.iparam[4] = localSolverList[i] ;

    numerics_solver_options.dparam[0] = tolerance ;
    numerics_solver_options.dparam[1] = 0.0 ;
    numerics_solver_options.dparam[2] = localtolerance ;
    numerics_solver_options.dparam[3] = 0.0 ;
    numerics_solver_options.dparam[4] = 0.0 ;

    info = fc3d_driver(problem,
		       reaction , velocity,
		       &numerics_solver_options);

    printf("\n");
    for (int k = 0 ; k < 3 * NC; k++)
    {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
    }
    printf("\n");
    printf("info =%i\n", info);
    /*  if (info) */
    /*      { */
    /*   char c; */
    /*   printf("Press a key to continue:"); */
    /*     scanf("%c",&c); */

    /*      } */

  }

  free(reaction);
  free(velocity);
  free(numerics_solver_options.iparam);
  free(numerics_solver_options.dparam);


  return info;
}

int main(void)
{

  printf("========= Starts Numerics tests for FrictionalContact_problem ========= \n");



  FrictionContactProblem * testproblem = NULL;

  int ntestfile = 4, i;
  char ** testfile = (char **)malloc(ntestfile * sizeof(char *));
  int ** localSolverList = (int **)malloc(ntestfile * sizeof(int *));
  int * nLocalSolver = (int *)malloc(ntestfile * sizeof(int));
  /*   for (int i =0 ; i< ntestfile; i++)  */
  /*       { */
  /*    testfile[i] = (char *)malloc(20*sizeof(char)); */
  /*       } */

  double*  tol = (double*)malloc(ntestfile * sizeof(double));
  double*  localtol = (double*)malloc(ntestfile * sizeof(double));
  int ilocal;
  ilocal = 0;
  testfile[0] = "./data/Example1_Fc3D_SBM.dat";
  tol[0] = 1e-16;
  localtol[0] = 1e-16;
  nLocalSolver[0] = 4;
  localSolverList[0] = (int *)malloc(nLocalSolver[0] * sizeof(int));
  for (i = 0; i < nLocalSolver[0] ; i++) localSolverList[0][i] = i;



  testfile[1] = "./data/Confeti-ex13-4contact-Fc3D-SBM.dat";
  tol[1] = 1e-5;
  localtol[1] = 1e-6;

  nLocalSolver[1] = 3;
  localSolverList[1] = (int *)malloc(nLocalSolver[1] * sizeof(int));
  for (i = 0; i < nLocalSolver[1] ; i++) localSolverList[1][i] = i;

  testfile[2] = "./data/Confeti-ex03-Fc3D-SBM.dat";
  tol[2] = 1e-5;
  localtol[2] = 1e-6;
  nLocalSolver[2] = 4;
  localSolverList[2] = (int *)malloc(nLocalSolver[2] * sizeof(int));
  for (i = 0; i < nLocalSolver[2] ; i++) localSolverList[2][i] = i;

  testfile[3] = "./data/Confeti-ex13-Fc3D-SBM.dat";
  tol[3] = 1e-5;
  localtol[3] = 1e-6;
  nLocalSolver[3] = 1;
  localSolverList[3] = (int *)malloc(nLocalSolver[3] * sizeof(int));
  /*   for (i=0; i<nLocalSolver[3] ; i++) localSolverList[3][i]=i;   */
  localSolverList[3][0] = 1;

  int info = -1;
  FILE * file;

  for (i = 0 ; i < ntestfile; i++)
  {
    printf("Test file number %i -- %s\n", i, testfile[i]);
    testproblem = (FrictionContactProblem*)malloc(sizeof(FrictionContactProblem));
    file = fopen(testfile[i], "r");
    if (file)
    {
      fc3d_newFromFile(testproblem, file);
    }
    fclose(file);


    /*    file = fopen("test.dat","w");     */
    /*    fc3d_printInFile(testproblem,file ); */
    /*    fclose(file); */



    info = test_Series_nsgs(testproblem, tol[i], localtol[i], nLocalSolver[i], localSolverList[i]);

    freeFrictionContactProblem(testproblem);

  }
  /*   for (int i =0 ; i< ntestfile; i++) free(testfile[i]); */

  free(testfile);
  free(tol);
  free(localtol);
  for (i = 0 ; i < ntestfile; i++) free(localSolverList[i]);
  free(localSolverList);
  free(nLocalSolver);



  printf("========= End Numerics tests for FrictionalContact_problem ========= \n");
  return info;
}

