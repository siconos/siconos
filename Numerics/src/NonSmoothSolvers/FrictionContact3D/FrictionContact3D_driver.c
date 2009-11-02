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
#include <time.h>
#include <float.h>
#include "LA.h"
#include "Numerics_Options.h"
#include "FrictionContact3D_Solvers.h"
#include "NonSmoothDrivers.h"

int frictionContact3D_driver(FrictionContact_Problem* problem, double *reaction , double *velocity, Solver_Options* options, Numerics_Options* global_options)
{

  /* Set global options */
  setNumericsOptions(global_options);

  /* If the options for solver have not been set, read default values in .opt file */
  int NoDefaultOptions = options->isSet; /* true(1) if the Solver_Options structure has been filled in else false(0) */

  if (!NoDefaultOptions)
    readSolverOptions(3, options);

  if (verbose > 0)
    printSolverOptions(options);

  /* Solver name */
  char * name = options->solverName;


  int info = -1 ;

  /* Non Smooth Gauss Seidel (NSGS) */
  if (strcmp(name, "NSGS") == 0)
  {
    if (verbose == 1)
      printf(" ========================== Call NSGS solver for Friction-Contact 3D problem ==========================\n");
    frictionContact3D_nsgs(problem, reaction , velocity , &info , options);
  }
  /* Proximal point algorithm */
  else if (strcmp(name, "PROX") == 0)
  {
    if (verbose == 1)
      printf(" ========================== Call PROX solver for Friction-Contact 3D problem ==========================\n");
    frictionContact3D_proximal(problem, reaction , velocity , &info , options);
  }
  /* Projected Gradient algorithm */
  else if (strcmp(name, "PG") == 0)
  {
    if (verbose == 1)
      printf(" ========================== Call Projected Gradinet (PG) solver for Friction-Contact 3D problem ==========================\n");
    frictionContact3D_projectedgradient(problem, reaction , velocity , &info , options);
  }
  else
  {
    fprintf(stderr, "Numerics, FrictionContact3D_driver failed. Unknown solver.\n");
    exit(EXIT_FAILURE);

  }

  return info;

}

int checkTrivialCase(int n, double* q, double* velocity, double* reaction, int* iparam, double* dparam)
{
  /* norm of vector q */
  double qs = DNRM2(n , q , 1);
  int i;
  int info = -1;
  if (qs <= DBL_EPSILON)
  {
    // q norm equal to zero (less than DBL_EPSILON)
    // -> trivial solution: reaction = 0 and velocity = q
    for (i = 0 ; i < n ; ++i)
    {
      velocity[i] = q[i];
      reaction[i] = 0.;
    }
    iparam[2] = 0;
    iparam[4] = 0;
    dparam[1] = 0.0;
    dparam[3] = 0.0;
    info = 0;
    if (verbose == 1)
      printf("FrictionContact3D driver, norm(q) = 0, trivial solution reaction = 0, velocity = q.\n");
  }
  return info;
}

int frictionContact3D_printInFile(FrictionContact_Problem*  problem, FILE* file)
{
  if (! problem)
  {
    fprintf(stderr, "Numerics, FrictionContact_Problem printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int nc = problem->numberOfContacts;
  fprintf(file, "%d\n", nc);
  printInFile(problem->M, file);
  for (int i = 0; i < problem->M->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->q[i]);
  }
  fprintf(file, "\n");
  for (int i = 0; i < nc; i++)
  {
    fprintf(file, "%32.24e ", problem->mu[i]);
  }
  fprintf(file, "\n");
  fprintf(file, "%d\n", problem->isComplete);
  return 0;
}

int frictionContact3D_newFromFile(FrictionContact_Problem* problem, FILE* file)
{
  int nc = 0;
  fscanf(file, "%d\n", &nc);
  problem->numberOfContacts = nc;
  problem->M = (NumericsMatrix *)malloc(sizeof(NumericsMatrix));

  readInFile(problem->M, file);

  problem->q = (double *) malloc(problem->M->size1 * sizeof(double));
  for (int i = 0; i < problem->M->size1; i++)
  {
    fscanf(file, "%lf ", &(problem->q[i]));
  }

  fscanf(file, "\n");
  problem->mu = (double *) malloc(nc * sizeof(double));
  for (int i = 0; i < nc; i++)
  {
    fscanf(file, "%lf ", &(problem->mu[i]));
  }
  fscanf(file, "\n");
  fscanf(file, "%d\n", &(problem->isComplete));
  return 0;
}
void freeFrictionContact_problem(FrictionContact_Problem* problem)
{

  freeNumericsMatrix(problem->M);
  free(problem->M);
  free(problem->mu);
  free(problem->q);
  free(problem);
  problem = NULL;

}
