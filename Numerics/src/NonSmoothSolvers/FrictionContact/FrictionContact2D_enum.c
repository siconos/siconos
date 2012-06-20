/* Siconos-Numerics, Copyright INRIA 2005-2011.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "LA.h"
#include "FrictionContact2D_Solvers.h"
#include "FrictionContact2D_compute_error.h"
#include "LCP_Solvers.h"
#include <assert.h>

void frictionContact2D_enum(FrictionContactProblem* problem, double *reaction, double *velocity, int *info, SolverOptions* options,  NumericsOptions* global_options)
{
  int i;
  // conversion into LCP
  LinearComplementarityProblem* lcp_problem = (LinearComplementarityProblem*)malloc(sizeof(LinearComplementarityProblem));

  FrictionContact2D_tolcp(problem, lcp_problem);
  /* frictionContact_display(problem); */
  /* LinearComplementarity_display(lcp_problem); */



  double *zlcp = (double*)malloc(lcp_problem->size * sizeof(double));
  double *wlcp = (double*)malloc(lcp_problem->size * sizeof(double));

  for (i = 0; i < lcp_problem->size; i++)
  {
    zlcp[i] = 0.0;
    wlcp[i] = 0.0;
  }


  /*  FILE * fcheck = fopen("lcp_relay.dat","w"); */
  /*  info = linearComplementarity_printInFile(lcp_problem,fcheck); */

  // Call the lcp_solver

  SolverOptions * lcp_options = options->internalSolvers;

  lcp_enum_init(lcp_problem, lcp_options, 1);
  //}
  * info = linearComplementarity_driver(lcp_problem, zlcp , wlcp, lcp_options, global_options);
  if (options->filterOn > 0)
    lcp_compute_error(lcp_problem, zlcp, wlcp, lcp_options->dparam[0], &(lcp_options->dparam[1]));
  lcp_enum_reset(lcp_problem, lcp_options, 1);

  /*       printf("\n"); */
  int nc = problem->numberOfContacts;
  // Conversion of result
  for (i = 0; i < nc; i++)
  {

    /* printf("Contact number = %i\n",i); */
    reaction[2 * i] = zlcp[3 * i];
    reaction[2 * i + 1] = 1.0 / 2.0 * (zlcp[3 * i + 1] - wlcp[3 * i + 2]);
    /* printf("reaction[ %i]=%12.10e\n", 2*i, reaction[2*i]); */
    /* printf("reaction[ %i]=%12.10e\n", 2*i+1, reaction[2*i+1]); */

    velocity[2 * i] = wlcp[3 * i];
    velocity[2 * i + 1] = wlcp[3 * i + 1] - zlcp[3 * i + 2];
    /* printf("velocity[ %i]=%12.10e\n", 2*i, velocity[2*i]);   */
    /* printf("velocity[ %i]=%12.10e\n", 2*i+1, velocity[2*i+1]); */
  }

  /*        for (i=0; i< lcp_problem->size; i++){  */
  /*     printf("zlcp[ %i]=%12.10e,\t wlcp[ %i]=%12.10e \n", i, zlcp[i],i, wlcp[i]); */
  /*        } */
  /*        printf("\n"); */

  /*        for (i=0; i< problem->size; i++){  */
  /*     printf("z[ %i]=%12.10e,\t w[ %i]=%12.10e\n", i, z[i],i, w[i]); */
  /*        } */


  /*        printf("\n"); */
  double error;
  *info = FrictionContact2D_compute_error(problem, reaction , velocity, options->dparam[0], &error);
  free(zlcp);
  free(wlcp);
  freeLinearComplementarityProblem(lcp_problem);
}


int frictionContact2D_enum_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the Lemke Solver for FrictionContact2D\n");
  }
  /*  strcpy(options->solverName,"Lemke");*/
  options->solverId = SICONOS_FRICTION_2D_ENUM;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->dparam[0] = 1e-6;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  LinearComplementarityProblem* problem = NULL;
  linearComplementarity_enum_setDefaultSolverOptions(problem, options->internalSolvers);

  return 0;
}


