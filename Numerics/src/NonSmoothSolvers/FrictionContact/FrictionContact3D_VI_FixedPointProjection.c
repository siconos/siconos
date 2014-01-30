/* Siconos-Numerics, Copyright INRIA 2005-2012.
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

#include "SiconosBlas.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FrictionContactProblem_as_VI.h"
#include "VariationalInequality_Solvers.h"
#include "FrictionContact3D_Solvers.h"
#include "FrictionContact3D_compute_error.h"

#include "SolverOptions.h"




void frictionContact3D_VI_FixedPointProjection(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  /* Dimension of the problem */
  int n = 3 * nc;


  VariationalInequality *vi = (VariationalInequality *)malloc(sizeof(VariationalInequality));

  //vi.self = &vi;
  vi->F = &Function_VI_FC3D;
  vi->ProjectionOnX = &Projection_VI_FC3D;

  int iter=0;
  double error=1e24;

  FrictionContactProblem_as_VI *fc3d_as_vi= (FrictionContactProblem_as_VI*)malloc(sizeof(FrictionContactProblem_as_VI));
  vi->env =fc3d_as_vi ;
  vi->size =  n;
  fc3d_as_vi->vi = vi;
  fc3d_as_vi->fc3d = problem;
  /* frictionContact_display(fc3d_as_vi->fc3d); */

  SolverOptions * visolver_options = (SolverOptions *) malloc(sizeof(SolverOptions));
  variationalInequality_setDefaultSolverOptions(visolver_options,
                                                SICONOS_VI_FPP);

  int isize = options->iSize;
  int dsize = options->dSize;
  int vi_isize = visolver_options->iSize;
  int vi_dsize = visolver_options->dSize;
  if (isize != vi_isize )
  {
    printf("size prolem in frictionContact3D_VI_FixedPointProjection\n");
  }
  if (dsize != vi_dsize )
  {
    printf("size prolem in frictionContact3D_VI_FixedPointProjection\n");
  }
  int i;
  for (i = 0; i < isize; i++)
  {
    visolver_options->iparam[i] = options->iparam[i] ;
  }
  for (i = 0; i < dsize; i++)
  {
    visolver_options->dparam[i] = options->dparam[i] ;
  }

  variationalInequality_FixedPointProjection(vi, reaction, velocity , info , visolver_options);



  /* **** Criterium convergence **** */
  FrictionContact3D_compute_error(problem, reaction , velocity, options->dparam[0], options, &error);

  /* for (i =0; i< n ; i++) */
  /* { */
  /*   printf("reaction[%i]=%f\t",i,reaction[i]);    printf("velocity[%i]=F[%i]=%f\n",i,i,velocity[i]); */
  /* } */

  error = visolver_options->dparam[1];
  iter = visolver_options->iparam[7];

  options->dparam[1] = error;
  options->iparam[7] = iter;


  if (verbose > 0)
  {
    printf("----------------------------------- FC3D - VI Fixed Point Projection (VI_FPP) - #Iteration %i Final Error = %14.7e\n", iter, error);
  }
  free(vi);

  deleteSolverOptions(visolver_options);
  free(visolver_options);
  free(fc3d_as_vi);



}


int frictionContact3D_VI_FixedPointProjection_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the FixedPointProjection Solver\n");
  }

  /*strcpy(options->solverName,"DSFP");*/
  options->solverId = SICONOS_FRICTION_3D_VI_FPP;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 8;
  options->dSize = 8;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;   options->callback = NULL; options->numericsOptions = NULL;
  for (i = 0; i < 8; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 20000;
  options->dparam[0] = 1e-3;
  options->dparam[3] = 1e-3;
  options->dparam[3] = -1.0; // rho is variable by default

  options->internalSolvers = NULL;

  return 0;
}
