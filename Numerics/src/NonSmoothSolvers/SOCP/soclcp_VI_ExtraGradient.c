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

#include "SecondOrderConeLinearComplementarityProblem_as_VI.h"
#include "VariationalInequality_Solvers.h"
#include "SOCLCP_Solvers.h"
#include "soclcp_compute_error.h"

#include "SolverOptions.h"




void soclcp_VI_ExtraGradient(SecondOrderConeLinearComplementarityProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{

  /* Dimension of the problem */
  int n = problem->n;


  VariationalInequality *vi = (VariationalInequality *)malloc(sizeof(VariationalInequality));

  //vi.self = &vi;
  vi->F = &Function_VI_SOCLCP;
  vi->ProjectionOnX = &Projection_VI_SOCLCP;

  int iter=0;
  double error=1e24;

  SecondOrderConeLinearComplementarityProblem_as_VI *soclcp_as_vi= (SecondOrderConeLinearComplementarityProblem_as_VI*)malloc(sizeof(SecondOrderConeLinearComplementarityProblem_as_VI));
  vi->env =soclcp_as_vi ;
  vi->size =  n;


  /*Set the norm of the VI to the norm of problem->q  */
  vi->normVI= cblas_dnrm2(n , problem->q , 1);
  vi->istheNormVIset=1;

  soclcp_as_vi->vi = vi;
  soclcp_as_vi->soclcp = problem;
  /* frictionContact_display(fc3d_as_vi->fc3d); */

  SolverOptions * visolver_options = (SolverOptions *) malloc(sizeof(SolverOptions));
  variationalInequality_setDefaultSolverOptions(visolver_options,
                                                SICONOS_VI_EG);

  int isize = options->iSize;
  int dsize = options->dSize;
  int vi_isize = visolver_options->iSize;
  int vi_dsize = visolver_options->dSize;

  if (isize != vi_isize )
  {
    printf("size prolem in frictionContact3D_VI_ExtraGradient\n");
  }
  if (dsize != vi_dsize )
  {
    printf("size prolem in frictionContact3D_VI_ExtraGradient\n");
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

  variationalInequality_ExtraGradient(vi, reaction, velocity , info , visolver_options);



  /* **** Criterium convergence **** */
  soclcp_compute_error(problem, reaction , velocity, options->dparam[0], options, &error);

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
    printf("----------------------------------- SOCLCP - VI Extra Gradient (VI_EG) - #Iteration %i Final Error = %14.7e\n", iter, error);
  }
  free(vi);

  deleteSolverOptions(visolver_options);
  free(visolver_options);
  visolver_options=NULL;
  free(soclcp_as_vi);



}


int soclcp_VI_ExtraGradient_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the ExtraGradient Solver\n");
  }

  /*strcpy(options->solverName,"DSFP");*/
  options->solverId = SICONOS_SOCLCP_VI_EG;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 8;
  options->dSize = 8;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  null_SolverOptions(options);
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
