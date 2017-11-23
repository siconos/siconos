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

#include "SiconosBlas.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "LinearComplementarityProblem.h"
#include "LinearComplementarityProblem_as_ConvexQP.h"
#include "ConvexQP_Solvers.h"
#include "SiconosCompat.h"
#include "SiconosFortran.h"
#include "LCP_Solvers.h"
#include "lcp_cst.h"

#include "SolverOptions.h"
#include "numerics_verbose.h"





void lcp_ConvexQP_ProjectedGradient(LinearComplementarityProblem* problem, double *z, double *w, int* info, SolverOptions* options)
{
  /* verbose=1; */
  /* Dimension of the problem */
  int n = problem->size;

  ConvexQP *cqp = (ConvexQP *)malloc(sizeof(ConvexQP));

  cqp->M = problem->M;
  cqp->q = problem->q;

  cqp->ProjectionOnC = &Projection_ConvexQP_LCP;

  int iter=0;
  double error=1e24;

  LinearComplementarityProblem_as_ConvexQP *lcp_as_cqp= (LinearComplementarityProblem_as_ConvexQP*)malloc(sizeof(LinearComplementarityProblem_as_ConvexQP));
  cqp->env = lcp_as_cqp ;
  cqp->size = n;

  /*set the norm of the ConvexQP to the norm of problem->q  */
  double norm_q = cblas_dnrm2(n , problem->q , 1);
  cqp->normConvexQP= norm_q;
  cqp->istheNormConvexQPset=1;

  lcp_as_cqp->cqp = cqp;
  lcp_as_cqp->lcp = problem;
  lcp_as_cqp->options = options;

  SolverOptions * cqpsolver_options = (SolverOptions *) malloc(sizeof(SolverOptions));
  convexQP_ProjectedGradient_setDefaultSolverOptions(cqpsolver_options);

  int isize = options->iSize;
  int dsize = options->dSize;
  int cqp_isize = cqpsolver_options->iSize;
  int cqp_dsize = cqpsolver_options->dSize;
  if (isize != cqp_isize )
  {
    printf("Warning: options->iSize in lcp_ConvexQP_ProjectedGradient is not consitent with options->iSize in ConvexQP_PG\n");
  }
  if (dsize != cqp_dsize )
  {
    printf("Warning: options->iSize in lcp_ConvexQP_ProjectedGradient is not consitent with options->iSize in ConvexQP_PG\n");
  }
  int i;
  for (i = 0; i < min(isize,cqp_isize); i++)
  {
    if (options->iparam[i] != 0 )
      cqpsolver_options->iparam[i] = options->iparam[i] ;
  }
  for (i = 0; i <  min(dsize,cqp_dsize); i++)
  {
    if (fabs(options->dparam[i]) >= 1e-24 )
      cqpsolver_options->dparam[i] = options->dparam[i] ;
  }

  convexQP_ProjectedGradient(cqp, z, w , info , cqpsolver_options);

  /* **** Criterium convergence **** */

  lcp_compute_error(problem, z , w, options->dparam[0], &error);

  error = cqpsolver_options->dparam[1];
  iter = cqpsolver_options->iparam[7];

  options->dparam[SICONOS_DPARAM_RESIDU] = error;
  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;


  if (verbose > 0)
  {
    printf("--------------- LCP - ConvexQP PG  - #Iteration %i Final Residual = %14.7e\n", iter, error);
  }
  free(cqp);

  solver_options_delete(cqpsolver_options);
  free(cqpsolver_options);
  free(lcp_as_cqp);



}


int linearComplementarity_ConvexQP_ProjectedGradient_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the CONVEX QP Solver\n");
  }
  convexQP_ProjectedGradient_setDefaultSolverOptions(options);

  options->solverId = SICONOS_LCP_CONVEXQP_PG;
  
  return 0;
}


