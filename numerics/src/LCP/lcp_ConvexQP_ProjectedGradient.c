/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include <math.h>                                      // for fabs
#include <stdio.h>                                     // for printf
#include <stdlib.h>                                    // for free, malloc
#include "SiconosBlas.h"     // for cblas_dnrm2
#include "ConvexQP.h"                                  // for ConvexQP
#include "ConvexQP_Solvers.h"                          // for convexQP_Proje...
#include "LCP_Solvers.h"                               // for lcp_compute_error
#include "LinearComplementarityProblem.h"              // for LinearCompleme...
#include "LinearComplementarityProblem_as_ConvexQP.h"  // for LinearCompleme...
#include "NumericsFwd.h"                               // for SolverOptions
#include "NSSTools.h"   // for min
#include "SolverOptions.h"                             // for SolverOptions
#include "lcp_cst.h"                                   // for SICONOS_LCP_CO...
#include "numerics_verbose.h"                          // for verbose


void lcp_ConvexQP_ProjectedGradient(LinearComplementarityProblem* problem, double *z, double *w, int* info, SolverOptions* options)
{
  /* verbose=1; */
  /* Dimension of the problem */
  int n = problem->size;

  ConvexQP *cqp = (ConvexQP *)malloc(sizeof(ConvexQP));

  cqp->M = problem->M;
  cqp->q = problem->q;

  cqp->ProjectionOnC = &Projection_ConvexQP_LCP;

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
  convexQP_ProjectedGradient(cqp, z, w , info , options);

  /* **** Criterium convergence **** */

  double error;
  lcp_compute_error(problem, z , w, options->dparam[SICONOS_DPARAM_TOL], &error);
  if (verbose > 0)
  {
    printf("--------------- LCP - ConvexQP PG  - #Iteration %i Final Residual = %14.7e\n",
           options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
  }
  free(cqp);
  free(lcp_as_cqp);


}
