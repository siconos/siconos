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
#include <math.h>                    // for sqrt
#include <stdio.h>                   // for printf, NULL
#include <stdlib.h>                  // for calloc, free
#include "FrictionContactProblem.h"  // for FrictionContactProblem
#include "Friction_cst.h"            // for SICONOS_FRICTION_3D_DSFP
#include "NumericsFwd.h"             // for SolverOptions, FrictionContactPr...
#include "NumericsMatrix.h"          // for NM_gemv
#include "SolverOptions.h"           // for SolverOptions, solver_options_nu...
#include "fc3d_Solvers.h"            // for fc3d_DeSaxceFixedPoint, fc3d_DeS...
#include "fc3d_compute_error.h"      // for fc3d_compute_error
#include "numerics_verbose.h"        // for verbose, numerics_error
#include "projectionOnCone.h"        // for projectionOnCone
#include "SiconosBlas.h"                   // for cblas_dcopy, cblas_dnrm2

void fc3d_DeSaxceFixedPoint(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  double* q = problem->q;
  NumericsMatrix* M = problem->M;
  double* mu = problem->mu;
  /* Dimension of the problem */
  int n = 3 * nc;
  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];
  double norm_q = cblas_dnrm2(nc*3 , problem->q , 1);

  /*****  Fixed point iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  int contact; /* Number of the current row of blocks in M */
  int nLocal = 3;
  double * velocitytmp = (double *)calloc(n, sizeof(double));

  double rho = 0.0;

  if (dparam[SICONOS_FRICTION_3D_NSN_RHO] > 0.0)
  {
    rho = dparam[SICONOS_FRICTION_3D_NSN_RHO];
    if (verbose > 0)
    {
      printf("--------------- FC3D - DeSaxce Fixed Point (DSFP) - Fixed stepsize with  rho = %14.7e \n", rho);
    }

  }
  else
  {
    numerics_error("fc3d_DeSaxceFixedPoint", "The De Saxce fixed point is implemented with a fixed time--step. Use FixedPointProjection (VI_FPP) method for a variable time--step");
  }

  double alpha = 1.0;
  double beta = 1.0;

    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      /* velocitytmp <- q  */
      cblas_dcopy(n , q , 1 , velocitytmp, 1);

      /* velocitytmp <- q + M * reaction  */
      beta = 1.0;
      NM_gemv(alpha, M, reaction, beta, velocitytmp);

      /* projection for each contact */
      for (contact = 0 ; contact < nc ; ++contact)
      {
        int pos = contact * nLocal;
        double  normUT = sqrt(velocitytmp[pos + 1] * velocitytmp[pos + 1] + velocitytmp[pos + 2] * velocitytmp[pos + 2]);
        reaction[pos] -= rho * (velocitytmp[pos] + mu[contact] * normUT);
        reaction[pos + 1] -= rho * velocitytmp[pos + 1];
        reaction[pos + 2] -= rho * velocitytmp[pos + 2];
        projectionOnCone(&reaction[pos], mu[contact]);
      }

      /* **** Criterium convergence **** */
      fc3d_compute_error(problem, reaction , velocity, tolerance, options, norm_q, &error);

      if (options->callback)
      {
        options->callback->collectStatsIteration(options->callback->env,
                                        nc * 3, reaction, velocity,
                                        error, NULL);
    }

      if (verbose > 0)
        printf("--------------- FC3D - DeSaxce Fixed Point (DSFP) - Iteration %i rho = %14.7e \tError = %14.7e\n", iter, rho, error);

      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }



  if (verbose > 0)
    printf("--------------- FC3D - DeSaxce Fixed point (DSFP) - #Iteration %i Final Residual = %14.7e\n", iter, error);
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  dparam[SICONOS_DPARAM_RESIDU] = error;
  free(velocitytmp);

}


void fc3d_dsfp_set_options(SolverOptions* options)
{
  options->dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1.0;
}
