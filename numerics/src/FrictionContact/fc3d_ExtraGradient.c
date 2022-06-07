/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include <math.h>                    // for sqrt, pow
#include <stdio.h>                   // for printf, NULL
#include <stdlib.h>                  // for calloc, free
#include "FrictionContactProblem.h"  // for FrictionContactProblem
#include "Friction_cst.h"            // for SICONOS_FRICTION_3D_EG
#include "NumericsFwd.h"             // for SolverOptions, FrictionContactPr...
#include "NumericsMatrix.h"          // for NM_gemv
#include "SolverOptions.h"           // for SolverOptions, solver_options_nu...
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"                   // for DEBUG_PRINTF, DEBUG_EXPR_WE
#include "fc3d_Solvers.h"            // for fc3d_ExtraGradient, fc3d_ExtraGr...
#include "fc3d_compute_error.h"      // for fc3d_compute_error
#include "numerics_verbose.h"        // for verbose
#include "projectionOnCone.h"        // for projectionOnCone
#include "SiconosBlas.h"                   // for cblas_dcopy, cblas_dnrm2, cblas_...

void fc3d_ExtraGradient(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
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
  double norm_q = cblas_dnrm2(nc*3, problem->q, 1);




  /*****  Fixed point iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  int contact; /* Number of the current row of blocks in M */
  int nLocal = 3;

  double * velocitytmp = (double *)calloc(n,sizeof(double));
  double * reactiontmp = (double *)calloc(n,sizeof(double));

  double rho = 0.0, rho_k =0.0;
  int isVariable = 0;

  if(dparam[SICONOS_FRICTION_3D_NSN_RHO] > 0.0)
  {
    rho = dparam[SICONOS_FRICTION_3D_NSN_RHO];
    if(verbose > 0)
    {
      printf("--------------- FC3D - Extra Gradient (EG) - Fixed stepsize with  rho = %14.7e \n", rho);
    }
  }
  else
  {
    /* Variable step in iterations*/
    isVariable = 1;
    rho = -dparam[SICONOS_FRICTION_3D_NSN_RHO];
    if(verbose > 0)
    {
      printf("--------------- FC3D - Extra Gradient (EG) - Variable stepsize with starting rho = %14.7e \n", rho);
    }

  }

  double alpha = 1.0;
  double beta = 1.0;

  /* Variable for Line_search */
  int success =0;
  double error_k;
  int ls_iter = 0;
  int ls_itermax = 10;
  double tau=0.6, L= 0.9, Lmin =0.3, taumin=0.7;
  double a1=0.0, a2=0.0;
  double * reaction_k =0;
  double * velocity_k =0;

  if(isVariable)
  {
    reaction_k = (double *)calloc(n,sizeof(double));
    velocity_k = (double *)calloc(n,sizeof(double));
  }

  if(!isVariable)
  {
    /*   double minusrho  = -1.0*rho; */
    while((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;

      /* reactiontmp <- reaction  */
      cblas_dcopy(n, reaction, 1, reactiontmp, 1);

      /* velocitytmp <- q  */
      cblas_dcopy(n, q, 1, velocitytmp, 1);


      NM_gemv(alpha, M, reactiontmp, beta, velocitytmp);
      // projection for each contact
      for(contact = 0 ; contact < nc ; ++contact)
      {
        int pos = contact * nLocal;
        double  normUT = sqrt(velocitytmp[pos + 1] * velocitytmp[pos + 1] + velocitytmp[pos + 2] * velocitytmp[pos + 2]);
        reactiontmp[pos] -= rho * (velocitytmp[pos] + mu[contact] * normUT);
        reactiontmp[pos + 1] -= rho * velocitytmp[pos + 1];
        reactiontmp[pos + 2] -= rho * velocitytmp[pos + 2];
        projectionOnCone(&reactiontmp[pos], mu[contact]);
      }
      cblas_dcopy(n, q, 1, velocitytmp, 1);
      NM_gemv(alpha, M, reactiontmp, beta, velocitytmp);
      // projection for each contact
      for(contact = 0 ; contact < nc ; ++contact)
      {
        int pos = contact * nLocal;
        double  normUT = sqrt(velocitytmp[pos + 1] * velocitytmp[pos + 1] + velocitytmp[pos + 2] * velocitytmp[pos + 2]);
        reaction[pos] -= rho * (velocitytmp[pos] + mu[contact] * normUT);
        reaction[pos + 1] -= rho * velocitytmp[pos + 1];
        reaction[pos + 2] -= rho * velocitytmp[pos + 2];
        projectionOnCone(&reaction[pos], mu[contact]);
      }

      /* **** Criterium convergence **** */
      fc3d_compute_error(problem, reaction, velocity, tolerance, options, norm_q, &error);

      if(options->callback)
      {
        options->callback->collectStatsIteration(options->callback->env, nc * 3,
            reaction, velocity,
            error, NULL);
      }

      if(verbose > 0)
      {
        printf("--------------- FC3D - Extra Gradient (EG) - Iteration %i rho = %14.7e \tError = %14.7e\n", iter, rho, error);
      }
      if(error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }

  if(isVariable)
  {
    while((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;


      /* Store the error */
      error_k = error;
      /* store the reaction at the beginning of the iteration */
      cblas_dcopy(n, reaction, 1, reaction_k, 1);

      /* velocity_k <- q  */
      cblas_dcopy(n, q, 1, velocity_k, 1);

      /* velocity_k <- q + M * reaction_k  */
      beta = 1.0;
      NM_gemv(alpha, M, reaction_k, beta, velocity_k);

      ls_iter = 0 ;
      success =0;

      while(!success && (ls_iter < ls_itermax))
      {

        rho_k = rho * pow(tau,ls_iter);

        /* projection for each contact */
        for(contact = 0 ; contact < nc ; ++contact)
        {
          int pos = contact * nLocal;
          double  normUT = sqrt(velocity_k[pos + 1] * velocity_k[pos + 1] + velocity_k[pos + 2] * velocity_k[pos + 2]);
          /* reaction[pos] = reaction_k[pos] -  rho_k * (velocity_k[pos] + mu[contact] * normUT); */
          /* reaction[pos + 1] = reaction_k[pos+1] - rho_k * velocity_k[pos + 1]; */
          /* reaction[pos + 2] = reaction_k[pos+2] - rho_k * velocity_k[pos + 2]; */

          reaction[pos] -= rho_k * (velocity_k[pos] + mu[contact] * normUT);
          reaction[pos + 1] -= rho_k * velocity_k[pos + 1];
          reaction[pos + 2] -= rho_k * velocity_k[pos + 2];
          projectionOnCone(&reaction[pos], mu[contact]);
        }


        /* velocity <- q + M * reaction  */
        beta = 1.0;
        cblas_dcopy(n, q, 1, velocity, 1);
        NM_gemv(alpha, M, reaction, beta, velocity);


        /* velocitytmp <- velocity */

        DEBUG_EXPR_WE(for(int i =0; i< 5 ; i++)
      {
        printf("reaction[%i]=%12.8e\t",i,reaction[i]);
          printf("velocity[%i]=F[%i]=%12.8e\n",i,i,velocity[i]);
        }
                     );
        cblas_dcopy(n, velocity, 1, velocitytmp, 1) ;
        /* velocitytmp <- modified velocity - velocity_k   */
        for(contact = 0 ; contact < nc ; ++contact)
        {
          int pos = contact * nLocal;
          double  normUT = sqrt(velocitytmp[pos + 1] * velocitytmp[pos + 1]
                                + velocitytmp[pos + 2] * velocitytmp[pos + 2]);
          double  normUT_k = sqrt(velocity_k[pos + 1] * velocity_k[pos + 1] + velocity_k[pos + 2] * velocity_k[pos + 2]);
          velocitytmp[pos] += mu[contact] * (normUT -normUT_k)  ;
        }

        /* for (contact = 0 ; contact < nc ; ++contact) */
        /* { */
        /*   int pos = contact * nLocal; */
        /*   double  normUT = sqrt(velocity_k[pos + 1] * velocity_k[pos + 1] + velocity_k[pos + 2] * velocity_k[pos + 2]); */
        /*   velocity_k[pos] += mu[contact] * normUT; */
        /* } */

        cblas_daxpy(n, -1.0, velocity_k, 1, velocitytmp, 1) ;

        a1 = cblas_dnrm2(n, velocitytmp, 1);
        DEBUG_PRINTF("a1 = %12.8e\n", a1);

        /* reactiontmp <- reaction */
        cblas_dcopy(n, reaction, 1, reactiontmp, 1) ;

        /* reactiontmp <- reaction - reaction_k   */
        cblas_daxpy(n, -1.0, reaction_k, 1, reactiontmp, 1) ;

        a2 = cblas_dnrm2(n, reactiontmp, 1) ;
        DEBUG_PRINTF("a2 = %12.8e\n", a2);

        success = (rho_k*a1 < L * a2)?1:0;

        /* printf("rho_k = %12.8e\t", rho_k); */
        /* printf("a1 = %12.8e\t", a1); */
        /* printf("a2 = %12.8e\t", a2); */
        /* printf("norm reaction = %12.8e\t",cblas_dnrm2(n, reaction, 1) ); */
        /* printf("success = %i\n", success); */

        ls_iter++;
      }
      /* velocitytmp <- M* reaction* q  */
      cblas_dcopy(n, q, 1, velocitytmp, 1);
      NM_gemv(alpha, M, reaction, beta, velocitytmp);

      // projection for each contact
      for(contact = 0 ; contact < nc ; ++contact)
      {
        int pos = contact * nLocal;
        double  normUT = sqrt(velocitytmp[pos + 1] * velocitytmp[pos + 1] +
                              velocitytmp[pos + 2] * velocitytmp[pos + 2]);
        reaction[pos] -= rho_k * (velocitytmp[pos] + mu[contact] * normUT);
        reaction[pos + 1] -= rho_k * velocitytmp[pos + 1];
        reaction[pos + 2] -= rho_k * velocitytmp[pos + 2];
        /* reaction[pos] = reaction_k[pos] -  rho_k * (velocitytmp[pos] + mu[contact] * normUT); */
        /* reaction[pos + 1] = reaction_k[pos+1] - rho_k * velocitytmp[pos + 1]; */
        /* reaction[pos + 2] = reaction_k[pos+2] - rho_k * velocitytmp[pos + 2]; */
        projectionOnCone(&reaction[pos], mu[contact]);
      }
      DEBUG_EXPR_WE(for(int i =0; i< 5 ; i++)
    {
      printf("reaction[%i]=%12.8e\t",i,reaction[i]);
        printf("velocity[%i]=F[%i]=%12.8e\n",i,i,velocity[i]);
      }
                   );

      /* **** Criterium convergence **** */
      fc3d_compute_error(problem, reaction, velocity, tolerance, options, norm_q, &error);
      DEBUG_PRINTF("error = %12.8e\t error_k = %12.8e\n",error,error_k);
      /*Update rho*/
      if((rho_k*a1 < Lmin * a2) && (error < error_k))
      {
        rho =rho_k/taumin;
      }
      else
        rho =rho_k;





      if(verbose > 0)
      {
        printf("--------------- FC3D - Extra Gradient (EG) - Iteration %i rho = %14.7e \tError = %14.7e\n", iter, rho, error);
      }
      if(error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }


  if(verbose > 0)
  {
    printf("--------------- FC3D -  Extra Gradient (EG) - #Iteration %i Final Residual = %14.7e\n", iter, error);
  }
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  free(velocitytmp);
  free(reactiontmp);
  if(isVariable)
  {
    free(reaction_k);
    free(velocity_k);
  }
}


void fc3d_eg_set_default(SolverOptions* options)
{
  options->dparam[SICONOS_FRICTION_3D_NSN_RHO] = -1.0; // rho is variable by default
}
