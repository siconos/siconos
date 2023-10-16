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
#include <assert.h>                              // for assert
#include <stdio.h>                               // for printf, fclose, fopen
#include <stdlib.h>                              // for malloc, free, exit
#include "FrictionContactProblem.h"              // for FrictionContactProblem
#include "GlobalFrictionContactProblem.h"        // for GlobalFrictionContac...
#include "NumericsFwd.h"                         // for NumericsMatrix, Fric...
#include "NumericsMatrix.h"                      // for NumericsMatrix, NM_gemv
#include "SiconosBlas.h"                         // for cblas_dcopy, cblas_d...
#include "fc3d_Solvers.h"                        // for fc3d_DeSaxceFixedPoint
#include "fc3d_nonsmooth_Newton_AlartCurnier.h"  // for fc3d_nonsmooth_Newto...
#include "gfc3d_Solvers.h"                       // for gfc3d_DeSaxceFixedPo...
#include "numerics_verbose.h"                    // for verbose, numerics_pr...
#include "gfc3d_compute_error.h"
#include "SolverOptions.h"                       // for SICONOS_DPARAM_TOL

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"                                // for DEBUG_EXPR, DEBUG_P...
#include "gfc3d_ipm.h"                           // for primalResidual, dualResidual ...


#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void  gfc3d_nsgs_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{

  /* verbose=1; */
  DEBUG_BEGIN("gfc3d_nsgs_wr\n");
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ...\n");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););
    if(verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    // call nsgs solver for the local problem
    fc3d_nsgs(localproblem, reaction, velocity, info, options);

    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    /* Number of contacts */
    int nc = problem->numberOfContacts;
    /* Dimension of the problem */
    int m = 3 * nc;
    int n = problem->M->size0;
    double norm_q = cblas_dnrm2(n, problem->q, 1);
    double norm_b = cblas_dnrm2(m, problem->b, 1);
    double error;
    gfc3d_compute_error(problem,  reaction, velocity, globalVelocity,  options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error);



    // Printing in the same style as in IPM solver
    if (options->iparam[SICONOS_FRICTION_3D_NSGS_PRINTING_LIKE_IPM] == SICONOS_FRICTION_3D_NSGS_PRINTING_LIKE_IPM_TRUE)
    {
      // Get data
      unsigned int nd = problem->H->size1;
      unsigned int d = problem->dimension;
      unsigned int n = problem->numberOfContacts;
      unsigned int m = problem->M->size0;

      NumericsMatrix *P_mu = NM_create(NM_SPARSE, nd, nd);
      NumericsMatrix *P_mu_inv = NM_create(NM_SPARSE, nd, nd);
      NM_triplet_alloc(P_mu, nd);
      NM_triplet_alloc(P_mu_inv, nd);
      P_mu->matrix2->origin = NSM_TRIPLET;
      P_mu_inv->matrix2->origin = NSM_TRIPLET;
      for(unsigned int i = 0; i < nd; ++i)
        if(i % d == 0)
        {
          NM_entry(P_mu, i, i, 1.);
          NM_entry(P_mu_inv, i, i, 1.);
        }
        else
        {
          NM_entry(P_mu, i, i, problem->mu[(int)(i/d)]);
          NM_entry(P_mu_inv, i, i, 1.0/problem->mu[(int)(i/d)]);
        }

      NumericsMatrix *M = problem->M;
      NumericsMatrix *H_tilde = NM_transpose(problem->H);
      double *w_tilde = problem->b;
      double *w = (double*)calloc(nd, sizeof(double));
      double *f = problem->q;

      double *velocity_tmp = (double*)calloc(nd, sizeof(double));
      double *reaction_tmp = (double*)calloc(nd, sizeof(double));
      double *primalConstraint = (double*)calloc(nd, sizeof(double));
      double *dualConstraint = (double*)calloc(m, sizeof(double));

      double pinfeas, dinfeas, complem, udotr, nub;

      // Change of variable
      NumericsMatrix *H = NM_multiply(P_mu, H_tilde);
      NM_gemv(1.0, P_mu, w_tilde, 0.0, w);

      NM_gemv(1.0, P_mu, velocity, 0.0, velocity_tmp);
      cblas_dcopy(nd, velocity_tmp, 1, velocity, 1);

      NM_gemv(1.0, P_mu_inv, reaction, 0.0, reaction_tmp);
      cblas_dcopy(nd, reaction_tmp, 1, reaction, 1);

      // Compute residuals
      // ATTENTION: Current velocity = [u0 + |ub|; ub]
      complem = complemResidualNorm(velocity, reaction, nd, n);
      udotr = cblas_ddot(nd, velocity, 1, reaction, 1)/n;

      // Transform current velocity into [u0; ub]
      for (unsigned int i = 0; i<nd; i+=d)
      {
        nub = cblas_dnrm2(2, velocity+i+1, 1);
        velocity[i] -= nub;
      }
      primalResidual(velocity, H, globalVelocity, w, primalConstraint, &pinfeas, 1e-8);
      dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, 1e-8);

      printf("\n============ Printing in the same style as in IPM solver ============\n");
      numerics_printf_verbose(-1, "| it  | pinfeas | dinfeas | |u o r| |  u'r/n  | prj err |");
      numerics_printf_verbose(-1, "---------------------------------------------------------");
      numerics_printf_verbose(-1, "| %3i | %.1e | %.1e | %.1e | %.1e | %.1e |",
                            options->iparam[SICONOS_IPARAM_ITER_DONE], pinfeas, dinfeas, complem, udotr, options->dparam[SICONOS_DPARAM_RESIDU]);
      printf("=====================================================================\n\n");


      if(H_tilde) {H_tilde = NM_free(H_tilde); H_tilde = NULL;}
      if(H) {H = NM_free(H); H = NULL;}
      if (w) {free(w); w = NULL;}
      if (velocity_tmp) {free(velocity_tmp); velocity_tmp = NULL;}
      if (reaction_tmp) {free(reaction_tmp); reaction_tmp = NULL;}
      if (primalConstraint) {free(primalConstraint); primalConstraint = NULL;}
      if (dualConstraint) {free(dualConstraint); dualConstraint = NULL;}
    }



    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
  DEBUG_END("gfc3d_nsgs_wr\n");
}


void  gfc3d_admm_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  DEBUG_BEGIN("gfc3d_admm_wr\n");
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ...\n");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););
    
    if(verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_admm(localproblem, reaction, velocity, info, options);
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
  DEBUG_END("gfc3d_admm_wr\n");
}

void  gfc3d_nonsmooth_Newton_AlartCurnier_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  DEBUG_BEGIN("gfc3d_nonsmooth_Newton_AlartCurnier_wr(...)\n");
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ...\n");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););
    
    
    numerics_printf("gfc3d_nonsmooth_Newton_AlartCurnier_wr - Call to the fc3d solver ...\n");

    fc3d_nonsmooth_Newton_AlartCurnier(localproblem, reaction, velocity, info, options);

    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }

  DEBUG_END("gfc3d_nonsmooth_Newton_AlartCurnier_wr(...)\n")


}

void  gfc3d_nsgs_velocity_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ...\n");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););
    if(verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_nsgs_velocity(localproblem, reaction, velocity, info, options);

    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
}

void  gfc3d_proximal_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ...\n");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););
    DEBUG_EXPR(frictionContact_display(localproblem););
    if(verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_proximal(localproblem, reaction, velocity, info, options);

    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
}

void  gfc3d_DeSaxceFixedPoint_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ...\n");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););

    if(verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_DeSaxceFixedPoint(localproblem, reaction, velocity, info, options);
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
}

void  gfc3d_TrescaFixedPoint_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ...\n");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););

    if(verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_TrescaFixedPoint(localproblem, reaction, velocity, info, options);
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }

}


void gfc3d_ipm_snm_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  verbose = 1;
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ...\n");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););

    if(verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_IPM_SNM(localproblem, reaction, velocity, info, options);
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
}


