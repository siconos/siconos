/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include <assert.h>                        // for assert
#include <float.h>                         // for DBL_EPSILON
#include <math.h>                          // for fabs, sqrt, INFINITY
#include <stdio.h>                         // for NULL, fprintf, printf, stderr
#include <stdlib.h>                        // for free, malloc, calloc
#include "FrictionContactProblem.h"        // for FrictionContactProblem
#include "Friction_cst.h"                  // for SICONOS_FRICTION_3D_IPARAM...
#include "LCP_Solvers.h"                   // for lcp_nsgs_SBM_buildLocalPro...
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NumericsFwd.h"                   // for SolverOptions, LinearCompl...
#include "NumericsMatrix.h"                // for NumericsMatrix, RawNumeric...
#include "SolverOptions.h"                 // for SolverOptions, SICONOS_DPA...
#include "SparseBlockMatrix.h"             // for SparseBlockStructuredMatrix
#include "fc2d_Solvers.h"                  // for fc2d_nsgs_sbm, fc2d_spa...
#include "fc2d_compute_error.h"            // for fc2d_compute_error
#include "numerics_verbose.h"              // for numerics_printf, verbose
#include "SiconosBlas.h"                         // for cblas_dnrm2

/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES 1 */
#include "siconos_debug.h"                         // for DEBUG_BEGIN, DEBUG_END
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif


#define SGN(x) ((x) < 0 ? -1 : (x) > 0 ? 1 : 0)

static void fc2d_nsgs_buildLocalProblem(int contact,
                                        FrictionContactProblem* problem,
                                        LinearComplementarityProblem* local_problem,
                                        double* reaction)
{

  NM_extract_diag_block2(problem->M, contact, &local_problem->M->matrix0);

  local_problem->M->size0 = 2; // Necessary ?
  local_problem->M->size1 = 2;

  local_problem->q[0]= problem->q[contact*2];
  local_problem->q[1]= problem->q[contact*2+1];
  NM_row_prod_no_diag2(2 * problem->numberOfContacts, contact, 2*contact,
                       problem->M, reaction, local_problem->q, false);

  DEBUG_EXPR(NM_display(local_problem->M););
  DEBUG_EXPR(NV_display(local_problem->q,2););
}
static void shuffle(int size, int * randnum) //size is the given range
{
  int i;
  int swap, randindex;
  /* for(i=0;i<size;i++) */
  /* { */
  /*  printf("Array before shuffling is : %d\n",randnum[i]); */
  /* } */
  for(i = 0; i < size; ++i)
  {
    swap = randnum[i];
    randindex = rand() % size;
    randnum[i] = randnum[randindex];
    randnum[randindex] = swap;
  }
  /* printf("\n\n\n"); */
  /* for(i=0;i<size;i++) */
  /* { */
  /*  printf("Array after shuffling is : %d\n",randnum[i]); */
  /* } */
}


static
double light_error_squared( double localreaction[2],
                            double *oldreaction)
{
  return (pow(oldreaction[0] - localreaction[0], 2) +
          pow(oldreaction[1] - localreaction[1], 2) );
}
static
double squared_norm(double localreaction[2])
{
  return (pow(localreaction[0], 2) +
          pow(localreaction[1], 2) );
}

static
void accumulateLightErrorSum(double *light_error_sum, double localreaction[3],
                             double *oldreaction)
{
  *light_error_sum += ((oldreaction[0] - localreaction[0])*(oldreaction[0] - localreaction[0]) +
                       (oldreaction[1] - localreaction[1])*(oldreaction[1] - localreaction[1])) ;
}
static
double calculateLightError(double light_error_sum, unsigned int nc, double *reaction, double *norm_r)
{
  double error = sqrt(light_error_sum);
  *norm_r = cblas_dnrm2(nc*2, reaction, 1);
  if(fabs(*norm_r) > DBL_EPSILON)
    error /= (*norm_r);
  return error;
}
static
int determine_convergence(double error, double tolerance, int iter,
                          SolverOptions *options)
{
  int has_not_converged = 1;
  if(error < tolerance)
  {
    has_not_converged = 0;
    numerics_printf("-- FC2D - NSGS - Iteration %i "
                    "Residual = %14.7e < %7.3e\n", iter, error, tolerance);
  }
  else
  {
    numerics_printf("-- FC2D - NSGS - Iteration %i "
                    "Residual = %14.7e > %7.3e\n", iter, error, tolerance);
  }
  return has_not_converged;
}

static
double calculateFullErrorFinal(FrictionContactProblem *problem, SolverOptions *options,
                               double *reaction, double *velocity, double tolerance,
                               double norm_q)
{
  double absolute_error;
  /* (*computeError)(problem, reaction , velocity, tolerance, */
  /*                 options, norm_q, &absolute_error); */

  fc2d_compute_error(problem, reaction, velocity, tolerance, norm_q, &absolute_error);


  if(verbose > 0)
  {
    if(absolute_error > options->dparam[SICONOS_DPARAM_TOL])
    {
      numerics_printf("-- FC2D - NSGS - Warning absolute "
                      "Residual = %14.7e is larger than required precision = %14.7e",
                      absolute_error, options->dparam[SICONOS_DPARAM_TOL]);
    }
    else
    {
      numerics_printf("-- FC2D - NSGS - absolute "
                      "Residual = %14.7e is smaller than required precision = %14.7e",
                      absolute_error, options->dparam[SICONOS_DPARAM_TOL]);
    }
  }
  return absolute_error;
}



static
int determine_convergence_with_full_final(FrictionContactProblem *problem, SolverOptions *options,
                                          double *reaction, double *velocity,
                                          double *tolerance, double norm_q, double error,
                                          int iter)
{
  int has_not_converged = 1;
  if(error < *tolerance)
  {
    has_not_converged = 0;
    numerics_printf("-- FC2D - NSGS - Iteration %i "
                    "Residual = %14.7e < %7.3e", iter, error, *tolerance);

    double absolute_error = calculateFullErrorFinal(problem, options,
                                                    reaction, velocity,
                                                    options->dparam[SICONOS_DPARAM_TOL],
                                                    norm_q);
    if(absolute_error > options->dparam[SICONOS_DPARAM_TOL])
    {
      *tolerance = error/absolute_error*options->dparam[SICONOS_DPARAM_TOL];
      assert(*tolerance > 0.0 && "tolerance has to be positive");
      /* if (*tolerance < DBL_EPSILON) */
      /* { */
      /*   numerics_warning("determine_convergence_with_full_fina", "We try to set a very smal tolerance"); */
      /*   *tolerance = DBL_EPSILON; */
      /* } */
      numerics_printf("-- FC2D - NSGS - We modify the required incremental precision to reach accuracy to %e", *tolerance);
      has_not_converged = 1;
    }
    else
    {
      numerics_printf("-- FC2D - NSGS - The incremental precision is sufficient to reach accuracy to %e", *tolerance);
    }




  }
  else
  {
    numerics_printf("-- FC2D - NSGS - Iteration %i "
                    "Residual = %14.7e > %7.3e", iter, error, *tolerance);
  }
  return has_not_converged;
}

static inline int fc2dLocalSolve(double *W, double *q, double mu, double *P, double *U)
{
  double D, muPn;

  /* | Wnn Wnt |
     | Wtn Wtt | */

#define Wnn W[0]
#define Wtn W[1]
#define Wnt W[2]
#define Wtt W[3]

  if(q[0] > 0)
  {
    P[0] = 0;
    P[1] = 0;
  }
  else
  {
    /* solve WP + q = 0  */
    D = Wnn * Wtt - Wnt * Wtn;
    if(D < DBL_EPSILON) return(1);

    P[0] = - (Wtt * q[0] - Wnt * q[1]) / D;
    P[1] = - (-Wtn * q[0] + Wnn * q[1]) / D;

    muPn = mu * P[0];

    if(fabs(P[1]) > muPn)
      /* outside cone */
    {

      if(P[1] + muPn < 0)
      {

        P[0] = - q[0] / (Wnn - mu * Wnt);
        P[1] = - mu * P[0];
      }
      else
      {

        P[0] = - q[0] / (Wnn + mu * Wnt);
        P[1] = mu * P[0];

      }
    }
  }


#undef Wnn
#undef Wnt
#undef Wtn
#undef Wtt

  return(0);
}
static
unsigned int* allocfreezingContacts(FrictionContactProblem *problem,
                                    SolverOptions *options)
{
  unsigned int *fcontacts = 0;
  unsigned int nc = problem->numberOfContacts;
  if(options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0)
  {
    fcontacts = (unsigned int *) malloc(nc * sizeof(unsigned int));
    for(unsigned int i = 0; i < nc ; ++i)
    {
      fcontacts[i] = 0;
    }
  }
  return fcontacts;
}

void fc2d_nsgs(FrictionContactProblem* problem, double *z, double *w,
               int *info, SolverOptions* options)
{
  /* Notes:
     - we suppose that the trivial solution case has been checked before,
     and that all inputs differs from NULL since this function is
     supposed to be called from lcp_driver_global().
  */

  /* Global Solver parameters*/
  int * iparam = options->iparam;
  double * dparam = options->dparam;

  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  double tolerance = dparam[SICONOS_DPARAM_TOL];

  /* Matrix M/vector q of the LCP */
  double * q = problem->q;

  int nc = problem->numberOfContacts;
  double norm_q = cblas_dnrm2(nc*2, problem->q, 1);
  double norm_r[] = {1e24};

  /* Local problem initialization */
  LinearComplementarityProblem * local_problem = (LinearComplementarityProblem *)
    malloc(sizeof(*local_problem));
  
  local_problem->M = NM_create(NM_DENSE,2,2);
  local_problem->q = (double*)malloc(2 * sizeof(double));

  double localreaction[2];

  /*****  Gauss-Seidel iterations *****/
  int iter = 0; /* Current iteration number */
  double error = INFINITY; /* Current error */
  int has_not_converged = 1;

  if(iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] >0)
  {
    unsigned int *freeze_contacts = NULL;
    freeze_contacts = allocfreezingContacts(problem, options);
    while((iter < itermax) && has_not_converged)
    {
      ++iter;

      double light_error_sum = 0.0;
      double light_error_2 = 0.0;
      /* Loop over the rows of blocks in blmat */
      /* contact: current row (of blocks) number */
      int number_of_freezed_contact=0;
      if(iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] >0)
      {
        for(unsigned int i = 0 ; i < nc ; ++i)
        {
          if (freeze_contacts[i] >0)
            number_of_freezed_contact++;
        }
      }
      for(unsigned int pos = 0, contact = 0; contact < nc; ++contact, ++pos, ++pos)
      {
        if (freeze_contacts[contact] >0)
        {
          /* we skip freeze contacts */
          freeze_contacts[contact] -=  1;
          continue;
        }

        /* store  old reaction */
        localreaction[0] = z[pos];
        localreaction[1] = z[pos+1];

        /* Local problem formalization */
        fc2d_nsgs_buildLocalProblem(contact, problem, local_problem, z);

        /* Solve local problem */
        int local_solver_info = fc2dLocalSolve(local_problem->M->matrix0,
                                               local_problem->q,
                                               problem->mu[contact],
                                               localreaction, &w[pos]);


        /* verbose if problem */
        if(local_solver_info)
        {
          /* Number of GS iterations */
          iparam[SICONOS_IPARAM_ITER_DONE] = iter;
          dparam[SICONOS_DPARAM_RESIDU] = error;
          fprintf(stderr,
                  "fc2d_nsgs error: local LCP solver failed at global iteration %d.\n", iter);
          fprintf(stderr,
                  "                for block-row number %d. Output info equal to %d.\n", contact, local_solver_info);
          *info = local_solver_info;
          goto free_and_return;
        }
        /* if(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT || */
        /*    iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL */
        /*   ) */
        /*   accumulateLightErrorSum(&light_error_sum, localreaction, &z[pos]); */

        light_error_2= light_error_squared(localreaction, &z[pos]);
        light_error_sum += light_error_2;

        if ((light_error_2*squared_norm(localreaction) <= tolerance*tolerance/(nc*nc*10)
             || squared_norm(localreaction) <=  (*norm_r* *norm_r/(nc*nc*1000)))
            && iter >=10)
        {
          /* we  freeze the contact for n iterations*/
          if (number_of_freezed_contact < nc-1)
          {
            number_of_freezed_contact++;
            freeze_contacts[contact] = iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT];
          }
          else
          {
            numerics_printf_verbose(2,"Number of freezed contacts too large w.r.t number of contact. we defreeze all contacts\n");
            for(unsigned int c = 0 ; c < nc ; ++c)  freeze_contacts[c] =0;
          }
          DEBUG_EXPR
            (printf("first criteria : light_error_2*squared_norm(localreaction) <= tolerance*tolerance/(nc*nc*10) ==> %e <= %e\n",
                    light_error_2*squared_norm(localreaction), tolerance*tolerance/(nc*nc*10));
             printf("second criteria :  squared_norm(localreaction) <=  (*norm_r* *norm_r/(nc*nc))/1000. ==> %e <= %e\n",
                    squared_norm(localreaction) ,  (*norm_r* *norm_r/(nc*nc))/1000.);
             printf("Contact % i is freezed for %i steps\n", contact,  iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT]);
              );
        }
        /* reaction update */
        z[pos]   = localreaction[0];
        z[pos+1] = localreaction[1];


      } //end for loop

      DEBUG_EXPR(
        int frozen_contact=0;
        for(unsigned int ii = 0 ; ii < nc ; ++ii) if (freeze_contacts[ii] >0) frozen_contact++;
        numerics_printf_verbose(2,"number of frozen contacts %i at iter : %i", frozen_contact, iter );
        );

      /* error evaluation */
      if(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
         SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT)
      {
        error = calculateLightError(light_error_sum, nc, z,  norm_r);
        has_not_converged = determine_convergence(error, tolerance, iter, options);
      }
      else if(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
              SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
      {
        error = calculateLightError(light_error_sum, nc, z, norm_r);
        has_not_converged = determine_convergence_with_full_final(problem,  options,
                                                                  z, w,
                                                                  &tolerance, norm_q, error, iter);
      }
    } //end while loop
    free(freeze_contacts);
  }
  else
  {
    while((iter < itermax) && has_not_converged)
    {
      ++iter;
      double light_error_sum = 0.0;
      /* Loop over the rows of blocks in blmat */
      for(unsigned int pos = 0, contact = 0; contact < nc; ++contact, ++pos, ++pos)
      {
        /* store  old reaction */
        localreaction[0] = z[pos];
        localreaction[1] = z[pos+1];

        /* Local problem formalization */
        fc2d_nsgs_buildLocalProblem(contact, problem, local_problem, z);

        /* Solve local problem */
        int local_solver_info = fc2dLocalSolve(local_problem->M->matrix0,
                                               local_problem->q,
                                               problem->mu[contact],
                                               localreaction, &w[pos]);


        /* verbose if problem */
        if(local_solver_info)
        {
          /* Number of GS iterations */
          iparam[SICONOS_IPARAM_ITER_DONE] = iter;
          dparam[SICONOS_DPARAM_RESIDU] = error;
          fprintf(stderr,
                  "fc2d_nsgs error: local LCP solver failed at global iteration %d.\n", iter);
          fprintf(stderr,
                  "                for block-row number %d. Output info equal to %d.\n", contact, local_solver_info);
          *info = local_solver_info;
          goto free_and_return;
        }
        if(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT ||
           iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL
          )
          accumulateLightErrorSum(&light_error_sum, localreaction, &z[pos]);

        z[pos]   = localreaction[0];
        z[pos+1] = localreaction[1];

      }// end for loop

      /*  error evaluation */
      if(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT)
      {
        error = calculateLightError(light_error_sum, nc, z, norm_r);
        has_not_converged = determine_convergence(error, tolerance, iter, options);
      }
      else if(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
      {
        error = calculateLightError(light_error_sum, nc, z, norm_r);
        has_not_converged = determine_convergence_with_full_final(problem,  options,
                                                                  z, w,
                                                                  &tolerance, norm_q, error,iter);
      }
    } //end while loop
  }
  /* Full criterium */
  if(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
  {
    error = calculateFullErrorFinal(problem, options, z, w,
                                    tolerance, norm_q);

    has_not_converged = determine_convergence(error,  dparam[SICONOS_DPARAM_TOL], iter, options);


  }


  // numerics_printf("Siconos Numerics : problem size=%d, nb iterations=%d, error=%g\n",
  //          blmat->blocknumber0,
  //          iter,
  //          error);

  *info = has_not_converged;

  /* Number of GS iterations */
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

  /* Resulting error */
  dparam[SICONOS_DPARAM_RESIDU] = error;

free_and_return:
  free(local_problem->q);
  free(local_problem->M);
  free(local_problem);
}



void fc2d_nsgs_dense(FrictionContactProblem* problem, double *reaction, double *velocity, int *info, SolverOptions* options)
{
  int nc = problem->numberOfContacts;
  double * vec = problem->M->matrix0;
  double *q = problem->q;
  double * mu = problem->mu;

  int i, j, k, kk, iter;
  int n = 2 * nc;
  int it_end = 0;
  int  incx, incy;

  double alpha, beta;
  double *y, res = INFINITY;
  double normr, avn, avt, det, gplus, gmoins;
  double apn, apt, zn, zt, den1, num1;
  double alm1;
  double aln1;
  int pivot;
  double factor1;
  int * randomContactList;

  int maxit = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double errmax = options->dparam[SICONOS_DPARAM_TOL];
  options->iparam[SICONOS_IPARAM_ITER_DONE]  = 0;
  options->dparam[SICONOS_DPARAM_RESIDU]  = 0.0;

  iter         = 0;

  y       = (double*) malloc(n  * sizeof(double));

  randomContactList = (int*) malloc(nc  * sizeof(int));

  for(i = 0; i < nc; i++)
  {
    randomContactList[i] = i;
  }


  for(i = 0; i < n; i++)
  {

    reaction[i]  = 0.0 ;
    velocity[i]  = 0.0 ;
  }

  normr    =   1.;

  while((iter < maxit) && (normr > errmax))
  {
    iter = iter + 1 ;

    if(options->iparam[SICONOS_IPARAM_NSGS_SHUFFLE] > 0)
    {
      shuffle(nc, randomContactList);
    }



    /*         Loop over contacts                */



    for(kk = 0; kk < nc; kk++)
    {


      i = randomContactList[kk];

      avn = 0.;
      avt = 0.;
      apn = 0.;
      apt = 0.;

      for(j = 0; j <= 2 * i - 1; j++)
      {

        avn = avn + vec[j * n + 2 * i] * reaction[j];
        avt = avt + vec[j * n + 2 * i + 1] * reaction[j];

      }

      for(k = 2 * i + 2; k < n; k++)
      {
        apn = apn + vec[k * n + 2 * i] * reaction[k];
        apt = apt + vec[k * n + 2 * i + 1] * reaction[k];
      }

      zn    = -q[2 * i] - avn - apn;
      zt    = -q[2 * i + 1] - avt - apt;








      if(-zn >= 0.0)
      {


        reaction[2 * i]   = 0.0; // PN
        velocity[2 * i]   = -zn; // UN
        reaction[2 * i + 1] = 0.0; // PT
        velocity[2 * i + 1] = -zt; // UT


      }
      else
      {

        velocity[2 * i]   = 0.0;
        velocity[2 * i + 1] = 0.0;


        det    = vec[2 * i + 2 * i * n] * vec[(2 * i + 1) + (2 * i + 1) * n]
          - vec[(2 * i + 1) + (2 * i) * n] * vec[(2 * i) + (2 * i + 1) * n];

        if(fabs(det) < 100* DBL_EPSILON)
        {
          if(verbose > 0)
          {
            printf("--------------- FC2D - NSGS -  Warning small denominator : %g . use of partial pivoting\n", fabs(det));
          }
          /* free(y); */
          /* free(randomContactList); */
          /* *info = 2; */
          /* return; */


          alm1 = fabs(vec[2 * i + (2 * i) * n]);
          aln1 = fabs(vec[(2 * i +1) + n * (2 * i)]);
          pivot = alm1 >= aln1 ? 0 : 1;
          switch(pivot)
          {
          case 0:
            if(alm1 < DBL_EPSILON)
            {
              *info = 1;
              return;
            }
            factor1 = vec[(2 * i +1) + n * (2 * i)]/vec[2 * i + (2 * i) * n];
            reaction[2 * i + 1]  = (zt - factor1*zn)/(vec[(2 * i + 1) + n * (2 * i + 1)] - factor1*vec[2 * i + (2 * i + 1) * n]);
            reaction[2 * i] = (zn - vec[2 * i + (2 * i + 1) * n]*reaction[2 * i + 1])/vec[2 * i + (2 * i) * n];
            break;
          case 1:
            if(aln1 < DBL_EPSILON)
            {
              *info = 1;
              return;
            }
            factor1 = vec[2 * i + (2 * i) * n]/vec[(2 * i +1) + n * (2 * i)];
            reaction[2 * i + 1]  = (zn - factor1*zt)/(vec[2 * i + (2 * i + 1) * n] - factor1*vec[(2 * i + 1) + n * (2 * i + 1)]);
            reaction[2 * i]= (zt - vec[(2 * i + 1) + n * (2 * i + 1)]*reaction[2 * i + 1])/vec[(2 * i +1) + n * (2 * i)];
            break;
          default:
            exit(EXIT_FAILURE);
          }
          DEBUG_PRINTF("contact %i , reaction[2 * i] = %g, reaction[2 * i + 1] = % g \n", i,  reaction[2 * i], reaction[2 * i + 1]);


        }
        else
        {
          reaction[2 * i]   = (zn * vec[(2 * i + 1) + n * (2 * i + 1)] - zt * vec[2 * i + (2 * i + 1) * n]) / det;
          reaction[2 * i + 1] = (-zn * vec[(2 * i +1) + n * (2 * i)]   + zt * vec[2 * i + (2 * i) * n]) / det;
          DEBUG_PRINTF("contact %i , reaction[2 * i] = %g, reaction[2 * i + 1] = % g \n", i, reaction[2 * i], reaction[2 * i + 1]);
        }

        if((reaction[2 * i] >= 0.0) && ((fabs(reaction[2 * i + 1]) - mu[i] * reaction[2 * i]) <= 0.0))
        {
          DEBUG_PRINTF("--------------- FC2D - NSGS - contact %i, Stick status \n", i);
        }
        else
        {


          velocity[2 * i]   = 0.0;


          gplus  = vec[2 * i + 2 * i * n] + mu[i] * vec[(2 * i) + (2 * i + 1) * n];


          if(fabs(gplus) < 1e-12)
          {
            if(verbose > 0)
              printf("--------------- FC2D - NSGS -  Warning small denominator (gplus) : %g \n", fabs(gplus));

            free(y);
            free(randomContactList);

            *info = 2;
            return;

          }
          else
          {

            velocity[2 * i + 1] = -zt + (zn / gplus) * (vec[2 * i + (2 * i + 1) * n] + mu[i] * vec[(2 * i + 1) + (2 * i + 1) * n]);

            reaction[2 * i]   = zn / gplus;
            reaction[2 * i + 1] = mu[i] * reaction[2 * i];

          }

          if((reaction[2 * i] >= 0.0) && (velocity[2 * i + 1] <= 0.0))
          {

            /*    printf("Slip+ status\n");*/

          }
          else
          {

            velocity[2 * i]   = 0.0;

            gmoins = vec[2 * i + 2 * i * n] - mu[i] * vec[(2 * i) + (2 * i + 1) * n];


            if(fabs(gmoins) < 1e-12)
            {
              if(verbose > 0)
                printf("--------------- FC2D - NSGS -  Warning small denominator (gmoins) : %g \n", fabs(gmoins));

              free(y);
              free(randomContactList);

              *info = 2;
              return;

            }
            else
            {


              velocity[2 * i + 1] = -zt + (zn / gmoins) * (vec[2 * i + (2 * i + 1) * n] - mu[i] * vec[(2 * i + 1) + (2 * i + 1) * n]);

              reaction[2 * i]   = zn / gmoins;
              reaction[2 * i + 1] = -mu[i] * reaction[2 * i];
            }

            /* printf("Slip- status\n");*/
          }
        }
      }

    }



    /*          Convergence criterium           */

    incx = 1;
    incy = 1;

    cblas_dcopy(n, q, incx, y, incy);

    alpha = 1.;
    beta  = 1.;
    cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, alpha, vec, n, reaction, incx, beta, y, incy);



    alpha = -1.;
    cblas_daxpy(n, alpha, velocity, incx, y, incy);


    num1 = cblas_ddot(n, y, incx, y, incy);
    den1 = cblas_ddot(n, q, incx, q, incy);


    normr = sqrt(num1 / den1);


    it_end = iter;
    res    = normr;
    if(verbose > 0)
      printf("--------------- FC2D - NSGS - Iteration %i "
             "Residual = %14.7e < %7.3e\n", iter, res, errmax);
  }


  options->iparam[SICONOS_IPARAM_ITER_DONE] = it_end;
  options->dparam[SICONOS_DPARAM_RESIDU] = res;



  if(normr > errmax)
  {

    if(verbose > 0)
      printf("--------------- FC2D - NSGS - No convergence after %i iterations"
             " residual = %14.7e < %7.3e\n", iter, res, errmax);

    *info = 1;
  }
  else
  {

    if(verbose > 0)
      printf("--------------- FC2D - NSGS - Convergence after %i iterations"
             " residual = %14.7e < %7.3e\n", iter, res, errmax);

    *info = 0;
  }




  free(y);
  free(randomContactList);



}





void fc2d_nsgs_set_default(SolverOptions *options)
{
  options->iparam[SICONOS_IPARAM_NSGS_SHUFFLE] = 0;
  options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] = SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL;
  //  useful only for the sparse nsgs case.
}

// options setup is done through fc2d_nsgs_set_default.
