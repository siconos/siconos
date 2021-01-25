/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#include <float.h>                   // for DBL_EPSILON
#include <math.h>                    // for fabs, sqrt, INFINITY
#include <stdio.h>                   // for printf
#include <stdlib.h>                  // for free, calloc, malloc, exit, rand
#include "FrictionContactProblem.h"  // for FrictionContactProblem
#include "Friction_cst.h"            // for SICONOS_FRICTION_2D_NSGS
#include "NumericsFwd.h"             // for SolverOptions, FrictionContactPr...
#include "NumericsMatrix.h"          // for NumericsMatrix, RawNumericsMatrix
#include "SolverOptions.h"           // for SolverOptions, solver_options_nu...
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"                   // for DEBUG_PRINTF
#include "fc2d_Solvers.h"            // for fc2d_nsgs, fc2d_nsgs_setDefaultS...
#include "numerics_verbose.h"        // for verbose
#include "SiconosBlas.h"                   // for cblas_ddot, cblas_daxpy, cblas_d...

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
static
int localContactolver(double * reaction, double * velocity,
                      double* vec, double *mu, double *q,
                      int i,
                      int n)
{
  int info =0;
  verbose=2;
  double  avn=0.0, avt=0.0, det=0.0, gplus=0.0, gmoins=0.0;
  double apn, apt, zn, zt;
  double alm1;
  double aln1;
  int pivot;
  double factor1;
  int i2 = 2*i;
  int i2_p_1 = i2+1;
  for(int j = 0; j <= i2 - 1; j++)
  {
    avn = avn + vec[j * n + i2] * reaction[j];
    avt = avt + vec[j * n + i2_p_1] * reaction[j];
  }
  for(int k = i2 + 2; k < n; k++)
  {
    apn = apn + vec[k * n + i2] * reaction[k];
    apt = apt + vec[k * n + i2_p_1] * reaction[k];
  }
  zn    = -q[i2] - avn - apn;
  zt    = -q[i2_p_1] - avt - apt;

  if(-zn >= 0.0)
  {
    reaction[i2]   = 0.0; // PN
    velocity[i2]   = -zn; // UN
    reaction[i2_p_1] = 0.0; // PT
    velocity[i2_p_1] = -zt; // UT
  }
  else
  {
    velocity[i2]   = 0.0;
    velocity[i2_p_1] = 0.0;
    det    = vec[i2 + i2 * n] * vec[(i2_p_1) + (i2_p_1) * n]
      - vec[(i2_p_1) + (i2) * n] * vec[(i2) + (i2_p_1) * n];
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

      alm1 = fabs(vec[i2 + (i2) * n]);
      aln1 = fabs(vec[(i2_p_1) + n * (i2)]);
      pivot = alm1 >= aln1 ? 0 : 1;
      switch(pivot)
      {
      case 0:
        if(alm1 < DBL_EPSILON)
        {
          info = 1;
          return info;
        }
        factor1 = vec[(i2_p_1) + n * (i2)]/vec[i2 + (i2) * n];
        reaction[i2_p_1]  = (zt - factor1*zn)/(vec[(i2_p_1) + n * (i2_p_1)] - factor1*vec[i2 + (i2_p_1) * n]);
        reaction[i2] = (zn - vec[i2 + (i2_p_1) * n]*reaction[i2_p_1])/vec[i2 + (i2) * n];
        break;
      case 1:
        if(aln1 < DBL_EPSILON)
        {
          info = 1;
          return 1;
        }
        factor1 = vec[i2 + (i2) * n]/vec[(i2_p_1) + n * (i2)];
        reaction[i2_p_1]  = (zn - factor1*zt)/(vec[i2 + (i2_p_1) * n] - factor1*vec[(i2_p_1) + n * (i2_p_1)]);
        reaction[i2]= (zt - vec[(i2_p_1) + n * (i2_p_1)]*reaction[i2_p_1])/vec[(i2_p_1) + n * (i2)];
        break;
      default:
        exit(EXIT_FAILURE);
      }
      DEBUG_PRINTF("contact %i , reaction[i2] = %g, reaction[i2_p_1] = % g \n", i,  reaction[i2], reaction[i2_p_1]);
    }
    else
    {
      reaction[i2]   = (zn * vec[(i2_p_1) + n * (i2_p_1)] - zt * vec[i2 + (i2_p_1) * n]) / det;
      reaction[i2_p_1] = (-zn * vec[(i2_p_1) + n * (i2)]   + zt * vec[i2 + (i2) * n]) / det;
      DEBUG_PRINTF("contact %i , reaction[i2] = %g, reaction[i2_p_1] = % g \n", i, reaction[i2], reaction[i2_p_1]);
    }

    if((reaction[i2] >= 0.0) && ((fabs(reaction[i2_p_1]) - mu[i] * reaction[i2]) <= 0.0))
    {
      DEBUG_PRINTF("--------------- FC2D - NSGS - contact %i, Stick status \n", i);
    }
    else
    {
      velocity[i2]   = 0.0;

      gplus  = vec[i2 + i2 * n] + mu[i] * vec[(i2) + (i2_p_1) * n];
      if(fabs(gplus) < 1e-12)
      {
        if(verbose > 0)
          printf("--------------- FC2D - NSGS -  Warning small denominator (gplus) : %g \n", fabs(gplus));
        info = 2;
        return info;

      }
      else
      {
        velocity[i2_p_1] = -zt + (zn / gplus) * (vec[i2 + (i2_p_1) * n] + mu[i] * vec[(i2_p_1) + (i2_p_1) * n]);

        reaction[i2]   = zn / gplus;
        reaction[i2_p_1] = mu[i] * reaction[i2];
      }

      if((reaction[i2] >= 0.0) && (velocity[i2_p_1] <= 0.0))
      {
        /*    printf("Slip+ status\n");*/
      }
      else
      {
        velocity[i2]   = 0.0;

        gmoins = vec[i2 + i2 * n] - mu[i] * vec[(i2) + (i2_p_1) * n];
        if(fabs(gmoins) < 1e-12)
        {
          if(verbose > 0)
            printf("--------------- FC2D - NSGS -  Warning small denominator (gmoins) : %g \n", fabs(gmoins));
          info = 2;
          return info;
        }
        else
        {
          velocity[i2_p_1] = -zt + (zn / gmoins) * (vec[i2 + (i2_p_1) * n] - mu[i] * vec[(i2_p_1) + (i2_p_1) * n]);

          reaction[i2]   = zn / gmoins;
          reaction[i2_p_1] = -mu[i] * reaction[i2];
        }

        /* printf("Slip- status\n");*/
      }
    }
  }
  return info;
}

static void fc2d_nsgs_new(FrictionContactProblem* problem, double *reaction, double *velocity, int *info, SolverOptions* options)
{
  verbose=3;
  int nc = problem->numberOfContacts;
  double * vec = problem->M->matrix0;
  double *q = problem->q;
  double * mu = problem->mu;

  int iter;
  int n = 2 * nc;
  int it_end = 0;
  int  incx, incy;

  double alpha, beta;
  double *y, res = INFINITY;
  double normr;
  double num1, den1;

  int maxit = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double tolerance = options->dparam[SICONOS_DPARAM_TOL];
  options->iparam[SICONOS_IPARAM_ITER_DONE]  = 0;
  options->dparam[SICONOS_DPARAM_RESIDU]  = 0.0;

  int * iparam = options->iparam;
  iter         = 0;

  y       = (double*) malloc(n  * sizeof(double));

  
  int * randomContactList =NULL;
  randomContactList = (int*) malloc(nc  * sizeof(int));
  for(int i = 0; i < nc; i++)
  {
    randomContactList[i] = i;
  }
  
  unsigned int *freeze_contacts = NULL;
 
  for(int i = 0; i < n; i++)
  {
    reaction[i]  = 0.0 ;
    velocity[i]  = 0.0 ;
  }

  normr    =   1.;
  double oldreaction[2];

  if(iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] >0)
  {
    freeze_contacts = allocfreezingContacts(problem, options);

    while((iter < maxit) && (normr > tolerance))
    {
      iter = iter + 1 ;

      if(iparam[SICONOS_IPARAM_NSGS_SHUFFLE] > 0)
      {
        shuffle(nc, randomContactList);
      }

      double light_error_sum = 0.0;
      double light_error_2 = 0.0;
      double norm_r=  cblas_dnrm2(nc*2, reaction, 1);
      /*         Loop over contacts                */

      for(int kk = 0; kk < nc; kk++)
      {


        int i = randomContactList[kk];
      
        if(iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] >0)
        {
          if(freeze_contacts[i] >0)
          {
            /* we skip freeze contacts */
            freeze_contacts[i] -=  1;
            continue;
          }
        }
      
        oldreaction[0]=reaction[2 * i];
        oldreaction[1]=reaction[2 * i+1];

        int info_local= localContactolver(reaction, velocity,
                                          vec, mu, q,
                                          i, n);

        if (info_local >0)
        {
          goto free_variable;
        }
      
        light_error_2= (oldreaction[0] - reaction[2*i])*(oldreaction[0] - reaction[2*i])
          + (oldreaction[1] - reaction[2*i+1])*(oldreaction[1] - reaction[2*i+1]);
        light_error_sum += light_error_2;
        double squared_norm=  reaction[2*i]*reaction[2*i]+ reaction[2*i+1]* reaction[2*i+1];
        if(iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] >0)
        {
          if((light_error_2*squared_norm <= tolerance*tolerance/(nc*nc*10)
              || squared_norm <= (norm_r* norm_r/(nc*nc*1000)))
             && iter >=10)
          {
            /* we  freeze the contact for n iterations*/
            numerics_printf_verbose(2,"Contact % i is freezed for %i steps", i,  iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT]);
            freeze_contacts[i] = iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] ;
          }
        }
      }
      if(iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] >0)
      {
        int frozen_contact=0;
        for(unsigned int ii = 0 ; ii < nc ; ++ii)
        {
          if (freeze_contacts[ii] >0)
          {
            frozen_contact++;
          }
        }
        printf("number of frozen contacts %i at iter : %i\n", frozen_contact, iter );
      }
      /*          Convergence criterium           */
      /*  Convergence criterium */
      cblas_dcopy(n, q, 1, y, 1);
      cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, 1., vec, n, reaction, 1, 1., y, 1);
      alpha = -1.;
      cblas_daxpy(n, -1., velocity, 1, y, 1);
      num1 = cblas_ddot(n, y, 1, y, 1);
      den1 = cblas_ddot(n, q, 1, q, 1);
      
      normr = sqrt(num1 / den1);
      it_end = iter;
      res    = normr;
      if(verbose > 0)
        printf("--------------- FC2D - NSGS - Iteration %i "
               "Residual = %14.7e < %7.3e\n", iter, res, tolerance);
    }
  }
  else
  {
    while((iter < maxit) && (normr > tolerance))
    {
      iter = iter + 1 ;

      if(iparam[SICONOS_IPARAM_NSGS_SHUFFLE] > 0)
      {
        shuffle(nc, randomContactList);
      }
      /*         Loop over contacts                */

      for(int kk = 0; kk < nc; kk++)
      {
        int i = randomContactList[kk];
        int info_local= localContactolver(reaction, velocity,
                                          vec, mu, q,
                                          i, n);
        
        if (info_local>0)
        {
          goto free_variable;
        }
      }
      /*  Convergence criterium */
      cblas_dcopy(n, q, 1, y, 1);
      cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, 1., vec, n, reaction, 1, 1., y, 1);
      alpha = -1.;
      cblas_daxpy(n, -1., velocity, 1, y, 1);
      num1 = cblas_ddot(n, y, 1, y, 1);
      den1 = cblas_ddot(n, q, 1, q, 1);
      
      normr = sqrt(num1 / den1);
      it_end = iter;
      res    = normr;
      if(verbose > 0)
        printf("--------------- FC2D - NSGS - Iteration %i "
               "Residual = %14.7e < %7.3e\n", iter, res, tolerance);
    }
  }
  
  iparam[SICONOS_IPARAM_ITER_DONE] = it_end;
  options->dparam[SICONOS_DPARAM_RESIDU] = res;



  if(normr > tolerance)
  {

    if(verbose > 0)
      printf("--------------- FC2D - NSGS - No convergence after %i iterations"
             " residual = %14.7e < %7.3e\n", iter, res, tolerance);

    *info = 1;
  }
  else
  {

    if(verbose > 0)
      printf("--------------- FC2D - NSGS - Convergence after %i iterations"
             " residual = %14.7e < %7.3e\n", iter, res, tolerance);

    *info = 0;
  }


free_variable:

  free(y);
  free(randomContactList);
  if(iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] >0)
    free(freeze_contacts);


}
void fc2d_nsgs_set_default(SolverOptions *options)
{
  options->iparam[SICONOS_IPARAM_NSGS_SHUFFLE] = 0;
  options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] = SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL;
  options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] =0;
  //  useful only for the sparse nsgs case.
}

void fc2d_nsgs(FrictionContactProblem* problem, double *reaction, double *velocity, int *info, SolverOptions* options)
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
