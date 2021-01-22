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
  unsigned int *freeze_contacts = NULL;
  freeze_contacts = allocfreezingContacts(problem, options);
  int maxit = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double tolerance = options->dparam[SICONOS_DPARAM_TOL];
  options->iparam[SICONOS_IPARAM_ITER_DONE]  = 0;
  options->dparam[SICONOS_DPARAM_RESIDU]  = 0.0;

  int * iparam = options->iparam;
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

  while((iter < maxit) && (normr > tolerance))
  {
    iter = iter + 1 ;

    if(options->iparam[SICONOS_IPARAM_NSGS_SHUFFLE] > 0)
    {
      shuffle(nc, randomContactList);
    }

    double light_error_sum = 0.0;
    double light_error_2 = 0.0;
    double norm_r=  cblas_dnrm2(nc*2, reaction, 1);
    /*         Loop over contacts                */



    for(kk = 0; kk < nc; kk++)
    {


      i = randomContactList[kk];
      if(iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] >0)
      {
        if (freeze_contacts[i] >0)
        {
          /* we skip freeze contacts */
          freeze_contacts[i] -=  1;
          continue;
        }
      }
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




      double oldreaction[2];
      oldreaction[0]=reaction[2 * i];
      oldreaction[1]=reaction[2 * i+1];


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
        light_error_2= (oldreaction[0] - reaction[2*i])*(oldreaction[0] - reaction[2*i])
          + (oldreaction[1] - reaction[2*i+1])*(oldreaction[1] - reaction[2*i+1]);
        light_error_sum += light_error_2;
        double squared_norm=          reaction[2*i]*reaction[2*i]+ reaction[2*i+1]* reaction[2*i+1];
        if(iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] >0)
        {
          if ((light_error_2*squared_norm <= tolerance*tolerance/(nc*nc*10)
               || squared_norm <=  (norm_r* norm_r/(nc*nc*1000)))
              && iter >=10)
          {
            /* we  freeze the contact for n iterations*/
            //printf("first criteria : light_error_2*squared_norm(localreaction) <= tolerance*tolerance/(nc*nc*10) ==> %e <= %e\n", light_error_2*squared_norm(localreaction), tolerance*tolerance/(nc*nc*10));
            //printf("second criteria :  squared_norm(localreaction) <=  (*norm_r* *norm_r/(nc*nc))/1000. ==> %e <= %e\n",  squared_norm(localreaction) ,  (*norm_r* *norm_r/(nc*nc))/1000.);
            //printf("Contact % i is freezed for %i steps\n", contact,  iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT]);
            numerics_printf_verbose(2,"Contact % i is freezed for %i steps", i,  iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT]);
            freeze_contacts[i] = iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] ;
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
             "Residual = %14.7e < %7.3e\n", iter, res, tolerance);
  }


  options->iparam[SICONOS_IPARAM_ITER_DONE] = it_end;
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




  free(y);
  free(randomContactList);



}
void fc2d_nsgs_set_default(SolverOptions *options)
{
  options->iparam[SICONOS_IPARAM_NSGS_SHUFFLE] = 0;
  options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] = SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL;
  options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] =0;
  //  useful only for the sparse nsgs case.
}
