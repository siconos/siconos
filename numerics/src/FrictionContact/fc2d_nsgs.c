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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SiconosBlas.h"
#include "fc2d_Solvers.h"
#include "numerics_verbose.h"
#include "NumericsMatrix.h"

void shuffle(int size, int * randnum);
void shuffle(int size, int * randnum) //size is the given range
{
  int i;
  int swap, randindex;
  /* for(i=0;i<size;i++) */
  /* { */
  /*  printf("Array before shuffling is : %d\n",randnum[i]); */
  /* } */
  for (i = 0; i < size; ++i)
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


void fc2d_nsgs(FrictionContactProblem* problem , double *reaction , double *velocity , int *info, SolverOptions* options)
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
  double apn, apt, zn , zt, den1, num1;

  int * randomContactList;

  int maxit = options->iparam[0];
  double errmax = options->dparam[0];
  options->iparam[1]  = 0;
  options->dparam[1]  = 0.0;

  iter         = 0;

  y       = (double*) malloc(n  * sizeof(double));

  randomContactList = (int*) malloc(nc  * sizeof(int));

  for (i = 0; i < nc; i++)
  {
    randomContactList[i] = i;
  }


  for (i = 0; i < n; i++)
  {

    reaction[i]  = 0.0 ;
    velocity[i]  = 0.0 ;
  }

  normr    =   1.;

  while ((iter < maxit) && (normr > errmax))
  {
    iter = iter + 1 ;

    if (options->iparam[2] > 0)
    {
      shuffle(nc, randomContactList);
    }



    /*         Loop over contacts                */



    for (kk = 0; kk < nc; kk++)
    {


      i = randomContactList[kk];

      avn = 0.;
      avt = 0.;
      apn = 0.;
      apt = 0.;

      for (j = 0; j <= 2 * i - 1; j++)
      {

        avn = avn + vec[j * n + 2 * i] * reaction[j];
        avt = avt + vec[j * n + 2 * i + 1] * reaction[j];

      }

      for (k = 2 * i + 2; k < n; k++)
      {
        apn = apn + vec[k * n + 2 * i] * reaction[k];
        apt = apt + vec[n * k + 2 * i + 1] * reaction[k];
      }

      zn    = -q[2 * i] - avn - apn;
      zt    = -q[2 * i + 1] - avt - apt;








      if (-zn >= 0.0)
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


        det    = vec[2 * i + 2 * i * n] * vec[(2 * i + 1) + (2 * i + 1) * n] - vec[(2 * i) + (2 * i + 1) * n] * vec[(2 * i) + (2 * i + 1) * n];

        if (fabs(det) < 1e-12)
        {

          if (verbose > 0)
            printf(" Warning denominator nul\n");

          free(y);
          free(randomContactList);
          *info = 2;
          return;

        }
        else
        {


          reaction[2 * i]   = (zn * vec[(2 * i + 1) + n * (2 * i + 1)] - zt * vec[2 * i + (2 * i + 1) * n]) / det;
          reaction[2 * i + 1] = (-zn * vec[(2 * i) + n * (2 * i + 1)] + zt * vec[2 * i + (2 * i) * n]) / det;

        }

        if ((reaction[2 * i] >= 0.0) && ((fabs(reaction[2 * i + 1]) - mu[i] * reaction[2 * i]) <= 0.0))
        {

          /*  printf("Stick status \n");*/
        }
        else
        {


          velocity[2 * i]   = 0.0;


          gplus  = vec[2 * i + 2 * i * n] + mu[i] * vec[(2 * i) + (2 * i + 1) * n];


          if (fabs(gplus) < 1e-12)
          {

            if (verbose > 0)
              printf(" Warning denominator nul\n");

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

          if ((reaction[2 * i] >= 0.0) && (velocity[2 * i + 1] <= 0.0))
          {

            /*    printf("Slip+ status\n");*/

          }
          else
          {

            velocity[2 * i]   = 0.0;

            gmoins = vec[2 * i + 2 * i * n] - mu[i] * vec[(2 * i) + (2 * i + 1) * n];


            if (fabs(gmoins) < 1e-12)
            {

              if (verbose > 0)
                printf(" Warning denominator nul\n");

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

  }


  options->iparam[1] = it_end;
  options->dparam[1] = res;



  if (normr > errmax)
  {
    if (verbose > 0)
      printf("No convergence after %d iterations, the residue is %g\n", iter, normr);

    *info = 1;
  }
  else
  {

    if (verbose > 0)
      printf("Convergence after %d iterations, the residue is %g \n", iter, normr);

    *info = 0;
  }




  free(y);
  free(randomContactList);



}
int fc2d_nsgs_setDefaultSolverOptions(SolverOptions *options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the 2D NSGS Solver\n");
  }

  /*  strcpy(options->solverName,"NLGS");*/
  options->solverId = SICONOS_FRICTION_2D_NSGS;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-4;
  options ->internalSolvers = NULL;
  return 0;
}

