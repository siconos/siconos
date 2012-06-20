/* Siconos-Numerics, Copyright INRIA 2005-2011.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"
#include "FrictionContact2D_Solvers.h"

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


void FrictionContact2D_nsgs(FrictionContactProblem* problem , double *reaction , double *velocity , int *info, SolverOptions* options)
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
  double normr, eps, avn, avt, det, gplus, gmoins;
  double apn, apt, zn , zt, den1, num1;

  int * randomContactList;

  int maxit = options->iparam[0];
  double errmax = options->dparam[0];
  options->iparam[1]  = 0;
  options->dparam[1]  = 0.0;

  iter         = 0;
  eps          = 1.e-08;

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

    DCOPY(n, q, incx, y, incy);

    alpha = 1.;
    beta  = 1.;
    DGEMV(LA_NOTRANS, n, n, alpha, vec, n, reaction, incx, beta, y, incy);



    alpha = -1.;
    DAXPY(n, alpha, velocity, incx, y, incy);


    num1 = DDOT(n, y, incx, y, incy);
    den1 = DDOT(n, q, incx, q, incy);


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
int frictionContact2D_nsgs_setDefaultSolverOptions(SolverOptions *options)
{
  int i;
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
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-4;
  options ->internalSolvers = NULL;
  return 0;
}

