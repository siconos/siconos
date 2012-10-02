/* Siconos-Numerics, Copyright INRIA 2005-2012.
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

#include "FrictionContact2D_Solvers.h"
#include "LA.h"
#include <math.h>

void FrictionContact2D_cpg(FrictionContactProblem* problem , double *reaction , double *velocity , int *info, SolverOptions* options)
{
  int nc = problem->numberOfContacts;
  double * vec = problem->M->matrix0;
  double * mu = problem->mu;

  int       n = 2 * nc, i, iter;
  int       incx = 1, incy = 1;
  int       *stat, *statusi, it_end;


  double    eps = 1.e-12;
  double    pAp, alpha, beta, wAp, rp, normr;
  double    alphaf, betaf, den, num, res;

  double    *p, *fric, *r;
  double    *fric1, *v, *w, *Ap, *xi, *z;


  int maxit = options->iparam[0];
  double tol = options->dparam[0];
  options->iparam[1]  = 0;
  options->dparam[1]  = 0.0;


  r       = (double*) malloc(n * sizeof(double));
  fric    = (double*) malloc(n * sizeof(double));
  p       = (double*) malloc(n * sizeof(double));
  v       = (double*) malloc(n * sizeof(double));
  w       = (double*) malloc(n * sizeof(double));
  Ap      = (double*) malloc(n * sizeof(double));
  xi      = (double*) malloc(n * sizeof(double));
  z       = (double*) malloc(n * sizeof(double));
  fric1   = (double*) malloc(n * sizeof(double));


  stat    = (int*)    malloc(nc * sizeof(int));
  statusi = (int*)    malloc(nc * sizeof(int));





  for (i = 0; i < n ; i++)
  {
    reaction[i]     = 0.0;
    xi[i]    = 0.0;
    r[i]     = 0.0;
    v[i]     = 0.0;
    p[i]     = 0.0;
    w[i]     = 0.0;
    Ap[i]    = 0.0;
    z[i]     = 0.0;
    fric1[i] = 1.0;
    fric[i]  = mu[i] * fric1[i];

    if (i < nc)
    {
      stat[i]    = 0;
      statusi[i] = 0;

    }

  }



  DCOPY(n, problem->q, incx, r, incy);

  alphaf = -1.;
  betaf  = -1.;

  DGEMV(LA_NOTRANS, n, n, alphaf, vec, n, reaction, incx, betaf, r, incy);





  /*             Check for initial status             */


  for (i = 0; i < nc; i++)
  {
    mu[i] = fric[i];
    if (reaction[2 * i] <= eps)
    {
      /*       No contact            */
      stat[i] = 0;
    }
    else if (reaction[2 * i + 1] <=  -mu[i]*reaction[2 * i])
    {
      /*     Slide backward         */
      stat[i] = 1;
    }
    else if (reaction[2 * i + 1] >=  mu[i]*reaction[2 * i])
    {
      /*   Slide forward          */
      stat[i] = 3;
    }
    else
    {
      /*     Stick contact        */
      stat[i] = 2;
    }
  }


  iter  = 0;
  normr = 1.0;



  while ((iter < maxit) && (normr > tol))
  {



    for (i = 0 ; i < nc ; i++)
      statusi[i] = stat[i];


    DCOPY(n, r, incx, v, incy);

    if (iter == 0)
    {
      DCOPY(n, r, incx, w, incy);

      DCOPY(n, w, incx, p, incy);
    }

    alphaf = 1.0;
    betaf  = 0.0;
    DGEMV(LA_NOTRANS, n, n, alphaf, vec, n, p, incx, betaf, Ap, incy);

    pAp    = DDOT(n, p, incx, Ap, incy);

    /*}
    else
    {
    alphaf = 1.0;
    betaf  = 0.0;
    dgemv_( &notrans, (integer *)&n, (integer *)&n, &alphaf, vec, (integer *)&n, p, &incx, &betaf, Ap, &incy );

    pAp    = ddot_( (integer *)&n, p, &incx, Ap, &incy );*/

    if (pAp == 0)
    {
      if (verbose > 0)
        printf("\n Operation non conform alpha at the iteration %d \n", iter);

      free(r);
      free(fric);
      free(p);
      free(v);
      free(w);
      free(Ap);
      free(xi);
      free(z);
      free(fric1);
      free(stat);
      free(statusi);

      *info = 2;

      return;
    }

    /*} */

    rp     = DDOT(n, r, incx, p, incy);

    alpha  = rp / pAp;

    DCOPY(n, reaction, incx, xi, incy);

    alphaf = alpha;
    DAXPY(n, alphaf, p, incx, xi, incy);

    FrictionContact2D_projc(xi, &n, statusi, p, fric, reaction, stat);


    /*         r(:)=b(:)-matmul(A,x)          */

    DCOPY(n, problem->q, incx, r, incy);

    alphaf = -1.;
    betaf  = -1.;
    DGEMV(LA_NOTRANS, n, n, alphaf, vec, n, reaction, incx, betaf, r, incy);

    FrictionContact2D_projf(statusi, &n, r, fric, w);

    FrictionContact2D_projf(statusi, &n, p, fric, z);


    wAp    = DDOT(n, w, incx, Ap, incy);

    beta   = - wAp / pAp;

    DCOPY(n, w, incx, p, incy);

    alphaf  = beta;
    DAXPY(n, alphaf, z, incx, p, incy);


    /*  alphaf  = 1.;
    betaf   = 0.;
    dgemv_( &notrans, (integer *)&n, (integer *)&n, &alphaf, vec , (integer *)&n, p, &incx, &betaf, Ap, &incy );

    pAp     = ddot_( (integer *)&n, p, &incx, Ap, &incy );*/

    DCOPY(n, r, incx, xi, incy);

    alphaf  = -1.;
    DAXPY(n, alphaf, v, incx, xi, incy);

    num     = DDOT(n, xi, incx, xi, incy);

    den     = DDOT(n, v, incx, v, incy);

    normr   = sqrt(num / den);

    it_end  = iter;
    res     = normr;


    options->iparam[1] = it_end;
    options->dparam[1] = res;


    iter = iter + 1;

  }




  if (normr < tol)
  {

    if (verbose > 0)
      printf("convergence after %d iterations with a residual %g\n", iter - 1, normr);

    *info = 0;


  }
  else
  {
    if (verbose > 0)
      printf("no convergence after %d iterations with a residual %g\n", iter - 1, normr);

    *info = 1;
  }


  alpha = -1.;
  DSCAL(n , alpha , r , incx);

  DCOPY(n, r, incx, velocity, incy);



  free(fric);
  free(p);
  free(v);
  free(w);
  free(Ap);
  free(xi);
  free(z);
  free(fric1);
  free(stat);
  free(statusi);
  free(r);




}
int frictionContact2D_cpg_setDefaultSolverOptions(SolverOptions *options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the CPG Solver\n");
  }

  /*  strcpy(options->solverName,"CPG");*/
  options->solverId =  SICONOS_FRICTION_2D_CPG;
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

