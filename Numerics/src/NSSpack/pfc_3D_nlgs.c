/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
/*!\file pfc_3D_nlgs.c
 *
 * This subroutine allows the primal resolution of contact problems with friction.\n
 *
 *   Try \f$(z,w)\f$ such that:\n
 *   \f$
 *    \left\lbrace
 *     \begin{array}{l}
 *      M z + q = w\\
 *      0 \le z_n \perp w_n \ge 0\\
 *      -w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
 *     \end{array}
 *    \right.
 *   \f$
 *
 * here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.
 *
 * \fn  pfc_3D_nlgs( int *nn , double *vec , double *q , double *z , double *w , int *info\n,
 *                   int *iparamLCP , double *dparamLCP )
 *
 * Generic pfc_3D parameters:\n
 *
 * \param nn      Unchanged parameter which represents the dimension of the system.
 * \param vec     Unchanged parameter which contains the components of the matrix with a fortran storage.
 * \param q       Unchanged parameter which contains the components of the right hand side vector.
 * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
 * \param w       Modified parameter which returns the solution of the problem.
 * \param info    Modified parameter which returns the termination value\n
 *                0 - convergence\n
 *                1 - iter = itermax\n
 *                2 - negative diagonal term\n
 *
 * Specific NLGS parameters:\n
 *
 * \param iparamLCP[0] = itermax Input unchanged parameter which represents the maximum number of iterations allowed.
 * \param iparamLCP[1] = ispeak  Input unchanged parameter which represents the output log identifiant\n
 *                       0 - no output\n
 *                       0 < active screen output\n
 * \param iparamLCP[2] = it_end  Output modified parameter which returns the number of iterations performed by the algorithm.
 *
 * \param dparamLCP[0] = mu      Input unchanged parameter which represents the friction coefficient.
 * \param dparamLCP[1] = tol     Input unchanged parameter which represents the tolerance required.
 * \param dparamLCP[2] = res     Output modified parameter which returns the final error value.
 *
 *
 * \author Mathieu Renouf & Nineb Sheherazade .
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"
#include <time.h>

void pfc_3D_nlgs(int *nn , double *vec , double *q , double *z , double *w , int *info,
                 int *iparamLCP , double *dparamLCP)
{


  FILE *f101;

  int n, in, it, is, ispeak, itermax, nc, i, iter;
  double err, zn , zt, zs, den, mrn, num, tol, mu, mu2;
  double qs, a1, b1;
  int incx, incy;
  double *diag, *ww;

  clock_t t1, t2;

  t1 = clock();

  ispeak = 0;
  nc     = *nn;
  incx   = 1;
  incy   = 1;
  n      = 3 * nc;

  /* Recup input */

  itermax = iparamLCP[0];
  ispeak  = iparamLCP[1];

  mu  = dparamLCP[0];
  tol = dparamLCP[1];
  mu2 = mu * mu;

  /* Initialize output */

  iparamLCP[2] = 0;
  dparamLCP[2] = 0.0;


  if (ispeak == 2) f101 = fopen("pfc_3D_nlgs.log" , "w+");

  iter = 0;

  /* Allocation */

  ww   = (double*)malloc(n * sizeof(double));
  diag = (double*)malloc(n * sizeof(double));

  /* Check for non trivial case */

  qs = DNRM2(n , q , incx);

  if (ispeak > 0) printf("\n ||q||= %g \n" , qs);

  if (qs > 1e-16) den = 1.0 / qs;
  else
  {
    for (i = 0 ; i < n ; ++i)
    {
      w[i] = 0.;
      z[i] = 0.;
    }

    free(ww);
    free(diag);
    *info = 0;
    return;
  }

  /* Intialization of w */

  for (i = 0 ; i < n ; ++i)
  {
    ww[i] = 0.;
    w[i]  = 0.;
  }



  /* Preparation of the diagonal of the inverse matrix */

  for (i = 0 ; i < nc ; ++i)
  {
    in = 3 * i;
    if (fabs(vec[in * n + in]) < 1e-16)
    {

      if (ispeak > 0)
      {
        printf(" Warning negative diagonal term \n");
        printf(" The local problem can be solved \n");
      }

      *info = 2;
      free(diag);
      free(ww);

      return;
    }
    else
    {
      diag[in  ] = 1.0 / vec[(in) * n + in  ];
      diag[in + 1] = 1.0 / vec[(in + 1) * n + in + 1];
      diag[in + 2] = 1.0 / vec[(in + 2) * n + in + 2];
    }
  }


  /*start iterations*/

  iter = 0;
  err  = 1.;

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    incx = 1;
    incy = 1;

    DCOPY(n , w , incx , ww , incy);
    DCOPY(n , q , incx ,  w , incy);

    for (i = 0 ; i < nc ; ++i)
    {

      in = 3 * i;
      it = 3 * i + 1;
      is = 3 * i + 2;

      incx = n;
      incy = 1;

      z[in] = 0.0;
      z[it] = 0.0;
      z[is] = 0.0;

      zn = q[in] + DDOT(n , &vec[in] , incx , z , incy);
      zt = q[it] + DDOT(n , &vec[it] , incx , z , incy);
      zs = q[is] + DDOT(n , &vec[is] , incx , z , incy);

      if (zn > 0.0)
      {
        z[3 * i  ] = 0.0;
        z[3 * i + 1] = 0.0;
        z[3 * i + 2] = 0.0;
      }
      else
      {
        z[in] = -zn * diag[in];
        z[it] = -zt * diag[it];
        z[is] = -zs * diag[is];

        mrn = z[it] * z[it] + z[is] * z[is];

        if (mrn > mu2 * z[in]*z[in])
        {
          num = mu * z[in] / sqrt(mrn);
          z[it] = z[it] * num;
          z[is] = z[is] * num;
        }
      }
    }

    /* **** Criterium convergence **** */

    incx =  1;
    incy =  1;

    a1 = 1.0;
    b1 = 1.0;

    DGEMV(LA_NOTRANS , n , n , a1 , vec , n , z , incx , b1 , w , incy);

    qs   = -1.0;
    DAXPY(n , qs , w , incx , ww , incy);

    num = DNRM2(n, ww , incx);
    err = num * den;

    /*  printf("Iteration %i Erreur = %24.8e\n",iter,err); */

    if (ispeak == 2) for (i = 0 ; i < n ; ++i) fprintf(f101, "%i  %i  %14.7e\n", iter - 1, i, z[i]);

  }

  iparamLCP[2] = iter;
  dparamLCP[2] = err;

  if (ispeak > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of NLGS after %i iterations\n" , iter);
      printf(" The residue is : %e \n", err);
      *info = 1;
    }
    else
    {
      printf(" Convergence of NLGS after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 0;
    }
  }
  else
  {
    if (err > tol) *info = 1;
    else *info = 0;
  }

  free(diag);
  free(ww);

  if (ispeak == 2) fclose(f101);
  t2 = clock();
  /*  printf("%.4lf seconds of processing\n", (t2-t1)/(double)CLOCKS_PER_SEC); */

}
