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
#include "blaslapack.h"

void pfc_3D_nlgs(int *nn , double *vec , double *q , double *z , double *w , int *info,
                 int *iparamLCP , double *dparamLCP)
{


  FILE *f101;

  int n, in, it, ispeak, itermax, incx, incy, nc, i, iter;
  double err, zn , zt, den, num, dft, dfn, tol, mu;
  double qs, a1, b1;

  double *det, *bfd, *ww;
  char NOTRANS = 'N';

  ispeak = 0;
  nc    = *nn;
  incx = 1;
  incy = 1;
  n    = 2 * nc;

  /* Recup input */

  itermax = iparamLCP[0];
  ispeak  = iparamLCP[1];

  mu  = dparamLCP[0];
  tol = dparamLCP[1];

  /* Initialize output */

  iparamLCP[2] = 0;
  dparamLCP[2] = 0.0;


  if (ispeak == 2) f101 = fopen("pfc_3D_nlgs.log" , "w+");

  iter = 0;

  /* Allocation */

  ww  = (double*)malloc(n * sizeof(double));
  det = (double*)malloc(nc * sizeof(double));
  bfd = (double*)malloc(nc * sizeof(double));

  /* Check for non trivial case */

  qs = dnrm2_(&n , q , &incx);

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
    free(det);
    free(bfd);
    *info = 0;
    return;
  }

  /* Intialization of w */

  for (i = 0 ; i < n ; ++i)
  {
    ww[i] = 0.;
    w[i]  = 0.;
  }

  incx = 1;
  incy = 1;

  dcopy_(&n , q , &incx , w , &incy);

  /* Preparation of the diagonal of the inverse matrix */

  for (i = 0 ; i < nc ; ++i)
  {
    in = 2 * i;
    if (fabs(vec[in * n + in]) < 1e-16)
    {

      if (ispeak > 0)
      {
        printf(" Warning negative diagonal term \n");
        printf(" The local problem can be solved \n");
      }

      *info = 2;
      free(det);
      free(bfd);
      free(ww);

      return;
    }
    else
    {
      it = 2 * i + 1;
      det[i] = vec[in * (n + 1)] * vec[it * (n + 1)] - vec[it * n + in] * vec[in * n + it];
      bfd[i] = mu * vec[it * n + in] / vec[in * (n + 1)];
    }
  }

  /*start iterations*/

  iter = 0;
  err  = 1.;

  incx = 1;
  incy = 1;

  dcopy_(&n , q , &incx , w , &incy);

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    incx = 1;
    incy = 1;

    dcopy_(&n , w , &incx , ww , &incy);
    dcopy_(&n , q , &incx ,  w , &incy);

    for (i = 0 ; i < nc ; ++i)
    {

      in = 2 * i;
      it = 2 * i + 1;

      incx = n;
      incy = 1;

      z[in] = 0.0;
      z[it] = 0.0;

      zn = q[in] + ddot_(&n , &vec[in] , &incx , z , &incy);
      zt = q[it] + ddot_(&n , &vec[it] , &incx , z , &incy);

      if (zn > 0.0)
      {
        z[2 * i  ] = 0.0;
        z[2 * i + 1] = 0.0;
      }
      else
      {
        dft = vec[in * (n + 1)] * zt - vec[in * n + it ] * zn;
        dfn = -vec[it * n + in ] * zt - vec[it * (n + 1)] * zn;

        a1 = dft + mu * dfn;
        b1 = dft - mu * dfn;

        if (a1 > 0.0)
        {
          z[2 * i  ] = -(zn / vec[in * (n + 1)]) / (1 - bfd[i]);
          z[2 * i + 1] = -mu * z[2 * i];
        }
        else if (b1 > 0.0)
        {
          z[2 * i  ] = -(zn / vec[in * (n + 1)]) / (1 + bfd[i]);
          z[2 * i + 1] =  mu * z[2 * i];
        }
        else
        {
          z[2 * i  ] = -dfn / det[i];
          z[2 * i + 1] = -dft / det[i];
        }
      }
    }

    /* **** Criterium convergence **** */

    incx =  1;
    incy =  1;

    a1 = 1.0;
    b1 = 1.0;

    dgemv_(&NOTRANS , &n , &n , &a1 , vec , &n , z , &incx , &b1 , w , &incy);

    qs   = -1.0;
    daxpy_(&n , &qs , w , &incx , ww , &incy);

    num = dnrm2_(&n, ww , &incx);
    err = num * den;

    if (ispeak == 2) for (i = 0 ; i < n ; ++i) fprintf(f101, "%d  %d  %14.7e\n", iter - 1, i, z[i]);

  }

  iparamLCP[2] = iter;
  dparamLCP[2] = err;

  if (ispeak > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of NLGS after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
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

  free(bfd);
  free(det);
  free(ww);

  if (ispeak == 2) fclose(f101);

}
