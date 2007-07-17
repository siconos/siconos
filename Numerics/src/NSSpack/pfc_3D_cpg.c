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
/*!\file pfc_3D_cpg.c
 *
 * This subroutine allows the primal resolution of contact problems with friction.\n
 *
 * Try \f$(z,w)\f$ such that:\n
 *  \f$
 *   \left\lbrace
 *    \begin{array}{l}
 *     M z + q = w\\
 *     0 \le z_n \perp w_n \ge 0\\
 *     -w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
 *    \end{array}
 *   \right.
 *  \f$
 *
 * here M is an (n x n) matrix, q, z and w n-vector.
 *
 * \fn  pfc_3D_cpg( int *nn , double *vec , double *q , double *z , double *w , int *info\n,
 *                  int *iparamLCP , double *dparamLCP )
 *
 * pfc_3D_cpg is a specific cpg (conjugated projected gradient) for primal contact problem with friction.\n
 * Ref: Renouf, M. and Alart, P. "" Comp. Method Appl. Mech. Engrg. (2004).
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
 *
 * Specific CPG parameters:\n
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
 * \author Nineb Sheherazade & Mathieu Renouf.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"
void pfc_3D_projf(int  , double * , double * , double * , double * , int *);
void pfc_3D_projc(int  , double , double *z , double * , int *);


void pfc_3D_cpg(int *nn , double *vec , double *q , double *z , double *w , int *info,
                int *iparamLCP , double *dparamLCP)
{

  FILE *f101;
  int n, nc;
  integer incx, incy;
  int i, iter, itermax, ispeak;

  double err, a1, b1, qs;

  double alpha, beta, rp, pMp;
  double den, num, tol, mu;

  char NOTRANS = 'N';

  int *status;
  double *zz , *pp , *rr , *ww , *Mp;

  printf("pfc_3D_* Algorithms are not reliable yet, report to siconos.gforge.inria.fr if you need it soon \n");
  return ;

  nc   = *nn;
  incx = 1;
  incy = 1;
  n    = 3 * nc;

  /* Recup input */

  itermax = iparamLCP[0];
  ispeak  = iparamLCP[1];

  mu  = dparamLCP[0];
  tol = dparamLCP[1];

  /* Initialize output */

  iparamLCP[2] = 0;
  dparamLCP[2] = 0.0;

  if (ispeak == 2) f101 = fopen("pfc_3D_cpg.log" , "w+");

  qs = dnrm2_((integer *)&n , q , &incx);

  if (ispeak > 0) printf("\n ||q||= %g \n" , qs);

  if (qs > 1e-16) den = 1.0 / qs;
  else
  {
    for (i = 0 ; i < n ; ++i)
    {
      w[i] = 0.;
      z[i] = 0.;
    }
    *info = 9;
    if (ispeak == 2) fclose(f101);
    return;
  }

  /* Allocations */

  status = (int*)malloc(nc * sizeof(int));

  ww = (double*)malloc(n * sizeof(double));
  rr = (double*)malloc(n * sizeof(double));
  pp = (double*)malloc(n * sizeof(double));
  zz = (double*)malloc(n * sizeof(double));

  Mp = (double*)malloc(n * sizeof(double));

  incx = 1;

  for (i = 0; i < n; ++i)
  {

    ww[i] = 0.;
    rr[i] = 0.;
    pp[i] = 0.;
    zz[i] = 0.;

    Mp[i] = 0.;

  }

  for (i = 0; i < nc ; ++i) status[i] = 0;


  /* rr = -Wz + q */

  incy = 1;

  dcopy_((integer *)&n , q , &incx , rr , &incy);

  a1 = -1.;
  b1 = -1.;

  dgemv_(&NOTRANS , (integer *)&n , (integer *)&n , &a1 , vec , (integer *)&n , z , &incx , &b1 , rr , &incy);

  /* Initialization of gradients */
  /* rr -> p and rr -> w */

  dcopy_((integer *)&n , rr , &incx , ww , &incy);
  dcopy_((integer *)&n , rr , &incx , pp , &incy);

  iter = 0;
  err  = 1.0 ;

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    /* Compute initial pMp */

    incx = 1;
    incy = 1;

    dcopy_((integer *)&n , pp , &incx , Mp , &incy);

    a1 = 1.0;
    b1 = 0.0;

    dgemv_(&NOTRANS , (integer *)&n , (integer *)&n , &a1 , vec , (integer *)&n , Mp , &incx , &b1 , w , &incy);

    pMp = ddot_((integer *)&n , pp , &incx , w  , &incy);
    //printf( " pWp = %10.4g  \n", pMp );
    if (fabs(pMp) < 1e-16)
    {

      if (ispeak > 0)
      {
        printf(" Operation no conform at the iteration %d \n", iter);
        printf(" Alpha can be obtained with pWp = %10.4g  \n", pMp);
      }

      free(Mp);
      free(ww);
      free(rr);
      free(pp);
      free(zz);
      free(status);

      iparamLCP[2] = iter;
      dparamLCP[2] = err;

      if (ispeak == 2) fclose(f101);
      *info = 3;
      return;
    }

    rp  = ddot_((integer *)&n , pp , &incx , rr , &incy);
    //printf( " rp = %10.4g  \n", rp );
    alpha = rp / pMp;

    /*
     * Iterate prediction:
     * z' = z + alpha*p
     *
     */

    daxpy_((integer *)&n , &alpha , pp , &incx , z , &incy);

    /* Iterate projection*/

    pfc_3D_projc(nc , mu , z , pp , status);

    /* rr = -Wz + q */

    dcopy_((integer *)&n , rr , &incx , w  , &incy);
    dcopy_((integer *)&n , q  , &incx , rr , &incy);

    a1 = -1.;
    b1 = -1.;

    dgemv_(&NOTRANS , (integer *)&n , (integer *)&n , &a1 , vec , (integer *)&n , z , &incx , &b1 , rr , &incy);

    /* Gradients projection
     * rr --> ww
     * pp --> zz
     */

    pfc_3D_projf(nc , ww , zz , rr , pp , status);

    /*   beta = -w.Mp / pMp  */

    rp = ddot_((integer *)&n , ww , &incx, w , &incy);

    beta = -rp / pMp;

    dcopy_((integer *)&n , ww , &incx , pp , &incy);
    daxpy_((integer *)&n, &beta , zz , &incx , pp , &incy);

    /* **** Criterium convergence **** */

    qs   = -1.0;
    daxpy_((integer *)&n , &qs , rr , &incx , w , &incy);

    /*
     * Note: the function dnrm2_ leads to the generation of core
     */
    num = sqrt(ddot_((integer *)&n , w , &incx, w , &incy));

    err = num * den;

    if (ispeak == 2) for (i = 0 ; i < n ; ++i) fprintf(f101, "%d  %d  %14.7e\n", iter , i , z[i]);

  }

  iparamLCP[2] = iter;
  dparamLCP[2] = err;

  dcopy_((integer *)&n , rr , &incx , w , &incy);

  qs   = -1.0;
  dscal_((integer *)&n , &qs , w , &incx);

  if (ispeak > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of CPG after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 1;
    }
    else
    {
      printf(" Convergence of CPG after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 0;
    }
  }
  else
  {
    if (err > tol) *info = 1;
    else *info = 0;
  }


  free(Mp);

  free(ww);
  free(rr);
  free(pp);
  free(zz);
  free(status);

  if (ispeak == 2) fclose(f101);

}
