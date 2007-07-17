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


/*!\file lcp_newton_FB.c
 *
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).
 * Try \f$(z,w)\f$ such that:
 *
 * \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *   0 \le z \perp M z + b = w \ge 0\\
 *   \end{array}
 *  \right.
 * \f$
 *
 * M is an (n x n)  matrix , q , w and z n-vectors.
 *
 *
 *!\fn lcp_newton_FB(  double *vec , double *q , int *nn , int *itermax , double *tol , int *ispeak , double *z , double *w , int *it_end , double *res , int *info )
 *
 * lcp_newton_min use a nonsmooth newton method based on the Fischer-Bursmeister convex function
 *
 * \f$
 *   0 \le z \perp w \ge 0 \Longrightarrow \phi(z,w)=\sqrt{z^2+w^2}-(z+w)=0
 * \f$

 * \f$
 *   \Phi(z) = \left[ \begin{array}{c}  \phi(z_1,w_1) \\ \phi(z_1,w_1) \\ \vdots \\  \phi(z_n,w_n)  \end{array}\right] =0\\
 * \f$
 *
 *
 * References: Alart & Curnier 1990, Pang 1990
 *
 * Generic lcp parameters:\n
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
 * \param dparamLCP[0] = tol     Input unchanged parameter which represents the tolerance required.
 * \param dparamLCP[1] = res     Output modified parameter which returns the final error value.
 *
 * \author Vincent Acary
 *
 * \todo Optimizing the memory allocation (Try to avoid the copy of JacH into A)
 * \todo Add rules for the computation of the penalization rho
 * \todo Add a globalization strategy based on a decrease of a merit function. (Nonmonotone LCP) Reference in Ferris Kanzow 2002
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"

void lcp_newton_FB(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP)
{


  int i, j, iter;
  int n = *nn, m, k;
  int itermax, ispeak;

  integer incx, incy;
  char NOTRANS = 'N';
  char TRANS = 'T';
  double err, tol, a1, b1;
  double alpha, normi;
  int infoDGESV;

  int *ipiv;
  double *beta, *mbeta;
  double *JacPhi, *JacPhi_copy, *Phi;

  printf("The Algorithm lcp_newton_FB is not reliable yet, report to siconos.gforge.inria.fr if you need it soon \n");
  return ;


  incx = 1;
  incy = 1;
  /*input*/

  itermax = iparamLCP[0];
  ispeak  = iparamLCP[1];

  tol   = dparamLCP[0];




  /*output*/

  iparamLCP[2] = 0;
  dparamLCP[1] = 0.0;

  for (i = 0; i < n; i++) z[i] = 0.0;


  // Creation of the gradient of the function H

  JacPhi = (double *)malloc(n * n * sizeof(double));
  JacPhi_copy   = (double *)malloc(n * n * sizeof(double));

  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++) JacPhi[j * n + i] = 0.0;


  // Creation of the RHS Phi,
  Phi = (double *)malloc(n * sizeof(double));
  for (i = 0; i < n; i++) Phi[i] = 0.0;

  beta = (double *)malloc(m * sizeof(double));
  mbeta = (double *)malloc(m * sizeof(double));

  iter = 0;
  err  = 1.;


  // Newton Iteration
  while ((iter < itermax) && (err > tol))
  {
    ++iter;

    // Construction of the directional derivatives of Phi, JacPhi
    // q --> w
    dcopy_((integer *)&n , q , &incx , w , &incy);
    // Mz+q --> w
    a1 = 1.;
    b1 = 1.;
    dgemv_(&TRANS , (integer *)&n , (integer *)&n , &a1 , vec , (integer *)&n , z , &incx , &b1 , w , &incy);
    for (i = 0; i < n; i++) printf("z[%i]=%e", i, z[i]);
    printf("\n");
    for (i = 0; i < n; i++) printf("w[%i]=%e", i, w[i]);
    printf("\n");



    for (i = 0; i < n; i++)
    {
      if ((z[i] == 0) && (w[i] == 0))
      {
        beta[i] = 1.0;
      }
      else
      {
        beta[i] = 0.0;
      }

    }
    for (i = 0; i < n; i++) printf("beta[%i]=%e", i, beta[i]);
    printf("\n");

    // M^T.beta --> mbeta
    a1 = 1.;
    b1 = 0.0;
    dgemv_(&NOTRANS , (integer *)&n , (integer *)&n , &a1 , vec , (integer *)&n , beta , &incx , &b1 , mbeta  , &incy);
    for (i = 0; i < n; i++) printf("mbeta[%i]=%e", i, mbeta[i]);
    printf("\n");



    for (i = 0; i < n; i++)
    {

      if ((z[i] == 0) && (w[i] == 0))
      {
        normi = sqrt(beta[i] * beta[i] + mbeta[i] * mbeta[i]);
        for (j = 0; j < n; j++)
        {
          JacPhi[j * n + i] = (mbeta[i] / normi - 1.0) * vec[j * n + i];
        }
        JacPhi[i * n + i] += (beta[i] / normi - 1.0);

      }
      else
      {
        normi = (z[i] * z[i] + w[i] * w[i]);
        printf("normi=%e", normi);
        for (j = 0; j < n; j++)
        {
          JacPhi[j * n + i] = (w[i] / normi - 1.0) * vec[j * n + i];

        }
        JacPhi[i * n + i] += (z[i] / normi - 1.0);


      }

    }
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++) printf("JacPhi[%i][%i]=%e", i, j, JacPhi[j * n + i]);
      printf("\n");
    }

    // Computation of the value Phi
    for (i = 0; i < n; i++)
    {
      Phi[i] = sqrt(z[i] * z[i] + w[i] * w[i]) - (z[i] + w[i]);

    }




    // Computation of the element of the subgradient.

    dcopy_((integer *)&n , JacPhi , &incx , JacPhi_copy , &incy);
    k = 1;
    F77NAME(dgesv)((integer *)&m, (integer *)&k, JacPhi_copy, (integer *)&m, (integer *)ipiv, beta, (integer *)&m, (integer *)&infoDGESV);

    if (infoDGESV)
    {
      if (ispeak > 0)
      {
        printf("Problem in DGESV\n");
      }
      iparamLCP[2] = iter;
      dparamLCP[1] = err;
      *info = 2;

      return ;

    }


    // iteration
    alpha = -1.0;
    daxpy_((integer *)&n , &alpha , beta , &incx , z , &incy);     //  z-beta --> z


    // Construction of the RHS for the next iterate and for the error evaluation
    // q --> w
    dcopy_((integer *)&n , q , &incx , w , &incy);
    // Mz+q --> w
    a1 = 1.;
    b1 = 1.;
    dgemv_(&TRANS , (integer *)&n , (integer *)&n , &a1 , vec , (integer *)&n , z , &incx , &b1 , w , &incy);

    for (i = 0; i < n; i++)
    {
      Phi[i] = sqrt(z[i] * z[i] + w[i] * w[i]) - (z[i] + w[i]);

    }


    // Error Evaluation



    err = dnrm2_((integer *)&n , Phi , &incx);
    err = 1 / 2 * err * err;

  }

  iparamLCP[2] = iter;
  dparamLCP[1] = err;

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

  free(JacPhi);
  free(JacPhi_copy);
  free(Phi);
  free(beta);
  free(mbeta);



}
