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
/*!\file lcp_psor.c

  This subroutine allows the resolution of MLCP (Mixed Linear Complementary Problem).\n
  Try \f$(u,v,w)\f$ such that:\n
  \f$
   \left\lbrace
    \begin{array}{l}
    A u + Cv +a =0\\
    D u + Bv +b = w
    0 \le v \perp  w \ge 0\\
    \end{array}
   \right.
  \f$
  where  A is an (\f$ n \times n\f$ ) matrix, B is an (\f$ m \times m\f$ ) matrix,  C is an (\f$ n \times m\f$ ) matrix, D is an (\f$ m \times n\f$ ) matrix,    a and u is an (\f$ n \f$ ) vectors b,v and w is an (\f$ m \f$ ) vectors.
*/
/*!\fn   void mlcp_pgs( int *nn , int* mm, double *A , double *B , double *C , double *D , double *a  double *b, double *u, double *v, double *w , int *info ,   int *iparamMLCP , double *dparamMLCP );


  mlcp_psor (projected successive overrelaxation method) is a solver for MLCP.\n

  \param A            On enter, a (\f$n \times n\f$)-vector of doubles which contains the components of the "A" MLCP matrix with a Fortran storage.
  \param B            On enter, a (\f$m \times m\f$)-vector of doubles which contains the components of the "B" MLCP matrix with a Fortran storage.
  \param C            On enter, a (\f$n \times m\f$)-vector of doubles which contains the components of the "C" MLCP matrix with a Fortran storage.
  \param D            On enter, a (\f$m \times n\f$)-vector of doubles which contains the components of the "D" MLCP matrix with a Fortran storage.
  \param a            On enter, a n-vector of doubles which contains the components of the constant right hand side vector.
  \param b            On enter, a m-vector of doubles which contains the components of the constant right hand side vector.
  \param nn           On enter, an integer which represents one the dimension of the MLCP problem.
  \param mm           On enter, an integer which represents one the dimension of the MLCP problem.
  \n \n
  \param u            On return, a n-vector of doubles which contains the solution of the problem.
  \param v            On return, a m-vector of doubles which contains the solution of the problem.
  \param w            On return, a m-vector of doubles which contains the complementary solution of the problem.
  \param info    On return, an integer which returns the termination value:\n
                 0 : convergence\n
                 1 : iter = itermax\n
                 2 : negative diagonal term

  \param iparamMLCP  On enter/return a vector of integers:\n
                - iparamMLCP[0] = itermax On enter, the maximum number of iterations allowed.
                - iparamMLCP[1] = verbose  On enter, the output log identifiant:\n
                        0 : no output\n
                        >0: active screen output\n
                - iparamMLCP[2] = it_end  On enter, the number of iterations performed by the algorithm.

  \param dparamMLCP  On enter/return a vector of doubles:\n
                - dparamMLCP[0] = tol     On enter, the tolerance required.
                - dparamMLCP[1] = omega   On enter, the relaxation parameter (not yet available).
                - dparamMLCP[2] = res     On return, the final error value.

  \author Vincent Acary

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "LA.h"
#include <math.h>

int mlcp_compute_error(int* n, int* mm,  double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v,  int chat, double *w, double * error);

void mlcp_psor(int *nn , int* mm, double *A , double *B , double *C , double *D , double *a, double *b, double *u, double *v, double *w , int *info ,  int *iparamMLCP , double *dparamMLCP)
{


  int n, m, incx, incy, incAx, incAy, incBx, incBy;
  int i, iter;
  int itermax, verbose;
  int incxn;
  double qs, err, den, vi;
  double tol, omega;
  double *wOld, *diagA, *diagB;

  n = *nn;
  m = *mm;
  incx = 1;
  incy = 1;
  incxn = n;
  /* Recup input */

  itermax = iparamMLCP[0];
  verbose  = iparamMLCP[1];

  tol   = dparamMLCP[0];
  omega = dparamMLCP[2];
  printf("omega %f\n", omega);

  /* Initialize output */

  iparamMLCP[2] = 0;
  dparamMLCP[2] = 0.0;

  /* Allocation */

  wOld   = (double*)malloc(m * sizeof(double));
  diagA = (double*)malloc(n * sizeof(double));
  diagB = (double*)malloc(m * sizeof(double));

  /* Check for non trivial case */
  /*   qs = DNRM2( n , q , incx ); */

  /*   if( verbose > 0 ) printf("\n ||q||= %g \n",qs); */

  /*   den = 1.0/qs; */

  /*   for( i = 0 ; i < n ; ++i ){ */
  /*     wOld[i] = 0.; */
  /*     w[i] = 0.; */
  /*   } */

  /* Intialization of w */

  incx = 1;
  incy = 1;
  DCOPY(m , b , incx , w , incy);   // b -> w

  /* Preparation of the diagonal of the inverse matrix */

  for (i = 0 ; i < n ; ++i)
  {
    if ((fabs(A[i * n + i]) < 1e-16))
    {

      if (verbose > 0)
      {
        printf(" Vanishing diagonal term \n");
        printf(" The local problem cannot be solved \n");
      }

      *info = 2;
      free(diagA);
      free(diagB);
      free(wOld);
      *info = 1;
      return;
    }
    else
    {
      diagA[i] = omega / A[i * n + i];

    }
  }
  for (i = 0 ; i < m ; ++i)
  {
    if ((fabs(B[i * m + i]) < 1e-16))
    {

      if (verbose > 0)
      {
        printf(" Vanishing diagonal term \n");
        printf(" The local problem cannot be solved \n");
      }

      *info = 2;
      free(diagA);
      free(diagB);
      free(wOld);

      return;
    }
    else
    {
      diagB[i] = omega / B[i * m + i];

    }
  }
  /*start iterations*/

  iter = 0;
  err  = 1.;

  incx = 1;
  incy = 1;
  incAx = n;
  incAy = 1;
  incBx = m;
  incBy = 1;

  DCOPY(m , b , incx , w , incy);       //  q --> w
  if (n >= 1)
  {
    mlcp_compute_error(nn, mm,  A , B , C , D , a , b, u, v, verbose, w,  &err);
  }
  else
  {
    lcp_compute_error(mm,   B , b, u, verbose, w ,  &err);
  }

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    incx = 1;
    incy = 1;

    DCOPY(n , w , incx , wOld , incy);   //  w --> wOld
    DCOPY(n , b , incx , w , incy);    //  b --> w

    for (i = 0 ; i < n ; ++i)
    {
      u[i] = 0.0;
      //uiprev = u[i];
      //zi = -( q[i] + DDOT( n , &vec[i] , incx , z , incy ))*diag[i];
      u[i] =  - (a[i] + DDOT(n , &A[i] , incAx , u , incy)   + DDOT(m , &C[i] , incAx , v , incy)) * diagA[i];
    }

    for (i = 0 ; i < m ; ++i)
    {
      //prevvi = v[i];
      v[i] = 0.0;
      //zi = -( q[i] + DDOT( n , &vec[i] , incx , z , incy ))*diag[i];
      vi = -(b[i] + DDOT(n , &D[i] , incBx , u , incy)   + DDOT(m , &B[i] , incBx , v , incy)) * diagB[i];

      if (vi < 0) v[i] = 0.0;
      else v[i] = vi;
    }



    /* **** Criterium convergence compliant with filter_result_MLCP **** */

    //mlcp_compute_error(n,vec,q,z,verbose,w, &err);

    if (n >= 1)
    {
      mlcp_compute_error(nn, mm,  A , B , C , D , a , b, u, v, verbose, w,  &err);
    }
    else
    {
      lcp_compute_error(mm,   B , b, u, verbose, w ,  &err);
    }
    //err = err ;



    if (verbose == 2)
    {
      printf(" # i%d -- %g : ", iter, err);
      for (i = 0 ; i < n ; ++i) printf(" %g", u[i]);
      for (i = 0 ; i < m ; ++i) printf(" %g", v[i]);
      for (i = 0 ; i < m ; ++i) printf(" %g", w[i]);
      printf("\n");
    }

    /* **** ********************* **** */

  }

  iparamMLCP[2] = iter;
  dparamMLCP[2] = err;

  if (err > tol)
  {
    printf("Siconos/Numerics: mlcp_psor: No convergence of PGS after %d iterations\n" , iter);
    printf("Siconos/Numerics: mlcp_psor: The residue is : %g \n", err);
    *info = 1;
  }
  else
  {
    if (verbose > 0)
    {
      printf("Siconos/Numerics: mlcp_psor: Convergence of PGS after %d iterations\n" , iter);
      printf("Siconos/Numerics: mlcp_psor: The residue is : %g \n", err);
    }
    *info = 0;
  }

  free(wOld);
  free(diagA);
  free(diagB);

  return;
}
