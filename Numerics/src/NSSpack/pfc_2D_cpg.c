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
/*!\file pfc_2D_cpg.c

  This subroutine allows the primal resolution of contact problems with friction in the 2D case (PFC_2D).

   Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
w - M z = q \\
0 \le z_n \perp w_n \ge 0\\
-w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
\end{array}
\right.
\f$

 here M is an (nn \f$\times\f$nn)-matrix, q an nn-dimensional vector, z an nn-dimensional  vector and w an nn-dimensional vector.

*/
/*!\fn void pfc_2D_cpg( int *nn , double *vec , double *b , double *x , double *rout , int *info,
     int *iparamPFC , double *dparamPFC )


   cfp_cpg is a specific cpg (conjugated projected gradient) solver for primal contact problems with friction in the 2D case.

   \param vec         On enter a (nn\f$\times\f$nn)-vector of doubles containing the components of the double matrix with a fortran90 allocation ( M ).
   \param b           On enter a nn-vector of doubles containing the components of the double vector ( q ).
   \param nn          On enter an integer, the dimension of the second member.


   \param iparamPFC   On enter/return a vector of integers:\n
                       - iparamPFC[0] = on enter, the maximum number of iterations allowed,
                       - iparamPFC[1] = on enter, the parameter which represents the output log identifiant,\n
                             0 - no output\n
           >0 - active screen output\n
           - iparamPFC[2] =  on return, the number of iterations performed by the algorithm.

  \param dparamPFC    On enter/return a vector of doubles:\n
                       - dparamPFC[0] = on enter, a positive double which represents the friction coefficient,
                       - dparamPFC[1] = on enter, a positive double which represents the tolerance required,
                       - dparamPFC[2] = on return, a positive double which represents the residu.



   \param x           On return a nn-vector of doubles, the solution of the problem ( z ).
   \param rout        On return a nn-vector of doubles, the solution of the problem ( w ).
   \param info        On return an integer, the termination reason:\n
                       0 = convergence,\n
           1 = no convergence,\n
           2 = Operation of alpha no conform.\n

   \author Nineb Sheherazade.

 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"

void pfc_2D_projc(double *, int *, int *  , double *, double *, double *, int *) ;
void pfc_2D_projf(int *, int *, double *, double *, double *);



void pfc_2D_cpg(int *nn , double *vec , double *b , double *x , double *rout , int *info,
                int *iparamPFC , double *dparamPFC)
{


  int       n = *nn, maxit, ispeak;
  int       nc = n / 2, i, iter;
  integer incx = 1, incy = 1;
  int       *stat, *statusi, it_end;


  double    mu, eps = 1.e-12;
  double    pAp, alpha, beta, wAp, rp, normr, tol;
  double    alphaf, betaf, den, num, res;

  double    *p, *fric, *r;
  double    *fric1, *v, *w, *Ap, *xi, *z;


  char      notrans = 'N';



  maxit         = iparamPFC[0];
  ispeak        = iparamPFC[1];

  mu            = dparamPFC[0];
  tol           = dparamPFC[1];

  iparamPFC[2]  = 0;

  dparamPFC[2]  = 0.0;


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
    x[i]     = 0.0;
    xi[i]    = 0.0;
    r[i]     = 0.0;
    v[i]     = 0.0;
    p[i]     = 0.0;
    w[i]     = 0.0;
    Ap[i]    = 0.0;
    z[i]     = 0.0;
    fric1[i] = 1.0;
    fric[i]  = mu * fric1[i];

    if (i < nc)
    {
      stat[i]    = 0;
      statusi[i] = 0;

    }

  }



  dcopy_((integer *)&n, b, &incx, r, &incy);

  alphaf = -1.;
  betaf  = -1.;

  dgemv_(&notrans, (integer *)&n, (integer *)&n, &alphaf, vec, (integer *)&n, x, &incx, &betaf, r, &incy);





  /*             Check for initial status             */


  for (i = 0; i < nc; i++)
  {
    mu = fric[i];
    if (x[2 * i] <= eps)
    {
      /*       No contact            */
      stat[i] = 0;
    }
    else if (x[2 * i + 1] <=  -mu * x[2 * i])
    {
      /*     Slide backward         */
      stat[i] = 1;
    }
    else if (x[2 * i + 1] >=  mu * x[2 * i])
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


    dcopy_((integer *)&n, r, &incx, v, &incy);

    if (iter == 0)
    {
      dcopy_((integer *)&n, r, &incx, w, &incy);

      dcopy_((integer *)&n, w, &incx, p, &incy);
    }

    alphaf = 1.0;
    betaf  = 0.0;
    dgemv_(&notrans, (integer *)&n, (integer *)&n, &alphaf, vec, (integer *)&n, p, &incx, &betaf, Ap, &incy);

    pAp    = ddot_((integer *)&n, p, &incx, Ap, &incy);

    /*}
      else
    {
    alphaf = 1.0;
    betaf  = 0.0;
    dgemv_( &notrans, (integer *)&n, (integer *)&n, &alphaf, vec, (integer *)&n, p, &incx, &betaf, Ap, &incy );

    pAp    = ddot_( (integer *)&n, p, &incx, Ap, &incy );*/

    if (pAp == 0)
    {
      if (ispeak > 0)
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

    rp     = ddot_((integer *)&n, r, &incx, p, &incy);

    alpha  = rp / pAp;

    dcopy_((integer *)&n, x, &incx, xi, &incy);

    alphaf = alpha;
    daxpy_((integer *)&n, &alphaf, p, &incx, xi, &incy);

    pfc_2D_projc(xi, &n, statusi, p, fric, x, stat);


    /*         r(:)=b(:)-matmul(A,x)          */

    dcopy_((integer *)&n, b, &incx, r, &incy);

    alphaf = -1.;
    betaf  = -1.;
    dgemv_(&notrans, (integer *)&n, (integer *)&n, &alphaf, vec, (integer *)&n, x, &incx, &betaf, r, &incy);

    pfc_2D_projf(statusi, &n, r, fric, w);

    pfc_2D_projf(statusi, &n, p, fric, z);


    wAp    = ddot_((integer *)&n, w, &incx, Ap, &incy);

    beta   = - wAp / pAp;

    dcopy_((integer *)&n, w, &incx, p, &incy);

    alphaf  = beta;
    daxpy_((integer *)&n, &alphaf, z, &incx, p, &incy);


    /*    alphaf  = 1.;
    betaf   = 0.;
    dgemv_( &notrans, (integer *)&n, (integer *)&n, &alphaf, vec , (integer *)&n, p, &incx, &betaf, Ap, &incy );

    pAp     = ddot_( (integer *)&n, p, &incx, Ap, &incy );*/

    dcopy_((integer *)&n, r, &incx, xi, &incy);

    alphaf  = -1.;
    daxpy_((integer *)&n, &alphaf, v, &incx, xi, &incy);

    num     = ddot_((integer *)&n, xi, &incx, xi, &incy);

    den     = ddot_((integer *)&n, v, &incx, v, &incy);

    normr   = sqrt(num / den);

    it_end  = iter;
    res     = normr;


    iparamPFC[2] = it_end;
    dparamPFC[2] = res;


    iter = iter + 1;

  }




  if (normr < tol)
  {

    if (ispeak > 0)
      printf("convergence after %d iterations with a residual %g\n", iter - 1, normr);

    *info = 0;


  }
  else
  {
    if (ispeak > 0)
      printf("no convergence after %d iterations with a residual %g\n", iter - 1, normr);

    *info = 1;
  }


  alpha = -1.;
  dscal_((integer *)&n , &alpha , r , &incx);

  dcopy_((integer *)&n, r, &incx, rout, &incy);



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
