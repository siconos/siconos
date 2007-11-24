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
/*!\file pfc_2D_nlgs.c

  This subroutine allows the primal resolution of contact problems with friction in the 2D case (PFC_2D).

   Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
w - M z = q\\
0 \le z_n \perp w_n \ge 0\\
-w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
\end{array}
\right.
\f$

 here M is an (nn \f$\times\f$nn)-matrix, q an nn-dimensional vector, z an nn-dimensional  vector and w an nn-dimensional vector.

*/

/*!\fn  void pfc_2D_nlgs( int *nn , double *vec , double *q , double *z , double *w , int *info,
      int *iparamPFC , double *dparamPFC )


   pfc_2D_gsnl is a specific nlgs (Non Linear Gauss Seidel ) solver for primal contact problem with friction in 2D case.

   \param vec         On enter a (nn\f$\times\f$nn)-vector of doubles containing the components of the double matrix with a fortran90 allocation.
   \param qq          On enter a nn-vector of doubles containing the components of the second member.
   \param nn          On enter an integer, the dimension of the second member.
   \param iparamPFC   On enter/return a vector of integers:\n
                       - iparamPFC[0] = on enter, the maximum number of iterations allowed,
                       - iparamPFC[1] = on enter, the parameter which represents the output log identifiant:\n
                             0 - no output\n
           >0 - active screen output\n
           - iparamPFC[2] =  on return, the number of iterations performed by the algorithm.\n

  \param dparamPFC    On enter/return a vector of doubles:\n
                       - dparamPFC[0] = on enter, a positive double which represents the friction coefficient,
                       - dparamPFC[1] = on enter, a positive double which represents the tolerance required,
                       - dparamPFC[2] = on return, a positive double which represents the residu.


   \param z           On return a nn-vector of doubles, the solution of the problem.
   \param w           On return a nn-vector of doubles, the solution of the problem.
   \param info        On return an integer, the termination reason:\n
                        0 = convergence,\n
            1 = no convergence,\n
            2 = nul term in denominator.



   \author Nineb Sheherazade.



 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"


void pfc_2D_nlgs(int *nn , double *vec , double *q , double *z , double *w , int *info, int *iparamPFC , double *dparamPFC)
{

  int i, j, k, iter, maxit;
  int n = *nn, nc = n / 2;
  int ispeak, it_end;
  int  incx, incy;

  double errmax, alpha, beta, mu;
  double *y, res;
  double normr, eps, avn, avt, det, gplus, gmoins;
  double apn, apt, zn , zt, den1, num1;

  maxit        = iparamPFC[0];
  ispeak       = iparamPFC[1];


  mu           = dparamPFC[0];
  errmax       = dparamPFC[1];


  iparamPFC[2] = 0;
  dparamPFC[2] = 0.0;


  iter         = 0;
  eps          = 1.e-08;



  y       = (double*) malloc(n  * sizeof(double));





  for (i = 0; i < n; i++)
  {

    z[i]  = 0.0 ;
    w[i]  = 0.0 ;
  }


  normr    =   1.;



  while ((iter < maxit) && (normr > errmax))
  {



    iter = iter + 1 ;


    /*         Loop over contacts                */



    for (i = 0; i < nc; i++)
    {
      avn = 0.;
      avt = 0.;
      apn = 0.;
      apt = 0.;

      for (j = 0; j <= 2 * i - 1; j++)
      {

        avn = avn + vec[j * n + 2 * i] * z[j];
        avt = avt + vec[j * n + 2 * i + 1] * z[j];

      }

      for (k = 2 * i + 2; k < n; k++)
      {

        apn = apn + vec[k * n + 2 * i] * z[k];
        apt = apt + vec[n * k + 2 * i + 1] * z[k];

      }



      zn    = -q[2 * i] - avn - apn;
      zt    = -q[2 * i + 1] - avt - apt;








      if (-zn >= 0.0)
      {


        z[2 * i]   = 0.0; // PN
        w[2 * i]   = -zn; // UN
        z[2 * i + 1] = 0.0; // PT
        w[2 * i + 1] = -zt; // UT


      }
      else
      {

        w[2 * i]   = 0.0;
        w[2 * i + 1] = 0.0;


        det    = vec[2 * i + 2 * i * n] * vec[(2 * i + 1) + (2 * i + 1) * n] - vec[(2 * i) + (2 * i + 1) * n] * vec[(2 * i) + (2 * i + 1) * n];

        if (fabs(det) < 1e-12)
        {

          if (ispeak > 0)
            printf(" Warning denominator nul\n");

          free(y);

          *info = 2;
          return;

        }
        else
        {


          z[2 * i]   = (zn * vec[(2 * i + 1) + n * (2 * i + 1)] - zt * vec[2 * i + (2 * i + 1) * n]) / det;
          z[2 * i + 1] = (-zn * vec[(2 * i) + n * (2 * i + 1)] + zt * vec[2 * i + (2 * i) * n]) / det;

        }

        if ((z[2 * i] >= 0.0) && ((fabs(z[2 * i + 1]) - mu * z[2 * i]) <= 0.0))
        {

          /*  printf("Stick status \n");*/
        }
        else
        {


          w[2 * i]   = 0.0;


          gplus  = vec[2 * i + 2 * i * n] + mu * vec[(2 * i) + (2 * i + 1) * n];


          if (fabs(gplus) < 1e-12)
          {

            if (ispeak > 0)
              printf(" Warning denominator nul\n");

            free(y);

            *info = 2;
            return;

          }
          else
          {

            w[2 * i + 1] = -zt + (zn / gplus) * (vec[2 * i + (2 * i + 1) * n] + mu * vec[(2 * i + 1) + (2 * i + 1) * n]);


            z[2 * i]   = zn / gplus;
            z[2 * i + 1] = mu * z[2 * i];

          }

          if ((z[2 * i] >= 0.0) && (w[2 * i + 1] <= 0.0))
          {

            /*    printf("Slip+ status\n");*/

          }
          else
          {

            w[2 * i]   = 0.0;

            gmoins = vec[2 * i + 2 * i * n] - mu * vec[(2 * i) + (2 * i + 1) * n];


            if (fabs(gmoins) < 1e-12)
            {

              if (ispeak > 0)
                printf(" Warning denominator nul\n");

              free(y);

              *info = 2;
              return;

            }
            else
            {


              w[2 * i + 1] = -zt + (zn / gmoins) * (vec[2 * i + (2 * i + 1) * n] - mu * vec[(2 * i + 1) + (2 * i + 1) * n]);

              z[2 * i]   = zn / gmoins;
              z[2 * i + 1] = -mu * z[2 * i];
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
    DGEMV(LA_NOTRANS, n, n, alpha, vec, n, z, incx, beta, y, incy);



    alpha = -1.;
    DAXPY(n, alpha, w, incx, y, incy);


    num1 = DDOT(n, y, incx, y, incy);
    den1 = DDOT(n, q, incx, q, incy);


    normr = sqrt(num1 / den1);


    it_end = iter;
    res    = normr;

  }


  iparamPFC[2] = it_end;
  dparamPFC[2] = res;



  if (normr > errmax)
  {
    if (ispeak > 0)
      printf("No convergence after %d iterations, the residue is %g\n", iter, normr);

    *info = 1;
  }
  else
  {

    if (ispeak > 0)
      printf("Convergence after %d iterations, the residue is %g \n", iter, normr);

    *info = 0;
  }




  free(y);



}
