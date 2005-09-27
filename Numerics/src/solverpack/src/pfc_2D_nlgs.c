/*!\file pfc_2D_nlgs.c

  This subroutine allows the primal resolution of contact problems with friction.

   Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
M z + q = w\\
0 \le z_n \perp w_n \ge 0\\
-w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
\end{array}
\right.
\f$

 here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.

*/
/*!\fn  void pfc_2D_nlgs( int *nn , double *vec , double *q , double *z , double *w , int *info,
      int *iparamPFC , double *dparamPFC )



   pfc_2D_gsnl is a specific nlgs (Non Linear Gauss Seidel ) solver for primal contact problem with friction.

   \param vec         On enter a double vector containing the components of the double matrix with a fortran90 allocation.
   \param qq          On enter a pointer over doubles containing the components of the double vector.
   \param nn          On enter a pointer over integers, the dimension of the second member.
   \param iparamPFC   On enter/return a vector of integers,

                       _ iparamPFC[0] = on enter, the maximum number of iterations allowed,
                       _ iparamPFC[1] = on enter, the parameter which represents the output log identifiant
                             0 - no output\n
           0 < active screen output\n
           _ iparamPFC[2] =  on return, the number of iterations performed by the algorithm

  \param dparamPFC    On enter/return a vector of doubles,

                       _ dparamPFC[0] = on enter, a positive double which represents the friction coefficient,
                       _ dparamPFC[1] = on enter, a positive double which represents the tolerance required,
                       _ dparamPFC[2] = on return, a positive double which represents the residu.


   \param z           On return double vector, the solution of the problem.
   \param w           On return double vector, the solution of the problem.
   \param info        On return a pointer over integers, the termination reason
                        0 = convergence
            1 = no convergence
            2 = nul term in diagonal of M.



   \author Nineb Sheherazade.



 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"


void pfc_2D_nlgs(int *nn , double *vec , double *q , double *z , double *w , int *info,
                 int *iparamPFC , double *dparamPFC)
{

  int i, j, k, iter, maxit;
  int n = *nn, incx = 1, incy = 1, nc = n / 2;
  int ispeak, it_end;


  double errmax, alpha, beta, mu;
  double *y, res;
  double *fric1, *fric;
  double normr, eps, avn, avt;
  double apn, apt, zn , zt, den1, num1;

  char  notrans = 'N';



  maxit        = iparamPFC[0];
  ispeak       = iparamPFC[1];


  mu           = dparamPFC[0];
  errmax       = dparamPFC[1];


  iparamPFC[2] = 0;
  dparamPFC[2] = 0.0;


  iter         = 0;
  eps          = 1.e-08;



  y       = (double*) malloc(n  * sizeof(double));
  fric1   = (double*) malloc(nc * sizeof(double));
  fric    = (double*) malloc(nc * sizeof(double));




  for (i = 0; i < n; i++)
  {

    z[i]  = 0.0 ;
    w[i]  = 0.0 ;

    if (i < nc)
    {

      fric1[i] = 1.0 ;
      fric[i]  = mu * fric1[i] ;

    }

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

        avn += vec[j * n + 2 * i] * z[j];
        avt += vec[j * n + 2 * i + 1] * z[j];

      }

      for (k = 2 * i + 2; k < n; k++)
      {

        apn = apn + vec[k * n + 2 * i] * z[k];
        apt = apt + vec[n * k + 2 * i + 1] * z[k];

      }


      zn    = -q[2 * i] - avn - apn;
      zt    = -q[2 * i + 1] - avt - apt;


      if (zn > eps)
      {

        if (fabs(vec[n * 2 * i + 2 * i]) < 1e-16)
        {
          if (ispeak > 0)
            printf("\n Warning nul diagonal term of M \n");

          free(fric);
          free(fric1);
          free(y);

          *info = 2;

          return;
        }
        else
        {

          z[2 * i]   = zn / vec[n * 2 * i + 2 * i] ;
          w[2 * i]   = 0.0 ;
        }

        if (fabs(vec[n * (2 * i + 1) + 2 * i + 1]) < 1e-16)
        {
          if (ispeak > 0)
            printf("\n Warning nul diagonal term of M \n");

          free(fric);
          free(fric1);
          free(y);

          *info = 2;

          return;

        }
        else
        {

          z[2 * i + 1] = zt / vec[n * (2 * i + 1) + 2 * i + 1] ;
          w[2 * i + 1] = 0.;

        }



        if (z[2 * i + 1] >  fric[i]*z[2 * i])
        {
          z[2 * i + 1] = fric[i] * z[2 * i];
          w[2 * i + 1] = -zt + vec[n * (2 * i + 1) + 2 * i + 1] * z[2 * i + 1];
        }
        else if (z[2 * i + 1] < -fric[i]*z[2 * i])
        {
          z[2 * i + 1] = -fric[i] * z[2 * i];
          w[2 * i + 1] = -zt + vec[n * (2 * i + 1) + 2 * i + 1] * z[2 * i + 1];
        }
      }
      else
      {
        z[2 * i]   = 0.0 ;
        w[2 * i]   = -zn ;
        z[2 * i + 1] = 0.0 ;
        w[2 * i + 1] = -zt ;

      }
    }


    dcopy_(&n, q, &incx, y, &incy);

    alpha = 1.;
    beta  = 1.;
    dgemv_(&notrans, &n, &n, &alpha, vec, &n, z, &incx, &beta, y, &incy);



    /*          Convergence criterium           */



    alpha = -1.;
    daxpy_(&n, &alpha, w, &incx, y, &incy);


    num1 = ddot_(&n, y, &incx, y, &incy);
    den1 = ddot_(&n, q, &incx, q, &incy);


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
  free(fric1);
  free(fric);


}
