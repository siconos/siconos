/*!\file pfc_2D_cgp.c

  This subroutine allows the primal resolution of contact problems with friction.

   Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
M z + q =  w\\
0 \le z_n \perp w_n \ge 0\\
-w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
\end{array}
\right.
\f$

 here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.

*/
/*!\fn void pfc_2D_cpg( int *nn , double *vec , double *b , double *x , double *rout , int *info,
     int *iparamPFC , double *dparamPFC )


   cfp_gcp is a specific gcp (gradient conjugated projected) solver for primal contact problem with friction.

   \param vec         On enter a double vector containing the components of the double matrix with a fortran90 allocation.
   \param b           On enter a pointer over doubles containing the components of the double vector ( q ).
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



   \param xout        On return double vector, the solution of the problem ( z ).
   \param rout        On return double vector, the solution of the problem ( w ).
   \param info        On return a pointer over integers, the termination reason
                       0 = convergenec
           1 = no convergence
           2 = Operation of alpha no conform.

   \author Nineb Sheherazade.

 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"

void pfc_2D_cpg(int *nn , double *vec , double *b , double *x , double *rout , int *info,
                int *iparamPFC , double *dparamPFC)
{


  int       n = *nn, maxit, ispeak;
  int       nc = n / 2, i, iter, incx = 1, incy = 1;
  int       *stat, *statusi, it_end;


  double    mu, eps = 1.e-16;
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



  dcopy_(&n, b, &incx, r, &incy);

  alphaf = -1.;
  betaf  = -1.;

  dgemv_(&notrans, &n, &n, &alphaf, vec, &n, x, &incx, &betaf, r, &incy);





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


    dcopy_(&n, r, &incx, v, &incy);

    if (iter == 0)
    {
      dcopy_(&n, r, &incx, w, &incy);

      dcopy_(&n, w, &incx, p, &incy);

      alphaf = 1.0;
      betaf  = 0.0;
      dgemv_(&notrans, &n, &n, &alphaf, vec, &n, p, &incx, &betaf, Ap, &incy);

      pAp    = ddot_(&n, p, &incx, Ap, &incy);
    }
    else
    {
      alphaf = 1.0;
      betaf  = 0.0;
      dgemv_(&notrans, &n, &n, &alphaf, vec, &n, p, &incx, &betaf, Ap, &incy);

      pAp    = ddot_(&n, p, &incx, Ap, &incy);

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

    }

    rp     = ddot_(&n, r, &incx, p, &incy);

    alpha  = rp / pAp;

    dcopy_(&n, x, &incx, xi, &incy);

    alphaf = alpha;
    daxpy_(&n, &alphaf, p, &incx, xi, &incy);

    pfc_2D_projc(xi, &n, statusi, p, fric, x, stat);


    /*         r(:)=b(:)-matmul(A,x)          */

    dcopy_(&n, b, &incx, r, &incy);

    alphaf = -1.;
    betaf  = -1.;
    dgemv_(&notrans, &n, &n, &alphaf, vec, &n, x, &incx, &betaf, r, &incy);

    pfc_2D_projf(statusi, &n, r, fric, w);

    pfc_2D_projf(statusi, &n, p, fric, z);


    wAp    = ddot_(&n, w, &incx, Ap, &incy);

    beta   = - wAp / pAp;

    dcopy_(&n, w, &incx, p, &incy);

    alphaf  = beta;
    daxpy_(&n, &alphaf, z, &incx, p, &incy);


    alphaf  = 1.;
    betaf   = 0.;
    dgemv_(&notrans, &n, &n, &alphaf, vec , &n, p, &incx, &betaf, Ap, &incy);

    pAp     = ddot_(&n, p, &incx, Ap, &incy);

    dcopy_(&n, r, &incx, xi, &incy);

    alphaf  = -1.;
    daxpy_(&n, &alphaf, v, &incx, xi, &incy);

    num     = ddot_(&n, xi, &incx, xi, &incy);

    den     = ddot_(&n, v, &incx, v, &incy);

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
  dscal_(&n , &alpha , r , &incx);

  dcopy_(&n, r, &incx, rout, &incy);



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
