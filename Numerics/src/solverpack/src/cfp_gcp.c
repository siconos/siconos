#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*!\file cfp_gcp.c

  This subroutine allows the primal resolution of contact problems with friction.

   Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
0 \le z_n \perp w_n \ge 0\\
-w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
\end{array}
\right.
\f$

 here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.

*/


double ddot_(int *, double [], int *, double [], int*);
//   !projc (double [],int [],double [],double [],double [],int []);


/*!\fn int cfp_gcp(double vec[],double *qq,int *nn,double *mu3,int * itermax, double * tol,double xout[],double rout[],int *it_end,double * res,int *info)

   cfp_gcp is a specific gcp (gradient conjugated projected) solver for primal contact problem with friction.

   \param vec On enter a double vector containing the components of the double matrix with a fortran90 allocation.
   \param qq On enter a pointer over doubles containing the components of the double vector.
   \param nn On enter a pointer over integers, the dimension of the second member.
   \param mu3 On enter a pointer over doubles, the friction coefficient.
   \param itermax On enter a pointer over integers, the maximum iterations required.
   \param tol On enter a pointer over doubles, the tolerance required.
   \param it_end On enter a pointer over integers, the number of iterations carried out.
   \param res On return a pointer over doubles, the error value.
   \param xout On return double vector, the solution of the problem.
   \param rout On return double vector, the solution of the problem.
   \param info On return a pointer over integers, the termination reason (0 is successful otherwise 1).
\return int : \todo tell whiwh result is good
   \author Nineb Sheherazade.

 */
int cfp_gcp(double vec[], double *qq, int *nn, double *mu3, int * itermax, double * tol, double xout[], double rout[], int *it_end, double * res, int *info)
{
  int n = *nn, maxit = *itermax;
  double mu = *mu3, eps = 1.e-08;
  int nc = n / 2, nr0, i, j, iter, eof1, eof2, k, jk, ik, ii, iii, incx, incy;
  double pAp, alpha, beta, rp0, wAp, wAp0, Ap0, pAp0, rp, Aikjk, Aij, bi, eta, normr;
  char trans;
  //  real(kind=8),dimension(:,:),allocatable::S,XX,dir
  double A[n][n], alphaf, betaf, temp, den, num;
  int *stat, *statusi, *projc2_out;
  double *resveccgp, *p, *fric, *projf_out, *projc1_out;
  double  *fric1, *v, *w, *Ap, *xi, *z, *y;
  double b[n], x[n], r[n], fa[nc];


  projf_out = (double*)malloc(n * sizeof(double));
  projc1_out = (double*)malloc(n * sizeof(double));
  projc2_out = (int*)malloc(nc * sizeof(int));
  fric = (double*)malloc(n * sizeof(double));
  p = (double*)malloc(n * sizeof(double));
  resveccgp = (double*)malloc((maxit + 1) * sizeof(double));
  v = (double*)malloc(n * sizeof(double));
  w = (double*)malloc(n * sizeof(double));
  Ap = (double*)malloc(n * sizeof(double));
  xi = (double*)malloc(n * sizeof(double));
  z = (double*)malloc(n * sizeof(double));
  y = (double*)malloc(n * sizeof(double));
  fric1 = (double*)malloc(n * sizeof(double));
  stat = (int*)malloc(nc * sizeof(int));
  statusi = (int*)malloc(nc * sizeof(int));


  for (i = 0; i < n; i++)
  {
    x[i] = 0.;
    xi[i] = 0.;
    projc1_out[i] = 0.;
    projf_out[i] = 0.;
    r[i] = 0.;
    v[i] = 0.;
    p[i] = 0.;
    w[i] = 0.;
    y[i] = 0.;
    Ap[i] = 0.;
    z[i] = 0.;
    fric1[i] = 1.;
    fric[i] = mu * fric1[i];
    b[i] = qq[i];
  }

  /*  incx=1;
  incy=1;
  dcopy_(&n,qq,&incx,b,&incy);*/

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A[i][j] = vec[i * n + j];


  for (i = 0; i < nc; i++)
  {
    stat[i] = 0;
    statusi[i] = 0;
    projc2_out[i] = 0.;
  }


  for (i = 0; i < (maxit + 1); i++)
  {
    resveccgp[i] = 0.;
  }


  //  r(:)=b(:)- matmul(A(:,:),x(:))
  incx = 1;
  incy = 1;
  dcopy_(&n, b, &incx, y, &incy);

  incx = 1;
  incy = 1;
  trans = 'T';
  alphaf = -1.;
  betaf = 1.;
  dgemv_(&trans, &n, &n, &alphaf, A, &n, x, &incx, &betaf, y, &incy);

  dcopy_(&n, y, &incx, r, &incy);


  iter = 1;
  //  nr0=sqrt(dot_product(r,r))

  temp = ddot_(&n, r, &incx, r, &incy);
  nr0 = sqrt(temp);
  resveccgp[1] = nr0;

  //  !Check for initial status
  for (i = 0; i < nc; i++)
  {
    mu = fric[i];
    if (x[2 * i] <= eps)
    {
      // !no contact
      stat[i] = 0;
    }
    else if (x[2 * i + 1] <= -mu * x[2 * i])
    {
      // !slide backward
      stat[i] = 1;
    }
    else if (x[2 * i + 1] >= mu * x[2 * i])
    {
      //  !slide forward
      stat[i] = 3;
    }
    else
    {
      // !stick contact
      stat[i] = 2;
    }
  }

  //  !loop over maxit iterations (unless convergence or failure)
  /*for(k=0;k<n;k++){
    printf("qq %g b %g  \n",qq[k],b[k]);
    }*/

  for (ii = 0; ii < maxit; ii++)
  {
    for (i = 0; i < nc; i++)
      statusi[i] = stat[i];

    //    dcopy_(&nc,stat,&incx,statusi,&incy);
    dcopy_(&n, r, &incx, v, &incy);

    if (ii == 0)
    {
      dcopy_(&n, r, &incx, w, &incy);
      dcopy_(&n, w, &incx, p, &incy);

      alphaf = 1.;
      betaf = 0.;
      dgemv_(&trans, &n, &n, &alphaf, A, &n, p, &incx, &betaf, Ap, &incy);
      //      Ap(:)=matmul(A,p)
      pAp = ddot_(&n, p, &incx, Ap, &incy);
      //      pAp=dot_product(p,Ap)
    }
    else
    {

      trans = 'T';
      alphaf = 1.;
      betaf = 0.;
      dgemv_(&trans, &n, &n, &alphaf, A, &n, p, &incx, &betaf, Ap, &incy);
      pAp = ddot_(&n, p, &incx, Ap, &incy);

      /*       for(k=0;k<n;k++){
        printf("k %d p %g \n",k,p[k]);
      }*/

      if (pAp == 0)
      {
        printf("operation non conform alpha at the iteration %d \n", iter);
        break;
      }

      rp = ddot_(&n, r, &incx, p, &incy);
      alpha = rp / pAp;


      dcopy_(&n, x, &incx, y, &incy);
      trans = 'T';
      alphaf = alpha;
      daxpy_(&n, &alphaf, p, &incx, y, &incy);
      dcopy_(&n, y, &incx, xi, &incy);
      //    xi(:)=x(:)+alpha*p(:)


      //   !    call projc(xi,statusi,p,fric)

      /*         for(k=0;k<n;k++){
        printf("alpha %g p %g \n",alpha,p[k]);

      }*/

      projc(xi, &n, statusi, p, fric, projc1_out, projc2_out);

      for (k = 0; k < nc; k++)
        printf("apres projc k %d statusi %d\n", k, statusi[k]);

      incx = 1;
      incy = 1;
      //           ICOPY_(&nc,statusi,&incx,fa,&incy);


      //    projc_out= projc(xi,statusi,p,fric)
      dcopy_(&n, projc1_out, &incx, x, &incy);

      for (i = 0; i < nc; i++)
        stat[i] = projc2_out[i];

      //      dcopy_(&nc,projc2_out,&incx,stat,&incy);
      //    x(:)= projc_out(:,1)
      //  status(:)=projc_out(1:nc,2)

      //   r(:)=b(:)-matmul(A,x)
      dcopy_(&n, b, &incx, y, &incy);
      alphaf = -1.;
      betaf = 1.;
      dgemv_(&trans, &n, &n, &alphaf, A, &n, x, &incx, &betaf, y, &incy);
      dcopy_(&n, y, &incx, r, &incy);


      for (k = 0; k < nc; k++)
        printf("avant projf k %d statusi %d \n", k, statusi[k]);

      // !   call projf(statusi,r,fric)
      projf(statusi, &n, r, fric, projf_out);
      //    dcopy_(&nc,projc2_out,&incx,stat,&incy);//en +

      /* for(k=0;k<n;k++){
        printf("k %d w %g     Ap %g\n",k,w[k],Ap[k]);
      }*/

      //    projf_out= projf(statusi,r,fric)
      dcopy_(&n, projf_out, &incx, w, &incy);
      //  w=projf_out(:)


      projf(statusi, &n, p, fric, projf_out);
      // z=projf(statusi,p,fric)
      dcopy_(&n, projf_out, &incx, z, &incy);

      ///////////

      wAp = ddot_(&n, w, &incx, Ap, &incy);
      beta = - wAp / pAp;

      dcopy_(&n, w, &incx, y, &incy);
      alphaf = beta;
      betaf = 1.;
      daxpy_(&n, &alphaf, z, &incx, y, &incy);


      dcopy_(&n, y, &incx, p, &incy);
      //    p(:)=w(:)+beta*z(:)


      // for(k=0;k<n;k++){
      //  printf("k %d w %g z %g\n",k,w[k],z[k]);
      //}


      alphaf = 1.;
      betaf = 0.;
      dgemv_(&trans, &n, &n, &alphaf, A, &n, p, &incx, &betaf, y, &incy);
      dcopy_(&n, y, &incx, Ap, &incy);
      /*
      for(k=0;k<n;k++){
        printf("k %d p %g     Ap %g\n",k,p[k],Ap[k]);

      }*/

      pAp = ddot_(&n, p, &incx, Ap, &incy);

      //    normr=sqrt(dot_product(r(:)-v(:),r(:)-v(:)))/sqrt(dot_product(v,v))
      /*            for(k=0;k<n;k++){
      printf("r %g v %g x %g\n",r[k],v[k],x[k]);
      }*/

      dcopy_(&n, r, &incx, y, &incy);
      trans = 'T';
      alphaf = -1.;
      daxpy_(&n, &alphaf, v, &incx, y, &incy);

      num = ddot_(&n, y, &incx, y, &incy);
      den = ddot_(&n, v, &incx, v, &incy);

      normr = sqrt(num / den);
      resveccgp[i + 1] = normr;

      //    S(:,i)=status(:)
      //  XX(:,iter)=x
      //  dir(:,i)=p

      if (normr < *tol)
      {
        iter = ii;
        printf("convergence after %d iterations with a residual %g\n", iter, normr);
        *info = 0;
        break;
      }
      iter = iter + 1;
      //  printf("iter %d\n",iter);
    }
  }



  /* for(k=0;k<n;k++){
        printf("k %d p %g     Ap %g\n",k,p[k],Ap[k]);

  }*/


  if (normr < *tol)
  {
    *info = 0;
  }
  else
  {
    *info = 1;
  }


  if (normr > *tol)
  {
    printf("no convergence after %d iterations with a residual %g\n", iter, normr);
  }

  it_end = &iter;
  res = &normr;
  dcopy_(&n, x, &incx, xout, &incy);


  //  xout=x(:)
  for (i = 0; i < n; i++)
    rout[i] = -r[i];

  free(projf_out);
  free(projc1_out);
  free(projc2_out);
  free(fric);
  free(p);
  free(resveccgp);
  free(v);
  free(w);
  free(Ap);
  free(xi);
  free(z);
  free(y);
  free(fric1);
  free(stat);
  free(statusi);
  return *info;
}
