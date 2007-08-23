/* Siconos-Numerics version 2.0.1, Copyright INRIA 2005-2007.
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
/*!\file pfc_3D_nlgsnewton.c
 *
 * This subroutine allows the primal resolution of contact problems with friction.\n
 *
 *
 *            | 1/rn*[wn - pos_part(wn - rn*zn)] |
 *       G =  |                                  | ;
 *            |   1/rt*[wT - proj(wT - rt*zT)]   |
 *
 * We try to solve G = 0 with Newton method.
 *
 * We use  z^{k+1} = z^k + alpha*[ -inv(JacG(z^k))*G(z^k)];
 *
 *
 * here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an  * n-dimensional vector.
 *
 * \fn  pfc_3D_nlgsnewton( int *nn , double *vec , double *q , double *z , double *w ,
 *                         int *info\n, int *iparamLCP , double *dparamLCP )
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
 * \author Houari Khenous.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"


//------ get matrix ------//
double *matrix_get(int m, int n)
{
  return (double *) malloc(sizeof(double) * m * n);
}


//------ copy one matrix into another------//
void matrix_cpy(int m, int n, double *A, double *B)
{
  int i;
  for (i = 0 ; i < m * n ; i++)
  {
    B[ i ] = A[ i ];
  }
}

//_/_/ multiplication mxn x nxl _/_//
void matrix_mult(int m, int n, int l, double *a , double *b, double *c)
{
  int i, j, k;
  double s;

  double *A = matrix_get(m, n);
  double *B = matrix_get(n, l);

  matrix_cpy(m, n, a, A);
  matrix_cpy(n, l, b, B);

  for (i = 0 ; i < m ; i++)
  {
    for (j = 0 ; j < l ; j++)
    {
      s = 0.0;
      for (k = 0 ; k < n ; k++)
        s += A[ i * n + k ] * B[ k * l + j ];
      c[ i * l + j ] = s;
    }
  }
  free(A);
  free(B);
}

//------ Addition of matrices A and B ------//
void matrix_add(int n, double *A , double *B, double *C)
{
  // n - number of ellements in A
  int i;
  for (i = 0 ; i < n ; i++)
  {
    C[i] =  A[i] + B[i];
  }
}

/* Compute function G */
void G_f(int m, double *G, double *Z , double *x , double *y , double rn, double rt, double coef)
{

  double zn , zt, zs, num, mrn, coef2, a2, b2, ab;
  coef2 = coef * coef;

  /* Projection on [0, +infty[ and on D(0, mu*zn) */
  zn = x[0] - rn * y[0];
  if (zn > 0)
  {
    G[0] = y[0];
  }
  else
  {
    G[0] = x[0] / rn;
  }
  zt = x[1] - rt * y[1];
  zs = x[2] - rt * y[2];

  a2 = zt * zt;
  b2 = zs * zs;
  ab = zt * zs;
  mrn = zt * zt + zs * zs;

  // if the radius is negative or the vector is null projection on the disk = 0
  if (coef = 0. || x[0] <= 0. || mrn < 1e-16)
  {
    G[1] = x[1] / rt;
    G[2] = x[2] / rt;
  }

  // if the radius is positive and the vector is non null, we compute projection on the disk
  else
  {
    if (mrn <= coef2 * x[0]*x[0])
    {
      G[1] = y[1];
      G[2] = y[2];
    }
    else
    {
      num  = coef / sqrt(mrn);
      G[1] = (x[1] - zt * x[0] * num) / rt;
      G[2] = (x[2] - zs * x[0] * num) / rt;
    }
  }
}

/* Compute Jacobian of function G */
void JacG_f(int m, double *A, double *B , double *x , double *y , double rn, double rt, double coef)
{
  double zn , zt, zs, mrn, mrn3, coef2, num, a2, b2, ab;

  coef2 = coef * coef;

  /* Projection on [0, +infty[ and on D(0, mu*zn) */
  zn = x[0] - rn * y[0];
  if (zn > 0)
  {
    A[0 * m + 0] = 1.;
  }
  else
  {
    B[0 * m + 0] = 1. / rn;
  }
  zt = x[1] - rt * y[1];
  zs = x[2] - rt * y[2];

  a2 = zt * zt;
  b2 = zs * zs;
  ab = zt * zs;
  mrn = zt * zt + zs * zs;

  // if the radius is negative or the vector is null projection on the disk = 0
  if (coef = 0. || x[0] < 0. || mrn < 1e-16)
    B[1 * m + 1] = B[2 * m + 2] = 1. / rt;
  // if the radius is positive and the vector is non null, we compute projection on the disk
  else
  {
    if (mrn <= coef2 * x[0]*x[0])
      A[1 * m + 1] = A[2 * m + 2] = 1.;
    else
    {
      num  = coef / sqrt(mrn);
      B[0 * m + 1] = - num * zt / rt;
      B[0 * m + 2] = - num * zs / rt;
      mrn3 = sqrt(mrn) * sqrt(mrn) * sqrt(mrn);
      A[1 * m + 1] =  coef * x[0] * b2 / mrn3;
      A[2 * m + 1] = -coef * x[0] * ab / mrn3;
      A[1 * m + 2] = -coef * x[0] * ab / mrn3;
      A[2 * m + 2] =  coef * x[0] * a2 / mrn3;
      B[1 * m + 1] = (1. - rt * coef * x[0] * b2 / mrn3) / rt;
      B[2 * m + 1] =  coef * x[0] * ab / mrn3;
      B[1 * m + 2] =  coef * x[0] * ab / mrn3;
      B[2 * m + 2] = (1. - rt * coef * x[0] * a2 / mrn3) / rt;
    }
  }
}

void pfc_3D_nlgsnewton(int *nn , double *vec , double *q , double *z , double *w , int *info, int *iparamLCP , double *dparamLCP)
{

  FILE *f101;

  int n, in, it, is, ispeak, itermax, nc, i, j, niter, iter, Gsize, mm;
  int nrhs = 1, infoDGESV;
  double err, nerr, nerr1, nerr2, tol, mu, mu2, an, at;
  double qs, a1, alpha, den, num; //,beta, det;
  integer incx, incy;
  double *ww, *www, *wwww, *G, *JacG, *A, *AA, *B, *zz, *W, *zzz, *zzzz, *C;
  int *ipiv;
  char NOTRANS = 'N';

  ispeak = 1;
  nc     = *nn;
  incx   = 1;
  incy   = 1;
  Gsize  = 3;
  n      = 3 * nc;
  mm = Gsize * Gsize;

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

  /* Allocation */

  JacG = (double*)malloc(Gsize * Gsize * sizeof(double));
  AA   = (double*)malloc(Gsize * Gsize * sizeof(double));
  A    = (double*)malloc(Gsize * Gsize * sizeof(double));
  B    = (double*)malloc(Gsize * Gsize * sizeof(double));
  C    = (double*)malloc(Gsize * Gsize * sizeof(double));

  G    = (double*)malloc(Gsize * sizeof(double));
  ww   = (double*)malloc(Gsize * sizeof(double));
  www  = (double*)malloc(Gsize * sizeof(double));
  wwww  = (double*)malloc(Gsize * sizeof(double));
  zz   = (double*)malloc(Gsize * sizeof(double));
  zzz  = (double*)malloc(Gsize * sizeof(double));
  zzzz = (double*)malloc(Gsize * sizeof(double));

  W    = (double*)malloc(n * sizeof(double));

  ipiv = (int *)malloc(Gsize * sizeof(int));

  /* Check for non trivial case */

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

    free(A);
    free(AA);
    free(G);
    free(B);
    free(C);
    free(zzz);
    free(JacG);
    free(ipiv);
    free(ww);
    free(W);
    free(www);
    free(wwww);
    free(zz);
    free(zzzz);

    *info = 0;
    return;
  }

  /* Intialization of w*/
  for (i = 0 ; i < n ; ++i)
  {
    W[i] = w[i] = 0;
  }

  /* Intialization of G, JacG, ww and www */
  for (i = 0 ; i < Gsize ; ++i)
  {
    G[i] = ww[i] = www[i] = zz[i] = zzz[i] = 0.;
    for (j = 0 ; j < Gsize ; ++j)
      JacG[j * Gsize + i] = AA[j * Gsize + i] = A[j * Gsize + i] = B[j * Gsize + i] = C[j * Gsize + i] = 0.;
  }

  /*start NlgsNewton iterations*/
  iter = 0;
  err  = 1.;

  while ((iter < itermax) && (err > tol))
  {

    ++iter;
    /*   printf("---------------- Iteration %i -----------------\n",iter); */

    incx = 1;
    incy = 1;

    dcopy_((integer *)&n , w , &incx , W , &incy);
    dcopy_((integer *)&n , q , &incx ,  w , &incy);

    for (i = 0 ; i < nc ; ++i)
    {

      /*     printf( " ---------------le point de contact %i-----------------\n" , i ); */

      in = 3 * i;
      it = 3 * i + 1;
      is = 3 * i + 2;

      /* rn and rt */
      an = 1. / vec[(in) * n + in];
      at = 1;
      /* alpha = vec[(it)*n+it] + vec[(is)*n+is]; */
      /*       det = vec[(it)*n+it]*vec[(is)*n+is] - vec[(it)*n+is]*vec[(is)*n+it]; */
      /*       beta = alpha*alpha - 4*det; */
      /*       if(beta > 0.) */
      /*  beta = sqrt(beta); */
      /*       else */
      /*  beta = 0.; */

      /*       at = 2*(alpha - beta)/((alpha + beta)*(alpha + beta)); */
      //an = at = 1.;

      C[0 * Gsize + 0] = vec[(in) * n + in];
      C[0 * Gsize + 1] = vec[(in) * n + it];
      C[0 * Gsize + 2] = vec[(in) * n + is];
      C[1 * Gsize + 0] = vec[(it) * n + in];
      C[1 * Gsize + 1] = vec[(it) * n + it];
      C[1 * Gsize + 2] = vec[(it) * n + is];
      C[2 * Gsize + 0] = vec[(is) * n + in];
      C[2 * Gsize + 1] = vec[(is) * n + it];
      C[2 * Gsize + 2] = vec[(is) * n + is];

      /*copy the unknown za  */
      zz[0] = z[in];
      zz[1] = z[it];
      zz[2] = z[is];


      incx = n;
      incy = 1;

      /*q + \sum_{b<a} Wab zb^{k+1} + \sum_{b>a} Wab zb^k */
      z[in] = 0.0;
      z[it] = 0.0;
      z[is] = 0.0;

      /*printf("la vitesse et la force de contact en %i est w[%i] = %8.5e et z[%i] = %8.5e\n",i,i,w[i],i, z[i]); */
      zzz[0] = q[in] + ddot_((integer *)&n , &vec[in] , &incx , z , &incy);
      zzz[1] = q[it] + ddot_((integer *)&n , &vec[it] , &incx , z , &incy);
      zzz[2] = q[is] + ddot_((integer *)&n , &vec[is] , &incx , z , &incy);

      if (zzz[0] > 0)
      {
        zz[0] = 0.;
        zz[1] = 0.;
        zz[2] = 0.;
      }
      else
      {


        niter = 0;
        nerr  = 1.;

        while ((niter < itermax) && (nerr > tol))
        {
          ++niter;

          ww[0] = vec[(in) * n + in] * zz[0] + vec[(it) * n + in] * zz[1] + vec[(is) * n + in] * zz[2] + zzz[0];
          ww[1] = vec[(in) * n + it] * zz[0] + vec[(it) * n + it] * zz[1] + vec[(is) * n + it] * zz[2] + zzz[1];
          ww[2] = vec[(in) * n + is] * zz[0] + vec[(it) * n + is] * zz[1] + vec[(is) * n + is] * zz[2] + zzz[2];

          G_f(Gsize , G , C , zz , ww , an , at , mu);
          JacG_f(Gsize , A , B , zz , ww , an , at , mu);

          nerr1 = dnrm2_((integer *)&Gsize, G , &incx);
          /*  printf("Iteration Newton %i Erreur = %14.7e\n",niter,nerr1); */

          matrix_mult(Gsize, Gsize, Gsize, A, C, C);
          matrix_add(mm, C, B, JacG);

          /* **** Criterium convergence **** */
          incx =  1;
          incy =  1;
          a1   = -1.0;

          /* compute the direction www s.t X^{k+1} = X^{k} + www, where X^{k} = za */
          dcopy_((integer *)&Gsize , G , &incx , www , &incy);
          dscal_((integer *)&Gsize , &a1 , www, &incx);
          dcopy_((integer *)&mm , JacG , &incx , AA , &incy);

          F77NAME(dgesv)((integer *)&Gsize, (integer *)&nrhs, AA, (integer *)&Gsize, (integer *)ipiv, www, (integer *)&Gsize, (integer *)&infoDGESV);

          if (infoDGESV)
          {
            if (ispeak > 0)
            {
              printf("Problem in DGESV\n");
            }
            iparamLCP[2] = iter;
            dparamLCP[1] = err;

            free(A);
            free(AA);
            free(G);
            free(B);
            free(C);
            free(zzz);
            free(JacG);
            free(ipiv);
            free(ww);
            free(W);
            free(www);
            free(wwww);
            free(zz);
            free(zzzz);
            *info = 2;
            return;
          }
          /* incx =  1; */
          /*    incy =  1; */
          /* qs   = 1.0; */
          /*    daxpy_((integer *)&Gsize , &qs , www , &incx , zz , &incy ); */

          if (nerr1 < tol)
          {
            incx =  1;
            incy =  1;
            qs = 1.;
            daxpy_((integer *)&Gsize , &qs , www , &incx , zz , &incy);
            nerr = nerr1;
          }
          else
          {
            alpha = 1.;
            while (alpha > 0.05)
            {
              incx = 1.;
              incy = 1.;
              dcopy_((integer *)&Gsize , zz , &incx , zzzz , &incy);
              daxpy_((integer *)&Gsize , &alpha , www , &incx , zzzz , &incy);

              wwww[0] = vec[(in) * n + in] * zzzz[0] + vec[(it) * n + in] * zzzz[1] + vec[(is) * n + in] * zzzz[2] + zzz[0];
              wwww[1] = vec[(in) * n + it] * zzzz[0] + vec[(it) * n + it] * zzzz[1] + vec[(is) * n + it] * zzzz[2] + zzz[1];
              wwww[2] = vec[(in) * n + is] * zzzz[0] + vec[(it) * n + is] * zzzz[1] + vec[(is) * n + is] * zzzz[2] + zzz[2];

              G_f(Gsize , G , C , zzzz , wwww , an , at , mu);
              nerr2 = dnrm2_((integer *)&Gsize, G , &incx);
              /* printf("Iteration %i Erreur = %14.7e\n",iter,nerr2); */
              if (nerr2 < nerr1 * nerr1) break;
              alpha = alpha * 0.5;
            }
            nerr = nerr2;
            /* printf("Iteration Newton %i Erreur = %14.7e\n",niter,nerr); */
            incx = 1.;
            incy = 1.;
            dcopy_((integer *)&Gsize , zzzz , &incx , zz , &incy);
            dcopy_((integer *)&Gsize , wwww , &incx , ww , &incy);
          }
        }
      }
      z[in] = zz[0] ;
      z[it] = zz[1] ;
      z[is] = zz[2] ;
    }

    /* **** Criterium convergence **** */

    incx =  1;
    incy =  1;
    a1 = 1.0;

    dgemv_(&NOTRANS , (integer *)&n , (integer *)&n , &a1 , vec , (integer *)&n , z , &incx , &a1 , w , &incy);

    qs   = -1.0;
    daxpy_((integer *)&n , &qs , w , &incx , W , &incy);
    num = dnrm2_((integer *)&n, W , &incx);
    err = num * den;

    /* printf("Iteration %i Erreur = %14.7e\n",iter,err); */

    /* for( i = 0 ; i < n ; ++i ) fprintf(f101,"%i  %i  %14.7e\n",iter-1,i,z[i]); */

  }

  iparamLCP[2] = iter;
  dparamLCP[2] = err;

  if (ispeak > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of NLGSNEWTON after %i iterations\n" , iter);
      printf(" The residue is : %e \n", err);
      *info = 1;
    }
    else
    {
      printf(" Convergence of NLGSNEWTON after %i iterations\n" , iter);
      printf(" The residue is : %e \n", err);
      *info = 0;
    }
  }
  else
  {
    if (err > tol) *info = 1;
    else *info = 0;
  }

  free(A);
  free(AA);
  free(G);
  free(B);
  free(C);
  free(zzz);
  free(zzzz);
  free(JacG);
  free(ipiv);
  free(ww);
  free(W);
  free(www);
  free(wwww);
  free(zz);

  if (ispeak == 2) fclose(f101);

}
