/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"
#include <time.h>
#include "blaslapack.h"
#include "pfc_3D_Alart_Curnier.h"



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

//_/_/ transposed matrix _/_//
void matrix_trans(int m, int n, double *a , double *b)
{
  int i, j;

  double *A = matrix_get(m, n);

  matrix_cpy(m, n, a, A);

  for (i = 0 ; i < m ; i++)
  {
    for (j = 0 ; j < n ; j++)
    {
      b[ j * m + i ] = A[ i * n + j ];
    }
  }

  free(A);
}

//------ formes a unit matrix A with dimension n------//
void matrix_I(int n, double *l)
{

  int i, j ;

  for (i = 0 ; i < n ; i++)
  {
    for (j = 0 ; j < n ; j++)
    {
      l[ i * n + j ] = 0;
    }
  }
  for (i = 0 ; i < n ; i++)
  {
    l[ i * n + i ] = 1;
  }

}

//_/_/   Inverse Matrix   _/_//
void matrix_inv(int n, double *a, double *b)
{
  int i, j, k;
  double p, q;

  double *A = matrix_get(n, n);
  matrix_cpy(n, n, a, A);
  matrix_I(n, b);

  for (k = 0 ; k < n ; ++k)
  {
    p = A[k * n + k];

    for (j = 0 ; j < n ; ++j)
    {
      b[k * n + j] /= p;
      A[k * n + j] /= p;
    }

    for (i = 0 ; i < n ; ++i)
    {
      if (i != k)
      {
        q = A[ i * n + k ];

        for (j = 0 ; j < n ; ++j)
        {
          A[ i * n + j ] -= q * A[ k * n + j ];
          b[ i * n + j ] -= q * b[ k * n + j ];
        }
      }
    }
  }
  free(A);
}



/* Compute function G */
void G_f(int m, double *G, double *x , double *y , double *C, double *b , double rn, double rt, double coef)
{

  double zn , zt, zs, num, mrn, coef2;
  coef2 = coef * coef;

  y[0] = C[0 * m + 0] * x[0] + C[1 * m + 0] * x[1] + C[2 * m + 0] * x[2] + b[0];
  y[1] = C[0 * m + 1] * x[0] + C[1 * m + 1] * x[1] + C[2 * m + 1] * x[2] + b[1];
  y[2] = C[0 * m + 2] * x[0] + C[1 * m + 2] * x[1] + C[2 * m + 2] * x[2] + b[2];

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

  mrn = zt * zt + zs * zs;

  // if the radius is negative or the vector is null projection on the disk = 0
  if (mrn <= 1e-16 || x[0] <= 1e-16)
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
void JacG_f(int m, double *JacG, double *A, double *B , double *x , double *y , double *C, double *b , double rn, double rt, double coef)
{
  double zn , zt, zs, mrn, mrn3, coef2, num, a2, b2, ab, a1;
  int mm, incx, incy;
  mm = m * m;

  coef2 = coef * coef;

  y[0] = C[0 * m + 0] * x[0] + C[1 * m + 0] * x[1] + C[2 * m + 0] * x[2] + b[0];
  y[1] = C[0 * m + 1] * x[0] + C[1 * m + 1] * x[1] + C[2 * m + 1] * x[2] + b[1];
  y[2] = C[0 * m + 2] * x[0] + C[1 * m + 2] * x[1] + C[2 * m + 2] * x[2] + b[2];


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
  if (mrn <= 1e-16 || x[0] <= 1e-16)
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

  incx =  3;
  incy =  1;
  a1   = 1.0;

  DCOPY(mm , B , incy ,  JacG , incy);
  DGEMM(LA_NOTRANS, LA_NOTRANS, m, m, m, a1 , A , m , C , incx , a1 , JacG , incx);
}


void pfc_3D_nlgsnewton(int nc , double *vec , double *q , double *z , double *w , double *mu, int *info,
                       int *iparamLCP , double *dparamLCP)
{


  FILE *f101;
  int n, in, it, is, ispeak, itermax, i, j, niter, iter, Gsize, mm;
  double err, nerr, nerr1, nerr2, tol, an, at;
  double qs, a1, b1, alpha, den, num, beta, det;
  int incx, incy;
  double *ww, *www, *wwww, *G, *JacG, *A, *AA, *B, *zz, *W, *zzz, *zzzz, *C, *AC, *diag;
  int *ipiv;

  clock_t t1, t2;

  t1 = clock();

  ispeak = 0;
  incx   = 1;
  incy   = 1;
  n      = 3 * nc;
  Gsize  = 3;
  mm = Gsize * Gsize;

  /* Recup input */

  itermax = iparamLCP[0];
  ispeak  = iparamLCP[1];

  tol = dparamLCP[0];

  /* Initialize output */

  iparamLCP[2] = 0;
  dparamLCP[1] = 0.0;


  if (ispeak == 2) f101 = fopen("pfc_3D_nlgs.log" , "w+");

  iter = 0;
  err  = 1.;

  /* Allocation */
  JacG   = (double*)malloc(Gsize * Gsize * sizeof(double));
  AA     = (double*)malloc(Gsize * Gsize * sizeof(double));
  A      = (double*)malloc(Gsize * Gsize * sizeof(double));
  B      = (double*)malloc(Gsize * Gsize * sizeof(double));
  C      = (double*)malloc(Gsize * Gsize * sizeof(double));
  AC     = (double*)malloc(Gsize * Gsize * sizeof(double));

  G    = (double*)malloc(Gsize * sizeof(double));
  ww   = (double*)malloc(Gsize * sizeof(double));
  www  = (double*)malloc(Gsize * sizeof(double));
  wwww = (double*)malloc(Gsize * sizeof(double));
  zz   = (double*)malloc(Gsize * sizeof(double));
  zzz  = (double*)malloc(Gsize * sizeof(double));
  zzzz = (double*)malloc(Gsize * sizeof(double));

  W    = (double*)malloc(n * sizeof(double));
  diag = (double*)malloc(n * sizeof(double));

  ipiv = (int *)malloc(Gsize * sizeof(int));

  /* Check for non trivial case */

  qs = DNRM2(n , q , incx);

  if (ispeak > 0) printf("\n ||q||= %g \n" , qs);

  if (qs > 1e-16) den = 1.0 / qs;
  else
  {
    for (i = 0 ; i < n ; ++i)
    {
      w[i] = 0.;
      z[i] = 0.;
    }

    free(diag);
    free(A);
    free(AC);
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

  /* Intialization of w */

  for (i = 0 ; i < n ; ++i)
  {
    W[i] = 0.;
    w[i]  = 0.;
  }

  /* Intialization of G, JacG, ww and www */
  for (i = 0 ; i < Gsize ; ++i)
  {
    G[i] = ww[i] = www[i] = zz[i] = zzz[i] = wwww[i] = zzzz[i] = 0.;
    for (j = 0 ; j < Gsize ; ++j)
      JacG[j * Gsize + i] = AA[j * Gsize + i] = AC[j * Gsize + i] = A[j * Gsize + i] = B[j * Gsize + i] = C[j * Gsize + i] =  0.;
  }


  /* Preparation of the diagonal of the inverse matrix */

  for (i = 0 ; i < nc ; ++i)
  {
    in = 3 * i;
    if (fabs(vec[in * n + in]) < 1e-16)
    {

      if (ispeak > 0)
      {
        printf(" Warning negative diagonal term \n");
        printf(" The local problem can be solved \n");
      }

      *info = 2;
      free(diag);
      free(A);
      free(AC);
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

      return;
    }
    else
    {
      diag[in  ] = 1.0 / vec[(in) * n + in  ];
      diag[in + 1] = 1.0 / vec[(in + 1) * n + in + 1];
      diag[in + 2] = 1.0 / vec[(in + 2) * n + in + 2];
    }
  }


  /*start iterations*/

  iter = 0;
  err  = 1.;

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    incx = 1;
    incy = 1;

    DCOPY(n , w , incx , W , incy);
    DCOPY(n , q , incx ,  w , incy);

    for (i = 0 ; i < nc ; ++i)
    {

      in = 3 * i;
      it = 3 * i + 1;
      is = 3 * i + 2;

      C[0 * Gsize + 0] = vec[(in) * n + in];
      C[0 * Gsize + 1] = vec[(in) * n + it];
      C[0 * Gsize + 2] = vec[(in) * n + is];
      C[1 * Gsize + 0] = vec[(it) * n + in];
      C[1 * Gsize + 1] = vec[(it) * n + it];
      C[1 * Gsize + 2] = vec[(it) * n + is];
      C[2 * Gsize + 0] = vec[(is) * n + in];
      C[2 * Gsize + 1] = vec[(is) * n + it];
      C[2 * Gsize + 2] = vec[(is) * n + is];

      /*  for( ii = 0 ; ii < Gsize ; ++ii ){ */
      /*  for( j = 0 ; j < Gsize ; ++j ){ */
      /*    printf("C[%i,%i] = %14.7e\t",ii,j,C[j*Gsize+ii]); */
      /*  } */
      /*  printf("\n"); */
      /*       } */
      /* rn and rt */
      an = 1. / C[0 * Gsize + 0];
      alpha = C[1 * Gsize + 1] + C[2 * Gsize + 2];
      det = C[1 * Gsize + 1] * C[2 * Gsize + 2] - C[2 * Gsize + 1] + C[1 * Gsize + 2];
      beta = alpha * alpha - 4 * det;
      if (beta > 0.)
        beta = sqrt(beta);
      else
        beta = 0.;

      at = 2 * (alpha - beta) / ((alpha + beta) * (alpha + beta));


      /*copy the unknown za  */
      zz[0] = z[in];
      zz[1] = z[it];
      zz[2] = z[is];

      incx = n;
      incy = 1;

      z[in] = 0.0;
      z[it] = 0.0;
      z[is] = 0.0;

      zzz[0] = q[in] + DDOT(n , &vec[in] , incx , z , incy);
      zzz[1] = q[it] + DDOT(n , &vec[it] , incx , z , incy);
      zzz[2] = q[is] + DDOT(n , &vec[is] , incx , z , incy);

      /*start NlgsNewton iterations*/
      niter = 0;
      nerr  = 1.;

      while ((niter < itermax) && (nerr > tol))
      {
        ++niter;

        G_f(Gsize , G , zz , ww , C , zzz, an , at , mu[i]);
        JacG_f(Gsize , JacG, A , B , zz , ww , C, zzz, an , at , mu[i]);
        incx = 1;
        nerr1 = DNRM2(Gsize, G , incx);

        /***** Criterium convergence *****/
        incx =  1;
        incy =  1;
        a1   = -1.0;

        /* compute the direction www s.t X^{k+1} = X^{k} + www, where X^{k} = za */

        DSCAL(Gsize , a1 , G, incx);

        matrix_inv3(JacG, AA);
        //matrix_inv(Gsize,JacG, AA);

        matrix_mult(Gsize, Gsize, 1, AA, G, www);


        alpha = 1.;
        while (alpha > 0.05)
        {
          incx = 1.;
          incy = 1.;
          DCOPY(Gsize , zz , incx , zzzz , incy);
          DAXPY(Gsize , alpha , www , incx , zzzz , incy);

          G_f(Gsize , G , zzzz , wwww , C, zzz, an , at , mu[i]);
          nerr2 = DNRM2(Gsize, G , incx);

          if (nerr2 < nerr1) break;
          alpha = alpha * 0.5;
        }
        nerr = nerr2;

        /* printf("Iteration Newton %i Erreur = %14.7e\n",niter,nerr); */
        incx = 1.;
        incy = 1.;
        DCOPY(Gsize , zzzz , incx , zz , incy);
        DCOPY(Gsize , wwww , incx , ww , incy);
      }

      z[in] = zz[0] ;
      z[it] = zz[1] ;
      z[is] = zz[2] ;

    }

    /* **** Criterium convergence **** */

    incx =  1;
    incy =  1;

    a1 = 1.0;
    b1 = 1.0;

    DGEMV(LA_NOTRANS , n , n , a1 , vec , n , z , incx , b1 , w , incy);

    qs   = -1.0;
    DAXPY(n , qs , w , incx , W , incy);

    num = DNRM2(n, W , incx);
    err = num * den;

    /*  printf("Iteration %i Erreur = %24.8e\n",iter,err); */

    if (ispeak == 2) for (i = 0 ; i < n ; ++i) fprintf(f101, "%i  %i  %14.7e\n", iter - 1, i, z[i]);

  }

  iparamLCP[2] = iter;
  dparamLCP[1] = err;

  if (ispeak > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of NLGS after %i iterations\n" , iter);
      printf(" The residue is : %e \n", err);
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

  free(diag);
  free(ww);
  free(A);
  free(AA);
  free(AC);
  free(G);
  free(B);
  free(C);
  free(zzz);
  free(zzzz);
  free(JacG);
  free(ipiv);
  free(W);
  free(www);
  free(wwww);
  free(zz);

  if (ispeak == 2) fclose(f101);
  t2 = clock();
  /*  printf("%.4lf seconds of processing\n", (t2-t1)/(double)CLOCKS_PER_SEC); */

}
