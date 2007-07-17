/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2006.
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
 *            |          M*z + q - w             |
 *       G =  | 1/rn*[wn - pos_part(wn - rn*zn)] | ;
 *            |   1/rt*[wT - proj(wT - rt*zT)]   |
 *
 * We try to solve G = 0 with Newton method.
 *
 * We use  (z^{k+1},w^{k+1}) = (z^k,w^k) + alpha*[ -inv(JacG(z^k,w^k))*G(z^k,w^k)];
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


void pfc_3D_nlgsnewton(int *nn , double *vec , double *q , double *z , double *w , int *info, int *iparamLCP , double *dparamLCP)
{


  FILE *f101;

  int n, in, it, is, ispeak, itermax, nc, i, j, iter, mm;
  int nrhs = 1, infoDGESV;
  double err, zn , zt, zs, zn1, mrn, mrn3, num, tol, mu, mu2, rn, rt, a2, b2, ab;
  double qs, a1, alpha, beta, det;
  integer incx, incy, Gsize;
  double *ww, *www, *G, *G0, *JacG, *A, *zz;
  int *ipiv;

  ispeak = 1;
  nc     = *nn;
  incx   = 1;
  incy   = 1;
  Gsize  = 6;
  n      = 3 * nc;
  mm = Gsize * Gsize;

  /* Augmentation parameter for the newton method (coming from resolvent formulation)*/

  rn = 1.;
  rt = 1.;


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

  iter = 0;

  /* Allocation */

  G    = (double*)malloc(Gsize * sizeof(double));
  G0   = (double*)malloc(Gsize * sizeof(double));
  JacG = (double*)malloc(Gsize * Gsize * sizeof(double));
  ww   = (double*)malloc(Gsize * sizeof(double));
  www  = (double*)malloc(Gsize * sizeof(double));
  A    = (double*)malloc(Gsize * Gsize * sizeof(double));
  zz    = (double*)malloc(3 * sizeof(double));

  ipiv = (int *)malloc(Gsize * sizeof(int));

  /* Intialization of w*/

  /* for( i = 0 ; i < n ; ++i ){ */
  /*     w[i] = 0.; */
  /*   } */

  /* Intialization of G, JacG, ww and www */
  for (i = 0 ; i < Gsize ; ++i)
  {
    G[i] = ww[i] = www[i] = 0.;
    for (j = 0 ; j < Gsize ; ++j)
      JacG[j * Gsize + i] = 0.;
  }


  /*  printf( " Nombre de contact %i\n" , nc ); */

  for (i = 0 ; i < nc ; ++i)
  {

    printf(" ---------------le point de contact %i-----------------\n" , i);

    /*  for( j = 0 ; j < Gsize ; ++j ){ */
    /*       printf("la vitesse et la force de contact en %i est (w,z)[%i] = %8.5e\n",i,j,www[j]); */
    /*       } */
    in = 3 * i;
    it = 3 * i + 1;
    is = 3 * i + 2;

    /* rn and rt */

    rn = 1. / vec[(in) * n + in];
    alpha = vec[(it) * n + it] + vec[(is) * n + is];
    det = vec[(it) * n + it] * vec[(is) * n + is] - vec[(it) * n + is] * vec[(is) * n + it];
    beta = alpha * alpha - 4 * det;
    if (beta > 0.)
      beta = sqrt(beta);
    else
      beta = 0.;

    rt = 2 * (alpha - beta) / ((alpha + beta) * (alpha + beta));


    JacG[0 * Gsize + 0] = JacG[1 * Gsize + 1] = JacG[2 * Gsize + 2] = 1.;

    JacG[3 * Gsize + 0] = -vec[(in) * n + in];
    JacG[3 * Gsize + 1] = -vec[(in) * n + it];
    JacG[3 * Gsize + 2] = -vec[(in) * n + is];
    JacG[4 * Gsize + 0] = -vec[(it) * n + in];
    JacG[4 * Gsize + 1] = -vec[(it) * n + it];
    JacG[4 * Gsize + 2] = -vec[(it) * n + is];
    JacG[5 * Gsize + 0] = -vec[(is) * n + in];
    JacG[5 * Gsize + 1] = -vec[(is) * n + it];
    JacG[5 * Gsize + 2] = -vec[(is) * n + is];

    incx = n;
    incy = 1;

    /*  q + \sum_{b<a} Wab zb^{k+1} + \sum_{b>a} Wab zb^k */
    z[in] = 0.0;
    z[it] = 0.0;
    z[is] = 0.0;

    //   printf("la vitesse et la force de contact en %i est w[%i] = %8.5e et z[%i] = %8.5e\n",i,i,w[i],i, z[i]);

    zz[0] = q[in] + ddot_((integer *)&n , &vec[in] , &incx , z , &incy);
    zz[1] = q[it] + ddot_((integer *)&n , &vec[it] , &incx , z , &incy);
    zz[2] = q[is] + ddot_((integer *)&n , &vec[is] , &incx , z , &incy);

    /*start Newton iterations*/

    iter = 0;
    err  = 1.;

    /* copy the unknown (wa,za) in www */
    www[0] = w[in];
    www[1] = w[it];
    www[2] = w[is];
    www[3] = z[in];
    www[4] = z[it];
    www[5] = z[is];

    /*  printf("---------------- Itermax = %i\n",itermax); */

    while ((iter < itermax) && (err > tol))
    {

      ++iter;
      printf("---------------- Iteration %i -----------------\n", iter);


      /* the unknown is X^k = ww = (wa,za) at each node a*/
      incx =  1;
      incy =  1;
      dcopy_(&Gsize , www , &incx , ww , &incy);

      /* Computation of G and its Jacobian matrix JacG */
      /*                               [Da,a   Da,a+1   Da,a+2   ]  */
      /* at each contact point a, Da = [Da+1,a Da+1,a+1 Da+1,a+2 ] is the delassus operator */
      /*
       [Da+2,a Da+2,a+1 Da+2,a+1 ]  */

      /* G = wa - Daa*za - q - \sum_{b<a} Dab zb^{k+1} - \sum_{b>a} Dab zb^k   */

      G[0] = ww[0] - vec[(in) * n + in] * ww[3] - vec[(it) * n + in] * ww[4] - vec[(is) * n + in] * ww[5] - zz[0];
      G[1] = ww[1] - vec[(in) * n + it] * ww[3] - vec[(it) * n + it] * ww[4] - vec[(is) * n + it] * ww[5] - zz[1];
      G[2] = ww[2] - vec[(in) * n + is] * ww[3] - vec[(it) * n + is] * ww[4] - vec[(is) * n + is] * ww[5] - zz[2];

      // Projection on [0, +infty[ and on D(0, mu*zn)

      zn = ww[3] - rn * ww[0];
      if (zn > 0)
      {
        G[3] = ww[0];
        JacG[0 * Gsize + 3] = 1.;
      }
      else
      {
        G[3] = ww[3] / rn;
        JacG[3 * Gsize + 3] = 1. / rn;
      }
      zt = ww[4] - rt * ww[1];
      zs = ww[5] - rt * ww[2];

      a2 = zt * zt;
      b2 = zs * zs;
      ab = zt * zs;

      mrn = zt * zt + zs * zs;

      // if the radius is negative or the vector is null projection on the disk = 0

      if (ww[3] < 0. || mrn < 1e-16)
      {
        G[4] = ww[4] / rt;
        G[5] = ww[5] / rt;
        JacG[4 * Gsize + 4] = JacG[5 * Gsize + 5] = 1. / rt;
      }

      // if the radius is positive and the vector is non null, we compute projection on the disk

      else
      {

        if (mrn < mu2 * ww[3]*ww[3])
        {
          G[4] = ww[1];
          G[5] = ww[2];
          JacG[1 * Gsize + 4] = JacG[2 * Gsize + 5] = 1.;
        }
        else
        {
          num  = mu / sqrt(mrn);
          G[4] = (ww[4] - zt * ww[3] * num) / rt;
          G[5] = (ww[5] - zs * ww[3] * num) / rt;
          JacG[3 * Gsize + 4] = - num * zt / rt;
          JacG[3 * Gsize + 5] = - num * zs / rt;
          mrn3 = sqrt(mrn) * sqrt(mrn) * sqrt(mrn);
          JacG[1 * Gsize + 4] = mu * ww[3] * b2 / mrn3;
          JacG[2 * Gsize + 4] = -mu * ww[3] * ab / mrn3;
          JacG[1 * Gsize + 5] = -mu * ww[3] * ab / mrn3;
          JacG[2 * Gsize + 5] = mu * ww[3] * a2 / mrn3;
          JacG[4 * Gsize + 4] = (1. - rt * mu * ww[3] * b2 / mrn3) / rt;
          JacG[5 * Gsize + 4] = mu * ww[3] * ab / mrn3;
          JacG[4 * Gsize + 5] = mu * ww[3] * ab / mrn3;
          JacG[5 * Gsize + 5] = (1. - rt * mu * ww[3] * a2 / mrn3) / rt;
        }
      }

      /* compute the norm of G */

      err = dnrm2_(&Gsize, G , &incx);

      printf("Iteration %i Erreur = %14.7e\n", iter, err);

      /* for( i = 0 ; i < Gsize ; ++i ){ */
      /*       printf("G[%i] = %8.5e\n",i,G[i]); */
      /*       } */

      /*       for( i = 0 ; i < Gsize ; ++i ){ */
      /*       for( j = 0 ; j < Gsize ; ++j ){ */
      /*       printf("JacG[%i,%i] = %8.5e\t",i,j,JacG[j*Gsize+i]); */
      /*       } */
      /*       printf("\n"); */
      /*       } */

      /* **** Criterium convergence **** */

      incx =  1;
      incy =  1;
      a1   = -1.0;

      /* compute the direction www s.t X^{k+1} = X^{k} + www, where X^{k} = ww = (wa,za) */

      dcopy_(&Gsize , G , &incx , www , &incy);
      dscal_(&Gsize , &a1 , www, &incx);
      dcopy_((integer *)&mm , JacG , &incx , A , &incy);


      F77NAME(dgesv)((integer *)&Gsize, (integer *)&nrhs, A, (integer *)&Gsize, (integer *)ipiv, www, (integer *)&Gsize, (integer *)&infoDGESV);

      if (infoDGESV)
      {
        if (ispeak > 0)
        {
          printf("Problem in DGESV\n");
        }
        iparamLCP[2] = iter;
        dparamLCP[1] = err;

        free(A);
        free(G);
        free(G0);
        free(JacG);
        free(ipiv);
        free(ww);
        free(www);
        free(zz);
        *info = 2;
        return;

      }
      /*
      for( i = 0 ; i < Gsize ; ++i ){
      printf("Iteration %i -- Direction www[%i] = %8.5e\n",iter,i,www[i]);
      }
      */
      qs   = 1.0;
      daxpy_(&Gsize , &qs , ww , &incx , www , &incy);
      /*
      for( i = 0 ; i < Gsize ; ++i ){
      printf("Iteration %i -- xk+1[%i] = %8.5e\n",iter,i,www[i]);
      }
      printf("Iteration %i Erreur = %14.7e\n",iter,err);
      for( j = 0 ; j < Gsize ; ++j ){
      printf("G[%i] = %8.5e\n",j,G[j]);
      }
      */

    }
    w[in] = www[0] ;
    w[it] = www[1] ;
    w[is] = www[2] ;
    z[in] = www[3] ;
    z[it] = www[4] ;
    z[is] = www[5] ;
  }

  /*for( i = 0 ; i < nc ; ++i ){
    printf("---------------------la vitesse et la force de contact en %i est w[%i] = %8.5e et z[%i] = %8.5e\n",i,i,w[i],i, z[i]);
    }*/

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
  free(G);
  free(G0);
  free(JacG);
  free(ipiv);
  free(ww);
  free(www);
  free(zz);


  if (ispeak == 2) fclose(f101);

}
