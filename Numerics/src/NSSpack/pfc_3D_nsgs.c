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
/*!\file pfc_3D_nsgs.c
 *
 * This subroutine allows the primal resolution of contact problems with friction.\n
 *
 *   Try \f$(z,w)\f$ such that:\n
 *   \f$
 *    \left\lbrace
 *     \begin{array}{l}
 *      M z + q = w\\
 *      0 \le z_n \perp w_n \ge 0\\
 *      -w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
 *     \end{array}
 *    \right.
 *   \f$
 *
 * here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an  * n-dimensional vector.
 *
 * \fn  pfc_3D_nsgs( int *nn , double *vec , double *q , double *z , double *w ,
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
#include "LA.h"
#include <time.h>
#include <NSSpack.h>
#include <pfc_3D_Alart_Curnier.h>
#include <pfc_3D_Fischer_Burmeister.h>

void (*pfc_3D_local_solver)(int n , double *C , double *zzz , double *zz , double *ww , double mu , Compute_G_function(*Compute_G), Compute_JacG_function(*Compute_JacG), int *iparam_local , double *dparam_local) = NULL;

void pfc_3D_nsgs(int *nn , double *vec , double *q , double *z , double *w , int *info, int *iparamLCP , double *dparamLCP)
{

  FILE *f101;

  int local_formulation, local_solver, n, in, it, is, ispeak, itermax, nc, i, j, iter, Gsize, mm;
  double err, tol, mu;
  double qs, a1, den, num;
  int incx, incy;
  double *W, *C, *ww, *zz, *zzz;

  Compute_G_function Compute_G;
  Compute_JacG_function Compute_JacG;

  int nb = 5;
  int     iparam_local[nb];
  double  dparam_local[nb];

  //  clock_t t1,t2;

  for (i = 0 ; i < nb ; ++i) iparam_local[i] = 0;
  for (i = 0 ; i < nb ; ++i) dparam_local[i] = 0.0;




  ispeak = 1;
  nc     = *nn;
  incx   = 1;
  incy   = 1;
  n      = 3 * nc;

  /* Recup input */

  itermax = iparamLCP[0];
  ispeak  = iparamLCP[1];

  mu  = dparamLCP[0];
  tol = dparamLCP[1];

  /* Initialize output */

  iparamLCP[2] = 0;
  dparamLCP[2] = 0.0;

  dparam_local[0] = dparamLCP[3];
  dparam_local[1] = dparamLCP[4];

  iparam_local[0] = iparamLCP[3];
  iparam_local[1] = iparamLCP[4];

  local_formulation = 0;//iparamLCP[5];
  local_solver      = 1;//iparamLCP[6];


  if (ispeak == 2) f101 = fopen("pfc_3D_nlgs.log" , "w+");

  /* local_solver */
  /* 0 for projection */
  /* 1 for newton with AC formulation */

  /* local_formulation */
  /* 0 for Alart-Curnier formulation */
  /* 1 for Fischer-Burmeister formulation */

  if (local_solver == 0)
  {
    pfc_3D_local_solver = &pfc_3D_projection;
    Gsize  = 3;
  }
  else
  {
    if (local_formulation == 0)
    {
      Gsize  = 3;
      Compute_G    = &Compute_G_AC;
      Compute_JacG = &Compute_JacG_AC;
      pfc_3D_local_solver = &pfc_3D_newton;
    }

    /*  else{ */
    /*       Gsize  = 5; */
    /*       Compute_G    = &Compute_G_FB; */
    /*       Compute_JacG = &Compute_JacG_FB; */
    /*       pfc_3D_local_solver = &pfc_3D_newton; */
    /*     } */
  }

  mm = Gsize * Gsize;

  /* Allocation */

  C    = (double*)malloc(mm * sizeof(double));
  W    = (double*)malloc(n * sizeof(double));
  ww   = (double*)malloc(Gsize * sizeof(double));
  zz   = (double*)malloc(Gsize * sizeof(double));
  zzz  = (double*)malloc(Gsize * sizeof(double));

  /* Intialization */
  for (i = 0 ; i < Gsize ; ++i)
  {
    ww[i] =  zz[i] = zzz[i] = 0.;
    for (j = 0 ; j < Gsize ; ++j)
      C[j * Gsize + i] = 0.;
  }

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

    free(W);
    free(C);
    free(ww);
    free(zz);
    free(zzz);

    *info = 0;
    return;
  }

  /* Intialization of w*/
  for (i = 0 ; i < n ; ++i)
  {
    W[i] = 0;
    w[i]  = 0.;
  }

  /*start NSGS iterations*/
  iter = 0;
  err  = 1.;

  while ((iter < itermax) && (err > tol))
  {

    ++iter;
    /*    printf("---------------- Iteration %i -----------------\n",iter); */

    incx = 1;
    incy = 1;

    DCOPY(n , w , incx , W , incy);
    DCOPY(n , q , incx , w , incy);

    for (i = 0 ; i < nc ; ++i)
    {

      /*       printf("---------------- contact point %i -----------------\n",i); */

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


      for (j = 0 ; j < Gsize ; ++j)
        zz[j] = z[Gsize * i + j];

      incx = n;
      incy = 1;

      z[in] = 0.0;
      z[it] = 0.0;
      z[is] = 0.0;

      zzz[0] = q[in] + DDOT(n , &vec[in] , incx , z , incy);
      zzz[1] = q[it] + DDOT(n , &vec[it] , incx , z , incy);
      zzz[2] = q[is] + DDOT(n , &vec[is] , incx , z , incy);

      (*pfc_3D_local_solver)(Gsize , C , zzz , zz , ww , mu , &Compute_G, &Compute_JacG, iparam_local , dparam_local);


      z[in] = zz[0];
      z[it] = zz[1];
      z[is] = zz[2];


    }
    /* **** Criterium convergence **** */

    incx =  1;
    incy =  1;
    a1 = 1.0;

    DGEMV(LA_NOTRANS , n , n , a1 , vec , n , z , incx , a1 , w , incy);

    qs   = -1.0;
    DAXPY(n , qs , w , incx , W , incy);
    num = DNRM2(n, W , incx);
    err = num * den;
    /*  printf("-----------------------------------Iteration %i Erreur = %14.7e\n",iter,err); */
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

  free(W);
  free(C);
  free(ww);
  free(zz);
  free(zzz);

  if (ispeak == 2) fclose(f101);

}
