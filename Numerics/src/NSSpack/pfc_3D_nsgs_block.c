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
/*!\file pfc_3D_nsgs_block.c
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
 * here M is an n by n  matrix, q, z and w some n-dimensional vectors.
 *
 * \fn  pfc_3D_nsgs_new( int *nn , SparseBlockStructuredMatrix *M , double *q , double *z , double *w ,
 *                         int *info\n, int *iparamLCP , double *dparamLCP )
 *
 * Generic pfc_3D parameters:\n
 *
 * \param nn      Unchanged parameter which represents the dimension of the system.
 * \param M       Unchanged parameter which contains the components of the sparse block matrix M.
 * \param q       Unchanged parameter which contains the components of the right hand side vector.
 * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
 * \param w       Modified parameter which returns the solution of the problem.
 * \param info    Modified parameter which returns the termination value\n
 *                0 - convergence\n
 *                1 - iter = itermax, ie the simulation reached the maximum number of iterations allowed\n
 *                2 - negative diagonal term(s) in M.\n
 *
 * Specific NSGS parameters:\n
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
 * \author Houari Khenous and Franck Perignon - Creation: 12/11/2007 - Last modification 14/11/2007.
 *
 *
 */

#include "LA.h"
#include <time.h>
#include <NSSpack.h>
#include <pfc_3D_Alart_Curnier.h>
#include <pfc_3D_Fischer_Burmeister.h>

void pfc_3D_nsgs_block(int *nn , SparseBlockStructuredMatrix *M , double *q , double *z , double *w , int *info, int *iparamLCP , double *dparamLCP)
{

  FILE *f101;

  int ispeak = iparamLCP[1]; // verbose mode
  int nc = *nn;              // number of contacts
  int n = 3 * nc;            // System dimension

  int i, j;

  pfc3D_fPtr Compute_G = NULL;
  pfc3D_fPtr Compute_JacG = NULL;
  PFC3D_local_solver pfc3D_local_solver = NULL;
  double mu  = dparamLCP[0]; // Friction coefficient
  double tol = dparamLCP[1]; // Required tolerance
  int itermax = iparamLCP[0]; // Max. number of iterations

  // Local parameters
  int nb = 5;
  int     iparam_local[nb];
  iparam_local[0] = iparamLCP[0]; /* local itermax   */
  iparam_local[1] = iparamLCP[4]; /* local iteration */
  iparam_local[2] = 0;
  /* local formulation */
  /* 0 for Alart-Curnier formulation */
  /* 1 for Fischer-Burmeister formulation */
  iparam_local[3] = 1; //iparamLCP[5];
  /* local_solver */
  /* 0 for projection */
  /* 1 for newton with AC formulation */
  iparam_local[4] = 1; //iparamLCP[6];

  double  dparam_local[nb];
  dparam_local[0] = tol;  //dparamLCP[3]; /* local tolerance */
  dparam_local[1] = dparamLCP[4]; /* local error     */
  dparam_local[2] = 0.0;
  dparam_local[3] = 0.0;
  dparam_local[4] = 0.0;

  /* Initialize output */

  iparamLCP[2] = 0;
  dparamLCP[2] = 0.0;

  if (ispeak == 2) f101 = fopen("pfc_3D_nsgs_block.log" , "w+");

  /* Connection to local solver functions */
  int Gsize;
  // 1-Projection solver
  if (iparam_local[4] == 0)
  {
    pfc3D_local_solver = &pfc_3D_projection;
    Gsize  = 3;
  }
  // 2-Newton solver
  else if (iparam_local[4] == 1)
  {
    // Alart-Curnier formulation
    if (iparam_local[3] == 0)
    {
      Gsize  = 3;
      Compute_G = &Compute_G_AC;
      Compute_JacG = &Compute_JacG_AC;
      pfc3D_local_solver = &pfc_3D_newton;
    }
    // Fischer-Burmeister formulation
    else if (iparam_local[3] == 1)
    {
      Gsize  = 5;
      Compute_G = &Compute_G_FB;
      Compute_JacG = &Compute_JacG_FB;
      pfc3D_local_solver = &pfc_3D_newton;
    }
    else
    {
      fprintf(stderr, "Numerics, pfc_3D_nsgs_block failed. Unknown formulation type.\n");
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    fprintf(stderr, "Numerics, pfc_3D_nsgs_block failed. Unknown solver type.\n");
    exit(EXIT_FAILURE);
  }

  int incx = 1, incy = 1;

  /****** Check for trivial case ******/

  /* norm of vector q */
  double qs = DNRM2(n , q , incx);

  if (ispeak > 0) printf("\n ||q||= %g \n" , qs);

  if (qs <= 1e-16)
  {
    // q norm equal to zero (less than 1e-16)
    // -> trivial solution: z = 0 and w = q
    for (i = 0 ; i < n ; ++i)
      w[i] = q[i];
    z[i] = 0.;
    *info = 0;
    return;
  }

  /****** Non trivial case ******/

  double den = 1.0 / qs;
  int mm = Gsize * Gsize;

  /* Memory allocation */
  double *wOld    = malloc(n * sizeof(*wOld)); /* w value computed at previous iteration */
  if (wOld == NULL)
  {
    fprintf(stderr, "pfc_3D_nsgs_block, memory allocation for wOld failed.\n");
    exit(EXIT_FAILURE);
  }
  double *wLocal   = malloc(Gsize * sizeof(*wLocal));
  if (wLocal == NULL)
  {
    fprintf(stderr, "pfc_3D_nsgs_block, memory allocation for wLocal failed.\n");
    exit(EXIT_FAILURE);
  }
  double *zBlock   = malloc(Gsize * sizeof(*zBlock));
  if (zBlock == NULL)
  {
    fprintf(stderr, "pfc_3D_nsgs_block, memory allocation for zBlock failed.\n");
    exit(EXIT_FAILURE);
  }
  double *qLocal  = malloc(Gsize * sizeof(*qLocal));
  if (qLocal == NULL)
  {
    fprintf(stderr, "pfc_3D_nsgs_block, memory allocation for qLocal failed.\n");
    exit(EXIT_FAILURE);
  }

  /* Initializations */
  for (i = 0 ; i < Gsize ; ++i)
  {
    wLocal[i] =  zBlock[i] = qLocal[i] = 0.;
  }
  for (i = 0 ; i < n ; ++i)
  {
    wOld[i] = w[i] = 0.0;
  }

  /**** Start NSGS iterations ****/

  int iter = 0; /* Current iteration number */
  double err  = 1.; /* Current error */

  /* Parameters for Alart-Curnier formulation */
  double *param1  = malloc(mm * sizeof(*param1));
  double *param2  = malloc(mm * sizeof(*param2));
  double *param3  = malloc(mm * sizeof(*param3));
  /* Parameters for Fischer-Burmeister formulation */
  double *IP = malloc(2 * sizeof(*IP));
  double *Ip = malloc(2 * 2 * sizeof(*Ip));
  double *I3 = malloc(2 * sizeof(*I3));
  double *V  = malloc(2 * sizeof(*V));

  /*   /\* Initialize w with q *\/ */
  /*   DCOPY( n , q , incx , w , incy ); */

  int in, it, is;
  double a1 = 1.0, num;
  double *MBlock = NULL; /* current sub-block of M*/

  while ((iter < itermax) && (err > tol))
  {
    ++iter;

    /* Loop through the contact points */
    for (i = 0 ; i < nc ; ++i)
    {
      in = 3 * i;
      it = in + 1;
      is = it + 1;

      /* Note FP: pb with Fischer-Burmeister, Gsize = 5 ie MBlock.size = 25, not taken into account here. */

      /* The part of M which corresponds to the current block is linked to MBlock */
      MBlock = M->block[i];

      /* --- Alart-Curnier formulation --- */
      if (iparam_local[3] == 0)
      {

        /* zBlock is initialized with the block of z that corresponds to the current contact */
        for (j = 0 ; j < Gsize ; ++j)
          zBlock[j] = z[Gsize * i + j];

        /****  Computation of qLocal = qBlock + sum over a row of blocks in M of the products MBlock.zBlock,
         excluding the block corresponding to the current contact. ****/

        /* z current block set to zero, to exclude current contact block */
        z[in] = 0.0;
        z[it] = 0.0;
        z[is] = 0.0;
        /* qLocal computation*/
        //        incx = n;

        DGEMV(LA_NOTRANS, 3, n, 1.0, MBlock, n, z , incx , 1.0, qLocal, incy);
        /*        qLocal[0] = q[in] + DDOT( n , &M[in] , incx , z , incy ); */
        /*        qLocal[1] = q[it] + DDOT( n , &M[it] , incx , z , incy ); */
        /*        qLocal[2] = q[is] + DDOT( n , &M[is] , incx , z , incy ); */


        /**** Local solver call ...
        (wLocal, zBlock) solution of wLocal = MBlock.zBlock + qLocal ****/
        for (j = 0 ; j < mm ; ++j)
        {
          param1[j] = param2[j] = param3[j] = 0.;
        }

        (*pfc3D_local_solver)(Gsize , MBlock , qLocal , zBlock , wLocal , mu ,  &Compute_G, &Compute_JacG, param1, param2, param3, iparam_local , dparam_local);

        /* z current block set to zBlock */
        z[in] = zBlock[0];
        z[it] = zBlock[1];
        z[is] = zBlock[2];

      }
      /* --- Fischer-Burmeister formulation --- */
      else if (iparam_local[3] == 1)
      {
        /*  double PI = 3.14; */
        /*  double alpha1 = PI/6; */
        /*  double alpha2 = 5*PI/6; */
        /*  double alpha3 = 3*PI/2; */

        /*  double deter = cos(alpha1)*sin(alpha2) - sin(alpha1)*cos(alpha2); */

        /*  /\* Ip = matrix_inv2(Ip`)*\/ */
        /*  /\* IP = matrix_inv2(Ip`)*mup *\/ */
        /*  /\* I3 = matrix_inv2(Ip)*e3 *\/ */

        /*  IP[0] = mu*(sin(alpha2)-sin(alpha1))/deter; */
        /*  IP[1] = mu*(cos(alpha1)-cos(alpha2))/deter; */

        /*  I3[0] = (cos(alpha3)*sin(alpha2) - sin(alpha3)*cos(alpha2))/deter; */
        /*  I3[1] = (sin(alpha3)*cos(alpha1) - cos(alpha3)*sin(alpha1))/deter; */

        /*  Ip[0*2+0] =  sin(alpha2)/deter; */
        /*  Ip[1*2+0] = -sin(alpha1)/deter; */
        /*  Ip[0*2+1] = -cos(alpha2)/deter; */
        /*  Ip[1*2+1] =  cos(alpha1)/deter; */
        /* Initialize IP, Ip, I3, V components to zero */
        V[0] = V[1] = 0.0;
        IP[0] = 0.0;
        IP[1] = 2.*mu;
        I3[0] = -1.;
        I3[1] = -1.;
        Ip[0] =  1. / sqrt(3.);
        Ip[2] = -1. / sqrt(3.);
        Ip[1] =  1.;
        Ip[3] =  1.;

        /* Initialize zBlock components */
        zBlock[0] = z[3 * i];
        zBlock[1] = mu * z[3 * i] - sqrt(3) * z[3 * i + 1] / 2. - z[3 * i + 2] / 2.;
        zBlock[2] = mu * z[3 * i] + sqrt(3) * z[3 * i + 1] / 2. - z[3 * i + 2] / 2.;
        zBlock[3] = 0.;
        zBlock[4] = 0.;

        /* Initialize qLocal components */
        qLocal[0] = q[in];
        /*  qLocal[1] = -(Ip[0*2+0]*q[it] +Ip[0*2+1]*q[is]); */
        /*  qLocal[2] = -(Ip[1*2+0]*q[it] +Ip[1*2+1]*q[is]); */

        /*     for( ii = 0 ; ii < Gsize ; ++ii ) */
        /*  printf("b[%i] =  %14.7e\n",ii,qLocal[ii]); printf("\n"); */

        /* Local solver call */
        (*pfc3D_local_solver)(Gsize , MBlock , qLocal , zBlock , wLocal , mu , &Compute_G, &Compute_JacG, Ip, IP, I3, iparam_local , dparam_local);

        /* set z current block */
        z[in] = zBlock[0];
        z[it] = Ip[0] * (mu * zBlock[0] - zBlock[1]) + Ip[2] * (mu * zBlock[0] - zBlock[2]);
        z[is] = Ip[1] * (mu * zBlock[0] - zBlock[1]) + Ip[3] * (mu * zBlock[0] - zBlock[2]);


        /*  double nrm = sqrt(z[it]*z[it]+z[is]*z[is]); */
        /*  /\* printf("nrm =  %14.7e\n",nrm); *\/ */
        /*  if(nrm){ */
        /*    V[0] = z[it]/nrm; */
        /*    V[1] = z[is]/nrm; */
        /*  } */
        /*  else{ */
        /*    V[0] = 0.; */
        /*    V[1] = -1.; */
        /*  } */
        /*  /\*   for( ii = 0 ; ii < 2 ; ++ii ) *\/ */
        /*  /\*   printf("V[%i] =  %14.7e\n",ii,V[ii]); printf("\n"); *\/ */

        /*  w[in] = wLocal[0]; */
        /*  w[it] = -(sqrt(3)*wLocal[1]/2. - sqrt(3)*wLocal[2]/2. + V[0]*wLocal[4]); */
        /*  w[is] = -(wLocal[1]/2. + wLocal[2]/2. + V[1]*wLocal[4]); */


      }
      else
      {
        fprintf(stderr, "pfc_3D_nsgs_block, unknown formulation type.\n");
        exit(EXIT_FAILURE);
      }
    }

    /* **** Criterium convergence **** */

    incx =  1;

    /* Copy w_(iter-1) = w into wOld */
    DCOPY(n , w , incx , wOld , incy);

    /* Compute w_(iter) = M*z_(iter) + q */
    DCOPY(n , q , incx , w , incy);
    DGEMV(LA_NOTRANS , n , n , a1 , M , n , z , incx , a1 , w , incy);

    /* Compute err = norm2(w_(iter) - w_(iter-1)) / norm2(q) */
    qs   = -1.0;
    DAXPY(n , qs , w , incx , wOld , incy);
    num = DNRM2(n, wOld , incx);
    err = num * den;
    /*    printf("-----------------------------------Iteration %i Erreur = %14.7e\n",iter,err); */
  }
  /***** End of iter/err loop *****/

  /* free memory for parameters */
  free(param1);
  free(param2);
  free(param3);
  free(IP);
  free(Ip);
  free(I3);
  free(V);

  iparamLCP[2] = iter;
  dparamLCP[2] = err;

  if (ispeak > 0)
  {
    if (err > tol)
    {
      printf(" Numerics, pfc_3D_nsgs_block: no convergence after %i iterations\n" , iter);
      printf(" The residue is : %e \n", err);
      *info = 1;
    }
    else
    {
      printf(" Numerics, pfc_3D_nsgs_block: convergence after %i iterations\n" , iter);
      printf(" The residue is : %e \n", err);
      *info = 0;
    }
  }
  else
  {
    if (err > tol) *info = 1;
    else *info = 0;
  }

  /* Free memory */
  free(wOld);
  free(MBlock);
  free(wLocal);
  free(zBlock);
  free(qLocal);

  if (ispeak == 2) fclose(f101);

}
