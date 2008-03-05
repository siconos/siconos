/* Siconos-Numerics version 2.0.1, Copyright INRIA 2005-2008.
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
#include "pfc_3D_Solvers.h"
#include "NCP.h"

void pfc_3D_nsgs_new(int nc , double *M , double *q , double *z , double *w, double *mu , int *info, int *iparam , double *dparam)
{
  FILE *f101;

  int info2 = 1;
  int ispeak = iparam[1];    // verbose mode
  int n = 3 * nc;            // System dimension
  int i, j;
  pfc3D_fPtr Compute_G = NULL;
  pfc3D_fPtr Compute_JacG = NULL;
  PFC3D_local_solver pfc3D_local_solver = NULL;
  double tol = dparam[0]; // Required tolerance
  int itermax = iparam[0]; // Max. number of iterations
  int incx = 1, incy = 1;
  int in, it, is;

  /* vectors of parameters for local solver */
  int iparam_local[7];
  iparam_local[0] = iparam[3];// local maximum number of iteration
  iparam_local[1] = iparam[1]; // ispeak;
  iparam_local[5] = iparam[5]; // formulation;
  iparam_local[6] = iparam[6]; // formulation;
  /* Output for iparam_local: [2], number of iterations */
  /* [3] and [4] are useless for iparam_local */
  double dparam_local[4];
  dparam_local[0] = dparam[2]; // local tolerance
  /* dparam_local[1]: output, local error. */

  if (ispeak == 2) f101 = fopen("pfc_3D_nsgs_new.log" , "w+");

  /****** Check for trivial case ******/

  /* norm of vector q */
  double qs = DNRM2(n , q , incx);

  if (ispeak > 0)
  {
    printf("===================================== PFC 3D solver call =============================\n");
    printf(" - Formulation type: ");
    if (iparam[5] == 1) printf(" Fischer-Burmeister with Newton solver");
    else if (iparam[5] == 0)
    {
      if (iparam[6] == 1)
        printf(" Alart-Curnier with Newton solver.");
      else if (iparam[6] == 0)
        printf(" Alart-Curnier with projection.");
      else
        printf(" Alart-Curnier with unknown solver.");
    }
    else
      printf(" unknown.");
    printf("\n");
    printf("\n ||q||= %g \n" , qs);
  }

  if (qs <= 1e-16)
  {
    // q norm equal to zero (less than 1e-16)
    // -> trivial solution: z = 0 and w = q
    for (i = 0 ; i < n ; ++i)
    {
      w[i] = q[i];
      z[i] = 0.;
    }
    *info = 0;
    iparam[2] = 0;
    dparam[1] = 0.0;
    return;
  }

  /****** Non trivial cases ******/

  /******************* Connection to local solver functions *******************/

  int Gsize;
  /* Sub-block of vector w, corresponding to the current contact. */
  double *wLocal = NULL;
  /* Parameters for Alart-Curnier and Fischer-Burmeister formulations */
  double *IP = NULL, *Ip = NULL, *I3 = NULL;
  /*Note that the allocation of memory for those vectors depends on the formulation/solver type */
  /* At the time, they are never used in AC case => not allocate. */

  // Alart-Curnier formulation
  if (iparam[5] == 0)
  {
    Gsize  = 3;
    // 1-Projection solver
    if (iparam[6] == 0)
    {
      pfc3D_local_solver = &pfc_3D_projection;
    }
    // 2-Newton solver
    else if (iparam[6] == 1)
    {
      Compute_G = &Compute_G_AC;
      Compute_JacG = &Compute_JacG_AC;
      pfc3D_local_solver = &pfc_3D_newton;
    }
  }
  // Fischer-Burmeister formulation
  else if (iparam[5] == 1)
  {
    if (iparam[6] != 1)
    {
      fprintf(stderr, "Numerics, pfc_3D_nsgs_new failed. Fisher-Burmeister formulation not yet implemented with the required local solver.\n");
      exit(EXIT_FAILURE);
    }
    Gsize  = 5;
    Compute_G = &Compute_G_FB;
    Compute_JacG = &Compute_JacG_FB;
    pfc3D_local_solver = &pfc_3D_newton;
    wLocal = (double*)malloc(Gsize * sizeof(*wLocal));
    /* wLocal is useless with AC but may be required for FB? See with Houari.*/
    if (wLocal == NULL)
    {
      fprintf(stderr, "pfc_3D_nsgs_new, memory allocation for wLocal failed.\n");
      exit(EXIT_FAILURE);
    }
    IP = malloc(2 * sizeof(*IP));
    Ip = malloc(2 * 2 * sizeof(*Ip));
    I3 = malloc(2 * sizeof(*I3));
    if (IP == NULL ||  Ip == NULL || I3 == NULL)
    {
      fprintf(stderr, "pfc_3D_nsgs_new, memory allocation for param failed.\n");
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    fprintf(stderr, "Numerics, pfc_3D_nsgs_new failed. Unknown formulation type.\n");
    exit(EXIT_FAILURE);
  }

  /******************* Memory allocation *******************/

  /* "Augmented" vector q, used in local problem. */
  double *qLocal  = (double*)malloc(Gsize * sizeof(*qLocal));
  if (qLocal == NULL)
  {
    fprintf(stderr, "pfc_3D_nsgs_new, memory allocation for qLocal failed.\n");
    exit(EXIT_FAILURE);
  }

  /* Sub Block of matrix M, corresponding to one contact. */
  double *MBlock = (double*)malloc(3 * 3 * sizeof(*MBlock));
  if (MBlock == NULL)
  {
    fprintf(stderr, "pfc_3D_nsgs_new, memory allocation for MBlock failed.\n");
    exit(EXIT_FAILURE);
  }

  /* Sub-block of vector z, corresponding to one specific contact. */
  double *zBlock   = (double*)malloc(Gsize * sizeof(*zBlock));
  if (zBlock == NULL)
  {
    fprintf(stderr, "pfc_3D_nsgs_new, memory allocation for zBlock failed.\n");
    exit(EXIT_FAILURE);
  }

  /** Note: at the time the last computed value of local number of iterations and local error are saved in i/dparam[3].*/

  /* temporary work vector */
  /*   double **work  = (double**)malloc(6*sizeof(**work)); */
  /*   /\*  */
  /*      work[0] = qLocal */
  /*      work[1] = zBlock */
  /*      work[2] = temp copy of z */
  /*      work[3] = res */
  /*      work[4] = wLocal */
  /*      work[5] = currentBlock */
  /*   *\/ */
  /*   work[0] = (double*)malloc(Gsize*sizeof(double)); */
  /*   work[1] = (double*)malloc(Gsize*sizeof(double)); */
  /*   work[2] = (double*)malloc(n*sizeof(double)); */
  /*   work[3] = (double*)malloc(Gsize*nc*sizeof(double)); */
  /*   work[4] = (double*)malloc(Gsize*sizeof(double)); */
  /*   work[5] = (double*)malloc(Gsize*Gsize*sizeof(double)); */

  /*   if ( work == NULL ) */
  /*     { */
  /*       fprintf(stderr,"pfc_3D_nsgs_new, memory allocation for work failed.\n"); */
  /*       exit(EXIT_FAILURE); */
  /*     }  */

  /*******************  NSGS Iterations *******************/
  int iter = 0; /* Current iteration number */
  double err = 1.; /* Current error */
  double tmp1, tmp2;

  while ((iter < itermax) && (info2 > 0))
  {
    ++iter;
    /* Loop through the contact points */
    for (i = 0 ; i < nc ; ++i)
    {
      in = 3 * i;
      it = in + 1;
      is = it + 1;

      /* The part of M which corresponds to the current block is copied into MBlock */
      MBlock[0] = M[(in) * n + in];
      MBlock[1] = M[(in) * n + it];
      MBlock[2] = M[(in) * n + is];
      MBlock[3] = M[(it) * n + in];
      MBlock[4] = M[(it) * n + it];
      MBlock[5] = M[(it) * n + is];
      MBlock[6] = M[(is) * n + in];
      MBlock[7] = M[(is) * n + it];
      MBlock[8] = M[(is) * n + is];

      /* --- Alart-Curnier formulation --- */
      if (iparam[5] == 0)
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
        incx = n;
        qLocal[0] = q[in] + DDOT(n , &M[in] , incx , z , incy);
        qLocal[1] = q[it] + DDOT(n , &M[it] , incx , z , incy);
        qLocal[2] = q[is] + DDOT(n , &M[is] , incx , z , incy);

        /**** Local solver call ...
        (wLocal, zBlock) solution of wLocal = MBlock.zBlock + qLocal
        wLocal = w[in it is]
        ****/

        (*pfc3D_local_solver)(Gsize , MBlock , qLocal , zBlock , &w[in], mu[i] ,  &Compute_G, &Compute_JacG, IP, Ip, I3, iparam_local , dparam_local);

        /* z current block set to zBlock */
        z[in] = zBlock[0];
        z[it] = zBlock[1];
        z[is] = zBlock[2];
      }
      /* --- Fischer-Burmeister formulation --- */
      else if (iparam[5] == 1)
      {

        /*  /\* Ip = matrix_inv2(Ip`)*\/ */
        /*  /\* IP = matrix_inv2(Ip`)*mup *\/ */
        /*  /\* I3 = matrix_inv2(Ip)*e3 *\/ */
        IP[0] = 0.0;
        IP[1] = 2.*mu[i];
        I3[0] = -1.;
        I3[1] = -1.;
        Ip[0] =  1. / sqrt(3.);
        Ip[2] = -1. / sqrt(3.);
        Ip[1] =  1.;
        Ip[3] =  1.;

        /* Initialize zBlock components */
        zBlock[0] = z[3 * i];
        zBlock[1] = mu[i] * z[3 * i] - sqrt(3) * z[3 * i + 1] / 2. - z[3 * i + 2] / 2.;
        zBlock[2] = mu[i] * z[3 * i] + sqrt(3) * z[3 * i + 1] / 2. - z[3 * i + 2] / 2.;
        zBlock[3] = 0.;
        zBlock[4] = 0.;

        /* Initialize qLocal components */
        incx = n;

        z[in] = 0.0;
        z[it] = 0.0;
        z[is] = 0.0;

        qLocal[0] = q[in] + DDOT(n , &M[in] , incx , z , incy);
        tmp1 = q[it] + DDOT(n , &M[it] , incx , z , incy);
        tmp2 = q[is] + DDOT(n , &M[is] , incx , z , incy);
        qLocal[1] = -(Ip[0] * tmp1 + Ip[1] * tmp2);
        qLocal[2] = -(Ip[2] * tmp1 + Ip[3] * tmp2);

        /* Local solver call */

        (*pfc3D_local_solver)(Gsize , MBlock , qLocal , zBlock , wLocal, mu[i] , &Compute_G, &Compute_JacG, Ip, IP, I3, iparam_local, dparam_local);

        /* set z current block */
        z[in] = zBlock[0];
        z[it] = Ip[0] * (mu[i] * zBlock[0] - zBlock[1]) + Ip[2] * (mu[i] * zBlock[0] - zBlock[2]);
        z[is] = Ip[1] * (mu[i] * zBlock[0] - zBlock[1]) + Ip[3] * (mu[i] * zBlock[0] - zBlock[2]);

        /*  w[in] = wLocal[0]; */
        /*  w[it] = -(sqrt(3)*wLocal[1]/2. - sqrt(3)*wLocal[2]/2. + V[0]*wLocal[4]); */
        /*  w[is] = -(wLocal[1]/2. + wLocal[2]/2. + V[1]*wLocal[4]); */
      }
      else
      {
        fprintf(stderr, "pfc_3D_nsgs_new, unknown formulation type.\n");
        exit(EXIT_FAILURE);
      }
    }
    /* **** Criterium convergence **** */

    NCP_compute_error(n , M , q , z , ispeak , w , &err);
    if (err < tol) info2 = 0;

    //info2 = filter_result_pfc(nc,Gsize,w,M,z,q,&Compute_G,mu,iparam,dparam, work);

    if (ispeak > 0)
      printf(" Iteration number %i - Error = %14.7e\n", iter, err);
  }
  /***** End of iter/err loop *****/


  dparam[3] = err;//dparam_local[1];
  iparam[4] = iparam_local[2];
  *info = info2;
  iparam[2] = iter;
  dparam[1] = err;

  if (ispeak > 0)
  {
    if (err > tol)
    {
      printf(" WARNING: no convergence after %i iterations\n" , iter);
      printf(" The residue is : %e \n", err);
    }
    else
    {
      printf(" Convergence after %i iterations\n" , iter);
      printf(" The residue is : %e \n", err);
    }
  }

  /* Free memory */
  if (iparam[5] == 1)
  {
    free(wLocal);
    free(IP);
    free(Ip);
    free(I3);
  }
  free(MBlock);
  free(zBlock);
  free(qLocal);
  /*   for(i=0;i<5;++i) */
  /*     free(work[i]); */
  /*   free(work); */
  if (ispeak == 2) fclose(f101);

}
