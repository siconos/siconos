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
#include "NonSmoothDrivers.h"
#include "pfc_3D_Alart_Curnier.h"
#include "pfc_3D_Fischer_Burmeister.h"

void pfc_3D_nsgs_block(int nc , SparseBlockStructuredMatrix *M, double *q , double *z , double *w, double *mu , int *info, int *iparam , double *dparam)
{
  FILE *f101;

  int info2 = 1;
  int ispeak = iparam[1];    // verbose mode
  int n = 3 * nc;            // System dimension
  int i, j;
  pfc3D_fPtr Compute_G = NULL;
  pfc3D_fPtr Compute_JacG = NULL;
  PFC3D_local_solver pfc3D_local_solver = NULL;
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

  if (ispeak == 2) f101 = fopen("pfc_3D_nsgs_block.log" , "w+");

  /****** Check for trivial case ******/

  /* norm of vector q */
  double qs = DNRM2(n , q , incx);

  if (ispeak > 0)
  {
    printf("===================================== PFC 3D solver-block call =============================\n");
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
  /* Sub-block of vector z, corresponding to one specific contact. */
  double *zBlock = NULL;
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
    /* wLocal and zBlock are useless with AC but may be required for FB? See with Houari.*/
    if (wLocal == NULL)
    {
      fprintf(stderr, "pfc_3D_nsgs_new, memory allocation for wLocal failed.\n");
      exit(EXIT_FAILURE);
    }
    zBlock = (double*)malloc(Gsize * sizeof(*zBlock));
    if (zBlock == NULL)
    {
      fprintf(stderr, "pfc_3D_nsgs_block, memory allocation for zBlock failed.\n");
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
    fprintf(stderr, "pfc_3D_nsgs_block, memory allocation for qLocal failed.\n");
    exit(EXIT_FAILURE);
  }

  /* Sub Block of matrix M, corresponding to one contact. */
  double * MBlock = NULL;

  /** Note: at the time the last computed value of local number of iterations and local error are saved in i/dparam[3].*/

  /*******************  NSGS Iterations *******************/
  int blockNum = 0;
  int diagPos = 0;
  int iter = 0; /* Current iteration number */
  double err = 1.; /* Current error */
  double tmp1, tmp2;
  while ((iter < itermax) && (info2 > 0))
  {
    ++iter;
    diagPos = 0;
    /* Loop through the contact points */
    for (i = 0 ; i < nc ; ++i)
    {
      in = 3 * i;
      it = in + 1;
      is = it + 1;
      MBlock = M->block[diagPos];
      diagPos = diagPos + 1 + nc;
      /****  Computation of qLocal = qBlock + sum over a row of blocks in M of the products MBlock.zBlock,
       excluding the block corresponding to the current contact.
       zBlock = z[in is it]
      ****/

      qLocal[0] = q[in] ;
      qLocal[1] = q[it] ;
      qLocal[2] = q[is] ;
      /* Loop through the columns(blocks) of M to compute qLocal */
      blockNum = i * nc;
      for (j = 0; j < nc ; ++j)
      {
        if (j != i)
        {
          DGEMV(LA_NOTRANS, 3, 3, 1.0, M->block[blockNum], Gsize, &z[3 * j], incx , 1.0, qLocal, incy);
        }
        blockNum = blockNum + 1;
      }
      /* --- Alart-Curnier formulation --- */
      if (iparam[5] == 0)
      {
        /**** Local solver call ...
        (wLocal, zBlock) solution of wLocal = MBlock.zBlock + qLocal ****/
        /* The part of M which corresponds to the current block is linked to MBlock */
        (*pfc3D_local_solver)(Gsize , MBlock , qLocal , &z[in] , &w[in], mu[i] ,  &Compute_G, &Compute_JacG, IP, Ip, I3, iparam_local , dparam_local);

        /* Note that w will be updated in filter function, during the call to G */

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

        tmp1 = qLocal[1];
        tmp2 = qLocal[2];
        qLocal[1] = -(Ip[0] * tmp1 + Ip[1] * tmp2);
        qLocal[2] = -(Ip[2] * tmp1 + Ip[3] * tmp2);

        /* Initialize zBlock components */
        zBlock[0] = z[3 * i];
        zBlock[1] = mu[i] * z[3 * i] - sqrt(3) * z[3 * i + 1] / 2. - z[3 * i + 2] / 2.;
        zBlock[2] = mu[i] * z[3 * i] + sqrt(3) * z[3 * i + 1] / 2. - z[3 * i + 2] / 2.;
        zBlock[3] = 0.;
        zBlock[4] = 0.;

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
        fprintf(stderr, "pfc_3D_nsgs_block, unknown formulation type.\n");
        exit(EXIT_FAILURE);
      }
    }

    /* **** Criterium convergence **** */
    // info2 = filter_result_pfc_block(nc,Gsize,w,M,z,q,&Compute_G,mu,iparam,dparam);
    NCP_block_compute_error(n, M, q, z, ispeak, w, &err);
    if (err < dparam[0]) info2 = 0;

    //info2 = filter_result_pfc(nc,Gsize,w,M,z,q,&Compute_G,mu,iparam,dparam, work);

    if (ispeak > 0)
      printf(" Iteration number %i - Error = %14.7e\n", iter, dparam_local[1]);
  }
  /***** End of iter/err loop *****/

  dparam[3] = dparam_local[1];
  iparam[4] = iparam_local[2];
  err = dparam_local[1];
  *info = info2;
  iparam[2] = iter;
  dparam[1] = err;

  if (ispeak > 0)
  {
    if (info2 > 0)
    {
      printf(" WARNING: no convergence after %i iterations\n" , iter);
    }
    else
    {
      printf(" Convergence of solver (pfc_3D_nsgs_block) after %i iterations\n" , iter);
      printf(" The residue is : %e \n", err);
    }
  }

  /* Free memory */
  if (iparam[5] == 1)
  {
    free(wLocal);
    free(zBlock);
    free(IP);
    free(Ip);
    free(I3);
  }

  free(qLocal);

  if (ispeak == 2) fclose(f101);

}
