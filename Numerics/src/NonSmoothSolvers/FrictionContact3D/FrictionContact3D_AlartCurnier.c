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

#include "LA.h"
#include "FrictionContact3D_Solvers.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*Static variables */

/* The global problem of size n= 3*nc, nc being the number of contacts, is locally saved in MGlobal and qGlobal */
/* mu corresponds to the vector of friction coefficients */
/* note that either MGlobal or MBGlobal is used, depending on the chosen storage */
static int n = 0;
static const NumericsMatrix* MGlobal = NULL;
static const double* qGlobal = NULL;
static const double* mu = NULL;

/* Local problem operators */
static const int nLocal = 3;
static double* MLocal;
static int isMAllocatedIn = 0; /* True if a malloc is done for MLocal, else false */
static double velocityLocal[3];
static double qLocal[3];
static double mu_i = 0.0;
static double an;
static double at;
static double projN;
static double projT;
static double projS;

void AC_fillMLocal(int contact)
{
  // Dense storage
  int storageType = MGlobal->storageType;
  if (storageType == 0)
  {
    int in = 3 * contact, it = in + 1, is = it + 1;
    int inc = n * in;
    double * MM = MGlobal->matrix0;
    /* The part of MM which corresponds to the current block is copied into MLocal */
    MLocal[0] = MM[inc + in];
    MLocal[1] = MM[inc + it];
    MLocal[2] = MM[inc + is];
    inc += n;
    MLocal[3] = MM[inc + in];
    MLocal[4] = MM[inc + it];
    MLocal[5] = MM[inc + is];
    inc += n;
    MLocal[6] = MM[inc + in];
    MLocal[7] = MM[inc + it];
    MLocal[8] = MM[inc + is];
  }
  else if (storageType == 1)
  {
    int diagPos = getDiagonalBlockPos(MGlobal->matrix1, contact);
    //MLocal = MGlobal->matrix1->block[diagPos];
    DCOPY(9, MGlobal->matrix1->block[diagPos], 1, MLocal, 1);

  }
  else
    numericsError("FrictionContact3D_AlartCurnier:AC_fillMLocal() -", "unknown storage type for matrix M");

}

void frictionContact3D_AC_initialize(int n0, const NumericsMatrix*const M0, const double*const q0, const double*const mu0)
{
  /*
    INPUT: the global problem operators: n0 (size), M0, q0 and mu0, vector of friction coefficients.
    In initialize, these operators are "connected" to their corresponding static variables, that will be used to build local problem
    for each considered contact.
    Local problem is built during call to update (which depends on the storage type for M).
  */

  n = n0;
  MGlobal = M0;
  qGlobal = q0;
  mu = mu0;
  /*  if(MGlobal->storageType == 0) */
  /*     { */
  MLocal = (double*)malloc(nLocal * nLocal * sizeof(*MLocal));
  isMAllocatedIn = 1;
  /*     } */
  /*   else */
  /*     isMAllocatedIn = 0; */
}

void frictionContact3D_AC_update(int contact, double * reaction)
{
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */
  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  AC_fillMLocal(contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock,
   excluding the block corresponding to the current contact. ****/

  int in = 3 * contact, it = in + 1, is = it + 1;
  /* reaction current block set to zero, to exclude current contact block */
  double rin = reaction[in] ;
  double rit = reaction[it] ;
  double ris = reaction[is] ;
  reaction[in] = 0.0;
  reaction[it] = 0.0;
  reaction[is] = 0.0;
  /* qLocal computation*/
  qLocal[0] = qGlobal[in];
  qLocal[1] = qGlobal[it];
  qLocal[2] = qGlobal[is];

  if (MGlobal->storageType == 0)
  {
    double * MM = MGlobal->matrix0;
    int incx = n, incy = 1;
    qLocal[0] += DDOT(n , &MM[in] , incx , reaction , incy);
    qLocal[1] += DDOT(n , &MM[it] , incx , reaction , incy);
    qLocal[2] += DDOT(n , &MM[is] , incx , reaction , incy);
  }
  else if (MGlobal->storageType == 1)
  {
    /* qLocal += rowMB * reaction
    with rowMB the row of blocks of MGlobal which corresponds to the current contact
    */
    rowProdNoDiagSBM(n, 3, contact, MGlobal->matrix1, reaction, qLocal, 0);
  }
  reaction[in] = rin;
  reaction[it] = rit;
  reaction[is] = ris;

  /* Friction coefficient for current block*/
  mu_i = mu[contact];
}

/* Compute function F(Reaction) */
void F_AC(int Fsize, double *reaction , double *F, int up2Date)
{
  if (F == NULL)
  {
    fprintf(stderr, "FrictionContact3D_AlartCurnier::F_AC error:  F == NULL.\n");
    exit(EXIT_FAILURE);
  }
  if (Fsize != nLocal)
  {
    fprintf(stderr, "FrictionContact3D_AlartCurnier::F error, wrong block size.\n");
    exit(EXIT_FAILURE);
  }

  /* up2Date = 1 = true if jacobianF(n, reaction,jacobianF) has been called just before jacobianFMatrix(...). In that case the computation of
     velocityLocal is not required again.
  */
  if (up2Date == 0)
  {
    /* velocityLocal = M.reaction + qLocal */
    int incx = 1, incy = 1;
    DCOPY(Fsize, qLocal, incx, velocityLocal, incy);
    DGEMV(LA_NOTRANS, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, velocityLocal, incy);
    /*   velocityLocal[0] = MLocal[0]*reaction[0] + MLocal[Fsize]*reaction[1] + MLocal[2*Fsize]*reaction[2] + qLocal[0]; */
    /*   velocityLocal[1] = MLocal[1]*reaction[0] + MLocal[Fsize+1]*reaction[1] + MLocal[2*Fsize+1]*reaction[2] + qLocal[1]; */
    /*   velocityLocal[2] = MLocal[2]*reaction[0] + MLocal[Fsize+2]*reaction[1] + MLocal[2*Fsize+2]*reaction[2] + qLocal[2]; */

    an = 1. / MLocal[0];
    double alpha = MLocal[Fsize + 1] + MLocal[2 * Fsize + 2];
    double det = MLocal[1 * Fsize + 1] * MLocal[2 * Fsize + 2] - MLocal[2 * Fsize + 1] + MLocal[1 * Fsize + 2];
    double beta = alpha * alpha - 4 * det;
    at = 2 * (alpha - beta) / ((alpha + beta) * (alpha + beta));

    /* Projection on [0, +infty[ and on D(0, mu*reactionn) */
    projN = reaction[0] - an * velocityLocal[0];
    projT = reaction[1] - at * velocityLocal[1];
    projS = reaction[2] - at * velocityLocal[2];
  }

  double num;
  double coef2 = mu_i * mu_i;
  if (projN > 0)
  {
    F[0] = velocityLocal[0];
  }
  else
  {
    F[0] = reaction[0] / an;
  }

  double mrn = projT * projT + projS * projS;
  if (mrn <= coef2 * reaction[0]*reaction[0])
  {
    F[1] = velocityLocal[1];
    F[2] = velocityLocal[2];
  }
  else
  {
    num  = mu_i / sqrt(mrn);
    F[1] = (reaction[1] - projT * reaction[0] * num) / at;
    F[2] = (reaction[2] - projS * reaction[0] * num) / at;
  }
  /*   for(int l = 0; l < 3 ; ++l) */
  /*     printf(" %lf", F[l]); */
  /*   printf("\n"); */

}

/* Compute Jacobian of function F */
void jacobianF_AC(int Fsize, double *reaction, double *jacobianFMatrix, int up2Date)
{
  if (jacobianFMatrix == NULL)
  {
    fprintf(stderr, "FrictionContact3D_AlartCurnier::jacobianF_AC error: jacobianMatrix == NULL.\n");
    exit(EXIT_FAILURE);
  }
  if (Fsize != nLocal)
  {
    fprintf(stderr, "FrictionContact3D_AlartCurnier::jacobianF_AC error, wrong block size.\n");
    exit(EXIT_FAILURE);
  }

  /* Warning: input reaction is not necessary equal to the last computed value of reactionBlock */

  /* up2Date = 1 = true if F(n, reaction,F) has been called just before jacobianFMatrix(...). In that case the computation of
     velocityLocal is not required again.
  */
  if (up2Date == 0)
  {
    /* velocityLocal = M.reaction + qLocal */
    int incx = 1, incy = 1;
    DCOPY(Fsize, qLocal, incx, velocityLocal, incy);
    DGEMV(LA_NOTRANS, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, velocityLocal, incy);

    an = 1. / MLocal[0];
    double alpha = MLocal[Fsize + 1] + MLocal[2 * Fsize + 2];
    double det = MLocal[1 * Fsize + 1] * MLocal[2 * Fsize + 2] - MLocal[2 * Fsize + 1] + MLocal[1 * Fsize + 2];
    double beta = alpha * alpha - 4 * det;
    at = 2 * (alpha - beta) / ((alpha + beta) * (alpha + beta));

    /* Projection on [0, +infty[ and on D(0, mu*reactionn) */
    projN = reaction[0] - an * velocityLocal[0];
    projT = reaction[1] - at * velocityLocal[1];
    projS = reaction[2] - at * velocityLocal[2];
  }

  double coef2 = mu_i * mu_i;

  int i, j;
  if (projN > 0)
  {
    for (j = 0; j < Fsize; ++j)
      jacobianFMatrix[j * Fsize] = MLocal[j * Fsize];
  }
  else
  {
    jacobianFMatrix[0] = 1.0 / an;
  }

  double mrn = projT * projT + projS * projS;
  double num, rcof, mrn3;
  double tmp;
  if (mrn <= coef2 * reaction[0]*reaction[0])
    for (i = 1; i < Fsize; ++i)
      for (j = 0; j < Fsize; ++j)
        jacobianFMatrix[j * Fsize + i] = MLocal[j * Fsize + i];
  else
  {
    num  = 1. / sqrt(mrn);
    mrn3 = 1. / sqrt(mrn) * sqrt(mrn) * sqrt(mrn);
    rcof = mu_i / at;
    tmp = at * mrn3 * (MLocal[1] * projT + MLocal[2] * projS);
    jacobianFMatrix[1] = -rcof * (num * projT + reaction[0] * projT * tmp);
    jacobianFMatrix[2] = -rcof * (num * projS + reaction[0] * projS * tmp);

    tmp = mrn3 * ((1 - at * MLocal[Fsize + 1]) * projT - at * MLocal[Fsize + 2] * projS);
    jacobianFMatrix[1 * Fsize + 1] = (1 - mu_i * reaction[0] * (num * (1 - at * MLocal[Fsize + 1]) - projT * tmp)) / at;
    jacobianFMatrix[1 * Fsize + 2] =  - rcof * reaction[0] * ((-num * at * MLocal[Fsize + 2]) - projS * tmp);

    tmp = mrn3 * ((1 - at * MLocal[2 * Fsize + 2]) * projS - at * MLocal[2 * Fsize + 1] * projT);
    jacobianFMatrix[2 * Fsize + 1] =  - rcof * reaction[0] * ((-num * at * MLocal[2 * Fsize + 1]) - projT * tmp);
    jacobianFMatrix[2 * Fsize + 2] = (1 - mu_i * reaction[0] * (num * (1 - at * MLocal[2 * Fsize + 2]) - projS * tmp)) / at;
  }
  /*   for(int l = 0; l < 9 ; ++l) */
  /*     printf(" %lf", jacobianFMatrix[l]); */
  /*   printf("\n"); */
}

void frictionContact3D_AC_post(int contact, double* reaction)
{
  /* This function is required in the interface but useless in Alart-Curnier case */
}

void computeFGlobal_AC(double* reaction, double* FGlobal)
{
  int contact, numberOfContacts = n / 3, diagPos = 0, position;
  int in, it, is, inc, incx;
  double * reactionLocal;
  double alpha, det, beta, num, coef2, mrn;
  for (contact = 0; contact < numberOfContacts; ++contact)
  {
    position = 3 * contact;
    if (MGlobal->storageType == 1) /* Sparse storage */
    {
      /* The part of MGlobal which corresponds to the current block is copied into MLocal */
      diagPos = numberOfContacts * contact + contact;
      MLocal = MGlobal->matrix1->block[diagPos];
    }
    else if (MGlobal->storageType == 0)
    {
      in = 3 * contact;
      it = in + 1;
      is = it + 1;
      inc = n * in;
      double *MM = MGlobal->matrix0;
      /* The part of MM which corresponds to the current block is copied into MLocal */
      MLocal[0] = MM[inc + in];
      MLocal[1] = MM[inc + it];
      MLocal[2] = MM[inc + is];
      inc += n;
      MLocal[3] = MM[inc + in];
      MLocal[4] = MM[inc + it];
      MLocal[5] = MM[inc + is];
      inc += n;
      MLocal[6] = MM[inc + in];
      MLocal[7] = MM[inc + it];
      MLocal[8] = MM[inc + is];
    }

    reactionLocal = &reaction[3 * contact];
    incx = 3;
    velocityLocal[0] = DDOT(3 , MLocal , incx , reactionLocal , 1) + qLocal[0];
    velocityLocal[1] = DDOT(3 , MLocal , incx , reactionLocal , 1) + qLocal[1];
    velocityLocal[2] = DDOT(3 , MLocal , incx , reactionLocal , 1) + qLocal[2];
    an = 1. / MLocal[0];
    alpha = MLocal[4] + MLocal[8];
    det = MLocal[4] * MLocal[8] - MLocal[7] + MLocal[5];
    beta = alpha * alpha - 4 * det;
    at = 2 * (alpha - beta) / ((alpha + beta) * (alpha + beta));
    projN = reactionLocal[0] - an * velocityLocal[0];
    projT = reactionLocal[1] - at * velocityLocal[1];
    projS = reactionLocal[2] - at * velocityLocal[2];
    coef2 = mu[contact] * mu[contact];
    if (projN > 0)
    {
      FGlobal[position] = velocityLocal[0];
    }
    else
    {
      FGlobal[position] = reactionLocal[0] / an;
    }

    mrn = projT * projT + projS * projS;
    if (mrn <= coef2 * reactionLocal[0]*reactionLocal[0])
    {
      FGlobal[position + 1] = velocityLocal[1];
      FGlobal[position + 2] = velocityLocal[2];
    }
    else
    {
      num  = mu[contact] / sqrt(mrn);
      FGlobal[position + 1] = (reactionLocal[1] - projT * reactionLocal[0] * num) / at;
      FGlobal[position + 2] = (reactionLocal[2] - projS * reactionLocal[0] * num) / at;
    }
  }
}

void frictionContact3D_AC_free()
{
  MGlobal = NULL;
  qGlobal = NULL;
  mu = NULL;
  if (isMAllocatedIn == 1)
    free(MLocal);
  MLocal = NULL;
}
