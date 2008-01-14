/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
#include "NSSTools.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/**Static variables */
static int n = 0;
static const double* M = NULL;
static const SparseBlockStructuredMatrix* MB = NULL;
static const double* q = NULL;
static const double* mu = NULL;

static const int nBlock = 3;
static double* MBlock;
static double reactionBlock[3];
static double velocityBlock[3];
static double qLocal[3];
static double mu_i = 0.0;
static double an;
static double at;
static double projN;
static double projT;
static double projS;

/* update pointer to function, used to switch to the function adapted to the storage for M. */
static void (*update)(int, double*);

void updateWithSparse(int contact, double * reaction)
{
  int in = 3 * contact, it = in + 1, is = it + 1;
  int j;
  int numberOfContact = n / 3;
  /* The part of M which corresponds to the current block is copied into MBlock */
  int diagPos = numberOfContact * contact + contact;
  MBlock = MB->block[diagPos];

  /* reactionBlock */
  for (j = 0 ; j < nBlock ; ++j)
    reactionBlock[j] = reaction[in + j];

  /****  Computation of qLocal = qBlock + sum over a row of blocks in M of the products MBlock.reactionBlock,
   excluding the block corresponding to the current contact. ****/

  /* reaction current block set to zero, to exclude current contact block */
  reaction[in] = 0.0;
  reaction[it] = 0.0;
  reaction[is] = 0.0;
  /* qLocal computation*/
  int incx = n, incy = 1;
  qLocal[0] = q[in];
  qLocal[1] = q[it];
  qLocal[2] = q[is];
  /* Loop through the columns(blocks) of M to compute qLocal */
  int blockNum = contact * numberOfContact;
  for (j = 0; j < numberOfContact ; ++j)
  {
    if (j != contact)
    {
      DGEMV(LA_NOTRANS, 3, 3, 1.0, MB->block[blockNum], 3, &reaction[3 * j], incx , 1.0, qLocal, incy);
    }
    blockNum = blockNum + 1;
  }
}

void updateNoSparse(int contact, double * reaction)
{
  int in = 3 * contact, it = in + 1, is = it + 1;
  int j, inc = n * in;

  /* The part of M which corresponds to the current block is copied into MBlock */
  MBlock[0] = M[inc + in];
  MBlock[1] = M[inc + it];
  MBlock[2] = M[inc + is];
  inc += n;
  MBlock[3] = M[inc + in];
  MBlock[4] = M[inc + it];
  MBlock[5] = M[inc + is];
  inc += n;
  MBlock[6] = M[inc + in];
  MBlock[7] = M[inc + it];
  MBlock[8] = M[inc + is];

  /* reactionBlock */
  for (j = 0 ; j < nBlock ; ++j)
    reactionBlock[j] = reaction[in + j];

  /****  Computation of qLocal = qBlock + sum over a contact of blocks in M of the products MBlock.reactionBlock,
   excluding the block corresponding to the current contact. ****/

  /* reaction current block set to zero, to exclude current contact block */
  reaction[in] = 0.0;
  reaction[it] = 0.0;
  reaction[is] = 0.0;
  /* qLocal computation*/
  int incx = n, incy = 1;
  qLocal[0] = q[in] + DDOT(n , &M[in] , incx , reaction , incy);
  qLocal[1] = q[it] + DDOT(n , &M[it] , incx , reaction , incy);
  qLocal[2] = q[is] + DDOT(n , &M[is] , incx , reaction , incy);
}

void initializeSolver_AC(int n0, const double*const M0, const double*const q0, const double*const mu0)
{
  n = n0;
  M = M0;
  q = q0;
  mu = mu0;
  update = &updateNoSparse;
  MBlock = (double*)malloc(nBlock * nBlock * sizeof(*MBlock));
}

void initializeSolver_AC_SB(int n0, const SparseBlockStructuredMatrix*const M0, const double*const q0, const double*const mu0)
{
  n = n0;
  MB = M0;
  q = q0;
  mu = mu0;
  update = &updateWithSparse;
}

void updateSolver_AC(int contact, double * reaction)
{
  /* Call the update function which depends on the storage for M/MB */
  (*update)(contact, reaction);

  /* Friction coefficient for current block*/
  mu_i = mu[contact];
}

/* Compute function F(Reaction) */
void F_AC(int Fsize, double *reactionTmp , double *F, int up2Date)
{
  /* Warning: input reactionTmp is not necessary equal to the last computed value of reactionBlock */

  if (Fsize != nBlock)
  {
    fprintf(stderr, "FrictionContact3D_AlartCurnier::F error, wrong block size.\n");
    exit(EXIT_FAILURE);
  }

  velocityBlock[0] = MBlock[0] * reactionTmp[0] + MBlock[Fsize] * reactionTmp[1] + MBlock[2 * Fsize] * reactionTmp[2] + qLocal[0];
  velocityBlock[1] = MBlock[1] * reactionTmp[0] + MBlock[Fsize + 1] * reactionTmp[1] + MBlock[2 * Fsize + 1] * reactionTmp[2] + qLocal[1];
  velocityBlock[2] = MBlock[2] * reactionTmp[0] + MBlock[Fsize + 2] * reactionTmp[1] + MBlock[2 * Fsize + 2] * reactionTmp[2] + qLocal[2];

  an = 1. / MBlock[0];
  double alpha = MBlock[Fsize + 1] + MBlock[2 * Fsize + 2];
  double det = MBlock[1 * Fsize + 1] * MBlock[2 * Fsize + 2] - MBlock[2 * Fsize + 1] + MBlock[1 * Fsize + 2];
  double beta = alpha * alpha - 4 * det;
  at = 2 * (alpha - beta) / ((alpha + beta) * (alpha + beta));

  double num;
  double coef2 = mu_i * mu_i;
  /* Projection on [0, +infty[ and on D(0, mu*reactionn) */
  projN = reactionTmp[0] - an * velocityBlock[0];
  projT = reactionTmp[1] - at * velocityBlock[1];
  projS = reactionTmp[2] - at * velocityBlock[2];

  if (projN > 0)
  {
    F[0] = velocityBlock[0];
  }
  else
  {
    F[0] = reactionTmp[0] / an;
  }

  double mrn = projT * projT + projS * projS;
  if (mrn <= coef2 * reactionTmp[0]*reactionTmp[0])
  {
    F[1] = velocityBlock[1];
    F[2] = velocityBlock[2];
  }
  else
  {
    num  = mu_i / sqrt(mrn);
    F[1] = (reactionTmp[1] - projT * reactionTmp[0] * num) / at;
    F[2] = (reactionTmp[2] - projS * reactionTmp[0] * num) / at;
  }
}

/* Compute Jacobian of function F */
void jacobianF_AC(int Fsize, double *reactionTmp, double *jacobianFMatrix, int up2Date)
{
  if (jacobianFMatrix == NULL)
  {
    fprintf(stderr, "FrictionContact3D_AlartCurnier::jacobianF_AC error: jacobianMatrix == NULL.\n");
    exit(EXIT_FAILURE);
  }
  if (Fsize != nBlock)
  {
    fprintf(stderr, "FrictionContact3D_AlartCurnier::jacobianF_AC error, wrong block size.\n");
    exit(EXIT_FAILURE);
  }

  /* Warning: input reactionTmp is not necessary equal to the last computed value of reactionBlock */

  /* up2Date = 1 = true if F(n, reactionTmp,F) has been called just before jacobianFMatrix(...). In that case the computation of
     velocityBlock is not required again.
  */
  if (up2Date == 0)
  {
    velocityBlock[0] = MBlock[0] * reactionTmp[0] + MBlock[Fsize] * reactionTmp[1] + MBlock[2 * Fsize] * reactionTmp[2] + qLocal[0];
    velocityBlock[1] = MBlock[1] * reactionTmp[0] + MBlock[Fsize + 1] * reactionTmp[1] + MBlock[2 * Fsize + 1] * reactionTmp[2] + qLocal[1];
    velocityBlock[2] = MBlock[2] * reactionTmp[0] + MBlock[Fsize + 2] * reactionTmp[1] + MBlock[2 * Fsize + 2] * reactionTmp[2] + qLocal[2];

    an = 1. / MBlock[0];
    double alpha = MBlock[Fsize + 1] + MBlock[2 * Fsize + 2];
    double det = MBlock[1 * Fsize + 1] * MBlock[2 * Fsize + 2] - MBlock[2 * Fsize + 1] + MBlock[1 * Fsize + 2];
    double beta = alpha * alpha - 4 * det;
    at = 2 * (alpha - beta) / ((alpha + beta) * (alpha + beta));

    /* Projection on [0, +infty[ and on D(0, mu*reactionn) */
    projN = reactionTmp[0] - an * velocityBlock[0];
    projT = reactionTmp[1] - at * velocityBlock[1];
    projS = reactionTmp[2] - at * velocityBlock[2];
  }

  double coef2 = mu_i * mu_i;

  int i, j;
  if (projN > 0)
  {
    for (j = 0; j < Fsize; ++j)
      jacobianFMatrix[j * Fsize] = MBlock[j * Fsize];
  }
  else
  {
    jacobianFMatrix[0] = 1.0 / an;
  }

  double mrn = projT * projT + projS * projS;
  double num, rcof, mrn3;
  double tmp;
  if (mrn <= coef2 * reactionTmp[0]*reactionTmp[0])
    for (i = 1; i < Fsize; ++i)
      for (j = 0; j < Fsize; ++j)
        jacobianFMatrix[j * Fsize + i] = MBlock[j * Fsize + i];
  else
  {
    num  = 1. / sqrt(mrn);
    mrn3 = 1. / sqrt(mrn) * sqrt(mrn) * sqrt(mrn);
    rcof = mu_i / at;
    tmp = at * mrn3 * (MBlock[1] * projT + MBlock[2] * projS);
    jacobianFMatrix[1] = -rcof * (num * projT + reactionTmp[0] * projT * tmp);
    jacobianFMatrix[2] = -rcof * (num * projS + reactionTmp[0] * projS * tmp);

    tmp = mrn3 * ((1 - at * MBlock[Fsize + 1]) * projT - at * MBlock[Fsize + 2] * projS);
    jacobianFMatrix[1 * Fsize + 1] = (1 - mu_i * reactionTmp[0] * (num * (1 - at * MBlock[Fsize + 1]) - projT * tmp)) / at;
    jacobianFMatrix[1 * Fsize + 2] =  - rcof * reactionTmp[0] * ((-num * at * MBlock[Fsize + 2]) - projS * tmp);

    tmp = mrn3 * ((1 - at * MBlock[2 * Fsize + 2]) * projS - at * MBlock[2 * Fsize + 1] * projT);
    jacobianFMatrix[2 * Fsize + 1] =  - rcof * reactionTmp[0] * ((-num * at * MBlock[2 * Fsize + 1]) - projT * tmp);
    jacobianFMatrix[2 * Fsize + 2] = (1 - mu_i * reactionTmp[0] * (num * (1 - at * MBlock[2 * Fsize + 2]) - projS * tmp)) / at;
  }
}

void postSolver_AC(int contact, double* reaction)
{
  /* This function is required in the interface but useless in Alart-Curnier case */
}

void computeFGlobal_AC(double* reaction, double* FGlobal)
{
  int contact, numberOfContacts = n / 3, diagPos = 0, position;
  int in, it, is, inc, incx;
  double * reactionTmp;
  double alpha, det, beta, num, coef2, mrn;
  for (contact = 0; contact < numberOfContacts; ++contact)
  {
    position = 3 * contact;
    if (MB != NULL) /* Sparse storage */
    {
      /* The part of M which corresponds to the current block is copied into MBlock */
      diagPos = numberOfContacts * contact + contact;
      MBlock = MB->block[diagPos];
    }
    else if (M != NULL)
    {
      in = 3 * contact;
      it = in + 1;
      is = it + 1;
      inc = n * in;

      /* The part of M which corresponds to the current block is copied into MBlock */
      MBlock[0] = M[inc + in];
      MBlock[1] = M[inc + it];
      MBlock[2] = M[inc + is];
      inc += n;
      MBlock[3] = M[inc + in];
      MBlock[4] = M[inc + it];
      MBlock[5] = M[inc + is];
      inc += n;
      MBlock[6] = M[inc + in];
      MBlock[7] = M[inc + it];
      MBlock[8] = M[inc + is];
    }

    reactionTmp = &reaction[3 * contact];
    incx = 3;
    velocityBlock[0] = DDOT(3 , MBlock , incx , reactionTmp , 1) + qLocal[0];
    velocityBlock[1] = DDOT(3 , MBlock , incx , reactionTmp , 1) + qLocal[1];
    velocityBlock[2] = DDOT(3 , MBlock , incx , reactionTmp , 1) + qLocal[2];
    an = 1. / MBlock[0];
    alpha = MBlock[4] + MBlock[8];
    det = MBlock[4] * MBlock[8] - MBlock[7] + MBlock[5];
    beta = alpha * alpha - 4 * det;
    at = 2 * (alpha - beta) / ((alpha + beta) * (alpha + beta));
    projN = reactionTmp[0] - an * velocityBlock[0];
    projT = reactionTmp[1] - at * velocityBlock[1];
    projS = reactionTmp[2] - at * velocityBlock[2];
    coef2 = mu[contact] * mu[contact];
    if (projN > 0)
    {
      FGlobal[position] = velocityBlock[0];
    }
    else
    {
      FGlobal[position] = reactionTmp[0] / an;
    }

    mrn = projT * projT + projS * projS;
    if (mrn <= coef2 * reactionTmp[0]*reactionTmp[0])
    {
      FGlobal[position + 1] = velocityBlock[1];
      FGlobal[position + 2] = velocityBlock[2];
    }
    else
    {
      num  = mu[contact] / sqrt(mrn);
      FGlobal[position + 1] = (reactionTmp[1] - projT * reactionTmp[0] * num) / at;
      FGlobal[position + 2] = (reactionTmp[2] - projS * reactionTmp[0] * num) / at;
    }
  }
}
