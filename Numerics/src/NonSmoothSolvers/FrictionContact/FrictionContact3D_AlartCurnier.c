/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#include "LA.h"
#include "FrictionContact3D_Solvers.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#define VERBOSE_DEBUG

//#define OPTI_RHO

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
    double det = MLocal[1 * Fsize + 1] * MLocal[2 * Fsize + 2] - MLocal[2 * Fsize + 1] * MLocal[1 * Fsize + 2];
    double beta = alpha * alpha - 4 * det;
    if (beta > 0.)
      beta = sqrt(beta);
    else
      beta = 0.;
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
    double det = MLocal[1 * Fsize + 1] * MLocal[2 * Fsize + 2] - MLocal[2 * Fsize + 1] * MLocal[1 * Fsize + 2];
    double beta = alpha * alpha - 4 * det;
    if (beta > 0.)
      beta = sqrt(beta);
    else
      beta = 0.;
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


int OldAlartCurnierNewton(int Fsize, double * reactionBlock, int *iparam, double *dparam)
{

  int ITMAX = iparam[0];
  int i, j, inew, INEWTON;
  int STATUS;
  double MU = mu_i;
  double TOL = dparam[0];
  double RVN, RVT, RVS, RV, RV1, RV3;
  double PhiT4, PhiS4, Phi3;
  double Threshold, DET, DET1;
  double d1, d2, d3, AA, BB, CC;
  double ERW, WREF;

  double Rloci[3], WRloc[3], WRloci[3], Rerror[3], Verror[3], Rloc[3], S[3];
  double A[3][3], AWB[3][3], W[3][3];
  Rloc[0] = reactionBlock[1];
  S[0] = qLocal[1];
  Rloc[1] = reactionBlock[2];
  S[1] = qLocal[2];
  Rloc[2] = reactionBlock[0];
  S[2] = qLocal[0];
  W[2][2] =  MLocal[0 + 0 * 3];
  W[2][1] =  MLocal[0 + 2 * 3];
  W[2][0] =  MLocal[0 + 1 * 3];
  W[1][2] =  MLocal[2 + 0 * 3];
  W[1][1] =  MLocal[2 + 2 * 3];
  W[1][0] =  MLocal[2 + 1 * 3];
  W[0][2] =  MLocal[1 + 0 * 3];
  W[0][1] =  MLocal[1 + 2 * 3];
  W[0][0] =  MLocal[1 + 1 * 3];
#ifdef VERBOSE_DEBUG
  printf("Initial conditions \n");
  printf("RN = %le\n", Rloc[2]);
  printf("RT = %le\n", Rloc[0]);
  printf("RS = %le\n", Rloc[1]);

  printf("WNN = %le\n", W[2][2]);
  printf("WNT = %le\n", W[2][0]);
  printf("WNS = %le\n", W[2][1]);
  printf("WTT = %le\n", W[0][0]);
  printf("WTS = %le\n", W[0][1]);
  printf("WSS = %le\n", W[1][1]);
#endif

  double RhoN;
  double RhoT;
  double SW;
  double DW;

  SW = W[0][0] + W[1][1];

  DW = SW * SW - 4.0 * (W[0][0] * W[1][1] - W[1][0] * W[0][1]);

  if (DW > 0.0) DW = sqrt(DW);
  else DW = 0.0;

  RhoN = 1.0 / W[2][2];
  RhoT = 2.0 * (SW - DW) / ((SW + DW) * (SW + DW));

  for (i = 0 ; i < 3 ; ++i)  WRloc[i] = W[i][0] * Rloc[0] + W[i][1] * Rloc[1] + W[i][2] * Rloc[2];


  for (inew = 0 ; inew < ITMAX ; ++inew)
  {
#ifdef VERBOSE_DEBUG
    printf("-----------------------------------    AlartCurnierNewton iteration = %i \n", inew);
#endif
    for (i = 0 ; i < 3 ; ++i)
    {
      WRloci[i] =  WRloc[i];
      Rloci[i]  =  Rloc[i];
      for (j = 0 ; j < 3 ; ++j)
      {
        AWB[i][j] = 0.0;
        A[i][j] = 0.0;
      }
    }
    /* On cherche à exprimer la matrice gradient:
     * Dphi=| I  -W  |
     *      | Ap  Bp |
     * p étant l'indice inew.
     * On effectue tour à tour les projections sur le cone
     * positif et sur le cone de Coulomb pour exprimer les différentes
     * valeurs des matrices Ap et Bp.
     */

    RVN = Rloc[2] - RhoN * (WRloc[2] + S[2]);

    if (RVN >= 0.0)
    {

      /* Projection sur le cone positif */

      Phi3 = RhoN * (WRloc[2] + S[2]);

      for (i = 0 ; i < 3 ; ++i) AWB[2][i] = RhoN * W[2][i];
      STATUS = 1;
    }
    else
    {
      Phi3 = Rloc[2];
      AWB[2][2] = 1.0;
      AWB[2][1] = 0.0;
      AWB[2][0] = 0.0;
      STATUS = 0;
    }
#ifdef VERBOSE_DEBUG
    printf("STATUS =%i\n", STATUS);
    printf("Phi3 =%le\n", Phi3);
#endif
    RVT = Rloc[1] - RhoT * (WRloc[1] + S[1]);
    RVS = Rloc[0] - RhoT * (WRloc[0] + S[0]);

    RV = sqrt(RVT * RVT + RVS * RVS);

    Threshold = MU * Rloc[2];

    if ((RV < Threshold) || (RV < 1.e-20))
    {

      PhiT4 = RhoT * (WRloc[1] + S[1]);
      PhiS4 = RhoT * (WRloc[0] + S[0]);

      for (i = 0 ; i < 2 ; ++i)
        for (j = 0 ; j < 3 ; ++j)
          AWB[i][j] = RhoT * W[i][j];
      STATUS = 1;
    }
    else
    {

      RV1 = 1.0 / RV;
      PhiT4 = Rloc[1] - Threshold * RVT * RV1;
      PhiS4 = Rloc[0] - Threshold * RVS * RV1;

      AWB[0][2] = -MU * RVS * RV1;
      AWB[1][2] = -MU * RVT * RV1;

      RV3 = RV1 * RV1 * RV1;

      A[1][1] = Threshold * (RV - RVS * RVS) * RV3;
      A[1][0] = -Threshold * RV3 * RVT * RVS;
      A[0][1] = A[1][0];
      A[0][0] = Threshold * (RV - RVT * RVT) * RV3;

      for (i = 0 ; i < 2 ; ++i) AWB[i][i] = 1.0 - A[i][i];
      AWB[1][0] = - A[1][0];
      AWB[0][1] = - A[0][1];

      for (i = 0 ; i < 2 ; ++i)
        for (j = 0 ; j < 2 ; ++j)
          A[i][j] = RhoT * A[i][j];

      for (i = 0 ; i < 2 ; ++i)
        for (j = 0 ; j < 2 ; ++j) AWB[i][j] += A[i][0] * W[0][j] + A[i][1] * W[1][j];
      STATUS = 2;
    }
#ifdef VERBOSE_DEBUG
    printf("STATUS =%i\n", STATUS);
    printf("PhiT4 =%le\n", PhiT4);
    printf("PhiS4 =%le\n", PhiS4);
#endif

    //
    d1 = AWB[1][1] * AWB[2][2] - AWB[1][2] * AWB[2][1];
    d2 = AWB[1][0] * AWB[2][2] - AWB[2][0] * AWB[1][2];
    d3 = AWB[1][0] * AWB[2][1] - AWB[2][0] * AWB[1][1];

    DET = AWB[0][0] * d1 - AWB[0][1] * d2 + AWB[0][2] * d3;

    if (fabs(DET) < 1.e-20)
    {
      printf("DET NULL\n");
    }
#ifdef VERBOSE_DEBUG
    printf("AWB[%i][%i] = %le\t", 0, 0, AWB[2][2]);
    printf("AWB[%i][%i] = %le\t", 0, 1, AWB[2][0]);
    printf("AWB[%i][%i] = %le\t", 0, 2, AWB[2][1]);
    printf("\n");
    printf("AWB[%i][%i] = %le\t", 1, 0, AWB[0][2]);
    printf("AWB[%i][%i] = %le\t", 1, 1, AWB[0][0]);
    printf("AWB[%i][%i] = %le\t", 1, 2, AWB[0][1]);
    printf("\n");
    printf("AWB[%i][%i] = %le\t", 2, 0, AWB[1][2]);
    printf("AWB[%i][%i] = %le\t", 2, 1, AWB[1][0]);
    printf("AWB[%i][%i] = %le\t", 2, 2, AWB[1][1]);
    printf("\n");
#endif
    DET1 = 1.0 / DET;

    AA
      = (AWB[0][1] * AWB[1][2] - AWB[1][1] * AWB[0][2]) * PhiS4
        - (AWB[0][0] * AWB[1][2] - AWB[1][0] * AWB[0][2]) * PhiT4
        + (AWB[0][0] * AWB[1][1] - AWB[1][0] * AWB[0][1]) * Phi3;

    //  Rloc[2] = fmax(0.0 , Rloc[2] - AA*DET1);
    Rloc[2] =  Rloc[2] - AA * DET1;

    BB
      = -(AWB[0][1] * AWB[2][2] - AWB[2][1] * AWB[0][2]) * PhiS4
        + (AWB[0][0] * AWB[2][2] - AWB[2][0] * AWB[0][2]) * PhiT4
        - (AWB[0][0] * AWB[2][1] - AWB[2][0] * AWB[0][1]) * Phi3;

    Rloc[1] = Rloc[1] - BB * DET1;

    CC = d1 * PhiS4 - d2 * PhiT4 + d3 * Phi3;

    Rloc[0] = Rloc[0] - CC * DET1;
#ifdef VERBOSE_DEBUG
    printf("RN = %le\n", Rloc[2]);
    printf("RT = %le\n", Rloc[0]);
    printf("RS = %le\n", Rloc[1]);
#endif
    for (i = 0 ; i < 3 ; ++i) WRloc[i] = W[i][0] * Rloc[0] + W[i][1] * Rloc[1] + W[i][2] * Rloc[2];

    WREF = 0.0;
    ERW  = 0.0;

    for (i = 0 ; i < 3 ; ++i)
    {
      Rerror[i] = Rloci[i] - Rloc[i];
      Verror[i] = WRloci[i] - WRloc[i];
      ERW  += Rerror[i] * Verror[i];
      WREF += WRloc[i] * Rloci[i];
    }

    if (WREF < 1.e-20) WREF = 1.0;
    ERW = ERW / WREF;

    INEWTON = inew + 1;
    if (fabs(ERW) < TOL)  break;
  }
  if (verbose > 1)
  {
    printf("-----------------------------------    AlartCurnierNewton number of iteration = %i  error = %.10e \n", INEWTON, ERW);
  }

  dparam[1] = ERW;
  reactionBlock[0] = Rloc[2];
  reactionBlock[1] = Rloc[0];
  reactionBlock[2] = Rloc[1];
  if (fabs(ERW) > TOL * TOL)  return 1;
  else return 0;


}
int frictionContact3D_AlartCurnierNewton_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the NSGS Solver\n");
  }

  strcpy(options->solverName, "AlartCurnierNewton");
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 10;
  options->dparam[0] = 1e-16;
  return 0;
}

int TruncatedCylinderAlartCurnierNewton(int Fsize, double * R, int *iparam, double *dparam)
/* int AlartCurnierNewton(int Fsize, double * R, int *iparam, double *dparam) */
{

  double mu = mu_i;
  double Tol = dparam[0];
  double itermax = iparam[0];
  int i, j, k, inew;
  double velocity[3];
  double A[3][3];
  A[0][1] = 0.0;
  A[0][2] = 0.0;
  A[1][0] = 0.0;
  A[2][0] = 0.0;
  double B[3][3];
  B[0][1] = 0.0;
  B[0][2] = 0.0;


  double AWplusB[3][3];


  double RhoN = 1.0;
  double RhoT = 1.0;
#ifdef OPTI_RHO
  double sw = MLocal[1 + 1 * 3] + MLocal[2 + 2 * 3];

  double dw = sw * sw - 4.0 * (MLocal[1 + 1 * 3] + MLocal[2 + 2 * 3] -  MLocal[2 + 1 * 3] + MLocal[1 + 2 * 3]);

  if (dw > 0.0) dw = sqrt(dw);
  else dw = 0.0;

  RhoN = 1.0 / MLocal[0 + 0 * 3];
  RhoT = 2.0 * (sw - dw) / ((sw + dw) * (sw + dw));
#ifdef VERBOSE_DEBUG
  printf("sw=%le = ", sw);
  printf("dw=%le = ", dw);
  printf("RhoN=%le = ", RhoN);
  printf("RhoT=%le\n", RhoT);
  printf("\n");
#endif
#endif

  for (i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * R[0] + qLocal[i]
                                          + MLocal[i + 1 * 3] * R[1]
                                          + MLocal[i + 2 * 3] * R[2] ;
#ifdef VERBOSE_DEBUG
  printf("velocity =\t");
  for (i = 0; i < 3; i++) printf("velocity[%i]= %le\t", i, velocity[i]);
  printf("\n");
#endif

  double RVN, RVT, RVS, RV, RV1, RV3, DET, DET1;
  double ResN, ResT, ResS, Radius;
  double d1, d2, d3, AA, BB, CC;
  double GammaTT, GammaTS, GammaST, GammaSS;

  RVN = R[0] - RhoN * velocity[0];
  RVT = R[1] - RhoT * velocity[1];
  RVS = R[2] - RhoT * velocity[2];
  RV = sqrt(RVT * RVT + RVS * RVS);
  Radius = mu * R[0];
#ifdef VERBOSE_DEBUG
  printf("RN = %le\t", R[0]);
  printf("RT = %le\t", R[1]);
  printf("RS = %le\n", R[2]);
#endif
  d1 = MLocal[1 + 1 * 3] * MLocal[2 + 2 * 3] - MLocal[1 + 2 * 3] * MLocal[2 + 1 * 3];
  d2 = MLocal[1 + 0 * 3] * MLocal[2 + 2 * 3] - MLocal[2 + 0 * 3] * MLocal[1 + 2 * 3];
  d3 = MLocal[1 + 0 * 3] * MLocal[2 + 1 * 3] - MLocal[2 + 0 * 3] * MLocal[1 + 1 * 3];

  DET = MLocal[0 + 0 * 3] * d1 - MLocal[0 + 1 * 3] * d2 + MLocal[0 + 2 * 3] * d3;
#ifdef VERBOSE_DEBUG
  printf("DET = %le\n", DET);
#endif
  if (fabs(DET) < 1.e-20)
  {
    printf("DET NULL\n");
    exit(EXIT_FAILURE);
  }

  for (inew = 0 ; inew < itermax ; ++inew)
  {



    // Compute the value of the Alart--Curnier Function and its gradient for the normal part


    if (RVN >= 0.0)
    {
#ifdef VERBOSE_DEBUG
      printf("Normal part in the cone\n");
#endif
      ResN = RhoN * (velocity[0]);
      A[0][0] =  RhoN;
      B[0][0] = 0.0;
    }
    else
    {
#ifdef VERBOSE_DEBUG
      printf("Normal part out the cone\n");
#endif

      ResN = R[0];
      A[0][0] = 0.0;
      B[0][0] = 1.0;
    }

    // Compute the value of the Alart--Curnier Function and its gradient for the tangential part


#ifdef VERBOSE_DEBUG
    printf("Radius=%le\n", Radius);
    printf("RV=%le\n", RV);
#endif
    if (RV < Radius) // We are in the disk and Radius is positive
    {
#ifdef VERBOSE_DEBUG
      printf("We are in the disk\n");
#endif
      ResT = RhoT * (velocity[1]);
      ResS = RhoT * (velocity[2]);
      A[1][1] = RhoT;
      A[1][2] = 0.0;
      A[2][1] = 0.0;
      A[2][2] = RhoT;
      B[1][0] = 0.0;
      B[1][1] = 0.0;
      B[1][2] = 0.0;
      B[2][0] = 0.0;
      B[2][1] = 0.0;
      B[2][2] = 0.0;
    }
    else if (RV >= Radius && Radius > 0) // We are out the disk and Radius is postive
    {
#ifdef VERBOSE_DEBUG
      printf("We are out the disk and Radius is positive\n");
#endif
      RV1 = 1.0 / RV;
      ResT = R[1] - Radius * RVT * RV1;
      ResS = R[2] - Radius * RVS * RV1;

      RV3 = RV1 * RV1 * RV1;
      GammaTT = RV1 - RVT * RVT * RV3;
      GammaTS =  - RVT * RVS * RV3;
      GammaST =  GammaTS;
      GammaSS = RV1 - RVS * RVS * RV3;

      A[1][1] = GammaTT * RhoT * Radius;

      A[1][2] = GammaTS * RhoT * Radius;
      A[2][1] = GammaST * RhoT * Radius;

      A[2][2] = GammaSS * RhoT * Radius;

      B[1][0] = -mu * RVT * RV1;

      B[1][1] = 1.0 - GammaTT * Radius ;
      B[1][2] = - GammaTS  * Radius ;

      B[2][0] = -mu * RVS * RV1;

      B[2][1] = - GammaST  * Radius;
      B[2][2] = 1.0 - GammaSS * Radius;
    }
    else // We are out the disk and Radius is negative
    {
#ifdef VERBOSE_DEBUG
      printf("We are out the disk and Radius is negative\n");
#endif

      /*Version original */
      ResT = R[1] ;
      ResS = R[2] ;
      A[1][1] = 0.0;
      A[1][2] = 0.0;
      A[2][1] = 0.0;
      A[2][2] = 0.0;

      B[1][0] = 0.0;
      B[1][1] = 1.0;
      B[1][2] = 0.0;
      B[2][0] = 0.0;
      B[2][1] = 0.0;
      B[2][2] = 1.0;
    }

    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
      {
        AWplusB[i][j] = 0.0;
        for (k = 0; k < 3; k++)
        {
          AWplusB[i][j] += A[i][k] * MLocal[k + j * 3];
        }
        AWplusB[i][j] += B[i][j];
      }
    }
#ifdef VERBOSE_DEBUG
    for (int l = 0; l < 3; l++)
    {
      for (int k = 0; k < 3; k++)
      {
        printf("A[%i][%i] = %le\t", l, k, A[l][k]);
      }
      printf("\n");
    }
    for (int l = 0; l < 3; l++)
    {
      for (int k = 0; k < 3; k++)
      {
        printf("B[%i][%i] = %le\t", l, k, B[l][k]);
      }
      printf("\n");
    }
    double diago = 0.0;
    for (int l = 0; l < 3; l++)
    {
      for (int k = 0; k < 3; k++)
      {
        if (k == l)  diago = 1.0;
        else diago = 0.0;
        printf("I-B[%i][%i] = %le\t", l, k, diago - B[l][k]);
      }
      printf("\n");
    }
    for (int l = 0; l < 3; l++)
    {
      for (int k = 0; k < 3; k++)
      {
        printf("AWplusB[%i][%i] = %le\t", l, k, AWplusB[l][k]);
      }
      printf("\n");
    }
    printf("ResN = %le\t", ResN);
    printf("ResT = %le\t", ResT);
    printf("ResS = %le\n", ResS);
#endif
    d1 = AWplusB[1][1] * AWplusB[2][2] - AWplusB[1][2] * AWplusB[2][1];
    d2 = AWplusB[1][0] * AWplusB[2][2] - AWplusB[2][0] * AWplusB[1][2];
    d3 = AWplusB[1][0] * AWplusB[2][1] - AWplusB[2][0] * AWplusB[1][1];

    DET = AWplusB[0][0] * d1 - AWplusB[0][1] * d2 + AWplusB[0][2] * d3;
#ifdef VERBOSE_DEBUG
    printf("DET = %le\n", DET);
#endif
    if (fabs(DET) < 1.e-20)
    {
      printf("DET NULL\n");
      exit(EXIT_FAILURE);
    }

    DET1 = 1.0 / DET;
    AA = +(AWplusB[1][1] * AWplusB[2][2] - AWplusB[2][1] * AWplusB[1][2]) * ResN
         - (AWplusB[0][1] * AWplusB[2][2] - AWplusB[2][1] * AWplusB[0][2]) * ResT
         + (AWplusB[0][1] * AWplusB[1][2] - AWplusB[1][1] * AWplusB[0][2]) * ResS ;

    BB = -(AWplusB[1][0] * AWplusB[2][2] - AWplusB[2][0] * AWplusB[1][2]) * ResN
         + (AWplusB[0][0] * AWplusB[2][2] - AWplusB[2][0] * AWplusB[0][2]) * ResT
         - (AWplusB[0][0] * AWplusB[1][2] - AWplusB[1][0] * AWplusB[0][2]) * ResS ;


    CC = +(AWplusB[1][0] * AWplusB[2][1] - AWplusB[2][0] * AWplusB[1][1]) * ResN
         - (AWplusB[0][0] * AWplusB[2][1] - AWplusB[2][0] * AWplusB[0][1]) * ResT
         + (AWplusB[0][0] * AWplusB[1][1] - AWplusB[1][0] * AWplusB[0][1]) * ResS ;

    R[0] = R[0] - AA * DET1;
    R[1] = R[1] - BB * DET1;
    R[2] = R[2] - CC * DET1;
#ifdef VERBOSE_DEBUG
    printf("RN = %le\n", R[0]);
    printf("RT = %le\n", R[1]);
    printf("RS = %le\n", R[2]);
#endif
    // compute new residue

    for (i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * R[0] + qLocal[i]
                                            + MLocal[i + 1 * 3] * R[1] +
                                            + MLocal[i + 2 * 3] * R[2] ;


    RVN = R[0] - RhoN * velocity[0];
    RVT = R[1] - RhoT * velocity[1];
    RVS = R[2] - RhoT * velocity[2];
    RV = sqrt(RVT * RVT + RVS * RVS);
    Radius = mu * R[0];
    if (RVN >= 0.0)
    {
      ResN = RhoN * (velocity[0]);
    }
    else
    {
      ResN = R[0];
    }
    if (RV < Radius) // We are in the disk and Radius is postive
    {
      ResT = RhoT * (velocity[1]);
      ResS = RhoT * (velocity[2]);
    }
    else if (RV >= Radius && Radius > 0) // We are out the disk and Radius is postive
    {
      RV1 = 1.0 / RV;
      ResT = R[1] - Radius * RVT * RV1;
      ResS = R[2] - Radius * RVS * RV1;
    }
    else // We are out the disk and Radius is negative
    {
      ResT = R[1] ;
      ResS = R[2] ;
    }

    dparam[1] = 0.5 * (ResN * ResN + ResT * ResT + ResS * ResS);


    if (verbose > 1)
    {
      printf("-----------------------------------    AlartCurnierNewton number of iteration = %i  error = %.10e \n", inew, dparam[1]);
    }
    if (dparam[1] < Tol)
    {
      return 0;

    }


  }


  return 1;



}

/* int CompleteCylinderAlartCurnierNewton(int Fsize, double * R, int *iparam, double *dparam) */
int AlartCurnierNewton(int Fsize, double * R, int *iparam, double *dparam)
{

  double mu = mu_i;
  double Tol = dparam[0];
  double itermax = iparam[0];
  int i, j, k, inew;
  double velocity[3];
  double A[3][3];
  A[0][1] = 0.0;
  A[0][2] = 0.0;
  A[1][0] = 0.0;
  A[2][0] = 0.0;
  double B[3][3];
  B[0][1] = 0.0;
  B[0][2] = 0.0;

  double AWplusB[3][3];


  double RhoN = 1.0;
  double RhoT = 1.0;
#ifdef OPTI_RHO
  double sw = MLocal[1 + 1 * 3] + MLocal[2 + 2 * 3];

  double dw = sw * sw - 4.0 * (MLocal[1 + 1 * 3] + MLocal[2 + 2 * 3] -  MLocal[2 + 1 * 3] + MLocal[1 + 2 * 3]);

  if (dw > 0.0) dw = sqrt(dw);
  else dw = 0.0;

  RhoN = 1.0 / MLocal[0 + 0 * 3];
  RhoT = 2.0 * (sw - dw) / ((sw + dw) * (sw + dw));
#ifdef VERBOSE_DEBUG
  printf("sw=%le = ", sw);
  printf("dw=%le = ", dw);
  printf("RhoN=%le = ", RhoN);
  printf("RhoT=%le\n", RhoT);
  printf("\n");
#endif
#endif
  for (i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * R[0] + qLocal[i]
                                          + MLocal[i + 1 * 3] * R[1]
                                          + MLocal[i + 2 * 3] * R[2] ;


  double RVN, RVT, RVS, RV, RV1, RV3, DET, DET1;
  double ResN, ResT, ResS, Radius;
  double d1, d2, d3, AA, BB, CC;
  double GammaTT, GammaTS, GammaST, GammaSS;

  RVN = R[0] - RhoN * velocity[0];
  RVT = R[1] - RhoT * velocity[1];
  RVS = R[2] - RhoT * velocity[2];
  RV = sqrt(RVT * RVT + RVS * RVS);
  Radius = mu * R[0];
#ifdef VERBOSE_DEBUG
  printf("RN = %le\n", R[0]);
  printf("RT = %le\n", R[1]);
  printf("RS = %le\n", R[2]);
#endif
  d1 = MLocal[1 + 1 * 3] * MLocal[2 + 2 * 3] - MLocal[1 + 2 * 3] * MLocal[2 + 1 * 3];
  d2 = MLocal[1 + 0 * 3] * MLocal[2 + 2 * 3] - MLocal[2 + 0 * 3] * MLocal[1 + 2 * 3];
  d3 = MLocal[1 + 0 * 3] * MLocal[2 + 1 * 3] - MLocal[2 + 0 * 3] * MLocal[1 + 1 * 3];

  DET = MLocal[0 + 0 * 3] * d1 - MLocal[0 + 1 * 3] * d2 + MLocal[0 + 2 * 3] * d3;
#ifdef VERBOSE_DEBUG
  printf("DET = %le\n", DET);
#endif
  if (fabs(DET) < 1.e-20)
  {
    printf("DET NULL\n");
    exit(EXIT_FAILURE);
  }

  for (inew = 0 ; inew < itermax ; ++inew)
  {



    // Compute the value of the Alart--Curnier Function and its gradient for the normal part

    if (RVN >= 0.0)
    {
      ResN = RhoN * (velocity[0]);
      A[0][0] =  RhoN;
      B[0][0] = 0.0;
    }
    else
    {
      ResN = R[0];
      A[0][0] = 0.0;
      B[0][0] = 1.0;
    }

    // Compute the value of the Alart--Curnier Function and its gradient for the tangential part


#ifdef VERBOSE_DEBUG
    printf("Radius=%le\n", Radius);
    printf("RV=%le\n", RV);
#endif
    if (RV < Radius || RV < 1e-20)  // We are in the disk
    {
#ifdef VERBOSE_DEBUG
      printf("We are in the disk \n");
#endif
      ResT = RhoT * (velocity[1]);
      ResS = RhoT * (velocity[2]);
      A[1][1] = RhoT;
      A[1][2] = 0.0;
      A[2][1] = 0.0;
      A[2][2] = RhoT;
      B[1][0] = 0.0;
      B[1][1] = 0.0;
      B[1][2] = 0.0;
      B[2][0] = 0.0;
      B[2][1] = 0.0;
      B[2][2] = 0.0;
    }
    else  // We are out the disk
    {
#ifdef VERBOSE_DEBUG
      printf("We are out the disk\n");
#endif
      /*        RV1 = 1.0/RV; */
      /*        ResT = R[1] - Radius*RVT*RV1; */
      /*        ResS = R[2] - Radius*RVS*RV1; */


      /*        RV3 = RV1*RV1*RV1; */
      /*        GammaTT = (RV - RVT*RVT)*RV3; */
      /*        GammaTS =  - RVT*RVS*RV3; */
      /*        GammaST =  GammaTS; */
      /*        GammaSS = (RV - RVS*RVS)*RV3; */


      RV1 = 1.0 / RV;
      ResT = R[1] - Radius * RVT * RV1;
      ResS = R[2] - Radius * RVS * RV1;

      RV3 = RV1 * RV1 * RV1;
      GammaTT = RV1 - RVT * RVT * RV3;
      GammaTS =  - RVT * RVS * RV3;
      GammaST =  GammaTS;
      GammaSS = RV1 - RVS * RVS * RV3;


      A[1][1] = GammaTT * RhoT * Radius;

      A[1][2] = GammaTS * RhoT * Radius;

      A[2][1] = GammaST * RhoT * Radius;
      A[2][2] = GammaSS * RhoT * Radius;

      B[1][0] = -mu * RVT * RV1;

      B[1][1] = 1.0 - GammaTT * Radius ;
      B[1][2] = - GammaTS * Radius ;

      B[2][0] = -mu * RVS * RV1;

      B[2][1] = - GammaST * Radius;
      B[2][2] = 1.0 - GammaSS * Radius;
    }


    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
      {
        AWplusB[i][j] = 0.0;
        for (k = 0; k < 3; k++)
        {
          AWplusB[i][j] += A[i][k] * MLocal[k + j * 3];
        }
        AWplusB[i][j] += B[i][j];
      }
    }
#ifdef VERBOSE_DEBUG
    for (int l = 0; l < 3; l++)
    {
      for (int k = 0; k < 3; k++)
      {
        printf("A[%i][%i] = %le\t", l, k, A[l][k]);
      }
      printf("\n");
    }
    for (int l = 0; l < 3; l++)
    {
      for (int k = 0; k < 3; k++)
      {
        printf("B[%i][%i] = %le\t", l, k, B[l][k]);
      }
      printf("\n");
    }
    for (int l = 0; l < 3; l++)
    {
      for (int k = 0; k < 3; k++)
      {
        printf("AWplusB[%i][%i] = %le\t", l, k, AWplusB[l][k]);
      }
      printf("\n");
    }
    printf("ResN = %le\n", ResN);
    printf("ResT = %le\n", ResT);
    printf("ResS = %le\n", ResS);
#endif
    d1 = AWplusB[1][1] * AWplusB[2][2] - AWplusB[1][2] * AWplusB[2][1];
    d2 = AWplusB[1][0] * AWplusB[2][2] - AWplusB[2][0] * AWplusB[1][2];
    d3 = AWplusB[1][0] * AWplusB[2][1] - AWplusB[2][0] * AWplusB[1][1];

    DET = AWplusB[0][0] * d1 - AWplusB[0][1] * d2 + AWplusB[0][2] * d3;
#ifdef VERBOSE_DEBUG
    printf("DET = %le\n", DET);
#endif
    if (fabs(DET) < 1.e-20)
    {
      printf("DET NULL\n");
      exit(EXIT_FAILURE);
    }

    DET1 = 1.0 / DET;
    AA = +(AWplusB[1][1] * AWplusB[2][2] - AWplusB[2][1] * AWplusB[1][2]) * ResN
         - (AWplusB[0][1] * AWplusB[2][2] - AWplusB[2][1] * AWplusB[0][2]) * ResT
         + (AWplusB[0][1] * AWplusB[1][2] - AWplusB[1][1] * AWplusB[0][2]) * ResS ;

    BB = -(AWplusB[1][0] * AWplusB[2][2] - AWplusB[2][0] * AWplusB[1][2]) * ResN
         + (AWplusB[0][0] * AWplusB[2][2] - AWplusB[2][0] * AWplusB[0][2]) * ResT
         - (AWplusB[0][0] * AWplusB[1][2] - AWplusB[1][0] * AWplusB[0][2]) * ResS ;


    CC = +(AWplusB[1][0] * AWplusB[2][1] - AWplusB[2][0] * AWplusB[1][1]) * ResN
         - (AWplusB[0][0] * AWplusB[2][1] - AWplusB[2][0] * AWplusB[0][1]) * ResT
         + (AWplusB[0][0] * AWplusB[1][1] - AWplusB[1][0] * AWplusB[0][1]) * ResS ;

    R[0] = R[0] - AA * DET1;
    R[1] = R[1] - BB * DET1;
    R[2] = R[2] - CC * DET1;
#ifdef VERBOSE_DEBUG
    printf("RN = %le\n", R[0]);
    printf("RT = %le\n", R[1]);
    printf("RS = %le\n", R[2]);
#endif
    // compute new residue

    for (i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * R[0] + qLocal[i]
                                            + MLocal[i + 1 * 3] * R[1] +
                                            + MLocal[i + 2 * 3] * R[2] ;


    RVN = R[0] - RhoN * velocity[0];
    RVT = R[1] - RhoT * velocity[1];
    RVS = R[2] - RhoT * velocity[2];
    RV = sqrt(RVT * RVT + RVS * RVS);
    Radius = mu * R[0];
    if (RVN >= 0.0)
    {
      ResN = RhoN * (velocity[0]);
    }
    else
    {
      ResN = R[0];
    }
    if (RV < Radius) // We are in the disk and Radius is postive
    {
      ResT = RhoT * (velocity[1]);
      ResS = RhoT * (velocity[2]);
    }
    else if (RV >= Radius && Radius > 0) // We are out the disk and Radius is postive
    {
      RV1 = 1.0 / RV;
      ResT = R[1] - Radius * RVT * RV1;
      ResS = R[2] - Radius * RVS * RV1;
    }
    else // We are out the disk and Radius is negative
    {
      ResT = R[1] ;
      ResS = R[2] ;
    }

    dparam[1] = 0.5 * (ResN * ResN + ResT * ResT + ResS * ResS);


    if (verbose > 1)
    {
      printf("-----------------------------------    AlartCurnierNewton number of iteration = %i  error = %.10e \n", inew, dparam[1]);
    }
    if (dparam[1] < Tol)
    {
      return 0;

    }


  }


  return 1;



}
