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

/*
   Notes FP

   - Remove WLocal => use MGlobal[in + ...] directly?
   - Use e1, e2 ... to write everything depending on these vectors?

*/

#include "LA.h"
#include "FrictionContact3D_Solvers.h"
#include "Numerics_Options.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*Static variables */

/* The global problem of size n= 3*nc, nc being the number of contacts, is locally saved in MGlobal and qGlobal */
/* mu corresponds to the vector of friction coefficients */


/* Global problem */
static int n = 0;
static const NumericsMatrix* MGlobal = NULL;
static const double* qGlobal = NULL;
static const double* mu = NULL;

/* Local FC3D problem */
static double* MLocal;
static int isMAllocatedIn = 0; /* True if a malloc is done for MLocal, else false */
static double qLocal[3];

/* Local "Glocker" variables */
static const int Gsize = 5;
static double reactionGlocker[5];
static double MGlocker[25];
static double qGlocker[5];
static double gGlocker[5];

/* Output */
static double jacobianFGlocker[25];
static double FGlocker[5];

static double mu_i = 0.0;

static double e1[3], e2[3], e3[2];
static double IpInv[4];
static double IpInvTranspose[4];
static double I[4];
# define PI 3.14159265358979323846  /* pi */

void computeE(unsigned int i, double* e)
{
  double alpha = (4 * i - 3) * PI / 6.0;
  e[0] = cos(alpha);
  e[1] = sin(alpha);
}

/* Compute and save MGlocker */
void computeMGlocker()
{
  /* Local function used in update */

  /* Requires: MLocal and mu_i*/

  int i;
  /* MGlocker is used to save MGlocker */

  /* row 0 */
  MGlocker[0]      =   MLocal[0] + 2 * mu_i * MLocal[6];
  MGlocker[Gsize]   =  -MLocal[6] - sqrt(3.) / 3 * MLocal[3] ;
  MGlocker[2 * Gsize] = sqrt(3.) / 3 * MLocal[3] - MLocal[6];
  MGlocker[3 * Gsize] = 0.0;
  MGlocker[4 * Gsize] = 0.0;

  /* row 1 */
  MGlocker[1]        = - MLocal[2] - sqrt(3.) / 3 * MLocal[1] - 2 * mu_i * MLocal[8]  - 2 * sqrt(3.) / 3 * mu_i * MLocal[7];
  MGlocker[Gsize + 1]   =   MLocal[8] + 1. / 3 * MLocal[4] + sqrt(3.) / 3 * (MLocal[5] + MLocal[7]) ;
  MGlocker[2 * Gsize + 1] =   MLocal[8] - 1. / 3 * MLocal[4] - sqrt(3.) / 3 * (MLocal[5] - MLocal[7]) ;
  MGlocker[3 * Gsize + 1] = 1.0;
  MGlocker[4 * Gsize + 1] = 0.0;

  /* row 2 */
  MGlocker[2]        =  - MLocal[2] + sqrt(3.) / 3 * MLocal[1] - 2 * mu_i * MLocal[8] + 2 * sqrt(3.) / 3 * mu_i * MLocal[7];
  MGlocker[Gsize + 2]   =    MLocal[8] - 1. / 3 * MLocal[4] + sqrt(3.) / 3 * (MLocal[5] - MLocal[7]) ;
  MGlocker[2 * Gsize + 2] =    MLocal[8] + 1. / 3 * MLocal[4] - sqrt(3.) / 3 * (MLocal[5] + MLocal[7]) ;;
  MGlocker[3 * Gsize + 2] = 1.0;
  MGlocker[4 * Gsize + 2] = 0.0;

  /* row 3 */
  MGlocker[3]        = 3 * mu_i;
  MGlocker[Gsize + 3]   = -1.0;
  MGlocker[2 * Gsize + 3] = -1.0;
  MGlocker[3 * Gsize + 3] = 0.0;
  MGlocker[4 * Gsize + 3] = 0.0;

  /* row 4 */
  for (i = 0; i < Gsize; ++i)
    MGlocker[i * Gsize + 4] = 0.0;
}

void computeGGlocker()
{
  /* Local function used in update */

  /* Requires: reactionGlocker and mu_i*/

  /* computation of gGlocker (added to FGlocker) */

  // gGlocker[0] = 0.0;

  FGlocker[1] += 4. / 3 * reactionGlocker[4] * (2 * reactionGlocker[1] + reactionGlocker[2] - 3 * mu_i * reactionGlocker[0]);
  FGlocker[2] += 4. / 3 * reactionGlocker[4] * (reactionGlocker[1] + 2 * reactionGlocker[2] - 3 * mu_i * reactionGlocker[0]);
  // gGlocker[3] = 0.0;
  /* norm_I of row1*/
  double tmp = 4. / 3 * (reactionGlocker[1] * reactionGlocker[1] + reactionGlocker[1] * reactionGlocker[2] + reactionGlocker[2] * reactionGlocker[2]);
  FGlocker[4] += 4 * mu_i * reactionGlocker[0] * reactionGlocker[1] - 3 * mu_i * mu_i * reactionGlocker[0] * reactionGlocker[0] + 4 * mu_i * reactionGlocker[0] * reactionGlocker[2] - tmp;

}

void NCPGlocker_fillMLocal(int contact)
{
  int storageType = MGlobal->storageType;
  if (storageType == 0)
  {
    int in = 3 * contact, it = in + 1, is = it + 1;
    int inc = n * in;
    // === Fill MLocal(3,3) according to the current contact number ===
    /* The part of MGlobal which corresponds to the current block is copied into MLocal */
    double * MM = MGlobal->matrix0;
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
    MLocal = MGlobal->matrix1->block[diagPos];
  }
  else
    numericsError("FrictionContact3D2NCP_Glocker::NCPGlocker_fillMLocal() -", "unknown storage type for matrix M");

}

void NCPGlocker_initialize(int n0, const NumericsMatrix*const M0, const double*const q0, const double*const mu0)
{
  /*
    INPUT: the global problem operators: n0 (size), M0, q0 and mu0, vector of friction coefficients.
    In initialize, these operators are "connected" to their corresponding static variables, that will be used to build local problem
    for each considered contact.
    Local problem is built during call to update (which depends on the storage type for M).

    Fill vectors/matrices of parameters: Ip, Ipinv ...
  */

  /* Connect static var to input (global problem) */
  n = n0;
  MGlobal = M0;
  qGlobal = q0;
  mu = mu0;

  if (MGlobal->storageType == 0)
  {
    MLocal = (double*)malloc(3 * 3 * sizeof(*MLocal));
    isMAllocatedIn = 1;
  }
  else
    isMAllocatedIn = 0;

  /* ei = [cos((4i-3)Pi/6), sin-((4i-3)Pi/6)]
     Ip = [e1 e2]
  */
  //computeE(1,e1);  computeE(2,e2);
  computeE(3, e3);

  /* compute I and inverse of Ip */
  //  computeI(e1,e2,e3, IpInv, I)
  IpInv[0] =  sqrt(3.) / 3 ;
  IpInv[2] =  1.;
  IpInv[1] =  -sqrt(3.) / 3 ;
  IpInv[3] =  1.;

  /* transpose of IpInv */
  IpInvTranspose[0] =  sqrt(3.) / 3;
  IpInvTranspose[2] = -sqrt(3.) / 3;
  IpInvTranspose[1] = 1.;
  IpInvTranspose[3] =  1.;

  /* I = IpInv * IpInvTranspose */
  I[0] = 4.0 / 3;
  I[1] = 2.0 / 3;
  I[2] = 2.0 / 3;
  I[3] = 4.0 / 3;

}

void NCPGlocker_update(int contact, double * reaction)
{
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */
  /* Update glocker/local M,q,reaction ... (functions of global M,q ...) depending on the considered contact */

  int in = 3 * contact, it = in + 1, is = it + 1;

  // === Friction coefficient for current block ===
  mu_i = mu[contact]; /* required in computeMGlocker */

  // === Fill MLocal(3,3) according to the current contact number ===
  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  NCPGlocker_fillMLocal(contact);

  // === computation of MGlocker = function(MLocal, mu_i, I, IpInvTranspose, IpInv, e3) ===/
  computeMGlocker();

  // === computation of qGlocker = function(qLocal, IpInv)
  // saved in FGlocker which is also initialized here ===
  //  - step 1: computes qLocal = qGlobal[in] + sum over a row of blocks in MGlobal of the products MLocal.reaction,
  //            excluding the block corresponding to the current contact.
  //  - step 2: computes qGlocker using qLocal values

  /* reaction current block set to zero, to exclude current contact block */

  double rin = reaction[in] ;
  double rit = reaction[it] ;
  double ris = reaction[is] ;
  reaction[in] = 0.0;
  reaction[it] = 0.0;
  reaction[is] = 0.0;
  /* qLocal computation*/
  double qLocal[3];
  qLocal[0] = qGlobal[in];
  qLocal[1] = qGlobal[it];
  qLocal[2] = qGlobal[is];
  int storageType = MGlobal->storageType;
  int incx = n, incy = 1;
  if (storageType == 0)
  {
    double * MM = MGlobal->matrix0;
    qLocal[0] += DDOT(n , &MM[in] , incx , reaction , incy);
    qLocal[1] += DDOT(n , &MM[it] , incx , reaction , incy);
    qLocal[2] += DDOT(n , &MM[is] , incx , reaction , incy);
  }
  else if (storageType == 1)
  {
    /* qLocal += rowMB * reaction
    with rowMB the row of blocks of MGlobal which corresponds to the current contact
    */
    subRowProdSBM(n, 3, contact, MGlobal->matrix1, reaction, qLocal, 0);
  }

  reaction[in] = rin;
  reaction[it] = rit;
  reaction[is] = ris;
  /* qGlocker (saved in FGlocker) */
  FGlocker[0] = qLocal[0];
  FGlocker[1] = -sqrt(3.) / 3 * qLocal[1] - qLocal[2];
  FGlocker[2] =  sqrt(3.) / 3 * qLocal[1] - qLocal[2];
  FGlocker[3] = 0.0;
  FGlocker[4] = 0.0;

  // === initialize reactionGlocker with reaction[currentContact] ===
  // reactionGlocker = function(reaction, mu_i)
  reactionGlocker[0] = reaction[in]; /* Pn */
  reactionGlocker[1] = mu_i * reaction[in] - sqrt(3) * reaction[it] / 2. - reaction[is] / 2.; /* SigmaP_1 */
  reactionGlocker[2] = mu_i * reaction[in] + sqrt(3) * reaction[it] / 2. - reaction[is] / 2.; /* SigmaP_2 */
  reactionGlocker[3] = 0.; /* k3 */
  reactionGlocker[4] = 0.; /* kD */

  // === Computation of gGlocker = function(reactionGlocker, I, mu_i) (added to FGlocker) ===/
  computeGGlocker();

  // End of this function:
  // - reactionGlocker is up to date
  // - MGlocker is up to date
  // - FGlocker = qGlocker + gGlocker
}

void computeJacobianGGlocker()
{
  /* At this point, jacobianFGlocker = M */
  /* We compute jacobianFGlocker += jacobian of g */
  /* row 0 = 0 */
  /* row 1 */

  // === jacobianFGlocker = MGlocker + jacobian_r g(r) ===
  DCOPY(Gsize * Gsize, MGlocker, 1, jacobianFGlocker, 1);

  double muR0 = 4.*mu_i * reactionGlocker[0];
  jacobianFGlocker[1]  -= 4.*mu_i * reactionGlocker[4];
  jacobianFGlocker[6]  += 8. / 3 * reactionGlocker[4];
  jacobianFGlocker[11] += 4. / 3 * reactionGlocker[4];
  jacobianFGlocker[21] += 8. / 3 * reactionGlocker[1] + 4. / 3 * reactionGlocker[2] - muR0;

  /* row 2 */
  jacobianFGlocker[2]  -= 4.*mu_i * reactionGlocker[4];
  jacobianFGlocker[7]  += 4. / 3 * reactionGlocker[4];
  jacobianFGlocker[12] += 8. / 3 * reactionGlocker[4];
  jacobianFGlocker[22] += 4. / 3 * reactionGlocker[1] + 8. / 3 * reactionGlocker[2] - muR0;

  /* row 3  = 0 */
  /* row 4 */
  jacobianFGlocker[4]  += 2.*mu_i * (2.*reactionGlocker[1] + 2.*reactionGlocker[2] - 3.*mu_i * reactionGlocker[0]);
  jacobianFGlocker[9]  += muR0 - 4. / 3 * reactionGlocker[2] - 8. / 3 * reactionGlocker[1];
  jacobianFGlocker[14] += muR0 - 8. / 3 * reactionGlocker[2] - 4. / 3 * reactionGlocker[1];

}

void  NCPGlocker_post(int contact, double * reaction0)
{
  /* Retrieve original vector reaction from reactionGlocker formulation
     Rn = reaction[0] =  reactionGlocker[0]

     RT = (reaction[1]) = IpInvTranspose*( reactionGlocker[0]*mu_i - (reactionGlocker[1]) )
     (reaction[2])                                             (reactionGlocker[2])

  */

  int in = 3 * contact, it = in + 1, is = in + 2;
  reaction0[in] = reactionGlocker[0];
  double tmp1 = mu_i * reactionGlocker[0] - reactionGlocker[1];
  double tmp2 = mu_i * reactionGlocker[0] - reactionGlocker[2];
  reaction0[it] = IpInvTranspose[0] * tmp1 + IpInvTranspose[2] * tmp2;
  reaction0[is] = IpInvTranspose[1] * tmp1 + IpInvTranspose[3] * tmp2;
}



void computeFGlocker(double** FOut, int up2Date)
{
  /* updateNCPGlocker must have been called before !!! */
  /* At this point, FGlocker = gGlocker + qGlocker */
  /* and jacobianFGlocker contains MGlocker */
  /* F = M.reaction + g(reaction) + q */
  DGEMV(LA_NOTRANS, Gsize, Gsize, 1.0, MGlocker, Gsize, reactionGlocker, 1, 1.0, FGlocker, 1);
  *FOut = FGlocker; /* pointer link */
}

void computeJacobianFGlocker(double** jacobianFOut, int up2Date)
{
  /* Computation of M (if required) and of the jacobian of g */
  if (up2Date == 0)
  {
    computeMGlocker();
    /* MGlocker is saved in jacobianF */
  }
  computeJacobianGGlocker(); /* add jacobianG to previously computed M in jacobianF */
  *jacobianFOut = jacobianFGlocker; /* pointer link */
}

void NCPGlocker_free()
{
  MGlobal = NULL;
  qGlobal = NULL;
  mu = NULL;
  if (isMAllocatedIn == 1)
    free(MLocal);
  MLocal = NULL;
}

