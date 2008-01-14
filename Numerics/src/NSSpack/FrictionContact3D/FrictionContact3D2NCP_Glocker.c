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

/*
   Notes FP

   - Remove qLocal => all computation in qGlocker?
   - Remove WLocal => use MGlobal[in + ...] directly?
   - Use e1, e2 ... to write everything depending on these vectors?

*/

#include "LA.h"
#include "NSSTools.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/**Static variables */
static int n = 0;
static const double* MGlobal = NULL;
static const double* qGlobal = NULL;
static const double* mu = NULL;

static double e3[2];
static double IpInv[4];
static double IpInvTranspose[4];
static double I[4];
static double mu_i = 0.0;

static const int Gsize = 5;
static double WLocal[9];
static double reactionGlocker[5];

static double jacobianFGlocker[25];
static double FGlocker[5];


void computeMGlocker()
{
  int i;
  /* jacobianFGlocker is used to save MGlocker */

  /* row 0 */
  jacobianFGlocker[0]      =   WLocal[0] + 2 * mu_i * WLocal[6];
  jacobianFGlocker[Gsize]   =  -WLocal[6] - sqrt(3.) / 3 * WLocal[3] ;
  jacobianFGlocker[2 * Gsize] = sqrt(3.) / 3 * WLocal[3] - WLocal[6];
  jacobianFGlocker[3 * Gsize] = 0.0;
  jacobianFGlocker[4 * Gsize] = 0.0;

  /* row 1 */
  jacobianFGlocker[1]        = - WLocal[2] - sqrt(3.) / 3 * WLocal[1] - 2 * mu_i * WLocal[8]  - 2 * sqrt(3.) / 3 * mu_i * WLocal[7];
  jacobianFGlocker[Gsize + 1]   =   WLocal[8] + 1. / 3 * WLocal[4] + sqrt(3.) / 3 * (WLocal[5] + WLocal[7]) ;
  jacobianFGlocker[2 * Gsize + 1] =   WLocal[8] - 1. / 3 * WLocal[4] - sqrt(3.) / 3 * (WLocal[5] - WLocal[7]) ;
  jacobianFGlocker[3 * Gsize + 1] = 1.0;
  jacobianFGlocker[4 * Gsize + 1] = 0.0;

  /* row 2 */
  jacobianFGlocker[2]        =  - WLocal[2] + sqrt(3.) / 3 * WLocal[1] - 2 * mu_i * WLocal[8] + 2 * sqrt(3.) / 3 * mu_i * WLocal[7];
  jacobianFGlocker[Gsize + 2]   =    WLocal[8] - 1. / 3 * WLocal[4] + sqrt(3.) / 3 * (WLocal[5] - WLocal[7]) ;
  jacobianFGlocker[2 * Gsize + 2] =    WLocal[8] + 1. / 3 * WLocal[4] - sqrt(3.) / 3 * (WLocal[5] + WLocal[7]) ;;
  jacobianFGlocker[3 * Gsize + 2] = 1.0;
  jacobianFGlocker[4 * Gsize + 2] = 0.0;

  /* row 3 */
  jacobianFGlocker[3]        = 3 * mu_i;
  jacobianFGlocker[Gsize + 3]   = -1.0;
  jacobianFGlocker[2 * Gsize + 3] = -1.0;
  jacobianFGlocker[3 * Gsize + 3] = 0.0;
  jacobianFGlocker[4 * Gsize + 3] = 0.0;

  /* row 4 */
  for (i = 0; i < Gsize; ++i)
    jacobianFGlocker[i * Gsize + 4] = 0.0;
}

void computeQGlocker(int in, double* reaction)
{
  /*   - computes qLocal = qGlobal[in] + sum over a row of blocks in MGlobal of the products WLocal.reaction,
   excluding the block corresponding to the current contact.
   - computes qGlocker using qLocal values
  ****/
  int it = in + 1, is = in + 2;
  /* reaction current block set to zero, to exclude current contact block */
  reaction[in] = 0.0;
  reaction[it] = 0.0;
  reaction[is] = 0.0;
  /* qGlocker computation*/
  int incx = n, incy = 1;
  FGlocker[0] = qGlobal[in] + DDOT(n , &MGlobal[in] , incx , reaction , incy);
  double tmp[2];
  tmp[0] = qGlobal[it] + DDOT(n , &MGlobal[it] , incx , reaction , incy);
  tmp[1] = qGlobal[is] + DDOT(n , &MGlobal[is] , incx , reaction , incy);
  /* FGlocker[1 ..2] = -IpInv * tmp */
  /*   DGEMV(LA_NOTRANS,2,2,-1.0,IpInv,2,tmp,incx,0.0,&FGlocker[1],incy); */
  FGlocker[1] = -sqrt(3.) / 3 * tmp[0] - tmp[1];
  FGlocker[2] =  sqrt(3.) / 3 * tmp[0] - tmp[1];
  FGlocker[3] = 0.0;
  FGlocker[4] = 0.0;

}

void computeGGlocker()
{
  /* computation of gGlocker (added to FGlocker) */

  // gGlocker[0] = 0.0;

  FGlocker[1] += 4. / 3 * reactionGlocker[4] * (2 * reactionGlocker[1] + reactionGlocker[2] - 3 * mu_i * reactionGlocker[0]);
  FGlocker[2] += 4. / 3 * reactionGlocker[4] * (reactionGlocker[1] + 2 * reactionGlocker[2] - 3 * mu_i * reactionGlocker[0]);
  // gGlocker[3] = 0.0;
  /* norm_I of row1*/
  double tmp = 4. / 3 * (reactionGlocker[1] * reactionGlocker[1] + reactionGlocker[1] * reactionGlocker[2] + reactionGlocker[2] * reactionGlocker[2]);
  FGlocker[4] += 4 * mu_i * reactionGlocker[0] * reactionGlocker[1] - 3 * mu_i * mu_i * reactionGlocker[0] * reactionGlocker[0] + 4 * mu_i * reactionGlocker[0] * reactionGlocker[2] - tmp;

}

void computeJacobianGGlocker()
{
  /* At this point, jacobianFGlocker = M */
  /* We compute jacobianFGlocker += jacobian of g */
  /* row 0 = 0 */
  /* row 1 */
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

void initializeNCPGlocker_SB(int n0, const SparseBlockStructuredMatrix*const M0, const double*const q0, const double*const mu0)
{
}

void initializeNCPGlocker(int n0, const double*const M0, const double*const q0, const double*const mu0)
{
  /* 1 - Get and save global objects (M,q, ...)
     2 - Fill vectors/matrices of parameters: Ip, Ipinv ...
  */

  /* 1 */
  MGlobal = M0;
  qGlobal = q0;
  mu = mu0;
  n = n0;
  /* ei = [cos((4i-3)Pi/6), sin-((4i-3)Pi/6)]
     Ip = [e1 e2]
  */

  e3[0] = 0.;
  e3[1] = -1.;

  /* inverse of Ip */
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

void updateNCPGlocker(int contact, double* reaction)
{
  /* Update glocker/local M,q,reaction ... (functions of global M,q ...) depending on the considered contact */

  int in = 3 * contact, it = in + 1, is = it + 1;
  int inc = n * in;
  /***** Friction coefficient for current block *****/
  mu_i = mu[contact];

  /****** initialize rG ******/
  reactionGlocker[0] = reaction[in]; /* Pn */
  reactionGlocker[1] = mu_i * reaction[in] - sqrt(3) * reaction[it] / 2. - reaction[is] / 2.; /* SigmaP_1 */
  reactionGlocker[2] = mu_i * reaction[in] + sqrt(3) * reaction[it] / 2. - reaction[is] / 2.; /* SigmaP_2 */
  reactionGlocker[3] = 0.; /* k3 */
  reactionGlocker[4] = 0.; /* kD */

  /****** computation of M (saved into jacobianFGlocker) ******/
  /* The part of MGlobal which corresponds to the current block is copied into WLocal */
  WLocal[0] = MGlobal[(in) * n + in];
  WLocal[1] = MGlobal[(in) * n + it];
  WLocal[2] = MGlobal[(in) * n + is];
  inc += n;
  WLocal[3] = MGlobal[(it) * n + in];
  WLocal[4] = MGlobal[(it) * n + it];
  WLocal[5] = MGlobal[(it) * n + is];
  inc += n;
  WLocal[6] = MGlobal[(is) * n + in];
  WLocal[7] = MGlobal[(is) * n + it];
  WLocal[8] = MGlobal[(is) * n + is];

  computeMGlocker();

  /****** Computation of qGlocker (saved into FGlocker) ******/
  computeQGlocker(in, reaction);

  /****** Computation of gGlocker (added to FGlocker) ******/
  computeGGlocker();

}

void postSolverNCPGlocker(int contact, double * reaction0)
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
  DGEMV(LA_NOTRANS, Gsize, Gsize, 1.0, jacobianFGlocker, Gsize, reactionGlocker, 1, 1.0, FGlocker, 1);
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

