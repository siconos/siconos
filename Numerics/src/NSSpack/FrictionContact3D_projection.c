
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
/*!\file pfc_3D_nlgs.c
 *
 * This subroutine allows the primal resolution of contact problems with friction.\n
 *
 *   Try \f$(reaction,velocity)\f$ such that:\n
 *   \f$
 *    \left\lbrace
 *     \begin{array}{l}
 *      M reaction + q = velocity\\
 *      0 \le reaction_n \perp velocity_n \ge 0\\
 *      -velocity_t \in \partial\psi_{[-\mu reaction_n, \mu reaction_n]}(reaction_t)\\
 *     \end{array}
 *    \right.
 *   \f$
 *
 * here M is an n by n  matrix, q an n-dimensional vector, reaction an n-dimensional  vector and velocity an n-dimensional vector.
 *
 * \fn  pfc_3D_projection( int n , double *C , double *b , double *zz , double *ww , double coef ,\n
 *                         void (*Compute_G), void (*Compute_JacG), int *iparam_local , double *dparam_local ) \n
 *
 * Generic pfc_3D parameters:\n
 *
 * \param nn      Unchanged parameter which represents the dimension of the system.
 * \param vec     Unchanged parameter which contains the components of the matrix with a fortran storage.
 * \param q       Unchanged parameter which contains the components of the right hand side vector.
 * \param reaction       Modified parameter which contains the initial solution and returns the solution of the problem.
 * \param w       Modified parameter which returns the solution of the problem.
 * \param coef    Unchanged parameter which represents the friction coefficient
 *
 *
 * \author Houari Khenous and Franck Perignon last modification 05/12/2007 .
 *
 */

#include "LA.h"
#include "NSSTools.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/**Static variables */
static int n = 0;
static const double* MGlobal = NULL;
static const SparseBlockStructuredMatrix* MBGlobal = NULL;
static const double* qGlobal = NULL;
static const double* mu = NULL;

static const int nBlock = 3;
static double* MBlock;
static double reactionBlock[3];
static double qLocal[3];
static double mu_i = 0.0;

/* update pointer to function, used to switch to the function adapted to the storage for M. */
static void (*update)(int, double*);

void proj_updateWithSparse(int contact, double * reaction)
{
  int in = 3 * contact, it = in + 1, is = it + 1;
  int j;
  int numberOfContact = n / 3;
  /* The part of M which corresponds to the current block is copied into MBlock */
  int diagPos = numberOfContact * contact + contact;
  MBlock = MBGlobal->block[diagPos];

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
  qLocal[0] = qGlobal[in];
  qLocal[1] = qGlobal[it];
  qLocal[2] = qGlobal[is];
  /* Loop through the columns(blocks) of M to compute qLocal */
  int blockNum = contact * numberOfContact;
  for (j = 0; j < numberOfContact ; ++j)
  {
    if (j != contact)
    {
      DGEMV(LA_NOTRANS, 3, 3, 1.0, MBGlobal->block[blockNum], 3, &reaction[3 * j], incx , 1.0, qLocal, incy);
    }
    blockNum = blockNum + 1;
  }
}

void proj_updateNoSparse(int contact, double * reaction)
{
  int in = 3 * contact, it = in + 1, is = it + 1;
  int j, inc = n * in;

  /* The part of MGlobal which corresponds to the current block is copied into MBlock */
  MBlock[0] = MGlobal[inc + in];
  MBlock[1] = MGlobal[inc + it];
  MBlock[2] = MGlobal[inc + is];
  inc += n;
  MBlock[3] = MGlobal[inc + in];
  MBlock[4] = MGlobal[inc + it];
  MBlock[5] = MGlobal[inc + is];
  inc += n;
  MBlock[6] = MGlobal[inc + in];
  MBlock[7] = MGlobal[inc + it];
  MBlock[8] = MGlobal[inc + is];

  /* reactionBlock */
  for (j = 0 ; j < nBlock ; ++j)
    reactionBlock[j] = reaction[in + j];

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MBlock.reactionBlock,
   excluding the block corresponding to the current contact. ****/

  /* reaction current block set to zero, to exclude current contact block */
  reaction[in] = 0.0;
  reaction[it] = 0.0;
  reaction[is] = 0.0;
  /* qLocal computation*/
  int incx = n, incy = 1;
  qLocal[0] = qGlobal[in] + DDOT(n , &MGlobal[in] , incx , reaction , incy);
  qLocal[1] = qGlobal[it] + DDOT(n , &MGlobal[it] , incx , reaction , incy);
  qLocal[2] = qGlobal[is] + DDOT(n , &MGlobal[is] , incx , reaction , incy);
}

void frictionContact3D_projection_initialize(int n0, const double*const M0, const double*const q0, const double*const mu0)
{
  n = n0;
  MGlobal = M0;
  qGlobal = q0;
  mu = mu0;
  update = &proj_updateNoSparse;
  MBlock = (double*)malloc(nBlock * nBlock * sizeof(*MBlock));
}

void frictionContact3D_projection_initialize_SB(int n0, const SparseBlockStructuredMatrix*const M0, const double*const q0, const double*const mu0)
{
  n = n0;
  MBGlobal = M0;
  qGlobal = q0;
  mu = mu0;
  update = &proj_updateWithSparse;
}

void frictionContact3D_projection_update(int contact, double* reaction)
{
  /* Call the update function which depends on the storage for M/MB */
  (*update)(contact, reaction);
  /* Friction coefficient for current block*/
  mu_i = mu[contact];
}

void frictionContact3D_projection_solve(int contact, int dimReaction, double* reaction, int* iparam, double* dparam)
{
  /* Current block position */
  int pos = contact * nBlock;
  double mrn, num, mu2 = mu_i * mu_i;

  frictionContact3D_projection_update(contact, reaction);

  if (qLocal[0] > 0.)
  {
    reaction[pos] = 0.;
    reaction[pos + 1] = 0.;
    reaction[pos + 2] = 0.;
  }
  else
  {
    if (MBlock[0] < 1e-16 || MBlock[nBlock + 1] < 1e-16 || MBlock[2 * nBlock + 2] < 1e-16)
    {
      fprintf(stderr, "FrictionContact3D_projection error: null term on MBlock diagonal.\n");
      exit(EXIT_FAILURE);
    }

    reaction[pos] = -qLocal[0] / MBlock[0];
    reaction[pos + 1] = -qLocal[1] / MBlock[nBlock + 1];
    reaction[pos + 2] = -qLocal[2] / MBlock[2 * nBlock + 2];

    mrn = reaction[pos + 1] * reaction[pos + 1] + reaction[pos + 2] * reaction[pos + 2];

    if (mrn > mu2 * reaction[pos]*reaction[pos])
    {
      num = mu_i * reaction[0] / sqrt(mrn);
      reaction[pos + 1] = reaction[pos + 1] * num;
      reaction[pos + 2] = reaction[pos + 2] * num;
    }
  }
}

void frictionContact3D_projection_free()
{
  MGlobal = NULL;
  qGlobal = NULL;
  mu = NULL;
  free(MBlock);
}
