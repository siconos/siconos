/* Siconos-Numerics, Copyright INRIA 2005-2015
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
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#include <assert.h>

#include "GlobalFrictionContact3D_Gams.h"
#include "GlobalFrictionContactProblem.h"
#include "SolverOptions.h"


/* Alart & Curnier solver for sparse global problem */
void globalFrictionContact3D_FixedPointCadoux(
  GlobalFrictionContactProblem* problem,
  double *reaction,
  double *velocity,
  double *globalVelocity,
  int *info,
  SolverOptions* options)
{

  assert(problem);
  assert(reaction);
  assert(velocity);
  assert(info);
  assert(options);

  assert(problem->dimension == 3);

  assert(options->iparam);
  assert(options->dparam);

  assert(problem->q);
  assert(problem->mu);
  assert(problem->M);
  assert(problem->H);

  assert(!problem->M->matrix0);

  /* M is square */
  assert(problem->M->size0 == problem->M->size1);

  assert(problem->M->size0 == problem->H->size0);

  unsigned int iter = 0;
  unsigned int itermax = options->iparam[0];
  unsigned int erritermax = options->iparam[7];

  if (erritermax == 0)
  {
    /* output a warning here */
    erritermax = 1;
  }

  assert(itermax > 0);
  assert(options->iparam[3] > 0);

  double tolerance = options->dparam[0];
  assert(tolerance > 0);

  /* sparse triplet storage */
//  NM_setup(problem->M);
//  NM_setup(problem->H);

  //Mlu = (double*)
  /* save multiple data M, H */



}
