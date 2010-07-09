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

#include "NCP_Path.h"
#include "NCP_FixedP.h"
#include "NCP_FischerBurmeister.h"
#include "NonSmoothNewton.h"
#include "FrictionContact3D_Solvers.h"
#include "FrictionContact3D2NCP_Glocker.h"
#include "LA.h"
#include <stdlib.h>
#include <stdio.h>
/* Pointer to function used to update the solver, to formalize the local problem for example. */
typedef void (*UpdateSolverPtr)(int, double*);

static UpdateSolverPtr updateSolver = NULL;
static PostSolverPtr postSolver = NULL;
static FreeSolverPtr freeSolver = NULL;

/* size of a block */
static int Fsize;

/** writes \f$ F(z) \f$ using Glocker formulation
 */
void F_GlockerFixedP(int sizeF, double* reaction, double* FVector, int up2Date)
{
  /* Glocker formulation */
  double* FGlocker = NULL;
  computeFGlocker(&FGlocker, up2Date);
  /* Note that FGlocker is a static var. in FrictionContact3D2NCP_Glocker and thus there is no memory allocation in
     the present file.
  */

  /* TMP COPY: review memory management for FGlocker ...*/
  DCOPY(sizeF , FGlocker , 1, FVector , 1);
  FGlocker = NULL;
}

/** writes \f$ \nabla_z F(z) \f$  using Glocker formulation and the Fischer-Burmeister function.
 */


void frictionContact3D_FixedP_initialize(FrictionContactProblem* problem, FrictionContactProblem* localproblem, SolverOptions * localsolver_options)
{

  /*
     Initialize solver (Compute F) according to the chosen formulation.
  */

  /* Glocker formulation */
  if (localsolver_options->solverId == SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint)
  {
    Fsize = 5;
    NCPGlocker_initialize(problem, localproblem);
    /*     updateSolver = &NCPGlocker_update; */
    postSolver = &NCPGlocker_post;
    freeSolver = &NCPGlocker_free;
  }
  else
  {
    fprintf(stderr, "Numerics, FrictionContact3D_nsgs failed. Unknown formulation type.\n");
    exit(EXIT_FAILURE);
  }
}

void frictionContact3D_FixedP_solve(FrictionContactProblem * localproblem , double* reaction, SolverOptions * options)
{
  int * iparam = options->iparam;
  double * dparam = options->dparam;

  double * reactionBlock = reaction;

  int info = Fixe(Fsize, reactionBlock, iparam, dparam);

  if (info > 0)
  {
    fprintf(stderr, "Numerics, FrictionContact3D_FixedP failed, reached max. number of iterations without convergence. Error = %f\n", dparam[1]);
    exit(EXIT_FAILURE);
  }

  /*   (*postSolver)(contact,reaction); */
}

void frictionContact3D_FixedP_free()
{
  updateSolver = NULL;
  postSolver = NULL;
  (*freeSolver)();
}

double frictionContact3D_FixedP_computeError(int contact, int dimReaction, double* reaction, double * error)
{
  return 0.0;
}
