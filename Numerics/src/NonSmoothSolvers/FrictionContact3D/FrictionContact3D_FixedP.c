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

#include "NCP_Path.h"
#include "NCP_FixedP.h"
#include "NonSmoothNewton.h"
#include "FrictionContact3D_Solvers.h"
#include "LA.h"
#include <stdlib.h>
#include <stdio.h>

static FuncEvalPtr F = NULL;
static JacEvalPtr jacobianF = NULL;
static UpdateSolverPtr updateSolver = NULL;
static PostSolverPtr postSolver = NULL;
static FreeSolverPtr freeSolver = NULL;

/* size of a block */
static int Fsize;

/** writes \f$ F(z) \f$ using Glocker formulation
 */
int F_GlockerFixedP(int sizeF, double* reaction, double* FVector)
{
  /* Glocker formulation */
  int up2Date = 0;
  double* FGlocker = NULL;
  computeFGlocker(&FGlocker, up2Date);
  /* Note that FGlocker is a static var. in FrictionContact3D2NCP_Glocker and thus there is no memory allocation in
     the present file.
  */

  /* TMP COPY: review memory management for FGlocker ...*/
  DCOPY(sizeF , FGlocker , 1, FVector , 1);

  FGlocker = NULL;
  return 1;
}

/** writes \f$ \nabla_z F(z) \f$  using Glocker formulation and the Fischer-Burmeister function.
 */

void frictionContact3D_FixedP_initialize(int n0, const double*const M0, const double*const q0, const double*const mu0, int* iparam)
{

  /*
     Initialize solver (Compute F) according to the chosen formulation.
  */

  /* Glocker formulation */
  if (iparam[4] == 4)
  {
    Fsize = 5;
    NCPGlocker_initialize(n0, M0, q0, mu0);
    F = &F_GlockerFixedP;
    postSolver = &NCPGlocker_post;
    freeSolver = &NCPGlocker_free;
  }
  else
  {
    fprintf(stderr, "Numerics, FrictionContact3D_nsgs failed. Unknown formulation type.\n");
    exit(EXIT_FAILURE);
  }
}

void frictionContact3D_FixedP_initialize_SBS(int n0, const SparseBlockStructuredMatrix*const M0, const double*const q0, const double*const mu0, int* iparam)
{
  /*
     Initialize solver (Compute F) according to the chosen formulation.
  */

  /* Glocker formulation */
  if (iparam[4] == 4)
  {
    Fsize = 5;
    NCPGlocker_initialize_SBS(n0, M0, q0, mu0);
    F = &F_GlockerFixedP;
    updateSolver = &NCPGlocker_update;
    postSolver = &NCPGlocker_post;
  }
  else
  {
    fprintf(stderr, "Numerics, FrictionContact3D_nsgs failed. Unknown formulation type.\n");
    exit(EXIT_FAILURE);
  }
}

void frictionContact3D_FixedP_solve(int contact, int dimReaction, double* reaction, int* iparam, double* dparam)
{
  (*updateSolver)(contact, reaction);
  int pos = Fsize * contact; /* Current block position */
  double * reactionBlock = &reaction[pos];

  int info = Fixe(Fsize, reactionBlock, F, iparam, dparam);

  if (info > 0)
  {
    fprintf(stderr, "Numerics, FrictionContact3D_Newton failed, reached max. number of iterations without convergence. Error = %f\n", dparam[1]);
    exit(EXIT_FAILURE);
  }

  (*postSolver)(contact, reaction);
}

void frictionContact3D_FixedP_free()
{
  F = NULL;
  jacobianF = NULL;
  updateSolver = NULL;
  postSolver = NULL;
  (*freeSolver)();
}

void frictionContact3D_FixedP_computeError(int n, double* velocity, double*reaction, double * error)
{
  /*   int numberOfContacts = n/3; */
  /*   int sizeGlobal = numberOfContacts*FSize; */
  /*   //  double * FGlobal = (double*)malloc(sizeGlobal*sizeof(*FGlobal));  */
  /*   (*computeFGlobal)(reaction,velocity); */
  /*   int i; */
  /*   double Fz; */
  /*   *error = 0; */
  /*   for(i=0;i<sizeGlobal;++i) */
  /*     { */
  /*       Fz = velocity[i]*reaction[i]; */
  /*       if(Fz>0) */
  /*  *error+=Fz; */
  /*       if(reaction[i]<0) */
  /*  *error+=reaction[i]; */
  /*       if(velocity[i]<0) */
  /*  *error+=velocity[i]; */
  /*     } */

  /*   // (*computeVelocity)(FGlobal); */

  /*   free(FGlobal); */

}
