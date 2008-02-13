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

#include "NonSmoothNewton.h"
#include "NCP.h"
#include "FrictionContact3D_Solvers.h"
#include <stdlib.h>
#include <stdio.h>

static NewtonFunctionPtr F = NULL;
static NewtonFunctionPtr jacobianF = NULL;
static UpdateSolverPtr updateSolver = NULL;
static PostSolverPtr postSolver = NULL;
static FreeSolverPtr freeSolver = NULL;
/* size of a block */
static int Fsize;

/** writes \f$ F(z) \f$ using Glocker formulation and the Fischer-Burmeister function.
 */
void F_GlockerFischerBurmeister(int sizeF, double* reaction, double* FVector, int up2Date)
{
  /* Glocker formulation */
  double* FGlocker = NULL;
  computeFGlocker(&FGlocker, up2Date);
  /* Note that FGlocker is a static var. in FrictionContact3D2NCP_Glocker and thus there is no memory allocation in
   the present file.
  */

  /* Call Fisher-Burmeister function => fill FVector */
  phi_FB(sizeF, reaction, FGlocker, FVector);
  FGlocker = NULL;
}

/** writes \f$ \nabla_z F(z) \f$  using Glocker formulation and the Fischer-Burmeister function.
 */
void jacobianF_GlockerFischerBurmeister(int sizeF, double* reaction, double* jacobianFMatrix, int up2Date)
{
  /* Glocker formulation */
  double* FGlocker = NULL, *jacobianFGlocker = NULL;
  computeFGlocker(&FGlocker, up2Date);
  computeJacobianFGlocker(&jacobianFGlocker, up2Date);
  /* Note that FGlocker and jacobianFGlocker are static var. in FrictionContact3D2NCP_Glocker and thus there is no memory allocation in
   the present file.
  */

  /* Call Fisher-Burmeister function => fill jacobianFMatrix */
  jacobianPhi_FB(sizeF, reaction, FGlocker, jacobianFGlocker, jacobianFMatrix);
  FGlocker = NULL;
  jacobianFGlocker = NULL;
}


void frictionContact3D_Newton_initialize(int n0, const double*const M0, const double*const q0, const double*const mu0, int* iparam)
{

  /*
     Initialize solver (Connect F and its jacobian, set local size ...) according to the chosen formulation.
  */

  /* Alart-Curnier formulation */
  if (iparam[4] == 1)
  {
    Fsize = 3;
    frictionContact3D_AC_initialize(n0, M0, q0, mu0);
    F = &F_AC;
    jacobianF = &jacobianF_AC;
    updateSolver = &frictionContact3D_AC_update;
    postSolver = &frictionContact3D_AC_post;
    freeSolver = &frictionContact3D_AC_free;

  }
  /* Glocker formulation - Fischer-Burmeister function used in Newton */
  else if (iparam[4] == 2)
  {
    Fsize = 5;
    NCPGlocker_initialize(n0, M0, q0, mu0);
    F = &F_GlockerFischerBurmeister;
    jacobianF = &jacobianF_GlockerFischerBurmeister;
    updateSolver = &NCPGlocker_update;
    postSolver = &NCPGlocker_post;
    freeSolver = &NCPGlocker_free;
  }
  else
  {
    fprintf(stderr, "Numerics, FrictionContact3D_nsgs failed. Unknown formulation type.\n");
    exit(EXIT_FAILURE);
  }
}

void frictionContact3D_Newton_initialize_SBS(int n0, const SparseBlockStructuredMatrix*const M0, const double*const q0, const double*const mu0, int* iparam)
{
  /*
     Initialize solver (Connect F and its jacobian, set local size ...) according to the chosen formulation.
  */

  /* Alart-Curnier formulation */
  if (iparam[4] == 1)
  {
    Fsize = 3;
    F = &F_AC;
    jacobianF = &jacobianF_AC;
    frictionContact3D_AC_initialize_SBS(n0, M0, q0, mu0);
    updateSolver = &frictionContact3D_AC_update;
    postSolver = &frictionContact3D_AC_post;

  }
  /* Glocker formulation - Fischer-Burmeister function used in Newton */
  else if (iparam[4] == 2)
  {
    Fsize = 5;
    NCPGlocker_initialize_SBS(n0, M0, q0, mu0);
    F = &F_GlockerFischerBurmeister;
    jacobianF = &jacobianF_GlockerFischerBurmeister;
    updateSolver = &NCPGlocker_update;
    postSolver = &NCPGlocker_post;
  }
  else
  {
    fprintf(stderr, "Numerics, FrictionContact3D_nsgs failed. Unknown formulation type.\n");
    exit(EXIT_FAILURE);
  }
}

void frictionContact3D_Newton_solve(int contact, int dimReaction, double* reaction, int* iparam, double* dparam)
{
  (*updateSolver)(contact, reaction);
  int pos = Fsize * contact; /* Current block position */
  double * reactionBlock = &reaction[pos];
  int info = nonSmoothNewton(Fsize, reactionBlock, &F, &jacobianF, iparam, dparam);
  if (info > 0)
  {
    fprintf(stderr, "Numerics, FrictionContact3D_Newton failed, reached max. number of iterations without convergence. Error = %f\n", dparam[1]);
    exit(EXIT_FAILURE);
  }

  (*postSolver)(contact, reaction);
}

void frictionContact3D_Newton_free()
{
  F = NULL;
  jacobianF = NULL;
  updateSolver = NULL;
  postSolver = NULL;
  (*freeSolver)();
}

void frictionContact3D_Newton_computeError(int n, double* velocity, double*reaction, double * error)
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
