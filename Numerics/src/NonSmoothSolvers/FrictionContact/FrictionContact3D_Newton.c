/* Siconos-Numerics, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
2 * the Free Software Foundation; either version 2 of the License, or
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

#include "NonSmoothNewton.h"
#include "NCP_Solvers.h"
#include "FrictionContact3D_Solvers.h"
#include <stdlib.h>
#include <stdio.h>
#include "FischerBurmeister.h"

/* Pointer to function used to update the solver, to formalize the local problem for example. */
typedef void (*UpdateSolverPtr)(int, double*);

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


void frictionContact3D_Newton_initialize(FrictionContactProblem* problem, FrictionContactProblem* localproblem,    SolverOptions * localsolver_options)
{

  /*
     Initialize solver (Connect F and its jacobian, set local size ...) according to the chosen formulation.
  */




  /* Alart-Curnier formulation */
  if (localsolver_options->solverId == SICONOS_FRICTION_3D_AlartCurnierNewton)
  {
    Fsize = 3;
    frictionContact3D_AC_initialize(problem, localproblem);
    F = &F_AC;
    jacobianF = &jacobianF_AC;
    /*     updateSolver = &frictionContact3D_AC_update; */
    postSolver = &frictionContact3D_AC_post;
    freeSolver = (FreeSolverPtr)&frictionContact3D_AC_free;

  }
  else if (localsolver_options->solverId == SICONOS_FRICTION_3D_DampedAlartCurnierNewton)
  {
    Fsize = 3;
    frictionContact3D_AC_initialize(problem, localproblem);
    F = &F_AC;
    jacobianF = &jacobianF_AC;
    /*     updateSolver = &frictionContact3D_AC_update; */
    postSolver = &frictionContact3D_AC_post;
    freeSolver = (FreeSolverPtr)&frictionContact3D_AC_free;

  }




  /* Glocker formulation - Fischer-Burmeister function used in Newton */
  else if (localsolver_options->solverId == SICONOS_FRICTION_3D_NCPGlockerFBNewton)
  {
    Fsize = 5;
    NCPGlocker_initialize(problem, localproblem);
    F = &F_GlockerFischerBurmeister;
    jacobianF = &jacobianF_GlockerFischerBurmeister;
    /*     updateSolver = &NCPGlocker_update; */
    postSolver = &NCPGlocker_post;
    freeSolver = (FreeSolverPtr)&NCPGlocker_free;
  }
  else
  {
    fprintf(stderr, "Numerics, FrictionContact3D_nsgs failed. Unknown formulation type.\n");
    exit(EXIT_FAILURE);
  }
}

int frictionContact3D_Newton_solve(FrictionContactProblem* localproblem, double* reaction, SolverOptions * options)
{


  /*  (*updateSolver)(contact, reaction); */


  double * reactionBlock = reaction;

  int * iparam = options->iparam;
  double * dparam = options->dparam;

  int info;
  if (options->solverId == SICONOS_FRICTION_3D_AlartCurnierNewton)
  {
    info = LocalNonsmoothNewtonSolver(localproblem, reactionBlock, iparam, dparam);
  }
  else if (options->solverId == SICONOS_FRICTION_3D_DampedAlartCurnierNewton)
  {
    info = DampedLocalNonsmoothNewtonSolver(localproblem, reactionBlock, iparam, dparam);
  }




  else
  {
    info = nonSmoothDirectNewton(Fsize, reactionBlock, &F, &jacobianF, iparam, dparam);
  }
  if (info > 0)
  {
    if (verbose > 0)
    {
      printf("Numerics, FrictionContact3D_Newton, warning. reached max. number of iterations without convergence. Error = %12.8e\n", dparam[1]);
      /* note : exit on failure should be done in DefaultCheckSolverOutput */
    }
  }
  return info;
  /*  (*postSolver)(contact,reaction); */
}

void frictionContact3D_Newton_free(FrictionContactProblem* localproblem)
{
  free(localproblem->M->matrix0);
  localproblem->M->matrix0 = NULL;
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
