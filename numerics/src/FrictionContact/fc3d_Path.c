/* Siconos-Numerics, Copyright INRIA 2005-2012.
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

#include "NonSmoothNewton.h"
#include "fc3d_Solvers.h"
#include "SiconosBlas.h"

#include <stdlib.h>
#include <stdio.h>
#include "Friction_cst.h"

#include "PathAlgebra.h"
#include "NonlinearComplementarityProblem.h"

#include "NCP_Solvers.h"

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

/* Pointer to function used to update the solver, to formalize the local problem for example. */
typedef void (*UpdateSolverPtr)(int, double*);


/* size of a block */

/** writes \f$ F(z) \f$ using Glocker formulation
 */
void F_GlockerPath(void* env, int sizeF, double* reaction, double* FVector)
{
  /* Glocker formulation */
  int up2Date = 0;
  double* FGlocker = NULL;
  computeFGlocker(&FGlocker, up2Date);
  /* Note that FGlocker is a static var. in fc3d2NCP_Glocker and thus there is no memory allocation in
     the present file.
  */

  /* TMP COPY: review memory management for FGlocker ...*/
  cblas_dcopy(sizeF , FGlocker , 1, FVector , 1);
  FGlocker = NULL;
}

/** writes \f$ \nabla_z F(z) \f$  using Glocker formulation and the Fischer-Burmeister function.
 */
void jacobianF_GlockerPath(void* env, int sizeF, double* reaction, NumericsMatrix* jacobianFMatrix)
{
  int up2Date = 0;
  /* Glocker formulation */
  double* FGlocker = NULL, *jacobianFGlocker = NULL;
  computeFGlocker(&FGlocker, up2Date);
  computeJacobianFGlocker(&jacobianFGlocker, up2Date);
  /* Note that FGlocker and jacobianFGlocker are static var. in fc3d2NCP_Glocker and thus there is no memory allocation in
   the present file.
  */

  FGlocker = NULL;
  jacobianFGlocker = NULL;
}


void fc3d_Path_initialize(FrictionContactProblem* problem, FrictionContactProblem* localproblem, SolverOptions * localsolver_options)
{

  /*
     Initialize solver (Connect F and its jacobian, set local size ...) according to the chosen formulation.
  */

  /* Glocker formulation */
  if (localsolver_options->solverId == SICONOS_FRICTION_3D_NCPGlockerFBPATH)
  {
    NCPGlocker_initialize(problem, localproblem);
  }
  else
  {
    fprintf(stderr, "Numerics, fc3d_Path failed. Unknown formulation type.\n");
    exit(EXIT_FAILURE);
  }
}

int fc3d_Path_solve(FrictionContactProblem * localproblem , double* reaction, SolverOptions * options)
{

  NonlinearComplementarityProblem NCP_struct = {
    5,
    &F_GlockerPath,
    &jacobianF_GlockerPath,
    NULL,
    NULL
  };

  double Fvec[5];
  int info;
  ncp_path(&NCP_struct, reaction, Fvec, &info, options);
  if (info > 0)
  {
    fprintf(stderr, "Numerics, fc3d_Path failed");
    exit(EXIT_FAILURE);
  }
  return info;
  /*   (*postSolver)(contact,reaction); */
}

void fc3d_Path_free()
{
  NCPGlocker_free();
}

void fc3d_Path_computeError(int n, double* velocity, double* reaction, double * error)
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
