/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "fc3d_GlockerFischerBurmeister_functions.h"

#include <stddef.h>  // for NULL

#include "FischerBurmeister.h"  // for jacobianPhi_FB, phi_FB
#include "fc3d_2NCP_Glocker.h"  // for computeFGlocker, computeJacobianFGlocker
#include "numerics_errors.h"

/* size of a block */
/* static int Fsize; */

/* static void F_GlockerFischerBurmeister(int sizeF, double* reaction, double* FVector, int
 * up2Date); */
/* static void jacobianF_GlockerFischerBurmeister(int sizeF, double* reaction, double*
 * jacobianFMatrix, int up2Date); */

/** writes \f$ F(z) \f$ using Glocker formulation and the Fischer-Burmeister function.
 */
void F_GlockerFischerBurmeister(int sizeF, double* reaction, double* FVector, int up2Date) {
  /* Glocker formulation */
  double* FGlocker = NULL;
  computeFGlocker(&FGlocker, up2Date);
  /* Note that FGlocker is a static var. in fc3d2NCP_Glocker and thus there is no memory
   allocation in the present file.
  */

  /* Call Fisher-Burmeister function => fill FVector */
  phi_FB(sizeF, reaction, FGlocker, FVector);
  FGlocker = NULL;
}

/** writes \f$ \nabla_z F(z) \f$  using Glocker formulation and the Fischer-Burmeister
 * function.
 */
void jacobianF_GlockerFischerBurmeister(int sizeF, double* reaction, double* jacobianFMatrix,
                                        int up2Date) {
  /* Glocker formulation */
  double *FGlocker = NULL, *jacobianFGlocker = NULL;
  computeFGlocker(&FGlocker, up2Date);
  computeJacobianFGlocker(&jacobianFGlocker, up2Date);
  /* Note that FGlocker and jacobianFGlocker are static var. in fc3d2NCP_Glocker and thus there
   is no memory allocation in the present file.
  */

  /* Call Fisher-Burmeister function => fill jacobianFMatrix */
  jacobianPhi_FB(sizeF, reaction, FGlocker, jacobianFGlocker, jacobianFMatrix);
  FGlocker = NULL;
  jacobianFGlocker = NULL;
}
