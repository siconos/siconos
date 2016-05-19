/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/


#include "fc3d_GlockerFischerBurmeister_functions.h"
#include "NonSmoothNewton.h"
#include "fc3d_Solvers.h"
#include "FischerBurmeister.h"
/* size of a block */
/* static int Fsize; */

/* static void F_GlockerFischerBurmeister(int sizeF, double* reaction, double* FVector, int up2Date); */
/* static void jacobianF_GlockerFischerBurmeister(int sizeF, double* reaction, double* jacobianFMatrix, int up2Date); */



/** writes \f$ F(z) \f$ using Glocker formulation and the Fischer-Burmeister function.
 */
void F_GlockerFischerBurmeister(int sizeF, double* reaction, double* FVector, int up2Date)
{
  /* Glocker formulation */
  double* FGlocker = NULL;
  computeFGlocker(&FGlocker, up2Date);
  /* Note that FGlocker is a static var. in fc3d2NCP_Glocker and thus there is no memory allocation in
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
  /* Note that FGlocker and jacobianFGlocker are static var. in fc3d2NCP_Glocker and thus there is no memory allocation in
   the present file.
  */

  /* Call Fisher-Burmeister function => fill jacobianFMatrix */
  jacobianPhi_FB(sizeF, reaction, FGlocker, jacobianFGlocker, jacobianFMatrix);
  FGlocker = NULL;
  jacobianFGlocker = NULL;
}
