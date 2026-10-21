/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2025 INRIA.
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
#ifndef FRICTIONCONTACT3D_compute_rho_H
#define FRICTIONCONTACT3D_compute_rho_H

/*!\file fc3d_compute_rho.h

  Functions for computing the rho parameter used in Alart-Curnier
  formulations for 3D frictional contact problems.

  The rho vector contains, for each contact, the diagonal entries of
  the preconditioner built from the local 3x3 diagonal block of the
  problem matrix.

 */
#include "FrictionContactProblem.h"   // for FrictionContactProblem
#include "FrictionContact_options.h"  // for SICONOS_FRICTION_3D_NSN_RHO_STRATEGY values
#include "NumericsFwd.h"              // for FrictionContactProblem

#if defined(__cplusplus)
extern "C" {
#endif

/** Compute rho for one contact using the split spectral norm
 *  conditioning strategy.
 *
 *  \param[in] localproblem local 3D friction contact problem (one contact)
 *  \param[out] rho rho values for the contact (size 3)
 */
void compute_rho_split_spectral_norm_cond(FrictionContactProblem* localproblem, double* rho);

/** Compute rho for one contact using the split spectral norm strategy.
 *
 *  \param[in] localproblem local 3D friction contact problem (one contact)
 *  \param[out] rho rho values for the contact (size 3)
 */
void compute_rho_split_spectral_norm(FrictionContactProblem* localproblem, double* rho);

/** Compute rho for one contact using the spectral norm strategy.
 *
 *  \param[in] localproblem local 3D friction contact problem (one contact)
 *  \param[out] rho rho values for the contact (size 3)
 */
void compute_rho_spectral_norm(FrictionContactProblem* localproblem, double* rho);

/** Compute rho values for all contacts of a 3D friction contact problem.
 *
 *  The output array \c rho must be allocated with size
 *  <tt>3 * problem->numberOfContacts</tt>. The values are stored
 *  contiguously per contact:
 *  <tt>[rho_N^0, rho_T^0, rho_S^0, rho_N^1, rho_T^1, rho_S^1, ...]</tt>.
 *
 *  \param[in] problem the global 3D friction contact problem
 *  \param[in] rho_strategy strategy for computing rho, see
 *             SICONOS_FRICTION_3D_NSN_RHO_STRATEGY_ENUM
 *  \param[in] rho_value constant rho value used when rho_strategy is
 *             SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_CONSTANT
 *  \param[out] rho vector of rho values (size 3 * numberOfContacts)
 *  \return 0 on success, a numerics error code otherwise
 */
int fc3d_compute_rho(FrictionContactProblem* problem, int rho_strategy, double rho_value,
                     double* rho);

#if defined(__cplusplus)
}
#endif

#endif
