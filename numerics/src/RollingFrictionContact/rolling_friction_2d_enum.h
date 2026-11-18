/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef ROLLINGFRICTIONCONTACT2DEnum_H
#define ROLLINGFRICTIONCONTACT2DEnum_H

/*!\file rolling_friction_2d_enum.h
  \brief The one-contact enumerative solver of 2D rolling friction
  (SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ENUM).

  The local problem u = W r + q, W 3 x 3, with r = (r_n, r_t, r_r) in the
  cone |r_t| <= mu r_n, |r_r| <= mu_r r_n, is piecewise linear: the cone is
  polyhedral.  Its cases: the opening (r = 0, q_n >= 0), then a closed
  contact (u_n = 0) where the tangent row sticks (u_t = 0) or slides
  (r_t = +-mu r_n) and the rolling row sticks (u_r = 0) or rolls
  (r_r = +-mu_r r_n).  Each case is a linear system of size 1 to 3.

  The stick-stick case comes from the inverse of the block, then the cases
  predicted by its trial (each row outside its bound slides or rolls with
  the sign of its trial reaction) in closed form, each with the conditions
  it leaves open (r_n >= 0, the bound of a sticking row, the velocity of a
  sliding or rolling row opposed to its reaction).  If none is accepted,
  every case is checked against the whole law in a fixed order: the first
  accepted is the solution; if round-off rejects every case, the case of
  smallest violation is kept and the violation is the local residual.

  A block without rolling row (zero third row and column: a contact whose
  relation has no rolling row) is a 2D Coulomb contact: r_r = 0, the
  inverse is the one of its normal and tangent block, and the same cases
  are solved on the normal and tangent rows.

  In rolling_friction_2d_nsgs, the diagonal blocks of W and their inverses
  are extracted once per solve (rolling_friction_2d_enum_initialize, into
  the dWork of the local solver options) and the local problem of a
  contact is built from them (rolling_friction_2d_enum_update); without
  them, solve enumerates every case.
*/

#include "NumericsFwd.h"  // for RollingFrictionContactProblem, SolverOptions

#if defined(__cplusplus)
extern "C" {
#endif

/** solve the local problem of one contact
 *
 *  \param localproblem the local problem (dimension 3, one contact)
 *  \param[in,out] reaction the local reaction, 3 values
 *  \param options the local solver options: the violation of the chosen
 *         case goes to dparam[SICONOS_DPARAM_RESIDU]
 *  \return 0 if a case satisfies the law to round-off, 1 otherwise
 */
int rolling_friction_2d_enum_solve(RollingFrictionContactProblem* localproblem,
                                   double* reaction, SolverOptions* options);

void rolling_friction_2d_enum_set_default(SolverOptions* options);

/** extract the diagonal blocks of the problem and their inverses into
 *  options->dWork (18 values per contact) */
void rolling_friction_2d_enum_initialize(RollingFrictionContactProblem* problem,
                                         RollingFrictionContactProblem* localproblem,
                                         SolverOptions* options);

/** the local problem of a contact: its diagonal block from options->dWork,
 *  its q from the reactions of the other contacts, mu and mu_r */
void rolling_friction_2d_enum_update(int contact, RollingFrictionContactProblem* problem,
                                     RollingFrictionContactProblem* localproblem,
                                     double* reaction, SolverOptions* options);

void rolling_friction_2d_enum_free(RollingFrictionContactProblem* problem,
                                   RollingFrictionContactProblem* localproblem,
                                   SolverOptions* options);

#if defined(__cplusplus)
}
#endif

#endif
