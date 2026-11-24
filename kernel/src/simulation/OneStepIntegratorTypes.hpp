/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/*! \file OneStepIntegratorTypes.hpp
  \brief enum of the available types for one-step time integrators.
*/

#ifndef OSITYPES_HPP
#define OSITYPES_HPP

namespace siconos::integrators {
/** List of OneStepIntegrator types/ids */
enum class IntegratorType {
  /** Euler-Moreau scheme*/
  EULERMOREAUOSI,
  /** Moreau-Jean scheme */
  MOREAUJEANOSI,
  /** ?? */
  MOREAUJEANGOSI,
  /** Lsodar (ode solver with rootfinding process) */
  LSODAROSI,
  /** odepack HEM5 (Hairer) solver */
  HEM5OSI,
  /** Moreau-Jean with direct projection */
  MOREAUDIRECTPROJECTIONOSI,
  /** Moreau-Jean with combined projection */
  MOREAUCOMBINEDPROJECTIONOSI,
  /** Intregrator based on discontinuous Galerkin methods*/
  D1MINUSLINEAROSI,
  /** Schatzman-Paoli scheme */
  SCHATZMANPAOLIOSI,
  /** Zero-order osi (?) */
  ZOHOSI,
  /** Newmark-like scheme*/
  NEWMARKALPHAOSI,
  /** Moreau-Jean-Bilbao */
  MOREAUJEANBILBAOOSI
};

}  // namespace siconos::integrators
#endif
