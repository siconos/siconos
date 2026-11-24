/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef MBTB_TSPROJ_H
#define MBTB_TSPROJ_H

#include "TimeSteppingDirectProjection.hpp"  // Base class

namespace siconos::mechanisms {

/**
 * \brief This class implements the time stepping with projection of a
 * multi-bodies system. It inherits from Siconos::TimeSteppingDirectProjection.
 * It consists in update the CAD word during the simulation.
 */
class MBTB_TimeSteppingProj : public siconos::simulation::TimeSteppingDirectProjection {
 public:
  /** Constructor with the time-discretisation.
   *  \param td pointer to a timeDiscretisation used in the integration
   *          (linked to the model that owns this simulation)
   *  \param osi one step integrator (default none)
   *  \param osnspb_velo one step non smooth problem (default none)
   *  \param osnspb_pos one step non smooth problem (default none)
   *  \param level
   */
  MBTB_TimeSteppingProj(
      std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds,
      std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
      std::shared_ptr<siconos::integrators::OneStepIntegrator> osi,
      std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem> osnspb_velo,
      std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem> osnspb_pos,
      unsigned int level)
      : TimeSteppingDirectProjection(nsds, td, osi, osnspb_velo, osnspb_pos, level) {};

  virtual ~MBTB_TimeSteppingProj() noexcept = default;

  /**  Update CAD model */
  virtual void updateWorldFromDS();
};
}  // namespace siconos::mechanisms
#endif
