/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef MBTB_TS_H
#define MBTB_TS_H

#include "TimeStepping.hpp"  // Base class

namespace siconos::mechanisms {
/**
 * \brief This class implements the time stepping of a multi-bodies system.
 * It inherits from siconos::TimeStepping.
 * It consists in update the CAD word during the simulation.
 */
class MBTB_TimeStepping : public siconos::simulation::TimeStepping {
 public:
  /** Standard constructor
   *
   *  \param nsds NonSmoothDynamicalSystem to be simulated
   *  \param td pointer to a timeDiscretisation used in the integration
   *  \param osi one step integrator (default none)
   *  \param osnspb one step non smooth problem (default none)
   */
  MBTB_TimeStepping(std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds,
                    std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
                    std::shared_ptr<siconos::integrators::OneStepIntegrator> osi,
                    std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem> osnspb)
      : siconos::simulation::TimeStepping(nsds, td, osi, osnspb) {}

  virtual ~MBTB_TimeStepping() noexcept = default;

  /**  Update CAD model */
  virtual void updateWorldFromDS();
};
}  // namespace siconos::mechanisms
#endif
