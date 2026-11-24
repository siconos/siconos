/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef OccTimeStepping_hpp
#define OccTimeStepping_hpp

#include "NonSmoothDynamicalSystem.hpp"
#include "TimeStepping.hpp"

namespace siconos::modeling {
class NonSmoothDynamicalSystem;
}

namespace siconos::mechanics::occ {
class OccTimeStepping : public siconos::simulation::TimeStepping {
 public:
  OccTimeStepping(std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds,
                  std::shared_ptr<siconos::simulation::TimeDiscretisation> td)
      : TimeStepping(nsds, td) {};

  virtual void updateWorldFromDS();
};
}  // namespace siconos::mechanics::occ
#endif
