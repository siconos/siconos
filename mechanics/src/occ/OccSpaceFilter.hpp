// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef OccSpaceFilter_hpp
#define OccSpaceFilter_hpp

#include "SpaceFilter.hpp"

namespace siconos::mechanics::occ {
class OccSpaceFilter : public siconos::collision::native::SpaceFilter {
 public:
  virtual void updateInteractions(std::shared_ptr<siconos::simulation::Simulation>){};
};
}  // namespace siconos::mechanics::occ
#endif
