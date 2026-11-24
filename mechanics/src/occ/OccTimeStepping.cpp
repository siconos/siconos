// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "OccTimeStepping.hpp"

#include <memory>

#include "NonSmoothDynamicalSystem.hpp"
#include "OccBody.hpp"
#include "SiconosVector.hpp"
#include "SimulationGraphs.hpp"

// namespace siconos::mechanics::occ::internal {
// struct UpdateShapes {
//   template <typename T>
//   void operator()(const T& ds) {
//     const_cast<T&>(ds).updateShapes();
//     const_cast<T&>(ds).updateContactShapes();
//   }
// };
// }  // namespace siconos::mechanics::occ::internal

void siconos::mechanics::occ::OccTimeStepping::updateWorldFromDS() {
  auto& dsg = *nonSmoothDynamicalSystem()->dynamicalSystems();
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsiend;
  std::tie(dsi, dsiend) = dsg.vertices();

  for (; dsi != dsiend; ++dsi) {
    if (auto ds = std::dynamic_pointer_cast<OccBody>(dsg.bundle(*dsi))) {
      ds->updateShapes();
      ds->updateContactShapes();
    }
  }
}
