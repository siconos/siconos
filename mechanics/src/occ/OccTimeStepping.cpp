/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
  auto& dsg = *_nsds->dynamicalSystems();
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsiend;
  std::tie(dsi, dsiend) = dsg.vertices();

  for (; dsi != dsiend; ++dsi) {
    if (auto ds = std::dynamic_pointer_cast<OccBody>(dsg.bundle(*dsi))) {
      ds->updateShapes();
      ds->updateContactShapes();
    }
  }
}
