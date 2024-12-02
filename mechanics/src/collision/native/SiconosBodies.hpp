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

/*! \file SiconosBodies.hpp
  \brief SiconosBodies class - model + plans + space filter
*/
#ifndef SiconosBodies_hpp
#define SiconosBodies_hpp

#include <memory>

#include "FMatrix.hpp"  // For FMatrix
#include "SiconosSerialization.hpp"

namespace siconos::algebra {
class SiconosMatrix;
}

namespace siconos::simulation {
class Simulation;
}

namespace siconos::collision::native {

class SpaceFilter;

/** SiconosBodies : a Siconos Model, some plans and space filtering capabilities
 */

class SiconosBodies {
 protected:
  ACCEPT_SERIALIZATION(SiconosBodies);

  std::shared_ptr<siconos::collision::native::FMatrix> _moving_plans{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _plans{nullptr};
  std::shared_ptr<siconos::simulation::Simulation> _sim{nullptr};
  std::shared_ptr<siconos::collision::native::SpaceFilter> _playground{nullptr};

 public:
  /** destructor */
  virtual ~SiconosBodies() noexcept = default;

  virtual void init() = 0;

  virtual void compute();

  std::shared_ptr<siconos::simulation::Simulation> simulation() { return _sim; }

  std::shared_ptr<siconos::collision::native::FMatrix> movingPlans() { return _moving_plans; }
  std::shared_ptr<siconos::algebra::SiconosMatrix> plans() { return _plans; }

  std::shared_ptr<siconos::collision::native::SpaceFilter> spaceFilter()
  {
    return _playground;
  };
};
}  // namespace siconos::collision::native
#endif  // SiconosBodies_hpp
