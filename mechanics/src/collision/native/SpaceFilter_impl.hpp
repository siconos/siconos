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

/*! \file SpaceFilter_impl.hpp
 *  \brief implementation details for moving plans
 */

#ifndef SpaceFilter_impl_hpp
#define SpaceFilter_impl_hpp
#include <memory>

#include "SiconosSerialization.hpp"

namespace siconos::modeling {
class DynamicalSystem;
}  // namespace siconos::modeling

// Internal stuff. Not supposed to be used/viewed in users'interface.
namespace siconos::collision::native::internal {

class Hashed : public std::enable_shared_from_this<Hashed> {
 protected:
  ACCEPT_SERIALIZATION(Hashed);

  Hashed() = default;
  Hashed(const Hashed&) = delete;
  Hashed(Hashed&&) = delete;
  Hashed operator=(const Hashed&) = delete;
  Hashed operator=(Hashed&&) = delete;

 public:
  std::shared_ptr<siconos::modeling::DynamicalSystem> body{nullptr};
  int i{0};
  int j{0};
  int k{0};

  Hashed(std::shared_ptr<siconos::modeling::DynamicalSystem> body, int i, int j, int k = 0)
      : body(body), i(i), j(j), k(k){};

  Hashed(int i, int j, int k = 0) : i(i), j(j), k(k){};

  ~Hashed() noexcept = default;
};

bool operator==(std::shared_ptr<Hashed> const& a, std::shared_ptr<Hashed> const& b);

std::size_t hash_value(std::shared_ptr<Hashed> const& h);

}  // namespace siconos::collision::native::internal
#endif
