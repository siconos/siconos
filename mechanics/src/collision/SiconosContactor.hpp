/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

/*! \file SiconosContactor.hpp
  \brief Definition of an abstract contactor and a set of such objects.
*/

#ifndef SiconosContactor_h
#define SiconosContactor_h

#include <memory>
#include <vector>

#include "SiconosSerialization.hpp"

namespace siconos::algebra {
class SiconosVector;

}  // namespace siconos::algebra

namespace siconos::collision {

class SiconosShape;

/** Class to hold the shape assigned to a body, and to associate each
 *  shape with an offset and collision group. */

class SiconosContactor {
 private:
  // Rule of five
  SiconosContactor() = delete;
  SiconosContactor(const SiconosContactor&) = delete;
  SiconosContactor(SiconosContactor&&) = delete;
  SiconosContactor& operator=(const SiconosContactor&) = delete;
  SiconosContactor& operator=(SiconosContactor&&) = delete;

 protected:
  ACCEPT_SERIALIZATION(SiconosContactor);

 public:
  SiconosContactor(std::shared_ptr<SiconosShape> shape,
                   std::shared_ptr<siconos::algebra::SiconosVector> offset = nullptr,
                   int collision_group = 0);

  std::shared_ptr<SiconosShape> shape{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> offset{nullptr};
  int collision_group{0};

  virtual ~SiconosContactor() noexcept = default;
};

class SiconosContactorSet : public std::vector<std::shared_ptr<SiconosContactor>> {
 protected:
  ACCEPT_SERIALIZATION(SiconosContactorSet);
  using iterator = std::vector<std::shared_ptr<SiconosContactor>>::iterator;

 public:
  void append(std::shared_ptr<SiconosContactor> b) { push_back(b); }
  void append(std::vector<std::shared_ptr<SiconosContactor>> b)
  {
    insert(end(), b.begin(), b.end());
  }
  void append(const SiconosContactorSet& b) { insert(end(), b.begin(), b.end()); }
  void append(const std::shared_ptr<SiconosContactorSet>& b)
  {
    insert(end(), b->begin(), b->end());
  }
};
}  // namespace siconos::collision
#endif /* SiconosContactor_h */
