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

/*! \file StaticBody.hpp
  \brief Definition of an abstract 3D rigid body above NewtonEulerDS
*/

#ifndef StaticBody_h
#define StaticBody_h

#include <memory>

#include "SiconosContactor.hpp"
#include "SiconosVector.hpp"

namespace siconos::collision {

class StaticBody : public std::enable_shared_from_this<StaticBody> {
 private:
  // Rule of five
  StaticBody(const StaticBody&) = delete;
  StaticBody(StaticBody&&) = delete;
  StaticBody& operator=(const StaticBody&) = delete;
  StaticBody& operator=(StaticBody&&) = delete;

  /** the contactor set associated with this body */
  std::shared_ptr<const SiconosContactorSet> contactors_{nullptr};

 public:
  StaticBody() = default;
  std::shared_ptr<siconos::algebra::SiconosVector> base{nullptr};
  int number{0};
  virtual ~StaticBody() noexcept = default;

  /** \return the contactor set associated with this body. */
  std::shared_ptr<const siconos::collision::SiconosContactorSet> contactors() const {
    return contactors_;
  }

  /** Provide a set of contactors to the body.
   *
   *  \param c A std::shared_ptr<SiconosContactorSet> */
  void setContactors(std::shared_ptr<const siconos::collision::SiconosContactorSet> c) {
    contactors_ = std::move(c);
  }
};

}  // namespace siconos::collision
#endif /* StaticBody_h */
