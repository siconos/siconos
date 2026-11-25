/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
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
