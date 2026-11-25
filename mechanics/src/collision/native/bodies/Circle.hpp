// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file Circle.hpp
 */
/** \class Circle
 * \brief Definition of a 2D Circle - Inherits from CircularDS
 */

#ifndef Circle_H
#define Circle_H

#include "CircularDS.hpp"

namespace siconos::collision::native::bodies {

class Circle : public CircularDS, public std::enable_shared_from_this<Circle> {
 private:
  ACCEPT_SERIALIZATION(Circle);

 public:
  /** constructor from initial state and velocity
   *
   *  \param R radius
   *  \param m mass
   *  \param position initial coordinates
   *  \param velocity initial velocity
   */
  Circle(double radius, double mass, const siconos::algebra::SiconosVector3& qinit,
         const siconos::algebra::SiconosVector3& vinit);

  /** destructor */
  virtual ~Circle() noexcept = default;

  virtual void acceptSP(modeling::dynamical_systems::Visitor& tourist) override {
    tourist.visit(shared_from_this());
  }
  virtual void accept(modeling::dynamical_systems::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::collision::native::bodies
#endif /* Circle_H */
