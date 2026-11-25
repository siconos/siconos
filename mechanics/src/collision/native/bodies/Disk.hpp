// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/** \file Disk.hpp
 */

#ifndef Disk_H
#define Disk_H

#include "CircularDS.hpp"

namespace siconos::collision::native::bodies {
/** \class Disk
 *  \brief Definition of a 2D disk - Inherits from LagrangianDS
 */

class Disk : public CircularDS, public std::enable_shared_from_this<Disk> {
 private:
  ACCEPT_SERIALIZATION(Disk);

 public:
  /** constructor from initial state and velocity
   *
   *  \param R radius
   *  \param m mass
   *  \param position initial coordinates
   *  \param velocity initial velocity
   */
  Disk(double radius, double mass, const siconos::algebra::SiconosVector3& position,
       const siconos::algebra::SiconosVector3& velocity);

  /** destructor */
  virtual ~Disk() noexcept = default;

  virtual void acceptSP(modeling::dynamical_systems::Visitor& tourist) override {
    tourist.visit(shared_from_this());
  }
  virtual void accept(modeling::dynamical_systems::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::collision::native::bodies
#endif /* Disk_H */
