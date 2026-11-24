// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file RigidBody2dDS.hpp
  \brief Definition of an abstract 3D rigid body above NewtonEulerDS
*/

#ifndef RigidBody2dDS_h
#define RigidBody2dDS_h

#include <LagrangianLinearTIDS.hpp>

#include "SiconosContactor.hpp"
// #include "StorageTools.hpp"

namespace siconos::collision {

class RigidBody2dDS : public siconos::modeling::LagrangianLinearTIDS,
                      public std::enable_shared_from_this<RigidBody2dDS> {
 protected:
  ACCEPT_SERIALIZATION(RigidBody2dDS);

  /** a scalar mass in the case of RigidBody2dDS */
  double scalarMass_{0.};

  std::shared_ptr<const siconos::collision::SiconosContactorSet> contactors_{nullptr};
  bool _useContactorInertia{true};

  /** If false, bodies connected to this body by a joint will not
   * collide. See also NewtonEulerJointR::_allowSelfCollide */
  bool _allowSelfCollide{true};

 public:
  RigidBody2dDS(const siconos::algebra::SiconosVector3& position,
                const siconos::algebra::SiconosVector3& velocity, double mass, double inertia);

  virtual ~RigidBody2dDS() noexcept = default;

  double scalarMass() const { return scalarMass_; };

  void setUseContactorInertia(bool use) { _useContactorInertia = use; }

  bool useContactorInertia() { return _useContactorInertia; }

  /** Return the value of the _allowSelfCollide flag. */
  bool allowSelfCollide() { return _allowSelfCollide; }

  /** Set the value of the _allowSelfCollide flag. */
  void setAllowSelfCollide(bool x) { _allowSelfCollide = x; }

  /** \return the contactor set associated with this body */
  std::shared_ptr<const siconos::collision::SiconosContactorSet> contactors() const {
    return contactors_;
  }

  /** Provide a set of contactors to the body.
   *
   *  \param c A std::shared_ptr<SiconosContactorSet> */
  void setContactors(std::shared_ptr<siconos::collision::SiconosContactorSet> c) {
    contactors_ = std::move(c);
  }

  /** Make the base position of the contactors equal to the DS q vector.
   *
   *  \return a std::shared_ptr<siconos::algebra::SiconosVector> */
  virtual std::shared_ptr<siconos::algebra::SiconosVector> base_position() { return q(); }

  virtual void accept(modeling::dynamical_systems::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::collision
#endif /* RigidBody2dDS_h */
