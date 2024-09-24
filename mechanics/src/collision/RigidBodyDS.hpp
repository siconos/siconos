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

/*! \file RigidBodyDS.hpp
  \brief Definition of an abstract 3D rigid body above NewtonEulerDS
*/

#ifndef RigidBodyDS_h
#define RigidBodyDS_h

#include "NewtonEulerDS.hpp"

namespace siconos::collision {
class SiconosContactorSet;
}

namespace siconos::collision {

class RigidBodyDS : public siconos::modeling::NewtonEulerDS,
                    public std::enable_shared_from_this<RigidBodyDS> {
 protected:
  ACCEPT_SERIALIZATION(RigidBodyDS);

  std::shared_ptr<siconos::collision::SiconosContactorSet> _contactors{nullptr};

  bool _useContactorInertia{true};

  /** If false, bodies connected to this body by a joint will not
   *  collide. See also NewtonEulerJointR::_allowSelfCollide */
  bool _allowSelfCollide = true;

  std::shared_ptr<siconos::algebra::SiconosVector> _qExtrapolated{nullptr};

 public:
  RigidBodyDS(Eigen::Ref<siconos::algebra::SiconosVector> position,
              Eigen::Ref<siconos::algebra::SiconosVector> velocity, double mass,
              Eigen::Ref<siconos::algebra::SiconosMatrix> inertia);

  virtual ~RigidBodyDS() noexcept = default;

  void setUseContactorInertia(bool use) { _useContactorInertia = use; }

  bool useContactorInertia() { return _useContactorInertia; }

  /** Return the value of the _allowSelfCollide flag. */
  bool allowSelfCollide() { return _allowSelfCollide; }

  /** Set the value of the _allowSelfCollide flag. */
  void setAllowSelfCollide(bool x) { _allowSelfCollide = x; }

  /** Access the contactor set associated with this body.
   *
   *  \return A std::shared_ptr<SiconosContactorSet> */
  std::shared_ptr<siconos::collision::SiconosContactorSet> contactors() const {
    return _contactors;
  }

  /** Provide a set of contactors to the body.
   *
   *  \param c A std::shared_ptr<SiconosContactorSet> */
  void setContactors(std::shared_ptr<siconos::collision::SiconosContactorSet> c) {
    _contactors = c;
  }

  /** Make the base position of the contactors equal to the DS q vector.
   *
   *  \return a std::shared_ptr<siconos::algebra::SiconosVector> */
  virtual std::shared_ptr<siconos::algebra::SiconosVector> base_position() { return q(); }

  virtual std::shared_ptr<siconos::algebra::SiconosVector> base_extrapolated_position() {
    return _qExtrapolated;
  };
  virtual void compute_extrapolated_position(double extrapolationCoefficient);

  void acceptSP(std::shared_ptr<siconos::internal::SiconosVisitor> tourist) const override;
};
}  // namespace siconos::collision
#endif /* RigidBodyDS_h */
