// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef Bullet1DR_hpp
#define Bullet1DR_hpp

#include "BulletDeclarations.h"
#include "NewtonEulerNormalVelocityR.hpp"

namespace siconos::collision::bullet {

class Bullet1DR : public siconos::modeling::NewtonEulerNormalVelocityR {
 private:
  ACCEPT_SERIALIZATION(Bullet1DR);

  std::shared_ptr<btManifoldPoint> _contactPoints{nullptr};

 public:
  Bullet1DR(std::shared_ptr<btManifoldPoint>);
  virtual ~Bullet1DR() noexcept = default;

  std::shared_ptr<btManifoldPoint> contactPoint() const { return _contactPoints; };

  /**
     to compute the output y = h(q) of the Relation

      \param[in] q1 generalized coordinates vector of the fist dynamical system involved
      in the relation
      \param[in] q2 generalized coordinates vector of the second dynamical system
      involved in the relation
      \param[in,out] y the resulting vector
  */
  virtual void computeh(
      const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
      Eigen::Ref<siconos::algebra::SiconosVector> y) override;
};
}  // namespace siconos::collision::bullet
#endif
