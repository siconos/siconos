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

#ifndef Bullet1DR_hpp
#define Bullet1DR_hpp

#include "BulletDeclarations.h"
#include "NewtonEuler1DR.hpp"

namespace siconos::collision::bullet {

class Bullet1DR : public siconos::modeling::NewtonEuler1DR {
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
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2,
      Eigen::Ref<siconos::algebra::SiconosVector> y) override;
};
}  // namespace siconos::collision::bullet
#endif
