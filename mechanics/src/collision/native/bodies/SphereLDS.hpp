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

/*! \file SphereLDS.hpp
  \brief Definition of a 3D sphere as a LagrangianDS (with Euler
         Angles)
*/

#ifndef SphereLDS_h
#define SphereLDS_h

#include "LagrangianDS.hpp"

namespace siconos::collision::native::bodies {

class SphereLDS : public siconos::modeling::LagrangianDS,
                  public std::enable_shared_from_this<SphereLDS> {
 protected:
  ACCEPT_SERIALIZATION(SphereLDS);

  double radius{0.};
  double massValue{0.};
  double I{0.};

 public:
  SphereLDS(double, double, std::shared_ptr<siconos::algebra::SiconosVector>,
            std::shared_ptr<siconos::algebra::SiconosVector>);

  ~SphereLDS() noexcept = default;

  double getQ(unsigned int pos);

  double getVelocity(unsigned int pos);

  inline double getMassValue() const { return massValue; };

  inline double getRadius() const { return radius; };

  void computeFGyr(std::shared_ptr<siconos::algebra::SiconosVector>,
                   std::shared_ptr<siconos::algebra::SiconosVector>) override;

  void computeFGyr() override;

  void computeJacobianFGyrq() override;
  void computeJacobianFGyrqDot() override;

  void computeJacobianFGyrq(std::shared_ptr<siconos::algebra::SiconosVector>,
                            std::shared_ptr<siconos::algebra::SiconosVector>) override;
  void computeJacobianFGyrqDot(std::shared_ptr<siconos::algebra::SiconosVector>,
                               std::shared_ptr<siconos::algebra::SiconosVector>) override;
};
}  // namespace siconos::collision::native::bodies
#endif /* SphereLDS_h */
