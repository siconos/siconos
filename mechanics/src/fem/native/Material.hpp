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

/*! \file Material.hpp */
#ifndef MATERIAL_H
#define MATERIAL_H
#include <numbers>  // for pi

namespace siconos::mechanics::fem {

/** Tag analysis types */
enum class AnalysisType2D { plane_strain, plane_stress, axysymmetric };

/** To describe the properties of a material (FE framework)
    Default = steel
 */
class Material {
 protected:
  /** mass density */
  double massDensity_{7850.};

  /** Young Modulus */
  double elasticYoungModulus_{2.1e11};

  /** Poisson coefficient */
  double poisson_sRatio_{0.3};

  /** Analysis type in 2D */
  AnalysisType2D analysisType2D_{AnalysisType2D::plane_strain};

  /** thickness in 2D plane stress analysis */
  double thickness_{1.0};

  /** Shear Modulus **/
  double shearModulus_{elasticYoungModulus_ / (2. * (1. + poisson_sRatio_))};

  /** Radius for circular Beams */
  double radius_{1.0};

  /** Moment of Inertia (assumed equal in x and y: Iy=Iz) */
  double momentOfInertia_{std::numbers::pi * radius_ * radius_ * radius_ * radius_ / 4.0};

  /** Cross section area  */
  double crossSectionArea_{std::numbers::pi * radius_ * radius_};

  /** Torsional constant (Polar moment = Ix) **/
  double torsionalConstant_{2 * momentOfInertia_};

  // Rule of five
  Material() = default;  // --> steel

 public:
  Material(const Material&) = default;
  Material(Material&&) = default;
  Material& operator=(const Material&) = default;
  Material& operator=(Material&&) = default;

  /** @brief constructor
   * @param massDensity
   * @param E elastic Young modulus
   * @param nu Poisson ratio
   */
  constexpr Material(double massDensity, double E, double nu, double radius = 1.)
      : massDensity_(massDensity), elasticYoungModulus_(E), poisson_sRatio_(nu) {}

  /** destructor */
  ~Material() noexcept = default;

  auto massDensity() const { return massDensity_; }

  auto elasticYoungModulus() const { return elasticYoungModulus_; }

  auto Poisson_s_ratio() const { return poisson_sRatio_; }

  auto analysisType2D() const { return analysisType2D_; }

  auto thickness() const { return thickness_; }

  auto shearModulus() const { return shearModulus_; }

  auto momentOfInertia() const { return momentOfInertia_; }

  auto crossSectionArea() const { return crossSectionArea_; }

  auto secondMomentOfArea() const { return torsionalConstant_; }
};

inline constexpr Material Steel{7850, 2.1e11, 0.3};

}  // namespace siconos::mechanics::fem
#endif
