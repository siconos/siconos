// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

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

  /** out-of-plane thickness in 2D analysis */
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
      : massDensity_{massDensity},
        elasticYoungModulus_{E},
        poisson_sRatio_{nu},
        radius_{radius} {}

  /** @brief constructor of a 2D continuum material
   * @param massDensity mass density (per unit volume)
   * @param E elastic Young modulus
   * @param nu Poisson ratio
   * @param analysisType plane strain (default of the other constructor) or plane stress
   * @param thickness out-of-plane thickness: the stiffness K, the mass M and the strain
   *        operator B (scaled by sqrt(area * thickness)) are integrated over it
   */
  constexpr Material(double massDensity, double E, double nu, AnalysisType2D analysisType,
                     double thickness = 1.0)
      : massDensity_{massDensity},
        elasticYoungModulus_{E},
        poisson_sRatio_{nu},
        analysisType2D_{analysisType},
        thickness_{thickness} {}

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
