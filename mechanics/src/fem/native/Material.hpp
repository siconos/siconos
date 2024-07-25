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

/*! \file Material.hpp

 */
#ifndef MATERIAL_H
#define MATERIAL_H

#include <iostream>
#include <vector>

namespace siconos::mechanics::fem {

/** Tag analysis types */
enum class AnalysisType2D { plane_strain, plane_stress, axysymmetric };

/** A simple class for material

    Default = steel
 */
class Material {
 protected:
  /** mass density */
  double _massDensity{7850.};

  /** Young Modulus */
  double _elasticYoungModulus{2.1e11};

  /** Poison coefficient */
  double _poissonCoefficient{0.3};

  /** Analysis type in 2D */
  AnalysisType2D _analysisType2D{AnalysisType2D::plane_strain};

  /** thickness in 2D plane stress analysis */
  double _thickness{1.0};

  /** default constructor */
  Material() = default;

 public:
  /** constructor */
  Material(double massDensity, double ElasticYoungModulus, double poissonCoefficient)
      : _massDensity(massDensity),
        _elasticYoungModulus(ElasticYoungModulus),
        _poissonCoefficient(poissonCoefficient) {}

  /** destructor */
  ~Material() noexcept = default;

  auto massDensity() { return _massDensity; }

  auto elasticYoungModulus() { return _elasticYoungModulus; }

  auto poissonCoefficient() { return _poissonCoefficient; }

  auto analysisType2D() { return _analysisType2D; }

  auto thickness() { return _thickness; }
};
}  // namespace siconos::mechanics::fem
#endif
