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

/*! \file SiconosVector.hpp
 */

#ifndef __SiconosVector__
#define __SiconosVector__

#include "EigenInclude.hpp"  // IWYU pragma: keep - Must be included before Eigen/Core
//
#include <Eigen/Core>
#include <cassert>
#include <memory>
#include <utility>  // in_range
namespace siconos::algebra {

using SiconosVector = Eigen::Matrix<double_t, Eigen::Dynamic, 1, Eigen::ColMajor>;
using SiconosVector2 = Eigen::Vector2d;
using SiconosVector3 = Eigen::Vector3d;
using SiconosVector4 = Eigen::Matrix<double_t, 4, 1, Eigen::ColMajor>;
using SiconosVector7 = Eigen::Matrix<double_t, 7, 1, Eigen::ColMajor>;
using SiconosVector6 = Eigen::Matrix<double_t, 6, 1, Eigen::ColMajor>;
using MapVectorType = Eigen::Map<SiconosVector>;
using MapVector3Type = Eigen::Map<SiconosVector3>;
using MapVector6Type = Eigen::Map<SiconosVector6>;
using MapVector7Type = Eigen::Map<SiconosVector7>;
using ConstMapVectorType = Eigen::Map<const SiconosVector>;
using ConstMapVector3Type = Eigen::Map<const SiconosVector3>;
using ConstMapVector6Type = Eigen::Map<const SiconosVector6>;
using ConstMapVector7Type = Eigen::Map<const SiconosVector7>;
using Index = Eigen::Index;

/** Vector of indices */
using IndicesVector = Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1, Eigen::ColMajor>;

namespace blocks {
/** Vector (std) of pointers to SiconosVector */
using SharedVector = std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>;

/** Vector (std) of SiconosVector */
using Vector = std::vector<siconos::algebra::SiconosVector>;

/** type of index used for vectors of vectors  */
using size_type = SharedVector::size_type;  // Same for Vector

}  // namespace blocks

void concatenateVectors(SiconosVector& target, const SiconosVector& a, const SiconosVector& b);
std::shared_ptr<SiconosVector> concatenateVectors(const SiconosVector& a,
                                                  const SiconosVector& b);

/** Utility function to ensure proper cast from size_t (or equivalent) to Eigen indices type */
template <std::integral T>  // only for  integer values
[[nodiscard]] Eigen::Index to_index(T value) {
  // if (!std::in_range<siconos::algebra::Index>(value))
  //   throw std::overflow_error("value too large for siconos::algebra::Index");
  assert(std::in_range<siconos::algebra::Index>(value));
  return static_cast<siconos::algebra::Index>(value);
}

/** Utility function to ensure proper cast from Eigen indices type to size_t (or equivalent) */
template <std::unsigned_integral T>  // Nodiscard to force a return value (else warning)
[[nodiscard]] T to_unsigned(Eigen::Index value) {
  // if (!std::in_range<T>(value))
  //   throw std::overflow_error(
  //       "Eigen::Index value cannot be represented as an unsigned integer");
  assert(std::in_range<T>(value));
  return static_cast<T>(value);
}

}  // namespace siconos::algebra

#endif