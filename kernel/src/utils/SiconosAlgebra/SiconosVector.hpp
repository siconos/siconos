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

#include "EigenInclude.hpp"  // Must be included before Eigen/Core
//
#include <Eigen/Core>

namespace siconos::algebra {

using SiconosVector = Eigen::Matrix<double_t, Eigen::Dynamic, 1, Eigen::ColMajor>;
using SiconosVector3 = Eigen::Matrix<double_t, 3, 1, Eigen::ColMajor>;
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

void concatenateVectors(SiconosVector& target, const SiconosVector& a, const SiconosVector& b);
std::shared_ptr<SiconosVector> concatenateVectors(const SiconosVector& a,
                                                  const SiconosVector& b);

}  // namespace siconos::algebra

#endif