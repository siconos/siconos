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
using SiconosVector7 = Eigen::Matrix<double_t, 7, 1, Eigen::ColMajor>;
using SiconosVector6 = Eigen::Matrix<double_t, 6, 1, Eigen::ColMajor>;

void concatenateVectors(SiconosVector& target, const SiconosVector& a, const SiconosVector& b);
std::shared_ptr<SiconosVector> concatenateVectors(const SiconosVector& a,
                                                  const SiconosVector& b);

}  // namespace siconos::algebra

#endif