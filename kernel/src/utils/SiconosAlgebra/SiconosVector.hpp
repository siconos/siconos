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

void concatenateVectors(SiconosVector& target, const SiconosVector& a, const SiconosVector& b);
std::shared_ptr<SiconosVector> concatenateVectors(const SiconosVector& a,
                                                  const SiconosVector& b);

/* function to compute an orthonormal basis form a given vector
 * adapted with eigen vectors from
 * Tom Duff, James Burgess, Per Christensen, Christophe Hery, Andrew Kensler, Max Liani, and
 * Ryusuke Villemin, Building an Orthonormal Basis, Revisited, Journal of Computer Graphics
 * Techniques (JCGT), vol. 6, no. 1, 1-8, 2017 Available online
 * http://jcgt.org/published/0006/01/01/ void branchlessONB(const Vec3f &n, Vec3f &b1, Vec3f
 * &b2)
 * {
 *   float sign = copysignf(1.0f, n.z);
 *   const float a = -1.0f / (sign + n.z);
 *   const float b = n.x * n.y * a;
 *   b1 = Vec3f(1.0f + sign * n.x * n.x * a, sign * b, -sign * n.x);
 *   b2 = Vec3f(b, sign + n.y * n.y * a, -n.y);
 * }
 */
/**  compute an orthonormal basis
 *   \param[in,out] A a reference vector (normalize after call)
 *   \param[out] A1 second base vector
 *   \param[out] A2 third base vector
 *   \return true if all went right else false (e.g. if A.norm =0)
 */
bool orthoBaseFromVector(SiconosVector3& A, SiconosVector3& A1, SiconosVector3& A2);

}  // namespace siconos::algebra

#endif