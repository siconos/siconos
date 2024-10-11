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

/*! \file FunctionTypes.hpp
  alias for functions prototypes used in dynamical systems to provide
  user-defined functions
*/

#ifndef FUNCTION_TYPES_H
#define FUNCTION_TYPES_H

#include <functional>

#include "SiconosMatrix.hpp"

namespace siconos::modeling {

namespace func_prototypes {

//  All function prototypes used to declare user-defined functions to be plugged
//  to compute dynamical systems or relations operators.
//
//  Naming convention: FunctionXYZ..._R, with:
//     - XYZ... the type of read-only in parameters, with S for scalar type, V for vector type
//     and M for matrix type
//     - R(=S, V or M) the read-write in-out parameter used to collect the result
//
// Warning FP:
//  - input parameters of vector of matrix type must be const Eigen::Ref<...> &
//  - output param must be Eigen::Ref<MapVectorType> or Eigen::Ref<MapType>
//
// This is the best (and maybe only) compromise to be able to
// - handle map, ref, vector ... as input parameters
// - ensure read-only properties on input parameters
// - deal with the pybind11 bindings and be able to get/set properly numpy arrays
//  without non expected copies or temps values.

/** fonction proto to compute f(vector,vector,t, result_vector) */
using FunctionVVS_V =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector> &,
                       const Eigen::Ref<const siconos::algebra::SiconosVector> &, double,
                       Eigen::Ref<siconos::algebra::MapVectorType>)>;

/** fonction proto to compute f(vector, result_vector) */
using FunctionV_V =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector> &,
                       Eigen::Ref<siconos::algebra::MapType>)>;

/** fonction proto to compute f(t, result_vector) */
using FunctionS_V = std::function<void(double, Eigen::Ref<siconos::algebra::MapVectorType>)>;

/** fonction proto to compute f(vector,vector,t, result_vector) */
using FunctionVV_V =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector> &,
                       const Eigen::Ref<const siconos::algebra::SiconosVector> &,
                       Eigen::Ref<siconos::algebra::MapVectorType>)>;

/** fonction proto to compute f(vector,vector,t, result_matrix) */
using FunctionVVS_M =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector> &,
                       const Eigen::Ref<const siconos::algebra::SiconosVector> &, double,
                       Eigen::Ref<siconos::algebra::MapType>)>;

/** fonction proto to compute f(vector,vector, result_matrix) */
using FunctionVV_M =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector> &,
                       const Eigen::Ref<const siconos::algebra::SiconosVector> &,
                       Eigen::Ref<siconos::algebra::MapType>)>;

/** fonction proto to compute f(matrix,vector, result_matrix) */
using FunctionMV_V =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosMatrix> &,
                       const Eigen::Ref<const siconos::algebra::SiconosVector> &,
                       Eigen::Ref<siconos::algebra::MapType>)>;
}  // namespace func_prototypes
}  // namespace siconos::modeling

#endif