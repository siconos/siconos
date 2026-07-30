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

namespace siconos::algebra {
class BlockVector;
}

namespace siconos::modeling {

namespace func_prototypes {

//  All function prototypes used to declare user-defined functions to be plugged
//  to compute dynamical systems or relations operators.
//
//  Naming convention: FunctionXYZ..._R or RFunctionXYZ with:
//     - XYZ... the type of read-only in parameters, with S for scalar type, V for vector type,
//     M for dense matrix type and Ms for sparse matrix type
//     - R(=S, V, M or Ms) the read-write in-out parameter used to collect the result
//     - RFunction if the return type is non void (and corresponds to the read-write in-out
//     parameter)
//
// Warning FP:
//  - input parameters of vector of matrix type must be const Eigen::Ref<...> &
//  - output param must be Eigen::Ref<MapVectorType> or Eigen::Ref<MapType>
//     or SiconosSparseMatrix& for the sparse case.
//  - RFunction are required for SparseMatrices to be handled properly in the pybind11 to scipy
//  interface
//
// This is the best (and maybe only) compromise to be able to
// - handle map, ref, vector ... as input parameters
// - ensure read-only properties on input parameters
// - deal with the pybind11 bindings and be able to get/set properly numpy arrays
//  without non expected copies or temps values.

/** function proto to compute f(vector,vector,t, result_vector) */
using FunctionS_S = std::function<double(double)>;

/** function proto to compute double = f(double) */
using FunctionVVS_V =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       const Eigen::Ref<const siconos::algebra::SiconosVector>&, double,
                       Eigen::Ref<siconos::algebra::MapVectorType>)>;

/** function proto to compute f(vector,t, result_vector) */
using FunctionVS_V =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector>&, double,
                       Eigen::Ref<siconos::algebra::MapVectorType>)>;

/** function proto to compute f(vector, result_vector) */
using FunctionV_V =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       Eigen::Ref<siconos::algebra::MapType>)>;

/** function proto to compute f(t, result_vector) */
using FunctionS_V = std::function<void(double, Eigen::Ref<siconos::algebra::MapVectorType>)>;

/** function proto to compute f(vector,vector,t, result_vector) */
using FunctionVV_V =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       Eigen::Ref<siconos::algebra::MapVectorType>)>;

/** function proto to compute f(block vector,t, result_vector) */
using FunctionBVS_V = std::function<void(const siconos::algebra::BlockVector&, double,
                                         Eigen::Ref<siconos::algebra::MapVectorType>)>;

/** function proto to compute f(block vector, result_vector) */
using FunctionBV_V = std::function<void(const siconos::algebra::BlockVector&,
                                        Eigen::Ref<siconos::algebra::MapVectorType>)>;

/** function proto to compute f(block vector, vector, result_vector) */
using FunctionBVV_V =
    std::function<void(const siconos::algebra::BlockVector&,
                       const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       Eigen::Ref<siconos::algebra::MapVectorType>)>;

/** function proto to compute f(block vector,t, vector, result_vector) */
using FunctionBVSV_V =
    std::function<void(const siconos::algebra::BlockVector&, double,
                       const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       Eigen::Ref<siconos::algebra::MapVectorType>)>;

/** function proto to compute f(block vector,t, vector, result_blockvector) */
using FunctionBVSV_BV = std::function<void(
    const siconos::algebra::BlockVector&, double,
    const Eigen::Ref<const siconos::algebra::SiconosVector>&, siconos::algebra::BlockVector&)>;

/** function proto to compute f(vector, result_blockvector) */
using FunctionV_BV = std::function<void(
    const Eigen::Ref<const siconos::algebra::SiconosVector>&, siconos::algebra::BlockVector&)>;

/** function proto to compute f(vector,t, result_matrix) */
using FunctionVS_M =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector>&, double,
                       Eigen::Ref<siconos::algebra::MapType>)>;

/** function proto to compute f(vector,vector,t, result_matrix) */
using FunctionVVS_M =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       const Eigen::Ref<const siconos::algebra::SiconosVector>&, double,
                       Eigen::Ref<siconos::algebra::MapType>)>;

/** function proto to compute f(vector,vector,t, result_matrix (sparse)) */
using FunctionVVS_Ms =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       const Eigen::Ref<const siconos::algebra::SiconosVector>&, double,
                       siconos::algebra::SiconosSparseMatrix&)>;

/** function proto to compute result_matrix = f(vector,vector,t, result_matrix (sparse)) */
using RFunctionVVS_Ms = std::function<siconos::algebra::SiconosSparseMatrix(
    const Eigen::Ref<const siconos::algebra::SiconosVector>&,
    const Eigen::Ref<const siconos::algebra::SiconosVector>&, double,
    siconos::algebra::SiconosSparseMatrix&)>;

/** function proto to compute f(vector,vector, result_matrix) */
using FunctionVV_M =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       Eigen::Ref<siconos::algebra::MapType>)>;

/*** function proto to compute f(vector,vector, result_matrix (sparse)) */
using FunctionVV_Ms =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       siconos::algebra::SiconosSparseMatrix&)>;
/*** function proto to compute result_matrix = f(vector,vector, result_matrix (sparse)) */
using RFunctionVV_Ms = std::function<siconos::algebra::SiconosSparseMatrix(
    const Eigen::Ref<const siconos::algebra::SiconosVector>&,
    const Eigen::Ref<const siconos::algebra::SiconosVector>&,
    siconos::algebra::SiconosSparseMatrix&)>;

/** function proto to compute f(matrix,vector, result_matrix) */
using FunctionMV_V =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosMatrix>&,
                       const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       Eigen::Ref<siconos::algebra::MapType>)>;

/** function proto to compute f(t, result_matrix) */
using FunctionS_M = std::function<void(double, Eigen::Ref<siconos::algebra::MapType>)>;

/** function proto to compute f(vector, result_matrix) */
using FunctionV_M =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       Eigen::Ref<siconos::algebra::MapType>)>;

/** function proto to compute f(vector, result_matrix (sparse)) */
using FunctionV_Ms =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       siconos::algebra::SiconosSparseMatrix&)>;

/** function proto to compute result_matrix = f(vector, result_matrix (sparse)) */
using RFunction_V_Ms = std::function<siconos::algebra::SiconosSparseMatrix(
    const Eigen::Ref<const siconos::algebra::SiconosVector>&,
    siconos::algebra::SiconosSparseMatrix&)>;

/** function proto to compute f(block vector, t, result_matrix) */
using FunctionBVS_M = std::function<void(const siconos::algebra::BlockVector&, double,
                                         Eigen::Ref<siconos::algebra::MapType>)>;

/** function proto to compute f(block vector, t, result_matrix_sparse) */
using FunctionBVS_Ms = std::function<void(const siconos::algebra::BlockVector&, double,
                                          siconos::algebra::SiconosSparseMatrix&)>;

/** function proto to compute result_matrix = f(block vector, t, result_matrix_sparse)*/
using RFunctionBVS_Ms = std::function<siconos::algebra::SiconosSparseMatrix(
    const siconos::algebra::BlockVector&, double, siconos::algebra::SiconosSparseMatrix&)>;

/** function proto to compute f(block vector, result_matrix) */
using FunctionBV_M = std::function<void(const siconos::algebra::BlockVector&,
                                        Eigen::Ref<siconos::algebra::MapType>)>;

/** function proto to compute f(block vector, vector, result_matrix) */
using FunctionBVV_M =
    std::function<void(const siconos::algebra::BlockVector&,
                       const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       Eigen::Ref<siconos::algebra::MapType>)>;

/** function proto to compute f(block vector, block vector, result_matrix) */
using FunctionBVBV_M = std::function<void(const siconos::algebra::BlockVector&,
                                          const siconos::algebra::BlockVector&,
                                          Eigen::Ref<siconos::algebra::MapType>)>;

/** function proto to compute f(block vector, t, vector, result_matrix) */
using FunctionBVSV_M =
    std::function<void(const siconos::algebra::BlockVector&, double,
                       const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       Eigen::Ref<siconos::algebra::MapType>)>;

}  // namespace func_prototypes
}  // namespace siconos::modeling

#endif