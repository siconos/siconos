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

#include "ReferenceClasses.hpp"

siconos::internal::devel_model::ClassA::ClassA(
    Eigen::Ref<siconos::algebra::SiconosVector> param1)
    : ndof_{param1.size()} {
  vector1_view_ = std::make_shared<siconos::algebra::MapVectorType>(param1.data(), ndof_);

  new (&vectorNameDirect_view_) siconos::algebra::MapVectorType(nullptr, 0);

  // Initialize var with vector1 content
  var_ = std::make_shared<siconos::algebra::SiconosVector>(*vector1_view_);

  (*var_) *= 2.;  // Test purpose. Just to be sure that memory is allocated (first touch policy
                  // and so on ...)
}

void siconos::internal::devel_model::ClassA::setConstantVector2(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for vector2_
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  vector2_internal_storage_ = nullptr;

  vector2_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), newValue.size());
  hasVector2_ = true;
  hasConstantVector2_ = true;
  computevector2_ = nullptr;
}

void siconos::internal::devel_model::ClassA::setComputeVector2Function(
    const siconos::modeling::func_prototypes::FunctionS_V& new_func) {
  // Ensure that memory is properly allocated for vector2_
  if (!vector2_internal_storage_) {
    vector2_internal_storage_ = std::make_unique<std::vector<double>>(ndof_);
  }
  vector2_view_ = std::make_shared<siconos::algebra::MapVectorType>(
      vector2_internal_storage_->data(), ndof_);

  hasVector2_ = true;
  hasConstantVector2_ = false;
  computevector2_ = new_func;
}

void siconos::internal::devel_model::ClassA::computeVector2(double time) {
  if (computevector2_)
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    computevector2_(time, *vector2_view_);
}

void siconos::internal::devel_model::ClassA::setConstantVectorNameSpan(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for vectorNameSpan_
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  vectorNameSpan_internal_storage_ = nullptr;

  vectorNameSpan_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), newValue.size());
  hasVectorNameSpan_ = true;
  hasConstantVectorNameSpan_ = true;
  computevectorNameSpan_ = nullptr;
}

void siconos::internal::devel_model::ClassA::setComputeVectorNameSpanFunction(
    const FunctionSpanT_V& new_func) {
  // Ensure that memory is properly allocated for vectorNameSpan_
  if (!vectorNameSpan_internal_storage_) {
    vectorNameSpan_internal_storage_ = std::make_unique<std::vector<double>>(ndof_);
  }
  vectorNameSpan_view_ = std::make_shared<siconos::algebra::MapVectorType>(
      vectorNameSpan_internal_storage_->data(), ndof_);

  hasVectorNameSpan_ = true;
  hasConstantVectorNameSpan_ = false;
  computevectorNameSpan_ = new_func;
}

void siconos::internal::devel_model::ClassA::computeVectorNameSpan(double time) {
  if (computevectorNameSpan_) {
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    // std::span<double> buff(vectorNameSpan_internal_storage_->data(), ndof_);
    // Or
    std::span<double> buff(vectorNameSpan_view_->data(), ndof_);
    computevectorNameSpan_(time, buff);
  }
}

void siconos::internal::devel_model::ClassA::setConstantVectorNameDirect(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for vectorNameDirect_
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  vectorNameDirect_internal_storage_ = nullptr;

  // new (&vectorNameDirect_view_)
  //     siconos::algebra::MapVectorType(newValue.data(), newValue.size());

  vectorNameDirect_view_.emplace(newValue.data(), newValue.size());
  hasVectorNameDirect_ = true;
  hasConstantVectorNameDirect_ = true;
  computevectorNameDirect_ = nullptr;
}

void siconos::internal::devel_model::ClassA::setComputeVectorNameDirectFunction(
    const siconos::modeling::func_prototypes::FunctionS_V& new_func) {
  // Ensure that memory is properly allocated for vectorNameDirect_
  if (!vectorNameDirect_internal_storage_) {
    vectorNameDirect_internal_storage_ = std::make_unique<std::vector<double>>(ndof_);
  }
  // new (&vectorNameDirect_view_)
  //     siconos::algebra::MapVectorType(vectorNameDirect_internal_storage_->data(), ndof_);
  vectorNameDirect_view_.emplace(vectorNameDirect_internal_storage_->data(), ndof_);

  hasVectorNameDirect_ = true;
  hasConstantVectorNameDirect_ = false;
  computevectorNameDirect_ = new_func;
}

void siconos::internal::devel_model::ClassA::computeVectorNameDirect(double time) {
  if (computevectorNameDirect_)
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    computevectorNameDirect_(time, vectorNameDirect_view_.value());
}

void siconos::internal::devel_model::ClassA::setConstantMatrix1(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for matrix1
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  matrix1_internal_storage_ = nullptr;

  matrix1_view_ = std::make_shared<siconos::algebra::MapType>(newValue.data(), ndof_, ndof_);
  hasMatrix1_ = true;
  hasConstantMatrix1_ = true;
  computematrix1_ = nullptr;
}

void siconos::internal::devel_model::ClassA::setComputeMatrix1Function(
    const FunctionVS_M& new_func) {
  // Ensure that memory is properly allocated for matrix1_
  if (!matrix1_internal_storage_) {
    matrix1_internal_storage_ = std::make_unique<std::vector<double>>(ndof_ * ndof_);
  }
  matrix1_view_ = std::make_shared<siconos::algebra::MapType>(
      matrix1_internal_storage_->data(), ndof_, ndof_);

  hasMatrix1_ = true;
  hasConstantMatrix1_ = false;
  computematrix1_ = new_func;
}

void siconos::internal::devel_model::ClassA::computeMatrix1(
    Eigen::Ref<siconos::algebra::SiconosVector> position, double time) {
  if (computematrix1_) {
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    // std::span<double> sppos(position.data(), ndof_);
    // std::span<double> spmat(matrix1_internal_storage_->data(),
    //                         matrix1_internal_storage_->size());
    computematrix1_(position, time, *matrix1_view_);
  }
}

// void siconos::internal::devel_model::ClassA::computeMatrix1(
//     double time, Eigen::Ref<siconos::algebra::SiconosVector> position) {
//   if (computematrix1_) {
//     // in that case, internal_storage must have been allocated by
//     // setCompute... call
//     std::span<double> sppos(position.data(), ndof_);
//     std::span<double> spmat(matrix1_internal_storage_->data(),
//                             matrix1_internal_storage_->size());
//     computematrix1_(time, sppos, spmat);
//   }
// }
void siconos::internal::devel_model::ClassA::display(bool brief) {
  std::cout << "----- Class A display -----\n";
  if (hasVector2_) {
    if (hasConstantVector2_)
      std::cout << "vector2 is ON and is constant\n";
    else
      std::cout << "vector2 is ON and is controlled by a plugin\n" << *vector2_view_ << "\n";
  }
  std::cout << *vector1_view_ << "\n";
  std::cout << *var_ << "\n";

  std::cout << "----- End of Class A display -----\n";
}
