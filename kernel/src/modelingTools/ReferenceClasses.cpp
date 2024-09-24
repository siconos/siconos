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
  vectorName1_view_ = std::make_shared<siconos::algebra::MapVectorType>(param1.data(), ndof_);

  new (&vectorNameDirect_view_) siconos::algebra::MapVectorType(nullptr, 0);

  // Initialize vectorName3 with vectorName1 content
  vectorName3_ = std::make_shared<siconos::algebra::SiconosVector>(*vectorName1_view_);

  (*vectorName3_) *=
      2.;  // just to be sure that memory is allocated (first touch policy and so on ...)
}



void siconos::internal::devel_model::ClassA::setConstantVectorName2(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for vectorName2_
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  vectorName2_internal_storage_ = nullptr;

  vectorName2_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), newValue.size());
  hasVectorName2_ = true;
  hasConstantVectorName2_ = true;
  computevectorName2_ = nullptr;
}

void siconos::internal::devel_model::ClassA::setComputeVectorName2Function(
    ExternalFunctionType new_func) {
  // Ensure that memory is properly allocated for vectorName2_
  if (!vectorName2_internal_storage_) {
    vectorName2_internal_storage_ = std::make_unique<std::vector<double>>(ndof_);
  }
  vectorName2_view_ = std::make_shared<siconos::algebra::MapVectorType>(
      vectorName2_internal_storage_->data(), ndof_);

  hasVectorName2_ = true;
  hasConstantVectorName2_ = false;
  computevectorName2_ = new_func;
}

void siconos::internal::devel_model::ClassA::computeVectorName2(double time) {
  if (computevectorName2_)
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    computevectorName2_(time, *vectorName2_view_);
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
    ExternalFunctionSpanType new_func) {
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
    ExternalFunctionType new_func) {
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

void siconos::internal::devel_model::ClassA::setConstantMatrixName(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for matrixName
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  matrixName_internal_storage_ = nullptr;

  matrixName_view_ =
      std::make_shared<siconos::algebra::MapType>(newValue.data(), ndof_, ndof_);
  hasMatrixName_ = true;
  hasConstantMatrixName_ = true;
  computematrixName_ = nullptr;
}

void siconos::internal::devel_model::ClassA::setComputeMatrixNameFunction(
    ExternalFunctionTypeMatrix new_func) {
  // Ensure that memory is properly allocated for matrixName_
  if (!matrixName_internal_storage_) {
    matrixName_internal_storage_ = std::make_unique<std::vector<double>>(ndof_ * ndof_);
  }
  matrixName_view_ = std::make_shared<siconos::algebra::MapType>(
      matrixName_internal_storage_->data(), ndof_, ndof_);

  hasMatrixName_ = true;
  hasConstantMatrixName_ = false;
  computematrixName_ = new_func;
}

void siconos::internal::devel_model::ClassA::computeMatrixName(
    double time, Eigen::Ref<siconos::algebra::SiconosVector> position) {
  if (computematrixName_) {
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    // std::span<double> sppos(position.data(), ndof_);
    // std::span<double> spmat(matrixName_internal_storage_->data(),
    //                         matrixName_internal_storage_->size());
    computematrixName_(position, time, *matrixName_view_);
  }
}

// void siconos::internal::devel_model::ClassA::computeMatrixName(
//     double time, Eigen::Ref<siconos::algebra::SiconosVector> position) {
//   if (computematrixName_) {
//     // in that case, internal_storage must have been allocated by
//     // setCompute... call
//     std::span<double> sppos(position.data(), ndof_);
//     std::span<double> spmat(matrixName_internal_storage_->data(),
//                             matrixName_internal_storage_->size());
//     computematrixName_(time, sppos, spmat);
//   }
// }
void siconos::internal::devel_model::ClassA::display(bool brief) {
  std::cout << "----- Class A display -----\n";
  if (hasVectorName2_) {
    if (hasConstantVectorName2_)
      std::cout << "vectorName2 is ON and is constant\n";
    else
      std::cout << "vectorName2 is ON and is controlled by a plugin\n";
    vectorName2_view_->display();
  }
  vectorName1_view_->display();
  vectorName3_->display();

  std::cout << "----- End of Class A display -----\n";
}
