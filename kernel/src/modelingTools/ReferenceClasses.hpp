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

/*! \file ReferenceClasses.hpp
  ReferenceClasses - Examples of how a Siconos class must be written
*/

#ifndef REFCLASS_H
#define REFCLASS_H

#include <functional>
#include <memory>
#include <optional>
#include <span>

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

namespace siconos::internal::devel_model {
/**

 Description of the class (aim, ...)

 All attributes and methods names MUST be self-sufficient to understand
 the purpose of the var or functions

Questions, remarks:

- move all private to public (to prepare 'data-oriented' design) ?

- convention: name all Map type attributes with xxx_view_ ?

*/
class ClassA {
 public:
  // A list of std::function types
  // Types that may be used to 'plug' an external function (lambda ...)
  // to a local method, used to compute attributes at runtime.
  using ExternalFunctionType =
      std::function<void(double, Eigen::Ref<siconos::algebra::MapVectorType>)>;
  using ExternalFunctionSpanType = std::function<void(double, std::span<double>)>;
  using ExternalFunctionTypeMatrix =
      std::function<void(Eigen::Ref<siconos::algebra::MapVectorType>, double,
                         Eigen::Ref<siconos::algebra::MapType>)>;
  // using ExternalFunctionTypeMatrix =
  //     std::function<void(double, std::span<double>, std::span<double>)>;

 protected:
  /** something which can be used to compute the size of all vectors, matrices attribute
   *  like number of degrees of freedom in a DS
   */
  siconos::algebra::SiconosVector::Index ndof_{0};

  /** A vector type attribute, that can be set by user
   *  - always required --> must be provided by all constructors
   *  - memory always allocated
   *  - can be accessed (--> read with a getter)
   *
   *  Example : q0 in LagrangianDS
   *
   *  better to copy? Keep map?
   */
  std::shared_ptr<siconos::algebra::MapVectorType> vectorName1_view_{nullptr};

  std::unique_ptr<std::vector<double>> vectorName1_view_internal_storage{nullptr};

  /** A vector type attribute, internal, can not be set by user
   *  - always required --> must be properly set by all constructors.
   *  - memory always allocated
   *  - can be accessed (--> read with a getter)
   *
   *  Example : q in LagrangianDS
   */
  std::shared_ptr<siconos::algebra::SiconosVector> vectorName3_{nullptr};

  /** A vector type attribute
   *  - optional -> need a "setXXX" function
   *  - no memory allocated by default
   *  - may be set by external user (set ...) or with a plugin
   *
   *  Example : fext in LagrangianDS
   *
   *  Just a Map (view) to some storage :
   *    - internal if required but not provided by user
   *       -> see setComputeVectorName2Function
   *       -> need internal_storage
   *
   *    - 'external' if required and provided by user
   *       -> see setVectorName2
   *       -> no internal storage, memory handled outside of this class
   *
   * Is it really necessary to keep a shared pointer ?
   *
   */
  std::shared_ptr<siconos::algebra::MapVectorType> vectorName2_view_{nullptr};

  /**  Internal storage for vectorName2
   *   Not nullptr only if setComputeVectorName2Function has been called
   *
   */
  std::unique_ptr<std::vector<double>> vectorName2_internal_storage_{nullptr};

  /** function wrapper used to compute vectorName2
   *
   *  Should we set it to lambda function returning nothing ? To nullptr ?
   */
  ExternalFunctionType computevectorName2_{nullptr};

  /** True if vectorName2 is required and constant */
  bool hasConstantVectorName2_{false};

  /** True if vectorName2 is required */
  bool hasVectorName2_{false};

  /* The same but with a span in the computeXX function.

    I keep it for test purpose but it seems to be better to use Eigen::Ref on our case?
  */
  std::shared_ptr<siconos::algebra::MapVectorType> vectorNameSpan_view_{nullptr};

  /**  Internal storage for vectorNameSpan
   *   Not nullptr only if setComputeVectorNameSpanFunction has been called
   *
   */
  std::unique_ptr<std::vector<double>> vectorNameSpan_internal_storage_{nullptr};

  /** function wrapper used to compute vectorNameSpan
   *
   *  Should we set it to lambda function returning nothing ? To nullptr ?
   */
  ExternalFunctionSpanType computevectorNameSpan_{nullptr};
  /** True if vectorName2 is required and constant */
  bool hasConstantVectorNameSpan_{false};

  /** True if vectorName2 is required */
  bool hasVectorNameSpan_{false};

  /* The same but without a shared pointer

   Test ...

 Can't be directly initialized at decl. It's possible to use
   std::optional<siconos::algebra::MapVectorType vectorNameSpan_view_> to avoid this
   but is it really necessary?
 */
  std::optional<siconos::algebra::MapVectorType> vectorNameDirect_view_;

  /**  Internal storage for vectorNameDirect
   *   Not nullptr only if setComputeVectorNameDirectFunction has been called
   *
   */
  std::unique_ptr<std::vector<double>> vectorNameDirect_internal_storage_{nullptr};

  /** function wrapper used to compute vectorNameDirect
   *
   *  Should we set it to lambda function returning nothing ? To nullptr ?
   */
  ExternalFunctionType computevectorNameDirect_{nullptr};

  /** True if vectorName2 is required and constant */
  bool hasConstantVectorNameDirect_{false};

  /** True if vectorName2 is required */
  bool hasVectorNameDirect_{false};

  /** A matrix type attribute
   *  - optional
   *  - may be set by external user as a constant
   *  - may be set by an external function
   *  Example : mass in Lagrangian system
   */
  /** mass of the system */
  std::shared_ptr<siconos::algebra::MapType> matrixName_view_{nullptr};

  /** mass internal storage */
  std::unique_ptr<std::vector<double>> matrixName_internal_storage_{nullptr};

  /** function wrapper used to compute matrixName
   *
   *  Should we set it to lambda function returning nothing ? To nullptr ?
   */
  ExternalFunctionTypeMatrix computematrixName_{nullptr};

  /** true if matrixName is constant */
  bool hasConstantMatrixName_{false};

  /** True if matrixName is required */
  bool hasMatrixName_{false};

  // Rule of five
  // Adapt it to each class case

  /** Default constructor */
  ClassA() = delete;
  ClassA(const ClassA &) = delete;
  ClassA(ClassA &&) = delete;
  ClassA &operator=(const ClassA &) = delete;
  ClassA &operator=(ClassA &&) = delete;

 public:
  // Rule 1: all required attributes must be properly set after constructor call
  // Rule 2: all attributes must have a defaul value (set with attr declaration)
  // Rule 3: constructors only with parameters for the required attributes.
  // Other attributes (e.g. vectorName2) should be set with setFunctions
  // Ok ?

  /** constructor with all required parameters
   *
   *  \param param1 decription used to set vectorName1
   */
  ClassA(Eigen::Ref<siconos::algebra::SiconosVector> param1);

  // ClassA(const Eigen::Map<siconos::algebra::SiconosVector>& param1);


  /** destructor */
  virtual ~ClassA() noexcept = default;

  /** \return describe what is returned (e.g. generalized coordinates, do not say 'returns q')
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> vectorName3() const {
    return vectorName3_;
  }

  /** \return describe what is returned (e.g. generalized coordinates, do not say 'returns a
   * pointer q')
   */
  inline siconos::algebra::SiconosVector &vectorName3_python() const {
    return *(vectorName3_);
  }

  // We do not provide setters for vectorName3-like variables. Should we?

  // Getters for attributes like vectorName1_

  /*    \return describe what is returned
   */
  inline std::shared_ptr<siconos::algebra::MapVectorType> vectorName1() const {
    return vectorName1_view_;
  }

  // Getters / Setters for attributes like vectorName2

  // Getters :
  // - standard (shared ptr) as it was formerly done in C++
  // - Map, as it must be done to be handled properly by pybind11
  // Should we keep both of them ? Don't know for the moment ...
  // Name convention: _view (or something else? _Map?) for the Map case

  // Setters: we keep only the Map signature (not the set(std::shared...)

  /** \return describe ...
   */
  inline std::shared_ptr<siconos::algebra::MapVectorType> vectorName2() const {
    return vectorName2_view_;
  }

  /** \return describe ...
   */
  inline siconos::algebra::MapVectorType &vectorName2_view() const {
    return *vectorName2_view_;
  }

  /** describe ... e.g. set a constant external forces vector for the system.
   *
   *  \param newValue external forces vector
   *
   *  Should we warn there that the memory is handled by newValue and that any change on
   * newValue impact vectorName2?
   */
  void setConstantVectorName2(Eigen::Ref<siconos::algebra::SiconosVector> newValue);

  /** describe ... e.g. \return True if external forces are taken into account in the system */
  bool hasVectorName2() const { return hasVectorName2_; }

  /** set a user-defined function to compute ...
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeVectorName2Function(ExternalFunctionType fext_func);

  /** default function to compute ...
   *
   *  \param time the current time
   */
  void computeVectorName2(double time);

  inline std::shared_ptr<siconos::algebra::MapVectorType> vectorNameSpan() const {
    return vectorNameSpan_view_;
  }

  /** \return describe ...
   */
  inline siconos::algebra::MapVectorType &vectorNameSpan_view() const {
    return *vectorNameSpan_view_;
  }

  /** describe ... e.g. set a constant external forces vector for the system.
   *
   *  \param newValue external forces vector
   *
   *  Should we warn there that the memory is handled by newValue and that any change on
   * newValue impact vectorNameSpan?
   */
  void setConstantVectorNameSpan(Eigen::Ref<siconos::algebra::SiconosVector> newValue);

  /** describe ... e.g. \return True if external forces are taken into account in the system */
  bool hasVectorNameSpan() const { return hasVectorNameSpan_; }

  /** set a user-defined function to compute ...
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeVectorNameSpanFunction(ExternalFunctionSpanType fext_func);

  /** default function to compute ...
   *
   *  \param time the current time
   */
  void computeVectorNameSpan(double time);

  /** \return describe ...
   */
  inline siconos::algebra::MapVectorType &vectorNameDirect_view() {
    return vectorNameDirect_view_.value();
  }

  /** describe ... e.g. set a constant external forces vector for the system.
   *
   *  \param newValue external forces vector
   *
   *  Should we warn there that the memory is handled by newValue and that any change on
   * newValue impact vectorNameDirect?
   */
  void setConstantVectorNameDirect(Eigen::Ref<siconos::algebra::SiconosVector> newValue);

  /** describe ... e.g. \return True if external forces are taken into account in the system */
  bool hasVectorNameDirect() const { return hasVectorNameDirect_; }

  /** set a user-defined function to compute ...
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeVectorNameDirectFunction(ExternalFunctionType fext_func);

  /** default function to compute ...
   *
   *  \param time the current time
   */
  void computeVectorNameDirect(double time);

  // Getters / Setters for attributes like matrixName
  // Setters: we keep only the Map signature (not the set(std::shared...)

  /*  \return the ... operator ... describe ... */
  inline std::shared_ptr<siconos::algebra::MapType> matrixName() const {
    return matrixName_view_;
  }

  /** \return describe ... */
  inline siconos::algebra::MapType &matrixName_view() const { return *matrixName_view_; }

  /** describe ... e.g. set a constant external forces vector for the system.
   *
   *  \param newValue external forces vector
   *
   *  Should we warn there that the memory is handled by newValue and that any change on
   * newValue impact matrixName?
   */
  void setConstantMatrixName(Eigen::Ref<siconos::algebra::SiconosMatrix> newValue);

  /** describe ... e.g. \return True if external forces are taken into account in the system */
  bool hasMatrixName() const { return hasMatrixName_; }

  /** set a user-defined function to compute ...
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeMatrixNameFunction(ExternalFunctionTypeMatrix fext_func);

  /** default function to compute ...
   *
   *  \param time the current time
   */
  void computeMatrixName(double time, Eigen::Ref<siconos::algebra::SiconosVector> position);

  /** print the data of the dynamical system on the standard output
   */
  void display(bool brief = true);
};

}  // namespace siconos::internal::devel_model
#endif  // REFCLASS_H
