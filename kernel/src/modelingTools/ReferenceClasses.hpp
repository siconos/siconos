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

#include "FunctionTypes.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

namespace siconos::internal::devel_model {

// A list of std::function types
// Types that may be used to 'plug' an external function (lambda ...)
// to a local method, used to compute attributes at runtime.

// Naming convention:
// FunctionXYZ_M: a function(x,y,z) which returns a matrix (V for vector)
// e.g for wrench(twist, q, t) in NewtonEulerDS : FunctionVQT_V
// Use:
//   - const Eigen::Ref<Map> & for [in] params
// //   - Eigen::Ref<Map> for [in,out]
// using siconos::modeling::func_prototypes::FunctionT_V =
//     std::function<void(double, Eigen::Ref<siconos::algebra::MapVectorType>)>;
// using FunctionVQT_V = std::function<void(const Eigen::Ref<siconos::algebra::MapVectorType>
// &,
//                                          const Eigen::Ref<siconos::algebra::MapVectorType>
//                                          &, double,
//                                          Eigen::Ref<siconos::algebra::MapVectorType>)>;

// // Extras, test purpose
using FunctionSpanT_V = std::function<void(double, std::span<double>)>;
using FunctionVS_M = std::function<void(Eigen::Ref<siconos::algebra::MapVectorType>, double,
                                        Eigen::Ref<siconos::algebra::MapType>)>;
// using FunctionVT_M =
//     std::function<void(double, std::span<double>, std::span<double>)>;

/**

 Description of the class (aim, ...)

 All attributes and methods names MUST be self-sufficient to understand
 the purpose of the var or functions

Questions, remarks:

- move all private to public (to prepare 'data-oriented' design) ?

- convention: name all Map type attributes with xxx_view_ ?

*/
class ClassA {
 protected:
  /** something which can be used to compute the size of all vectors, matrices attribute
   *  like number of degrees of freedom in a DS
   */
  siconos::algebra::SiconosVector::Index ndof_{0};

  // Attribute representing a variable of the problem

  //   - must/can be read-write accessible
  //   - memory handled by the class (for the moment. Later, external "pool"?)
  //   example: q for SecondOrderDS, twist for NewtonEuler ...
  //
  //  -> shared_ptr to SiconosVector
  //  -> a getter/setter returning a shared_ptr (--> write access)
  //     a Map would probably be enough, but because of DSlink
  //     in OSI, we keep shared_ptr for the time being
  //  -> a read-only getter (ConstMap)
  std::shared_ptr<siconos::algebra::SiconosVector> var_{nullptr};

  // Attribute representing a 'parameter' or operator of the system
  //  - memory might be externally handled (shared)
  //  - can not be modified after construction
  //  example: q0, v0 for LagrangianDS
  //
  //  -> attributes: a shared_ptr to a Map AND an internal_storage (unique_ptr)
  //  -> must be set by constructor
  //  -> has a read-only getter (ConstMap)

  /** a vector type attribute */
  std::shared_ptr<siconos::algebra::MapVectorType> vector1_view_{nullptr};

  std::unique_ptr<std::vector<double>> vector1_internal_storage_{nullptr};

  // Attribute representing a 'parameter' or operator of the system
  //  - memory might be externally handled (shared)
  //  - can be modified after construction: either with a constant input (matrix or vector)
  //    or with a used-defined function
  //  - can be accessed (read-only)
  //  example: mass, fext for LagrangianDS
  //
  //  -> attributes:
  //    - a shared_ptr to a Map and a unique_ptr to std::vector as internal_storage.
  //      internal_storage allocated only if a pluggin is used
  //    - bool hasConstantXXX_ and hasXXX_ + hasXXX() method
  //    - a std::function to be plugged
  //  -> not set by constructor
  //  -> a setConstantXX() function and a setComputeXXFunction
  //  -> a read-only getter (ConstMap)
  //  -> a computeXX() function
  /** A vector type attribute */
  std::shared_ptr<siconos::algebra::MapVectorType> vector2_view_{nullptr};
  std::unique_ptr<std::vector<double>> vector2_internal_storage_{nullptr};
  bool hasConstantVector2_{false};
  bool hasVector2_{false};
  siconos::modeling::func_prototypes::FunctionS_V computevector2_{nullptr};

  // Attribute representing a 'parameter' or operator of the system
  //  - internal use only (can not be/need not to be accessed)
  //  - no memory
  //  - can be computed with a user defined function
  //  example: fint for NewtonEulerDS
  //
  //  -> attributes:
  //    - bool hasXXX_ + hasXXX() method
  //    - a std::function to be plugged
  //  -> a setComputeXXFunction
  /** true if vector3 is taken into account in the model */
  bool hasVector3_{false};

  /** A function to compute vector3 */
  siconos::modeling::func_prototypes::FunctionS_V computeVector3_{nullptr};

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
  FunctionSpanT_V computevectorNameSpan_{nullptr};
  /** True if vector2 is required and constant */
  bool hasConstantVectorNameSpan_{false};

  /** True if vector2 is required */
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
  siconos::modeling::func_prototypes::FunctionS_V computevectorNameDirect_{nullptr};

  /** True if vector2 is required and constant */
  bool hasConstantVectorNameDirect_{false};

  /** True if vector2 is required */
  bool hasVectorNameDirect_{false};

  /** A matrix type attribute
   *  - optional
   *  - may be set by external user as a constant
   *  - may be set by an external function
   *  Example : mass in Lagrangian system
   */
  /** mass of the system */
  std::shared_ptr<siconos::algebra::MapType> matrix1_view_{nullptr};

  /** mass internal storage */
  std::unique_ptr<std::vector<double>> matrix1_internal_storage_{nullptr};

  /** function wrapper used to compute matrix1
   *
   *  Should we set it to lambda function returning nothing ? To nullptr ?
   */
  FunctionVS_M computematrix1_{nullptr};

  /** true if matrix1 is constant */
  bool hasConstantMatrix1_{false};

  /** True if matrix1 is required */
  bool hasMatrix1_{false};

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
  // Other attributes (e.g. vector2) should be set with setFunctions
  // Ok ?

  /** constructor ...
   *
   *  \param param1 decription used to set vector1
   */
  ClassA(Eigen::Ref<siconos::algebra::SiconosVector> param1);

  /** destructor */
  virtual ~ClassA() noexcept = default;

  // ------ var -------

  /** \return a read-only view on the generalized coordinates vector (size=dimension()) */
  inline const siconos::algebra::ConstMapVectorType var_read() const {
    return siconos::algebra::ConstMapVectorType(var_->data(), var_->size());
  }

  /** \return the generalized coordinates of the system (pointer link) */
  inline std::shared_ptr<siconos::algebra::SiconosVector> var() const { return var_; }

  // ----- vector1 ------
  /** \return a read-only view on ... */
  inline const siconos::algebra::ConstMapVectorType vector1() const {
    return siconos::algebra::ConstMapVectorType(vector1_view_->data(), vector1_view_->size());
  }

  // ----- vector2 ------
  /** \return  a read-only view on ... */
  inline auto vector2() const {
    return siconos::algebra::ConstMapVectorType(vector2_view_->data(), vector2_view_->size());
  }

  /** set a constant...
   *
   *  \param newVector2 ...
   */
  void setConstantVector2(Eigen::Ref<siconos::algebra::SiconosVector> newVector2);

  /** True if ... taken into account */
  bool hasVector2() const { return hasVector2_; }

  /** set a user-defined function to ...
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeVector2Function(const siconos::modeling::func_prototypes::FunctionS_V &fct);

  /** Update ...
   *
   *  \param time ...
   */
  void computeVector2(double time);

  // ----- vector3 ------

  /** set a user-defined function to compute ...
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeVector3Function(const siconos::modeling::func_prototypes::FunctionVVS_V &fct);

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
  void setComputeVectorNameSpanFunction(const FunctionSpanT_V &fext_func);

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
  void setComputeVectorNameDirectFunction(
      const siconos::modeling::func_prototypes::FunctionS_V &fext_func);

  /** default function to compute ...
   *
   *  \param time the current time
   */
  void computeVectorNameDirect(double time);

  // Getters / Setters for attributes like matrix1
  // Setters: we keep only the Map signature (not the set(std::shared...)

  /*  \return a read-only view on ... */
  inline auto matrix1() const {
    return siconos::algebra::ConstMapType(matrix1_view_->data(), matrix1_view_->rows(),
                                          matrix1_view_->cols());
  }

  /** describe ... e.g. set a constant external forces vector for the system.
   *
   *  \param newValue external forces vector
   *
   *  Should we warn there that the memory is handled by newValue and that any change on
   * newValue impact matrix1?
   */
  void setConstantMatrix1(Eigen::Ref<siconos::algebra::SiconosMatrix> newValue);

  /** describe ... e.g. \return True if external forces are taken into account in the system */
  bool hasMatrix1() const { return hasMatrix1_; }

  /** set a user-defined function to compute ...
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeMatrix1Function(const FunctionVS_M &fext_func);

  /** default function to compute ...
   *
   *  \param time the current time
   */
  void computeMatrix1(Eigen::Ref<siconos::algebra::SiconosVector> position, double time);

  /** print the data of the dynamical system on the standard output
   */
  void display(bool brief = true);
};

}  // namespace siconos::internal::devel_model
#endif  // REFCLASS_H
