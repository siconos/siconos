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
/*! \file NonSmoothDynamicalSystem.hpp
 * \brief container for DynamicalSystem and Interaction
 */
#ifndef LinearComplementaritySystemsNSDS_H
#define LinearComplementaritySystemsNSDS_H

#include "NonSmoothDynamicalSystem.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "StorageTools.hpp"

namespace siconos::algebra {}  // namespace siconos::algebra

namespace siconos::modeling {

class FirstOrderLinearDS;
class FirstOrderLinearTIR;
class ComplementarityConditionNSL;

/**
    The LinearComplementaritySystemsNSDS_H inherits frim NSDS
    for a direct instanciation of a LCS
*/
class LinearComplementaritySystemsNSDS : public NonSmoothDynamicalSystem {
 private:
  ACCEPT_SERIALIZATION(LinearComplementaritySystemsNSDS);

  /** a first order linear TI dynamical systems */
  std::shared_ptr<FirstOrderLinearDS> _ds{nullptr};
  /** a first order linear TI relation */
  std::shared_ptr<FirstOrderLinearTIR> _relation{nullptr};
  /** a complementarity condition */
  std::shared_ptr<ComplementarityConditionNSL> _nslaw{nullptr};
  /** an interaction*/
  std::shared_ptr<Interaction> _interaction{nullptr};

 public:
  /** @brief constructor with t0 and T
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means:
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   *
   *  @param t0 initial time
   *  @param T final time
   *  @param x0 initial state
   *  @param A ds matrix
   *  @param B relation matrix
   *  @param C relation matrix
   *  @param D relation matrix
   *  @param b ds vector
   *  @param e relation vector
   *  @param tag Pass siconos::algebra::alias_t to select this overload
   *       (rather than copy version)
   */
  LinearComplementaritySystemsNSDS(double t0, double T,
                                   Eigen::Ref<siconos::algebra::SiconosVector> x0,
                                   Eigen::Ref<siconos::algebra::SiconosMatrix> A,
                                   Eigen::Ref<siconos::algebra::SiconosMatrix> B,
                                   Eigen::Ref<siconos::algebra::SiconosMatrix> C,
                                   Eigen::Ref<siconos::algebra::SiconosMatrix> D,
                                   Eigen::Ref<siconos::algebra::SiconosVector> b,
                                   Eigen::Ref<siconos::algebra::SiconosVector> e,
                                   siconos::algebra::AliasTag tag);

  //   /** @brief constructor with t0 and T
  //    *
  //    * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
  //    * pointing directly to the memory provided by the argument.
  //    *
  //    * This means:
  //    *  - ownership stays external (!! only implemented for DS op.)
  //    *  - modifications to the original vector are reflected inside the class
  //    *
  //    *  @param t0 initial time
  //    *  @param T final time
  //    *  @param x0 initial state
  //    *  @param A ds matrix
  //    *  @param B relation matrix
  //    *  @param C relation matrix
  //    *  @param D relation matrix
  //    *  @param b ds vector
  //    *  @param e relation vector
  //    *  @param tag Pass siconos::algebra::alias_t to select this overload
  //    *       (rather than copy version)

  //    */
  //   LinearComplementaritySystemsNSDS(double t0, double T,
  //                                    const siconos::algebra::SiconosVector& x0,
  //                                    const siconos::algebra::SiconosDenseMatrix& A,
  //                                    const siconos::algebra::SiconosDenseMatrix& B,
  //                                    const siconos::algebra::SiconosDenseMatrix& C,
  //                                    const siconos::algebra::SiconosDenseMatrix& D,
  //                                    const siconos::algebra::SiconosVector& b,
  //                                    const siconos::algebra::SiconosVector& e,
  //                                    siconos::algebra::CopyTag tag);
  /** destructor
   */
  ~LinearComplementaritySystemsNSDS() noexcept = default;

  // --- GETTERS/SETTERS ---

  auto ds() { return _ds; };

  std::shared_ptr<FirstOrderLinearTIR> relation() { return _relation; };

  std::shared_ptr<ComplementarityConditionNSL> nslaw() { return _nslaw; };
  std::shared_ptr<Interaction> interaction() { return _interaction; };
};
}  // namespace siconos::modeling

#endif
