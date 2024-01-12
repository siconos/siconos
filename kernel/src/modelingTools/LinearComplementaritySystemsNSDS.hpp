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

namespace siconos::algebra {

class SiconosVector;
class SimpleMatrix;
}  // namespace siconos::algebra

namespace siconos::modeling {

class FirstOrderLinearTIDS;
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
  std::shared_ptr<FirstOrderLinearTIDS> _ds{nullptr};
  /** a first order linear TI relation */
  std::shared_ptr<FirstOrderLinearTIR> _relation{nullptr};
  /** a complementarity condition */
  std::shared_ptr<ComplementarityConditionNSL> _nslaw{nullptr};
  /** an interaction*/
  std::shared_ptr<Interaction> _interaction{nullptr};

 public:
  /** constructor with t0 and T
   *
   *  \param t0 initial time
   *  \param T final time
   */
  LinearComplementaritySystemsNSDS(double t0, double T,
                                   std::shared_ptr<siconos::algebra::SiconosVector> x0,
                                   std::shared_ptr<siconos::algebra::SimpleMatrix> A,
                                   std::shared_ptr<siconos::algebra::SimpleMatrix> B,
                                   std::shared_ptr<siconos::algebra::SimpleMatrix> C,
                                   std::shared_ptr<siconos::algebra::SimpleMatrix> D,
                                   std::shared_ptr<siconos::algebra::SiconosVector> a,
                                   std::shared_ptr<siconos::algebra::SiconosVector> b);

  /** destructor
   */
  ~LinearComplementaritySystemsNSDS() noexcept = default;

  // --- GETTERS/SETTERS ---

  std::shared_ptr<FirstOrderLinearTIDS> ds() { return _ds; };

  std::shared_ptr<FirstOrderLinearTIR> relation() { return _relation; };

  std::shared_ptr<ComplementarityConditionNSL> nslaw() { return _nslaw; };
  std::shared_ptr<Interaction> interaction() { return _interaction; };
};
}  // namespace siconos::modeling

#endif
