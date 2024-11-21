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
#include "LinearComplementaritySystemsNSDS.hpp"

#include "ComplementarityConditionNSL.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "Interaction.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include "siconos_debug.h"

//  constructor
siconos::modeling::LinearComplementaritySystemsNSDS::LinearComplementaritySystemsNSDS(
    double t0, double T, Eigen::Ref<siconos::algebra::SiconosVector> x0,
    Eigen::Ref<siconos::algebra::SiconosMatrix> A,
    Eigen::Ref<siconos::algebra::SiconosMatrix> B,
    Eigen::Ref<siconos::algebra::SiconosMatrix> C,
    Eigen::Ref<siconos::algebra::SiconosMatrix> D,
    Eigen::Ref<siconos::algebra::SiconosVector> b,
    Eigen::Ref<siconos::algebra::SiconosVector> e)
    : NonSmoothDynamicalSystem(t0, T) {
  _ds = std::make_shared<FirstOrderLinearTIDS>(x0, A);
  _ds->setbPtr(b);

  _relation = std::make_shared<FirstOrderLinearTIR>(C, B);

  // todo: check sizes --> done during rel->initialize()
  _relation->setConstantD(D);
  _relation->setConstanteVector(e);

  _nslaw = std::make_shared<ComplementarityConditionNSL>(C.rows());
  _interaction = std::make_shared<Interaction>(_nslaw, _relation);

  link(_interaction, _ds);

  display();
};
