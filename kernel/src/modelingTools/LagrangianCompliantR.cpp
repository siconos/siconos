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

// \todo : create a work vector for all tmp vectors used in computeg, computeh ...

#include "LagrangianCompliantR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosAlgebraAddons.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
#include "Tools.hpp"

void siconos::modeling::LagrangianCompliantR::initialize(Interaction& inter) {
  auto sizeY = inter.dimension();
  auto sizeDS = inter.getSizeOfDS();

  if (!jacobianhOver_q_internal_storage_) {
    jacobianhOver_q_internal_storage_ = std::make_unique<std::vector<double>>(sizeY * sizeDS);
  }
  jacobianhOver_q_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobianhOver_q_internal_storage_->data(), sizeY, sizeDS);

  if (!jacobianhOver_lambda_)
    jacobianhOver_lambda_ = std::make_shared<siconos::algebra::SiconosMatrix>(sizeY, sizeY);
}

void siconos::modeling::LagrangianCompliantR::setComputehFunction(
    const siconos::modeling::func_prototypes::FunctionBVV_V& h_func) {
  computeh_ = h_func;
}

void siconos::modeling::LagrangianCompliantR::computeh(
    const siconos::algebra::BlockVector& positions,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  if (computeh_) computeh_(positions, lambda, y);
}

void siconos::modeling::LagrangianCompliantR::setComputeJacobianhOver_qFunction(
    const siconos::modeling::func_prototypes::FunctionBVV_M& func) {
  computejacobianhOver_q_ = func;
}

void siconos::modeling::LagrangianCompliantR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector& positions,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda) {
  if (computejacobianhOver_q_)
    computejacobianhOver_q_(positions, lambda, *jacobianhOver_q_view_);
}

void siconos::modeling::LagrangianCompliantR::setComputeJacobianhOver_lambdaFunction(
    const siconos::modeling::func_prototypes::FunctionBVV_M& func) {
  computejacobianhOver_lambda_ = func;
}

void siconos::modeling::LagrangianCompliantR::computeJacobianhOver_lambda(
    const siconos::algebra::BlockVector& positions,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda) {
  if (computejacobianhOver_lambda_)
    computejacobianhOver_lambda_(positions, lambda, *jacobianhOver_lambda_);
}

void siconos::modeling::LagrangianCompliantR::computeOutput(double time, Interaction& inter,
                                                            unsigned int derivativeNumber) {
  auto& DSlink = inter.linkToDSVariables();
  if (derivativeNumber == 0) {
    auto& y = *inter.y(0);
    auto& lambda = *inter.lambda(0);
    computeh(*DSlink[tools::enum_to_index(WorkDS::q0)], lambda, y);
  } else {
    auto& y = *inter.y(derivativeNumber);
    auto& lambda = *inter.lambda(derivativeNumber);
    computeJacobianhOver_q(*DSlink[tools::enum_to_index(WorkDS::q0)], lambda);
    computeJacobianhOver_lambda(*DSlink[tools::enum_to_index(WorkDS::q0)], lambda);
    if (derivativeNumber == 1) {
      // y = Jach[0] q1 + Jach[1] lambda
      siconos::algebra::matrixBlockVector_prod(*jacobianhOver_q_view_,
                                               *DSlink[tools::enum_to_index(WorkDS::q1)], y);
      y += *jacobianhOver_lambda_ * lambda;
    } else if (derivativeNumber == 2)
      siconos::algebra::matrixBlockVector_prod(
          *jacobianhOver_q_view_, *DSlink[tools::enum_to_index(WorkDS::q2)],
          y);  // Approx: y[2] = Jach[0]q[2], other terms are neglected ...
    else
      THROW_EXCEPTION(
          "siconos::modeling::LagrangianCompliantR::computeOutput, index out of range or not "
          "yet implemented.");
  }
}

void siconos::modeling::LagrangianCompliantR::computeInput(double time, Interaction& inter,
                                                           unsigned int level) {
  // get lambda of the concerned interaction

  auto& lambda = *inter.lambda(level);
  auto& DSlink = inter.linkToDSVariables();
  computeJacobianhOver_q(*DSlink[tools::enum_to_index(WorkDS::q0)], lambda);
  // data[name] += trans(G) * lambda
  siconos::algebra::transposeMatrixVector_prod_toBlock(
      lambda, *jacobianhOver_q_view_, *DSlink[tools::enum_to_index(WorkDS::p0) + level],
      false);
}

void siconos::modeling::LagrangianCompliantR::computeJach(double time, Interaction& inter) {
  auto& DSlink = inter.linkToDSVariables();
  auto& lambda = *inter.lambda(0);
  computeJacobianhOver_q(*DSlink[tools::enum_to_index(WorkDS::q0)], lambda);
  computeJacobianhOver_lambda(*DSlink[tools::enum_to_index(WorkDS::q0)], lambda);
}
