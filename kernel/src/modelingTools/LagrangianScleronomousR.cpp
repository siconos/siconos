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

#include "LagrangianScleronomousR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
#include "Tools.hpp"

// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
#include "SiconosException.hpp"
#include "siconos_debug.h"

void siconos::modeling::LagrangianScleronomousR::initialize(Interaction& inter) {
  auto sizeY = inter.dimension();
  auto sizeDS = inter.getSizeOfDS();
  if (hasConstantJacobianhOver_q_) {
    assert(jacobianhOver_q_view_);
    assert(jacobianhOver_q_view_->rows() == sizeY);
    assert(jacobianhOver_q_view_->cols() == sizeDS);
  } else {
    if (!jacobianhOver_q_internal_storage_) {
      jacobianhOver_q_internal_storage_ =
          std::make_unique<std::vector<double>>(sizeY * sizeDS);
    }
    jacobianhOver_q_view_ = std::make_shared<siconos::algebra::MapType>(
        jacobianhOver_q_internal_storage_->data(), sizeY, sizeDS);
  }

  if (hasJacobianhOver_q_dot_) {
    // True only if setComputejacobianhOver_q_dotFunction has been called
    // Ensure that memory is properly allocated for hdot_
    if (!jacobianhOver_q_dot_) {
      jacobianhOver_q_dot_ = std::make_shared<siconos::algebra::SiconosMatrix>(sizeY, sizeDS);
    }
  }
}

void siconos::modeling::LagrangianScleronomousR::setComputehFunction(
    const siconos::modeling::func_prototypes::FunctionBV_V& h_func) {
  computeh_ = h_func;
}

void siconos::modeling::LagrangianScleronomousR::computeh(
    const siconos::algebra::BlockVector& positions,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  if (computeh_) computeh_(positions, y);
}

void siconos::modeling::LagrangianScleronomousR::setConstantJacobianhOver_q(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  jacobianhOver_q_internal_storage_ = nullptr;

  jacobianhOver_q_view_ = std::make_shared<siconos::algebra::MapType>(
      newValue.data(), newValue.rows(), newValue.cols());
  hasConstantJacobianhOver_q_ = true;
  computejacobianhOver_q_ = nullptr;
}

void siconos::modeling::LagrangianScleronomousR::setComputeJacobianhOver_qFunction(
    const siconos::modeling::func_prototypes::FunctionBV_M& func) {
  computejacobianhOver_q_ = func;
}

void siconos::modeling::LagrangianScleronomousR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector& positions) {
  if (computejacobianhOver_q_) computejacobianhOver_q_(positions, *jacobianhOver_q_view_);
}

void siconos::modeling::LagrangianScleronomousR::setComputejacobianhOver_q_dotFunction(
    const siconos::modeling::func_prototypes::FunctionBVBV_M& func) {
  hasJacobianhOver_q_dot_ = true;
  computejacobianhOver_q_dot_ = func;
}

void siconos::modeling::LagrangianScleronomousR::computejacobianhOver_q_dot(
    const siconos::algebra::BlockVector& q, const siconos::algebra::BlockVector& qdot) {
  if (computejacobianhOver_q_dot_) computejacobianhOver_q_dot_(q, qdot, *jacobianhOver_q_dot_);
}

void siconos::modeling::LagrangianScleronomousR::computeJacobianhOver_q_dot_X_qdot(
    double time, Interaction& inter) {
  auto& DSlink = inter.linkToDSVariables();
  DEBUG_PRINT("siconos::modeling::LagrangianScleronomousR::computedotjacqhXqdot starts");

  computejacobianhOver_q_dot(*DSlink[tools::enum_to_index(WorkDS::q0)],
                             *DSlink[tools::enum_to_index(WorkDS::q1)]);
  if (!jacobianhOver_q_dot_X_qdot_) {
    jacobianhOver_q_dot_X_qdot_ =
        std::make_shared<siconos::algebra::SiconosVector>(jacobianhOver_q_dot_->rows());
  }
  siconos::algebra::matrixBlockVector_prod(*jacobianhOver_q_dot_,
                                           *DSlink[tools::enum_to_index(WorkDS::q1)],
                                           *jacobianhOver_q_dot_X_qdot_);
  DEBUG_PRINT("siconos::modeling::LagrangianScleronomousR::computedotjacqhXqdot ends");
}

void siconos::modeling::LagrangianScleronomousR::computeOutput(double time, Interaction& inter,
                                                               unsigned int derivativeNumber) {
  DEBUG_PRINTF(
      "siconos::modeling::LagrangianScleronomousR::computeOutput(double time, Interaction& "
      "inter, InteractionProperties& interProp, unsigned int derivativeNumber) with time = "
      "%f "
      "and derivativeNumber = %i\n",
      time, derivativeNumber);
  auto& DSlink = inter.linkToDSVariables();
  auto& y = *inter.y(derivativeNumber);
  if (derivativeNumber == 0) {
    computeh(*DSlink[tools::enum_to_index(WorkDS::q0)], y);
  } else {
    computeJacobianhOver_q(*DSlink[tools::enum_to_index(WorkDS::q0)]);

    if (derivativeNumber == 1) {
      siconos::algebra::matrixBlockVector_prod(*jacobianhOver_q_view_,
                                               *DSlink[tools::enum_to_index(WorkDS::q1)], y);
    } else if (derivativeNumber == 2) {
      siconos::algebra::matrixBlockVector_prod(*jacobianhOver_q_view_,
                                               *DSlink[tools::enum_to_index(WorkDS::q2)], y);
      computejacobianhOver_q_dot(*DSlink[tools::enum_to_index(WorkDS::q0)],
                                 *DSlink[tools::enum_to_index(WorkDS::q1)]);
      siconos::algebra::matrixBlockVector_prod(
          *jacobianhOver_q_dot_, *DSlink[tools::enum_to_index(WorkDS::q1)], y, false);

    } else
      THROW_EXCEPTION(
          "siconos::modeling::LagrangianScleronomousR::computeOutput(double time, "
          "Interaction& inter, InteractionProperties& interProp, unsigned int "
          "derivativeNumber), index out of range");
  }
}

void siconos::modeling::LagrangianScleronomousR::computeInput(double time, Interaction& inter,
                                                              unsigned int level) {
  DEBUG_BEGIN(
      "void siconos::modeling::LagrangianScleronomousR::computeInput(double time, "
      "Interaction& inter, InteractionProperties& interProp, unsigned int level) \n");

  DEBUG_PRINTF("level = %i\n", level);
  auto& DSlink = inter.linkToDSVariables();
  computeJacobianhOver_q(*DSlink[tools::enum_to_index(WorkDS::q0)]);
  // get lambda of the concerned interaction
  auto& lambda = *inter.lambda(level);
  DEBUG_EXPR(siconos::algebra::print(lambda););
  // data[name] += trans(G) * lambda
  siconos::algebra::transposeMatrixVector_prod_toBlock(
      lambda, *jacobianhOver_q_view_, *DSlink[tools::enum_to_index(WorkDS::p0) + level],
      false);
  DEBUG_EXPR(siconos::algebra::print(*DSlink[tools::enum_to_index(WorkDS::p0) + level]););
  DEBUG_END(
      "void siconos::modeling::LagrangianScleronomousR::computeInput(double time, "
      "Interaction& inter, InteractionProperties& interProp, unsigned int level) \n");
}

void siconos::modeling::LagrangianScleronomousR::computeJach(double time, Interaction& inter) {
  DEBUG_BEGIN(
      "void siconos::modeling::LagrangianScleronomousR::computeJach(double time, "
      "Interaction& "
      "inter) \n");
  auto& DSlink = inter.linkToDSVariables();
  DEBUG_EXPR(siconos::algebra::print(inter););
  computeJacobianhOver_q(*DSlink[tools::enum_to_index(WorkDS::q0)]);
  computejacobianhOver_q_dot(*DSlink[tools::enum_to_index(WorkDS::q0)],
                             *DSlink[tools::enum_to_index(WorkDS::q1)]);

  DEBUG_END(
      "void siconos::modeling::LagrangianScleronomousR::computeJach(double time, "
      "Interaction& "
      "inter) \n");
}
