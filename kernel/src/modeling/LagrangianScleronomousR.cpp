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
#include "SiconosAlgebraAddons.hpp"  // for matrix-vector prod
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "Tools.hpp"
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
#include "siconos_debug.h"

void siconos::modeling::LagrangianScleronomousR::initialize(Interaction& inter) {
  auto sizeY = inter.dimension();
  auto sizeDS = inter.getSizeOfDS();
  if (hasConstantH_) {
    assert(H_view_);
    assert(H_view_->rows() == sizeY);
    assert(H_view_->cols() == sizeDS);
  } else {
    if (!H_storage_) {
      H_storage_ = std::make_unique<std::vector<double>>(sizeY * sizeDS);
    }
    H_view_ = std::make_shared<siconos::algebra::MapType>(H_storage_->data(), sizeY, sizeDS);
  }

  if (hasHdot_) {
    // True only if setComputeHdotFunction has been called
    // Ensure that memory is properly allocated for hdot_
    if (!Hdot_) {
      Hdot_ = std::make_shared<siconos::algebra::SiconosMatrix>(sizeY, sizeDS);
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

void siconos::modeling::LagrangianScleronomousR::setConstantH(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  H_storage_ = nullptr;

  H_view_ = std::make_shared<siconos::algebra::MapType>(newValue.data(), newValue.rows(),
                                                        newValue.cols());
  hasConstantH_ = true;
  computeH_ = nullptr;
}

void siconos::modeling::LagrangianScleronomousR::setComputeHFunction(
    const siconos::modeling::func_prototypes::FunctionBV_M& func) {
  computeH_ = func;
}

void siconos::modeling::LagrangianScleronomousR::computeH(
    const siconos::algebra::BlockVector& positions) {
  if (!computeHLegacy_) {  // dispatch to a class overriding the deprecated name
    detail::LegacyCallGuard guard(computeHLegacy_);
    SICONOS_IGNORE_DEPRECATED_BEGIN
    computeJacobianhOver_q(positions);
    SICONOS_IGNORE_DEPRECATED_END
    return;
  }
  if (computeH_) computeH_(positions, *H_view_);
}

void siconos::modeling::LagrangianScleronomousR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector& positions) {
  computeH(positions);
}

void siconos::modeling::LagrangianScleronomousR::setComputeHdotFunction(
    const siconos::modeling::func_prototypes::FunctionBVBV_M& func) {
  hasHdot_ = true;
  computeHdot_ = func;
}

void siconos::modeling::LagrangianScleronomousR::computeHdot(
    const siconos::algebra::BlockVector& q, const siconos::algebra::BlockVector& qdot) {
  if (computeHdot_) computeHdot_(q, qdot, *Hdot_);
}

void siconos::modeling::LagrangianScleronomousR::computeHdotv(double time,
                                                              Interaction& inter) {
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  DEBUG_PRINT("siconos::modeling::LagrangianScleronomousR::computedotjacqhXqdot starts");

  computeHdot(*ds_vars[tools::enum_to_index(ds_var::q0)],
              *ds_vars[tools::enum_to_index(ds_var::q1)]);
  if (!Hdotv_) {
    Hdotv_ = std::make_shared<siconos::algebra::SiconosVector>(Hdot_->rows());
  }
  siconos::algebra::matrixBlockVector_prod(*Hdot_, *ds_vars[tools::enum_to_index(ds_var::q1)],
                                           *Hdotv_);
  DEBUG_PRINT("siconos::modeling::LagrangianScleronomousR::computedotjacqhXqdot ends");
}

void siconos::modeling::LagrangianScleronomousR::computeOutput(
    double time, Interaction& inter, siconos::algebra::blocks::size_type derivativeNumber) {
  DEBUG_PRINTF(
      "siconos::modeling::LagrangianScleronomousR::computeOutput(double time, Interaction& "
      "inter, InteractionProperties& interProp, siconos::algebra::blocks::size_type "
      "derivativeNumber) with time = "
      "%f "
      "and derivativeNumber = %i\n",
      time, derivativeNumber);
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  auto& y = *inter.y(derivativeNumber);
  if (derivativeNumber == 0) {
    computeh(*ds_vars[tools::enum_to_index(ds_var::q0)], y);
  } else {
    computeH(*ds_vars[tools::enum_to_index(ds_var::q0)]);

    if (derivativeNumber == 1) {
      siconos::algebra::matrixBlockVector_prod(*H_view_,
                                               *ds_vars[tools::enum_to_index(ds_var::q1)], y);
    } else if (derivativeNumber == 2) {
      siconos::algebra::matrixBlockVector_prod(*H_view_,
                                               *ds_vars[tools::enum_to_index(ds_var::q2)], y);
      computeHdot(*ds_vars[tools::enum_to_index(ds_var::q0)],
                  *ds_vars[tools::enum_to_index(ds_var::q1)]);
      siconos::algebra::matrixBlockVector_prod(
          *Hdot_, *ds_vars[tools::enum_to_index(ds_var::q1)], y, false);

    } else
      THROW_EXCEPTION(
          "siconos::modeling::LagrangianScleronomousR::computeOutput(double time, "
          "Interaction& inter, InteractionProperties& interProp, "
          "siconos::algebra::blocks::size_type "
          "derivativeNumber), index out of range");
  }
}

void siconos::modeling::LagrangianScleronomousR::computeInput(
    double time, Interaction& inter, siconos::algebra::blocks::size_type level) {
  DEBUG_BEGIN(
      "void siconos::modeling::LagrangianScleronomousR::computeInput(double time, "
      "Interaction& inter, InteractionProperties& interProp, "
      "siconos::algebra::blocks::size_type level) \n");

  DEBUG_PRINTF("level = %i\n", level);
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  computeH(*ds_vars[tools::enum_to_index(ds_var::q0)]);
  // get lambda of the concerned interaction
  auto& lambda = *inter.lambda(level);
  DEBUG_EXPR(siconos::algebra::print(lambda););
  // data[name] += trans(G) * lambda
  siconos::algebra::transposeMatrixVector_prod_toBlock(
      lambda, *H_view_, *ds_vars[tools::enum_to_index(ds_var::p0) + level], false);
  DEBUG_EXPR(siconos::algebra::print(*ds_vars[tools::enum_to_index(ds_var::p0) + level]););
  DEBUG_END(
      "void siconos::modeling::LagrangianScleronomousR::computeInput(double time, "
      "Interaction& inter, InteractionProperties& interProp, "
      "siconos::algebra::blocks::size_type level) \n");
}

void siconos::modeling::LagrangianScleronomousR::computeJach(double time, Interaction& inter) {
  DEBUG_BEGIN(
      "void siconos::modeling::LagrangianScleronomousR::computeJach(double time, "
      "Interaction& "
      "inter) \n");
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  DEBUG_EXPR(siconos::algebra::print(inter););
  computeH(*ds_vars[tools::enum_to_index(ds_var::q0)]);
  computeHdot(*ds_vars[tools::enum_to_index(ds_var::q0)],
              *ds_vars[tools::enum_to_index(ds_var::q1)]);

  DEBUG_END(
      "void siconos::modeling::LagrangianScleronomousR::computeJach(double time, "
      "Interaction& "
      "inter) \n");
}
