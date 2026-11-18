
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

#include "LagrangianRheonomousR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosAlgebraAddons.hpp"  // for matrix-vector prod
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "Tools.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "SiconosException.hpp"
#include "siconos_debug.h"

void siconos::modeling::LagrangianRheonomousR::initialize(Interaction& inter) {
  auto sizeY = inter.dimension();
  auto sizeDS = inter.getSizeOfDS();
  if (!H_storage_) {
    H_storage_ = std::make_unique<std::vector<double>>(sizeY * sizeDS);
  }
  H_view_ = std::make_shared<siconos::algebra::MapType>(H_storage_->data(), sizeY, sizeDS);

  if (computedhdt_) {  // True only if setComputedhdtFunction has been called
    // Ensure that memory is properly allocated for dhdt_
    if (!dhdt_) {
      dhdt_ = std::make_shared<siconos::algebra::SiconosVector>(sizeY);
    }
  }
}

void siconos::modeling::LagrangianRheonomousR::setComputedhdtFunction(
    const siconos::modeling::func_prototypes::FunctionBVS_V& hdot_func) {
  // Warning: we can't allocate dhdt here, since we have no clue on the
  // required size
  computedhdt_ = hdot_func;
}

void siconos::modeling::LagrangianRheonomousR::computedhdt(
    const siconos::algebra::BlockVector& positions, double time) {
  if (!computedhdtLegacy_) {  // dispatch to a class overriding the deprecated name
    detail::LegacyCallGuard guard(computedhdtLegacy_);
    SICONOS_IGNORE_DEPRECATED_BEGIN
    computehdot(positions, time);
    SICONOS_IGNORE_DEPRECATED_END
    return;
  }
  if (computedhdt_)
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    computedhdt_(positions, time, *dhdt_);
}

void siconos::modeling::LagrangianRheonomousR::setComputehFunction(
    const siconos::modeling::func_prototypes::FunctionBVS_V& h_func) {
  computeh_ = h_func;
}

void siconos::modeling::LagrangianRheonomousR::computeh(
    const siconos::algebra::BlockVector& positions, double time,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  DEBUG_PRINT(" siconos::modeling::LagrangianRheonomousR::computeh()")
  if (computeh_) computeh_(positions, time, y);
}

void siconos::modeling::LagrangianRheonomousR::setComputeHFunction(
    const siconos::modeling::func_prototypes::FunctionBVS_M& func) {
  // Warning: we can't allocate dhdt here, since we have no clue on the
  // required size
  computeH_ = func;
}

void siconos::modeling::LagrangianRheonomousR::computeH(
    const siconos::algebra::BlockVector& positions, double time) {
  if (!computeHLegacy_) {  // dispatch to a class overriding the deprecated name
    detail::LegacyCallGuard guard(computeHLegacy_);
    SICONOS_IGNORE_DEPRECATED_BEGIN
    computeJacobianhOver_q(positions, time);
    SICONOS_IGNORE_DEPRECATED_END
    return;
  }
  if (computeH_) computeH_(positions, time, *H_view_);
}

void siconos::modeling::LagrangianRheonomousR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector& positions, double time) {
  computeH(positions, time);
}

void siconos::modeling::LagrangianRheonomousR::computehdot(
    const siconos::algebra::BlockVector& positions, double time) {
  computedhdt(positions, time);
}

void siconos::modeling::LagrangianRheonomousR::computeOutput(
    double time, Interaction& inter, siconos::algebra::blocks::size_type derivativeNumber) {
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  auto& y = *inter.y(derivativeNumber);
  if (derivativeNumber == 0)
    computeh(*ds_vars[tools::enum_to_index(ds_var::q0)], time, y);
  else {
    computeH(*ds_vars[tools::enum_to_index(ds_var::q0)], time);
    if (derivativeNumber == 1) {
      if (!dhdt_) {
        auto sizeY = inter.dimension();
        dhdt_ = std::make_shared<siconos::algebra::SiconosVector>(sizeY);
      }
      // Computation of the partial derivative w.r.t time of h(q,t)
      computedhdt(*ds_vars[tools::enum_to_index(ds_var::q0)], time);
      // Computation of the partial derivative w.r.t q of h(q,t) : \nabla_q h(q,t) \dot q
      siconos::algebra::matrixBlockVector_prod(*H_view_,
                                               *ds_vars[tools::enum_to_index(ds_var::q1)], y);
      // Sum of the terms
      y += *dhdt_;
    } else if (derivativeNumber == 2) {
      siconos::algebra::matrixBlockVector_prod(*H_view_,
                                               *ds_vars[tools::enum_to_index(ds_var::q2)],
                                               y);  // Approx:,  ...
      // \warning : the computation of y[2] (in event-driven
      // simulation for instance) is approximated by y[2] =
      // Jach[0]q[2]. For the moment, other terms are neglected
      // (especially, partial derivatives with respect to time).
    } else
      THROW_EXCEPTION(
          "siconos::modeling::LagrangianRheonomousR::computeOutput index >2 not yet "
          "implemented.");
  }
}

void siconos::modeling::LagrangianRheonomousR::computeInput(
    double time, Interaction& inter, siconos::algebra::blocks::size_type level) {
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  computeH(*ds_vars[tools::enum_to_index(ds_var::q0)], time);
  // get lambda of the concerned interaction
  auto& lambda = *inter.lambda(level);
  // data[name] += trans(G) * lambda
  siconos::algebra::transposeMatrixVector_prod_toBlock(
      lambda, *H_view_, *ds_vars[tools::enum_to_index(ds_var::p0) + level], false);
}

void siconos::modeling::LagrangianRheonomousR::computeJach(double time, Interaction& inter) {
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  computeH(*ds_vars[tools::enum_to_index(ds_var::q0)], time);
  // computeJachqDot(time, inter);
  //    computeDotJachq(time, q, z);
  // computeJachlambda(time, inter);
  computedhdt(*ds_vars[tools::enum_to_index(ds_var::q0)], time);
}
