
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
#include "PluggedObject.hpp"
#include "PluggedObject.hpp"  // getPluginFunctionname ...
#include "PluginTypes.hpp"    // FPtr4 ...
#include "SiconosMatrix.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
#include "Tools.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "SiconosException.hpp"
#include "siconos_debug.h"

void siconos::modeling::LagrangianRheonomousR::initialize(Interaction& inter) {
  auto sizeY = inter.dimension();
  auto sizeDS = inter.getSizeOfDS();
  if (!jacobianhOver_q_internal_storage_) {
    jacobianhOver_q_internal_storage_ = std::make_unique<std::vector<double>>(sizeY * sizeDS);
  }
  jacobianhOver_q_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobianhOver_q_internal_storage_->data(), sizeY, sizeDS);

  if (computehdot_) {  // True only if setComputehdotFunction has been called
    // Ensure that memory is properly allocated for hdot_
    if (!hdot_) {
      hdot_ = std::make_shared<siconos::algebra::SiconosVector>(sizeY);
    }
  }
}

void siconos::modeling::LagrangianRheonomousR::setComputehdotFunction(
    const siconos::modeling::func_prototypes::FunctionBVS_V& hdot_func) {
  // Warning: we can't allocate hdot here, since we have no clue on the
  // required size
  computehdot_ = hdot_func;
}

void siconos::modeling::LagrangianRheonomousR::computehdot(
    const siconos::algebra::BlockVector& positions, double time) {
  if (computehdot_)
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    computehdot_(positions, time, *hdot_);
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

void siconos::modeling::LagrangianRheonomousR::setComputeJacobianhOver_qFunction(
    const siconos::modeling::func_prototypes::FunctionBVS_M& func) {
  // Warning: we can't allocate hdot here, since we have no clue on the
  // required size
  computejacobianhOver_q_ = func;
}

void siconos::modeling::LagrangianRheonomousR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector& positions, double time) {
  if (computejacobianhOver_q_)
    computejacobianhOver_q_(positions, time, *jacobianhOver_q_view_);
}

void siconos::modeling::LagrangianRheonomousR::computeOutput(double time, Interaction& inter,
                                                             unsigned int derivativeNumber) {
  auto& DSlink = inter.linkToDSVariables();
  auto& y = *inter.y(derivativeNumber);
  if (derivativeNumber == 0)
    computeh(*DSlink[tools::enum_to_index(WorkDS::q0)], time, y);
  else {
    computeJacobianhOver_q(*DSlink[tools::enum_to_index(WorkDS::q0)], time);
    if (derivativeNumber == 1) {
      if (!hdot_) {
        auto sizeY = inter.dimension();
        hdot_ = std::make_shared<siconos::algebra::SiconosVector>(sizeY);
      }
      // Computation of the partial derivative w.r.t time of h(q,t)
      computehdot(*DSlink[tools::enum_to_index(WorkDS::q0)], time);
      // Computation of the partial derivative w.r.t q of h(q,t) : \nabla_q h(q,t) \dot q
      siconos::algebra::matrixBlockVector_prod(*jacobianhOver_q_view_,
                                               *DSlink[tools::enum_to_index(WorkDS::q1)], y);
      // Sum of the terms
      y += *hdot_;
    } else if (derivativeNumber == 2) {
      siconos::algebra::matrixBlockVector_prod(*jacobianhOver_q_view_,
                                               *DSlink[tools::enum_to_index(WorkDS::q2)],
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

void siconos::modeling::LagrangianRheonomousR::computeInput(double time, Interaction& inter,
                                                            unsigned int level) {
  auto& DSlink = inter.linkToDSVariables();
  computeJacobianhOver_q(*DSlink[tools::enum_to_index(WorkDS::q0)], time);
  // get lambda of the concerned interaction
  auto& lambda = *inter.lambda(level);
  // data[name] += trans(G) * lambda
  siconos::algebra::transposeMatrixVector_prod_toBlock(
      lambda, *jacobianhOver_q_view_, *DSlink[tools::enum_to_index(WorkDS::p0) + level], false);
}

void siconos::modeling::LagrangianRheonomousR::computeJach(double time, Interaction& inter) {
  auto& DSlink = inter.linkToDSVariables();
  computeJacobianhOver_q(*DSlink[tools::enum_to_index(WorkDS::q0)], time);
  // computeJachqDot(time, inter);
  //    computeDotJachq(time, q, z);
  // computeJachlambda(time, inter);
  computehdot(*DSlink[tools::enum_to_index(WorkDS::q0)], time);
}
