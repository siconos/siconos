
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

#include "LagrangianSparseRheonomousR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosAlgebraAddons.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
#include "Tools.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "SiconosException.hpp"
#include "siconos_debug.h"

void siconos::modeling::LagrangianSparseRheonomousR::initialize(Interaction& inter) {
  auto sizeY = inter.dimension();
  auto sizeDS = inter.getSizeOfDS();

  // -- Jacobian h over q --
  //  --> Always required
  // Case 1: user-defined function --> owned storage required and not set with
  // a call to setCompute...function
  if (computejacobianhOver_q_) {
    if (!std::holds_alternative<siconos::algebra::OwnedSparse>(jacobianhOver_q_storage_)) {
      // not owned storage? --> allocate
      jacobianhOver_q_storage_ =
          std::make_unique<siconos::algebra::SiconosSparseMatrix>(sizeY, sizeDS);
    } else {
      // storage done and owned --> check dimensions
      bool are_dim_complient = false;
      useJacobianhOver_q([&are_dim_complient, sizeY, sizeDS](const auto& M) {
        are_dim_complient = M.rows() != sizeY || M.cols() != sizeDS;
      });
      if (!are_dim_complient)
        jacobianhOver_q_storage_ =
            std::make_unique<siconos::algebra::SiconosSparseMatrix>(sizeY, sizeDS);
    }
  } else {
    // Other cases:
    // - previous call to setConstantJacobianhOver_q with alias or copy
    //   --> memory allocated in both cases (owned or not)
    //   --> just check if dimensions are compliant with the interaction
    // - no call to any set... --> error

    // Check if the variant is empty (i.e. monostate), if so error.
    if (std::holds_alternative<std::monostate>(jacobianhOver_q_storage_)) {
      // -> no previous call to setConstant or setCompute, this is not allowed
      THROW_EXCEPTION(
          "You must call setConstantJacobianhOver_q or setComputeJacobianhOver_qFunction "
          "before "
          "initialize.")
    }
    useJacobianhOver_q([&](const auto& M) {
      assert(M.rows() == sizeY);
      assert(M.cols() == sizeDS);
    });
  }

  if (computehdot_) {  // True only if setComputehdotFunction has been called
    // Ensure that memory is properly allocated for hdot_
    if (!hdot_) {
      hdot_ = std::make_shared<siconos::algebra::SiconosVector>(sizeY);
    }
  }
}

void siconos::modeling::LagrangianSparseRheonomousR::setComputehdotFunction(
    const siconos::modeling::func_prototypes::FunctionBVS_V& hdot_func) {
  // Warning: we can't allocate hdot here, since we have no clue on the
  // required size
  computehdot_ = hdot_func;
}

void siconos::modeling::LagrangianSparseRheonomousR::setConstantJacobianhOver_q(
    Eigen::Map<siconos::algebra::SiconosSparseMatrix>& newValue, siconos::algebra::AliasTag) {
  jacobianhOver_q_storage_ =
      std::make_shared<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>(newValue);

  // Warning: we can't check dimensions. This will be done during initialize.
  hasJacobianhOver_q_ = true;
  hasConstantJacobianhOver_q_ = true;
  computejacobianhOver_q_ = nullptr;
  computejacobianhOver_q_python_ = nullptr;
}
void siconos::modeling::LagrangianSparseRheonomousR::computehdot(
    const siconos::algebra::BlockVector& positions, double time) {
  if (computehdot_)
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    computehdot_(positions, time, *hdot_);
}

void siconos::modeling::LagrangianSparseRheonomousR::setComputehFunction(
    const siconos::modeling::func_prototypes::FunctionBVS_V& h_func) {
  computeh_ = h_func;
}

void siconos::modeling::LagrangianSparseRheonomousR::computeh(
    const siconos::algebra::BlockVector& positions, double time,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  DEBUG_PRINT(" siconos::modeling::LagrangianSparseRheonomousR::computeh()")
  if (computeh_) computeh_(positions, time, y);
}

void siconos::modeling::LagrangianSparseRheonomousR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector& positions, double time) {
  // Warning: the jacobian must be allocated and of the right size before a call to this
  // function
  auto isComputed = siconos::algebra::computeSparse(
      jacobianhOver_q_storage_, computejacobianhOver_q_, computejacobianhOver_q_python_,
      "jacobianhOver_q_storage_", positions, time);
  // isComputed == false => means jacobian is constant, functions not set
}

void siconos::modeling::LagrangianSparseRheonomousR::computeOutput(
    double time, Interaction& inter, siconos::algebra::blocks::size_type derivativeNumber) {
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  auto& y = *inter.y(derivativeNumber);
  if (derivativeNumber == 0)
    computeh(*ds_vars[tools::enum_to_index(ds_var::q0)], time, y);
  else {
    computeJacobianhOver_q(*ds_vars[tools::enum_to_index(ds_var::q0)], time);
    if (derivativeNumber == 1) {
      if (!hdot_) {
        auto sizeY = inter.dimension();
        hdot_ = std::make_shared<siconos::algebra::SiconosVector>(sizeY);
      }
      // Computation of the partial derivative w.r.t time of h(q,t)
      computehdot(*ds_vars[tools::enum_to_index(ds_var::q0)], time);
      // Computation of the partial derivative w.r.t q of h(q,t) : \nabla_q h(q,t) \dot q
      useJacobianhOver_q([&](const auto& M) {
        siconos::algebra::matrixBlockVector_prod(M, *ds_vars[tools::enum_to_index(ds_var::q1)],
                                                 y);
      });
      // Sum of the terms
      y += *hdot_;
    } else if (derivativeNumber == 2) {
      useJacobianhOver_q([&](const auto& M) {
        siconos::algebra::matrixBlockVector_prod(M, *ds_vars[tools::enum_to_index(ds_var::q2)],
                                                 y);
      });
      // \warning : the computation of y[2] (in event-driven
      // simulation for instance) is approximated by y[2] =
      // Jach[0]q[2]. For the moment, other terms are neglected
      // (especially, partial derivatives with respect to time).
    } else
      THROW_EXCEPTION(
          "siconos::modeling::LagrangianSparseRheonomousR::computeOutput index >2 not yet "
          "implemented.");
  }
}

void siconos::modeling::LagrangianSparseRheonomousR::computeInput(
    double time, Interaction& inter, siconos::algebra::blocks::size_type level) {
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  computeJacobianhOver_q(*ds_vars[tools::enum_to_index(ds_var::q0)], time);
  // get lambda of the concerned interaction
  auto& lambda = *inter.lambda(level);
  // data[name] += trans(G) * lambda
  useJacobianhOver_q([&](const auto& M) {
    siconos::algebra::transposeMatrixVector_prod_toBlock(
        lambda, M, *ds_vars[tools::enum_to_index(ds_var::p0) + level], false);
  });
}

void siconos::modeling::LagrangianSparseRheonomousR::computeJach(double time,
                                                                 Interaction& inter) {
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  computeJacobianhOver_q(*ds_vars[tools::enum_to_index(ds_var::q0)], time);
  // computeJachqDot(time, inter);
  //    computeDotJachq(time, q, z);
  // computeJachlambda(time, inter);
  computehdot(*ds_vars[tools::enum_to_index(ds_var::q0)], time);
}
