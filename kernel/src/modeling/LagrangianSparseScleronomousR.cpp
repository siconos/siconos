/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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

#include "LagrangianSparseScleronomousR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosAlgebraAddons.hpp"  // for matrix-vector prod
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "Tools.hpp"

// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
#include "SiconosException.hpp"
#include "siconos_debug.h"

void siconos::modeling::LagrangianSparseScleronomousR::initialize(Interaction& inter) {
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

  // -- \frac{\partial}{\partial t}(J^h_q) --
  //  --> Optional
  // Case 1: user-defined function --> owned storage required and not set with
  // a call to setCompute...function
  if (computejacobianhOver_q_dot_) {
    if (!std::holds_alternative<siconos::algebra::OwnedSparse>(jacobianhOver_q_dot_storage_)) {
      // not owned storage? --> allocate
      jacobianhOver_q_dot_storage_ =
          std::make_unique<siconos::algebra::SiconosSparseMatrix>(sizeY, sizeDS);
    } else {
      // storage done and owned --> check dimensions
      bool are_dim_complient = false;
      useJacobianhOver_q_dot([&are_dim_complient, sizeY, sizeDS](const auto& M) {
        are_dim_complient = M.rows() != sizeY || M.cols() != sizeDS;
      });
      if (!are_dim_complient)
        jacobianhOver_q_dot_storage_ =
            std::make_unique<siconos::algebra::SiconosSparseMatrix>(sizeY, sizeDS);
    }
  } else {
    // Other cases:
    // - previous call to setConstantJacobianhOver_q with alias or copy
    //   --> memory allocated in both cases (owned or not)
    //   --> just check if dimensions are compliant with the interaction

    // Check if the variant is not empty (i.e. monostate), if so, check dimensions
    if (!std::holds_alternative<std::monostate>(jacobianhOver_q_dot_storage_)) {
      // -> previous call to setConstant or setCompute
      useJacobianhOver_q_dot([&](const auto& M) {
        assert(M.rows() == sizeY);
        assert(M.cols() == sizeDS);
      });
    }
  }
}

void siconos::modeling::LagrangianSparseScleronomousR::setComputehFunction(
    const siconos::modeling::func_prototypes::FunctionBV_V& h_func) {
  computeh_ = h_func;
}

void siconos::modeling::LagrangianSparseScleronomousR::computeh(
    const siconos::algebra::BlockVector& positions,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  if (computeh_) computeh_(positions, y);
}

void siconos::modeling::LagrangianSparseScleronomousR::setConstantJacobianhOver_q(
    Eigen::Map<siconos::algebra::SiconosSparseMatrix>& newValue, siconos::algebra::AliasTag) {
  jacobianhOver_q_storage_ =
      std::make_shared<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>(newValue);

  // Warning: we can't check dimensions. This will be done during initialize.
  hasJacobianhOver_q_ = true;
  hasConstantJacobianhOver_q_ = true;
  computejacobianhOver_q_ = nullptr;
  computejacobianhOver_q_python_ = nullptr;
}

void siconos::modeling::LagrangianSparseScleronomousR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector& positions) {
  // Warning: the jacobian must be allocated and of the right size before a call to this
  // function
  auto isComputed = siconos::algebra::computeSparse(
      jacobianhOver_q_storage_, computejacobianhOver_q_, computejacobianhOver_q_python_,
      "jacobianhOver_q_storage_", positions);
  // isComputed == false => means jacobian is constant, functions not set
}

void siconos::modeling::LagrangianSparseScleronomousR::setConstantJacobianhOver_q_dot(
    Eigen::Map<siconos::algebra::SiconosSparseMatrix>& newValue, siconos::algebra::AliasTag) {
  jacobianhOver_q_dot_storage_ =
      std::make_shared<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>(newValue);

  // Warning: we can't check dimensions. This will be done during initialize.
  hasJacobianhOver_q_dot_ = true;
  hasConstantJacobianhOver_q_dot_ = true;
  computejacobianhOver_q_dot_ = nullptr;
  computejacobianhOver_q_dot_python_ = nullptr;
}

void siconos::modeling::LagrangianSparseScleronomousR::computeJacobianhOver_q_dot(
    const siconos::algebra::BlockVector& q, const siconos::algebra::BlockVector& qdot) {
  auto isComputed = siconos::algebra::computeSparse(
      jacobianhOver_q_dot_storage_, computejacobianhOver_q_dot_,
      computejacobianhOver_q_dot_python_, "jacobianhOver_q_storage_", q, qdot);
  // isComputed == false => means jacobian is constant, functions not set
}

void siconos::modeling::LagrangianSparseScleronomousR::computeJacobianhOver_q_dot_X_qdot(
    double time, Interaction& inter) {
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  DEBUG_PRINT("siconos::modeling::LagrangianSparseScleronomousR::computedotjacqhXqdot starts");

  computeJacobianhOver_q_dot(*ds_vars[tools::enum_to_index(ds_var::q0)],
                             *ds_vars[tools::enum_to_index(ds_var::q1)]);
  if (!jacobianhOver_q_dot_X_qdot_) {
    useJacobianhOver_q_dot([&](const auto& M) {
      jacobianhOver_q_dot_X_qdot_ =
          std::make_shared<siconos::algebra::SiconosVector>(M.rows());

      siconos::algebra::matrixBlockVector_prod(M, *ds_vars[tools::enum_to_index(ds_var::q1)],
                                               *jacobianhOver_q_dot_X_qdot_);
    });
  }
  DEBUG_PRINT("siconos::modeling::LagrangianSparseScleronomousR::computedotjacqhXqdot ends");
}

void siconos::modeling::LagrangianSparseScleronomousR::computeJach(double time,
                                                                   Interaction& inter) {
  DEBUG_BEGIN(
      "void siconos::modeling::LagrangianSparseScleronomousR::computeJach(double time, "
      "Interaction& "
      "inter) \n");
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  DEBUG_EXPR(siconos::algebra::print(inter););
  computeJacobianhOver_q(*ds_vars[tools::enum_to_index(ds_var::q0)]);
  computeJacobianhOver_q_dot(*ds_vars[tools::enum_to_index(ds_var::q0)],
                             *ds_vars[tools::enum_to_index(ds_var::q1)]);

  DEBUG_END(
      "void siconos::modeling::LagrangianSparseScleronomousR::computeJach(double time, "
      "Interaction& "
      "inter) \n");
}

void siconos::modeling::LagrangianSparseScleronomousR::computeOutput(
    double time, Interaction& inter, siconos::algebra::blocks::size_type derivativeNumber) {
  DEBUG_PRINTF(
      "siconos::modeling::LagrangianSparseScleronomousR::computeOutput(double time, "
      "Interaction& "
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
    computeJacobianhOver_q(*ds_vars[tools::enum_to_index(ds_var::q0)]);

    if (derivativeNumber == 1) {
      useJacobianhOver_q([&](const auto& M) {
        siconos::algebra::matrixBlockVector_prod(M, *ds_vars[tools::enum_to_index(ds_var::q1)],
                                                 y);
      });
    } else if (derivativeNumber == 2) {
      useJacobianhOver_q([&](const auto& M) {
        siconos::algebra::matrixBlockVector_prod(M, *ds_vars[tools::enum_to_index(ds_var::q2)],
                                                 y);
      });

      computeJacobianhOver_q_dot(*ds_vars[tools::enum_to_index(ds_var::q0)],
                                 *ds_vars[tools::enum_to_index(ds_var::q1)]);
      useJacobianhOver_q_dot([&](const auto& M) {
        siconos::algebra::matrixBlockVector_prod(M, *ds_vars[tools::enum_to_index(ds_var::q1)],
                                                 y);
      });
    } else
      THROW_EXCEPTION(
          "siconos::modeling::LagrangianSparseScleronomousR::computeOutput(double time, "
          "Interaction& inter, InteractionProperties& interProp, "
          "siconos::algebra::blocks::size_type "
          "derivativeNumber), index out of range");
  }
}

void siconos::modeling::LagrangianSparseScleronomousR::computeInput(
    double time, Interaction& inter, siconos::algebra::blocks::size_type level) {
  DEBUG_BEGIN(
      "void siconos::modeling::LagrangianSparseScleronomousR::computeInput(double time, "
      "Interaction& inter, InteractionProperties& interProp, "
      "siconos::algebra::blocks::size_type level) \n");

  DEBUG_PRINTF("level = %i\n", level);
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  computeJacobianhOver_q(*ds_vars[tools::enum_to_index(ds_var::q0)]);
  // get lambda of the concerned interaction
  auto& lambda = *inter.lambda(level);
  DEBUG_EXPR(siconos::algebra::print(lambda););
  // data[name] += trans(G) * lambda

  useJacobianhOver_q([&](const auto& M) {
    siconos::algebra::transposeMatrixVector_prod_toBlock(
        lambda, M, *ds_vars[tools::enum_to_index(ds_var::p0) + level], false);
  });

  DEBUG_EXPR(siconos::algebra::print(*ds_vars[tools::enum_to_index(ds_var::p0) + level]););
  DEBUG_END(
      "void siconos::modeling::LagrangianSparseScleronomousR::computeInput(double time, "
      "Interaction& inter, InteractionProperties& interProp, "
      "siconos::algebra::blocks::size_type level) \n");
}
