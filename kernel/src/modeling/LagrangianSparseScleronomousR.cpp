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
  if (computeH_) {
    if (!std::holds_alternative<siconos::algebra::OwnedSparse>(H_storage_)) {
      // not owned storage? --> allocate
      H_storage_ = std::make_unique<siconos::algebra::SiconosSparseMatrix>(sizeY, sizeDS);
    } else {
      // storage done and owned --> check dimensions
      bool are_dim_complient = false;
      useH([&are_dim_complient, sizeY, sizeDS](const auto& M) {
        are_dim_complient = M.rows() != sizeY || M.cols() != sizeDS;
      });
      if (!are_dim_complient)
        H_storage_ = std::make_unique<siconos::algebra::SiconosSparseMatrix>(sizeY, sizeDS);
    }
  } else {
    // Other cases:
    // - previous call to setConstantH with alias or copy
    //   --> memory allocated in both cases (owned or not)
    //   --> just check if dimensions are compliant with the interaction
    // - no call to any set... --> error

    // Check if the variant is empty (i.e. monostate), if so error.
    if (std::holds_alternative<std::monostate>(H_storage_)) {
      // -> no previous call to setConstant or setCompute, this is not allowed
      THROW_EXCEPTION(
          "You must call setConstantH or setComputeHFunction "
          "before "
          "initialize.")
    }
    useH([&](const auto& M) {
      assert(M.rows() == sizeY);
      assert(M.cols() == sizeDS);
    });
  }

  // -- \frac{\partial}{\partial t}(J^h_q) --
  //  --> Optional
  // Case 1: user-defined function --> owned storage required and not set with
  // a call to setCompute...function
  if (computeHdot_) {
    if (!std::holds_alternative<siconos::algebra::OwnedSparse>(Hdot_storage_)) {
      // not owned storage? --> allocate
      Hdot_storage_ = std::make_unique<siconos::algebra::SiconosSparseMatrix>(sizeY, sizeDS);
    } else {
      // storage done and owned --> check dimensions
      bool are_dim_complient = false;
      useHdot([&are_dim_complient, sizeY, sizeDS](const auto& M) {
        are_dim_complient = M.rows() != sizeY || M.cols() != sizeDS;
      });
      if (!are_dim_complient)
        Hdot_storage_ = std::make_unique<siconos::algebra::SiconosSparseMatrix>(sizeY, sizeDS);
    }
  } else {
    // Other cases:
    // - previous call to setConstantH with alias or copy
    //   --> memory allocated in both cases (owned or not)
    //   --> just check if dimensions are compliant with the interaction

    // Check if the variant is not empty (i.e. monostate), if so, check dimensions
    if (!std::holds_alternative<std::monostate>(Hdot_storage_)) {
      // -> previous call to setConstant or setCompute
      useHdot([&](const auto& M) {
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

void siconos::modeling::LagrangianSparseScleronomousR::setConstantH(
    Eigen::Map<siconos::algebra::SiconosSparseMatrix>& newValue, siconos::algebra::AliasTag) {
  H_storage_ = std::make_shared<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>(newValue);

  // Warning: we can't check dimensions. This will be done during initialize.
  hasH_ = true;
  hasConstantH_ = true;
  computeH_ = nullptr;
  computeH_python_ = nullptr;
}

void siconos::modeling::LagrangianSparseScleronomousR::computeH(
    const siconos::algebra::BlockVector& positions) {
  // Warning: the jacobian must be allocated and of the right size before a call to this
  // function
  siconos::algebra::computeSparse(H_storage_, computeH_, computeH_python_, "H_storage_",
                                  positions);
  // Return value == false => means jacobian is constant, functions not set
}

void siconos::modeling::LagrangianSparseScleronomousR::setConstantHdot(
    Eigen::Map<siconos::algebra::SiconosSparseMatrix>& newValue, siconos::algebra::AliasTag) {
  Hdot_storage_ =
      std::make_shared<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>(newValue);

  // Warning: we can't check dimensions. This will be done during initialize.
  hasHdot_ = true;
  hasConstantHdot_ = true;
  computeHdot_ = nullptr;
  computeHdot_python_ = nullptr;
}

void siconos::modeling::LagrangianSparseScleronomousR::computeHdot(
    const siconos::algebra::BlockVector& q, const siconos::algebra::BlockVector& qdot) {
  siconos::algebra::computeSparse(Hdot_storage_, computeHdot_, computeHdot_python_,
                                  "H_storage_", q, qdot);
  // Return value == false => means jacobian is constant, functions not set
}

void siconos::modeling::LagrangianSparseScleronomousR::computeHdotv(double time,
                                                                    Interaction& inter) {
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  DEBUG_PRINT("siconos::modeling::LagrangianSparseScleronomousR::computedotjacqhXqdot starts");

  computeHdot(*ds_vars[tools::enum_to_index(ds_var::q0)],
              *ds_vars[tools::enum_to_index(ds_var::q1)]);
  if (!Hdotv_) {
    useHdot([&](const auto& M) {
      Hdotv_ = std::make_shared<siconos::algebra::SiconosVector>(M.rows());

      siconos::algebra::matrixBlockVector_prod(M, *ds_vars[tools::enum_to_index(ds_var::q1)],
                                               *Hdotv_);
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
  computeH(*ds_vars[tools::enum_to_index(ds_var::q0)]);
  computeHdot(*ds_vars[tools::enum_to_index(ds_var::q0)],
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
    computeH(*ds_vars[tools::enum_to_index(ds_var::q0)]);

    if (derivativeNumber == 1) {
      useH([&](const auto& M) {
        siconos::algebra::matrixBlockVector_prod(M, *ds_vars[tools::enum_to_index(ds_var::q1)],
                                                 y);
      });
    } else if (derivativeNumber == 2) {
      useH([&](const auto& M) {
        siconos::algebra::matrixBlockVector_prod(M, *ds_vars[tools::enum_to_index(ds_var::q2)],
                                                 y);
      });

      computeHdot(*ds_vars[tools::enum_to_index(ds_var::q0)],
                  *ds_vars[tools::enum_to_index(ds_var::q1)]);
      useHdot([&](const auto& M) {
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
  computeH(*ds_vars[tools::enum_to_index(ds_var::q0)]);
  // get lambda of the concerned interaction
  auto& lambda = *inter.lambda(level);
  DEBUG_EXPR(siconos::algebra::print(lambda););
  // data[name] += trans(G) * lambda

  useH([&](const auto& M) {
    siconos::algebra::transposeMatrixVector_prod_toBlock(
        lambda, M, *ds_vars[tools::enum_to_index(ds_var::p0) + level], false);
  });

  DEBUG_EXPR(siconos::algebra::print(*ds_vars[tools::enum_to_index(ds_var::p0) + level]););
  DEBUG_END(
      "void siconos::modeling::LagrangianSparseScleronomousR::computeInput(double time, "
      "Interaction& inter, InteractionProperties& interProp, "
      "siconos::algebra::blocks::size_type level) \n");
}
