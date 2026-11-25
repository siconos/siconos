
/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
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
  if (computeContactMatrix_) {
    if (!std::holds_alternative<siconos::algebra::OwnedSparse>(contactMatrix_storage_)) {
      // not owned storage? --> allocate
      contactMatrix_storage_ =
          std::make_unique<siconos::algebra::SiconosSparseMatrix>(sizeY, sizeDS);
    } else {
      // storage done and owned --> check dimensions
      bool are_dim_complient = false;
      useContactMatrix([&are_dim_complient, sizeY, sizeDS](const auto& M) {
        are_dim_complient = M.rows() != sizeY || M.cols() != sizeDS;
      });
      if (!are_dim_complient)
        contactMatrix_storage_ =
            std::make_unique<siconos::algebra::SiconosSparseMatrix>(sizeY, sizeDS);
    }
  } else {
    // Other cases:
    // - previous call to setConstantContactMatrix with alias or copy
    //   --> memory allocated in both cases (owned or not)
    //   --> just check if dimensions are compliant with the interaction
    // - no call to any set... --> error

    // Check if the variant is empty (i.e. monostate), if so error.
    if (std::holds_alternative<std::monostate>(contactMatrix_storage_)) {
      // -> no previous call to setConstant or setCompute, this is not allowed
      THROW_EXCEPTION(
          "You must call setConstantContactMatrix or setComputeContactMatrixFunction "
          "before "
          "initialize.")
    }
    useContactMatrix([&](const auto& M) {
      assert(M.rows() == sizeY);
      assert(M.cols() == sizeDS);
    });
  }

  if (computedhdt_) {  // True only if setComputedhdtFunction has been called
    // Ensure that memory is properly allocated for dhdt_
    if (!dhdt_) {
      dhdt_ = std::make_shared<siconos::algebra::SiconosVector>(sizeY);
    }
  }
}

void siconos::modeling::LagrangianSparseRheonomousR::setComputedhdtFunction(
    const siconos::modeling::func_prototypes::FunctionBVS_V& hdot_func) {
  // Warning: we can't allocate dhdt here, since we have no clue on the
  // required size
  computedhdt_ = hdot_func;
}

void siconos::modeling::LagrangianSparseRheonomousR::setConstantContactMatrix(
    Eigen::Map<siconos::algebra::SiconosSparseMatrix>& newValue, siconos::algebra::AliasTag) {
  contactMatrix_storage_ =
      std::make_shared<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>(newValue);

  // Warning: we can't check dimensions. This will be done during initialize.
  hasContactMatrix_ = true;
  hasConstantContactMatrix_ = true;
  computeContactMatrix_ = nullptr;
  computeContactMatrix_python_ = nullptr;
}

void siconos::modeling::LagrangianSparseRheonomousR::computehdot(
    const siconos::algebra::BlockVector& positions, double time) {
  computedhdt(positions, time);
}

void siconos::modeling::LagrangianSparseRheonomousR::computedhdt(
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

void siconos::modeling::LagrangianSparseRheonomousR::computeContactMatrix(
    const siconos::algebra::BlockVector& positions, double time) {
  // Warning: the jacobian must be allocated and of the right size before a call to this
  // function
  siconos::algebra::computeSparse(contactMatrix_storage_, computeContactMatrix_,
                                  computeContactMatrix_python_, "contactMatrix_storage_",
                                  positions, time);
  // Return value == false => means jacobian is constant, functions not set
}

void siconos::modeling::LagrangianSparseRheonomousR::computeOutput(
    double time, Interaction& inter, siconos::algebra::blocks::size_type derivativeNumber) {
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  auto& y = *inter.y(derivativeNumber);
  if (derivativeNumber == 0)
    computeh(*ds_vars[tools::enum_to_index(ds_var::q0)], time, y);
  else {
    computeContactMatrix(*ds_vars[tools::enum_to_index(ds_var::q0)], time);
    if (derivativeNumber == 1) {
      if (!dhdt_) {
        auto sizeY = inter.dimension();
        dhdt_ = std::make_shared<siconos::algebra::SiconosVector>(sizeY);
      }
      // Computation of the partial derivative w.r.t time of h(q,t)
      computedhdt(*ds_vars[tools::enum_to_index(ds_var::q0)], time);
      // Computation of the partial derivative w.r.t q of h(q,t) : \nabla_q h(q,t) \dot q
      useContactMatrix([&](const auto& M) {
        siconos::algebra::matrixBlockVector_prod(M, *ds_vars[tools::enum_to_index(ds_var::q1)],
                                                 y);
      });
      // Sum of the terms
      y += *dhdt_;
    } else if (derivativeNumber == 2) {
      useContactMatrix([&](const auto& M) {
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
  computeContactMatrix(*ds_vars[tools::enum_to_index(ds_var::q0)], time);
  // get lambda of the concerned interaction
  auto& lambda = *inter.lambda(level);
  // data[name] += trans(G) * lambda
  useContactMatrix([&](const auto& M) {
    siconos::algebra::transposeMatrixVector_prod_toBlock(
        lambda, M, *ds_vars[tools::enum_to_index(ds_var::p0) + level], false);
  });
}

void siconos::modeling::LagrangianSparseRheonomousR::computeJach(double time,
                                                                 Interaction& inter) {
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  computeContactMatrix(*ds_vars[tools::enum_to_index(ds_var::q0)], time);
  // computeJachqDot(time, inter);
  //    computeDotJachq(time, q, z);
  // computeJachlambda(time, inter);
  computedhdt(*ds_vars[tools::enum_to_index(ds_var::q0)], time);
}
