// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

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
  if (hasConstantContactMatrix_) {
    assert(contactMatrix_view_);
    assert(contactMatrix_view_->rows() == sizeY);
    assert(contactMatrix_view_->cols() == sizeDS);
  } else {
    if (!contactMatrix_storage_) {
      contactMatrix_storage_ = std::make_unique<std::vector<double>>(sizeY * sizeDS);
    }
    contactMatrix_view_ = std::make_shared<siconos::algebra::MapType>(
        contactMatrix_storage_->data(), sizeY, sizeDS);
  }

  if (hasContactMatrixDot_) {
    // True only if setComputeContactMatrixDotFunction has been called
    // Ensure that memory is properly allocated for contactMatrixDot_
    if (!contactMatrixDot_) {
      contactMatrixDot_ = std::make_shared<siconos::algebra::SiconosMatrix>(sizeY, sizeDS);
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

void siconos::modeling::LagrangianScleronomousR::setConstantContactMatrix(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  contactMatrix_storage_ = nullptr;

  contactMatrix_view_ = std::make_shared<siconos::algebra::MapType>(
      newValue.data(), newValue.rows(), newValue.cols());
  hasConstantContactMatrix_ = true;
  computeContactMatrix_ = nullptr;
}

void siconos::modeling::LagrangianScleronomousR::setComputeContactMatrixFunction(
    const siconos::modeling::func_prototypes::FunctionBV_M& func) {
  computeContactMatrix_ = func;
}

void siconos::modeling::LagrangianScleronomousR::computeContactMatrix(
    const siconos::algebra::BlockVector& positions) {
  if (!computeContactMatrixLegacy_) {  // dispatch to a class overriding the deprecated name
    detail::LegacyCallGuard guard(computeContactMatrixLegacy_);
    SICONOS_IGNORE_DEPRECATED_BEGIN
    computeJacobianhOver_q(positions);
    SICONOS_IGNORE_DEPRECATED_END
    return;
  }
  if (computeContactMatrix_) computeContactMatrix_(positions, *contactMatrix_view_);
}

void siconos::modeling::LagrangianScleronomousR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector& positions) {
  computeContactMatrix(positions);
}

void siconos::modeling::LagrangianScleronomousR::setComputeContactMatrixDotFunction(
    const siconos::modeling::func_prototypes::FunctionBVBV_M& func) {
  hasContactMatrixDot_ = true;
  computeContactMatrixDot_ = func;
}

void siconos::modeling::LagrangianScleronomousR::computeContactMatrixDot(
    const siconos::algebra::BlockVector& q, const siconos::algebra::BlockVector& qdot) {
  if (computeContactMatrixDot_) computeContactMatrixDot_(q, qdot, *contactMatrixDot_);
}

void siconos::modeling::LagrangianScleronomousR::computeContactMatrixDotV(double time,
                                                                          Interaction& inter) {
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  DEBUG_PRINT("siconos::modeling::LagrangianScleronomousR::computedotjacqhXqdot starts");

  computeContactMatrixDot(*ds_vars[tools::enum_to_index(ds_var::q0)],
                          *ds_vars[tools::enum_to_index(ds_var::q1)]);
  if (!contactMatrixDotV_) {
    contactMatrixDotV_ =
        std::make_shared<siconos::algebra::SiconosVector>(contactMatrixDot_->rows());
  }
  siconos::algebra::matrixBlockVector_prod(
      *contactMatrixDot_, *ds_vars[tools::enum_to_index(ds_var::q1)], *contactMatrixDotV_);
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
    computeContactMatrix(*ds_vars[tools::enum_to_index(ds_var::q0)]);

    if (derivativeNumber == 1) {
      siconos::algebra::matrixBlockVector_prod(*contactMatrix_view_,
                                               *ds_vars[tools::enum_to_index(ds_var::q1)], y);
    } else if (derivativeNumber == 2) {
      siconos::algebra::matrixBlockVector_prod(*contactMatrix_view_,
                                               *ds_vars[tools::enum_to_index(ds_var::q2)], y);
      computeContactMatrixDot(*ds_vars[tools::enum_to_index(ds_var::q0)],
                              *ds_vars[tools::enum_to_index(ds_var::q1)]);
      siconos::algebra::matrixBlockVector_prod(
          *contactMatrixDot_, *ds_vars[tools::enum_to_index(ds_var::q1)], y, false);

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
  computeContactMatrix(*ds_vars[tools::enum_to_index(ds_var::q0)]);
  // get lambda of the concerned interaction
  auto& lambda = *inter.lambda(level);
  DEBUG_EXPR(siconos::algebra::print(lambda););
  // data[name] += trans(G) * lambda
  siconos::algebra::transposeMatrixVector_prod_toBlock(
      lambda, *contactMatrix_view_, *ds_vars[tools::enum_to_index(ds_var::p0) + level], false);
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
  computeContactMatrix(*ds_vars[tools::enum_to_index(ds_var::q0)]);
  computeContactMatrixDot(*ds_vars[tools::enum_to_index(ds_var::q0)],
                          *ds_vars[tools::enum_to_index(ds_var::q1)]);

  DEBUG_END(
      "void siconos::modeling::LagrangianScleronomousR::computeJach(double time, "
      "Interaction& "
      "inter) \n");
}
