// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

// \todo : create a work vector for all tmp vectors used in computeg, computeh ...

#include "LagrangianCompliantR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosAlgebraAddons.hpp"  // for matrix-vector prod
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "Tools.hpp"

void siconos::modeling::LagrangianCompliantR::initialize(Interaction& inter) {
  auto sizeY = inter.dimension();
  auto sizeDS = inter.getSizeOfDS();

  if (!contactMatrix_storage_) {
    contactMatrix_storage_ = std::make_unique<std::vector<double>>(sizeY * sizeDS);
  }
  contactMatrix_view_ = std::make_shared<siconos::algebra::MapType>(
      contactMatrix_storage_->data(), sizeY, sizeDS);

  if (!D_) D_ = std::make_shared<siconos::algebra::SiconosMatrix>(sizeY, sizeY);
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

void siconos::modeling::LagrangianCompliantR::setComputeContactMatrixFunction(
    const siconos::modeling::func_prototypes::FunctionBVV_M& func) {
  computeContactMatrix_ = func;
}

void siconos::modeling::LagrangianCompliantR::computeContactMatrix(
    const siconos::algebra::BlockVector& positions,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda) {
  if (!computeContactMatrixLegacy_) {  // dispatch to a class overriding the deprecated name
    detail::LegacyCallGuard guard(computeContactMatrixLegacy_);
    SICONOS_IGNORE_DEPRECATED_BEGIN
    computeJacobianhOver_q(positions, lambda);
    SICONOS_IGNORE_DEPRECATED_END
    return;
  }
  if (computeContactMatrix_) computeContactMatrix_(positions, lambda, *contactMatrix_view_);
}

void siconos::modeling::LagrangianCompliantR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector& positions,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda) {
  computeContactMatrix(positions, lambda);
}

void siconos::modeling::LagrangianCompliantR::setComputeDFunction(
    const siconos::modeling::func_prototypes::FunctionBVV_M& func) {
  computeD_ = func;
}

void siconos::modeling::LagrangianCompliantR::computeD(
    const siconos::algebra::BlockVector& positions,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda) {
  if (!computeDLegacy_) {  // dispatch to a class overriding the deprecated name
    detail::LegacyCallGuard guard(computeDLegacy_);
    SICONOS_IGNORE_DEPRECATED_BEGIN
    computeJacobianhOver_lambda(positions, lambda);
    SICONOS_IGNORE_DEPRECATED_END
    return;
  }
  if (computeD_) computeD_(positions, lambda, *D_);
}

void siconos::modeling::LagrangianCompliantR::computeJacobianhOver_lambda(
    const siconos::algebra::BlockVector& positions,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda) {
  computeD(positions, lambda);
}

void siconos::modeling::LagrangianCompliantR::computeOutput(
    double time, Interaction& inter, siconos::algebra::blocks::size_type derivativeNumber) {
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  if (derivativeNumber == 0) {
    auto& y = *inter.y(0);
    auto& lambda = *inter.lambda(0);
    computeh(*ds_vars[tools::enum_to_index(ds_var::q0)], lambda, y);
  } else {
    auto& y = *inter.y(derivativeNumber);
    auto& lambda = *inter.lambda(derivativeNumber);
    computeContactMatrix(*ds_vars[tools::enum_to_index(ds_var::q0)], lambda);
    computeD(*ds_vars[tools::enum_to_index(ds_var::q0)], lambda);
    if (derivativeNumber == 1) {
      // y = Jach[0] q1 + Jach[1] lambda
      siconos::algebra::matrixBlockVector_prod(*contactMatrix_view_,
                                               *ds_vars[tools::enum_to_index(ds_var::q1)], y);
      y += *D_ * lambda;
    } else if (derivativeNumber == 2)
      siconos::algebra::matrixBlockVector_prod(
          *contactMatrix_view_, *ds_vars[tools::enum_to_index(ds_var::q2)],
          y);  // Approx: y[2] = Jach[0]q[2], other terms are neglected ...
    else
      THROW_EXCEPTION(
          "siconos::modeling::LagrangianCompliantR::computeOutput, index out of range or not "
          "yet implemented.");
  }
}

void siconos::modeling::LagrangianCompliantR::computeInput(
    double time, Interaction& inter, siconos::algebra::blocks::size_type level) {
  // get lambda of the concerned interaction

  auto& lambda = *inter.lambda(level);
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  computeContactMatrix(*ds_vars[tools::enum_to_index(ds_var::q0)], lambda);
  // data[name] += trans(G) * lambda
  siconos::algebra::transposeMatrixVector_prod_toBlock(
      lambda, *contactMatrix_view_, *ds_vars[tools::enum_to_index(ds_var::p0) + level], false);
}

void siconos::modeling::LagrangianCompliantR::computeJach(double time, Interaction& inter) {
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  auto& lambda = *inter.lambda(0);
  computeContactMatrix(*ds_vars[tools::enum_to_index(ds_var::q0)], lambda);
  computeD(*ds_vars[tools::enum_to_index(ds_var::q0)], lambda);
}
