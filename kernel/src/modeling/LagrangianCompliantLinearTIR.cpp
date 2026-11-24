// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "LagrangianCompliantLinearTIR.hpp"

#include <iostream>

// #include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosAlgebraAddons.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "Tools.hpp"

// Minimum data (C as pointer) constructor
siconos::modeling::LagrangianCompliantLinearTIR::LagrangianCompliantLinearTIR(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newC,
    Eigen::Ref<siconos::algebra::SiconosMatrix> newD)
    : LagrangianR{RelationSubType::CompliantLinearTIR} {
  contactMatrix_view_ =  // C
      std::make_shared<siconos::algebra::MapType>(newC.data(), newC.rows(), newC.cols());

  D_view_ = std::make_shared<siconos::algebra::MapType>(newD.data(), newD.rows(), newD.cols());
}

siconos::modeling::LagrangianCompliantLinearTIR::LagrangianCompliantLinearTIR(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newC,
    Eigen::Ref<siconos::algebra::SiconosMatrix> newD,
    Eigen::Ref<siconos::algebra::SiconosVector> newe)
    : LagrangianR{RelationSubType::CompliantLinearTIR} {
  contactMatrix_view_ =  // C
      std::make_shared<siconos::algebra::MapType>(newC.data(), newC.rows(), newC.cols());

  D_view_ = std::make_shared<siconos::algebra::MapType>(newD.data(), newD.rows(), newD.cols());

  eVector_view_ = std::make_shared<siconos::algebra::MapVectorType>(newe.data(), newe.size());
}

void siconos::modeling::LagrangianCompliantLinearTIR::initialize(Interaction& inter) {
  checkSize(inter);
}

void siconos::modeling::LagrangianCompliantLinearTIR::checkSize(
    const Interaction& inter) const {
  auto sizeY = inter.dimension();

  if (!(contactMatrix_view_) || contactMatrix_view_->cols() != inter.getSizeOfDS() ||
      contactMatrix_view_->rows() != sizeY)
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianCompliantLinearTIR::checkSize inconsistent sizes "
        "between H matrix and the interaction.");

  if ((D_view_) && (D_view_->rows() != sizeY || D_view_->cols() != sizeY))
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianCompliantLinearTIR::checkSize inconsistent sizes "
        "between D matrix and the interaction.");

  if ((eVector_view_) && eVector_view_->size() != sizeY)
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianCompliantLinearTIR::checkSize inconsistent sizes "
        "between e vector and the dimension of the interaction.");
}

void siconos::modeling::LagrangianCompliantLinearTIR::computeInput(
    double time, Interaction& inter, siconos::algebra::blocks::size_type level) {
  // get lambda of the concerned interaction
  auto& lambda = *inter.lambda(level);
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  // computation of p = Ht lambda
  siconos::algebra::transposeMatrixVector_prod_toBlock(
      lambda, *contactMatrix_view_, *ds_vars[tools::enum_to_index(ds_var::p0) + level], false);
}
void siconos::modeling::LagrangianCompliantLinearTIR::computeOutput(
    double time, Interaction& inter, siconos::algebra::blocks::size_type derivativeNumber) {
  // get y and lambda of the interaction
  auto& y = *inter.y(derivativeNumber);
  auto& lambda = *inter.lambda(derivativeNumber);
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  siconos::algebra::matrixBlockVector_prod(
      *contactMatrix_view_, *ds_vars[tools::enum_to_index(ds_var::q0) + derivativeNumber], y);
  y += *D_view_ * lambda;
  if (derivativeNumber == 0) {
    if (eVector_view_) y += *eVector_view_;
  }
}

void siconos::modeling::LagrangianCompliantLinearTIR::display() const {
  LagrangianR::display();
  std::cout << "===== Lagrangian Linear Relation display ===== " << std::endl;
  std::cout << " C: \n" << *contactMatrix_view_ << "\n";
  std::cout << " e: " << std::endl;
  if (eVector_view_)
    std::cout << *eVector_view_ << "\n";
  else
    std::cout << " -> nullptr " << std::endl;
  std::cout << " D: \n" << *D_view_ << "\n";
  std::cout << "===================================== " << std::endl;
}
