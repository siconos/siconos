// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "SolidLinearTIDS.hpp"

#include <stdio.h>

#include <cmath>
#include <iostream>

#include "FENode.hpp"
#include "FiniteElementLinearTIDS.hpp"
#include "FiniteElementModel.hpp"
#include "Material.hpp"
#include "Mesh.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES

#include "siconos_debug.h"

siconos::mechanics::fem::SolidLinearTIDS::SolidLinearTIDS(
    std::shared_ptr<Mesh> mesh, const std::map<int, const Material>& materials,
    MassType massType)
    : FiniteElementLinearTIDS(mesh, materials, massType) {
  DEBUG_BEGIN("SolidLinearTIDS::SolidLinearTIDS(...)\n");

  auto nElements = FEModel_->elements().size();
  auto dim = FEModel_->mesh()->dim();

  // Computes the size of the stress tensor
  size_t dimStressInElem = (dim == 2) ? 3 : 6;  // Works also for Beam cases
  dimStress_ = siconos::algebra::to_index(dimStressInElem * nElements);

  x_size_ = ndof_ + dimStress_;
  // aka the true unknowns are v and sigma
  //    _q0.reset(new siconos::algebra::SiconosVector(ndof_,0.0));
  //    _velocity0.reset(new siconos::algebra::SiconosVector(ndof_,0.0));

  sigma0_ = std::make_unique<siconos::algebra::SiconosVector>(dimStress_);
  sigma0_->setZero();
  epsilon_p0_ = std::make_unique<siconos::algebra::SiconosVector>(dimStress_);
  epsilon_p0_->setZero();

  // -- Memory allocation for stress --
  //    _sigma = std::make_shared<siconos::algebra::SiconosVector>(*sigma0_);
  stress_[0] = std::make_shared<siconos::algebra::SiconosVector>(*sigma0_);
  stress_[1] = std::make_shared<siconos::algebra::SiconosVector>(*sigma0_);
  stress_[2] = std::make_shared<siconos::algebra::SiconosVector>(*sigma0_);
  epsilon_p_[0] = std::make_shared<siconos::algebra::SiconosVector>(*epsilon_p0_);
  epsilon_p_[1] = std::make_shared<siconos::algebra::SiconosVector>(*epsilon_p0_);

  // -- Initialize the elasticity matrix --
  elasticity_matrix_storage_ =
      std::make_unique<siconos::algebra::SiconosSparseMatrix>(dimStress_, dimStress_);
  siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
      elasticity_matrix_storage_,
      [&](auto& matrix) { FEModel_->computeElasticityMatrix(matrix, materials_); },
      "elasticity_matrix_storage_");

  // -- Initialize the B matrix --
  B_matrix_storage_ =
      std::make_unique<siconos::algebra::SiconosSparseMatrix>(dimStress_, ndof_);
  siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
      B_matrix_storage_, [&](auto& matrix) { FEModel_->computeBMatrix(matrix, materials_); },
      "B_matrix_storage_");

  DEBUG_END(
      "SolidLinearTIDS::SolidLinearTIDS(std::shared_ptr<Mesh> mesh, std::shared_ptr<Material> "
      "material\n");
}

// --- Functions for memory handling ---
void siconos::mechanics::fem::SolidLinearTIDS::initMemory(
    siconos::algebra::blocks::size_type steps) {
  LagrangianSparseDS::initMemory(steps);
  stressMemory_.setMemorySize(steps, dimStress_);
  plasticDeformationMemory_.setMemorySize(steps, dimStress_);
  // plasticRate_Memory_.setMemorySize(steps, dimStress_);
}

void siconos::mechanics::fem::SolidLinearTIDS::swapInMemory() {
  LagrangianSparseDS::swapInMemory();
  stressMemory_.swap(*stress_[0]);
  // stressRateMemory_.swap(*stress_[1]);
  plasticDeformationMemory_.swap(*epsilon_p_[0]);
  // plasticRate_Memory_.swap(*epsilon_p_[1]);
}

void siconos::mechanics::fem::SolidLinearTIDS::resetAllNonSmoothParts() {
  FiniteElementLinearTIDS::resetAllNonSmoothParts();
  if (epsilon_p_[1]) epsilon_p_[1]->setZero();
}

void siconos::mechanics::fem::SolidLinearTIDS::resetNonSmoothPart(
    siconos::algebra::blocks::size_type level) {
  FiniteElementLinearTIDS::resetNonSmoothPart(level);
  if (level == 1 && epsilon_p_[1]) epsilon_p_[1]->setZero();
}

void siconos::mechanics::fem::SolidLinearTIDS::display(bool brief) const {
  std::cout << "===== SolidLinearTIDS display ===== " << std::endl;
  FiniteElementLinearTIDS::display();
}

double siconos::mechanics::fem::SolidLinearTIDS::elementStressScaling(size_t element) const {
  const auto& elements = FEModel_->elements();
  if (element >= elements.size())
    THROW_EXCEPTION("SolidLinearTIDS::elementStressScaling: element index out of range");
  auto& fe = *elements[element];
  if (fe.type() != FiniteElementType::T3)
    THROW_EXCEPTION("SolidLinearTIDS::elementStressScaling: only implemented for T3 elements");
  const auto& nd = fe.nodes();
  double twoA = (nd[1]->x() - nd[0]->x()) * (nd[2]->y() - nd[0]->y()) -
                (nd[2]->x() - nd[0]->x()) * (nd[1]->y() - nd[0]->y());
  double thickness = materials_.at(fe.mElement()->tags(0)).thickness();
  return std::sqrt(0.5 * twoA * thickness);
}

siconos::algebra::SiconosMatrix
siconos::mechanics::fem::SolidLinearTIDS::elementStressSelection(size_t element) const {
  double s = elementStressScaling(element);
  siconos::algebra::SiconosMatrix C(3, stressDimension());
  C.setZero();
  auto e = static_cast<siconos::algebra::Index>(3 * element);
  for (siconos::algebra::Index i = 0; i < 3; i++) C(i, e + i) = 1.0 / s;
  return C;
}
