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
#include "FiniteElementLinearTIDS.hpp"

#include "BoundaryCondition.hpp"
#include "FENode.hpp"
#include "FiniteElementModel.hpp"
#include "Material.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
// #include "siconos_debug.h"

siconos::mechanics::fem::FiniteElementLinearTIDS::FiniteElementLinearTIDS(
    std::shared_ptr<Mesh> mesh, const std::map<int, const Material>& materials)
    : LagrangianSparseLinearTIDS(), mesh_(mesh), materials_(materials) {
  // Warning FP: the DS is built from default empty constructor in an unusual
  // way. Care must be taken to properly set all attributes in DS, SecondOrder,
  // Lagrangian ... It may be better to :
  // - build FEModel from mesh
  // - compute ndof
  // - build ds from ndof or q0, v0

  FEModel_ = std::make_shared<FiniteElementModel>(mesh);
  ndof_ = FEModel_->init();

  q0_storage_ = std::make_unique<siconos::algebra::SiconosVector>(ndof_);
  velocity0_storage_ = std::make_unique<siconos::algebra::SiconosVector>(ndof_);
  use_q0([&](auto& v) { v.setZero(); });
  use_velocity0([&](auto& v) { v.setZero(); });

  // -- Memory allocation for vector and matrix members --
  state_q_[0] = std::make_shared<siconos::algebra::SiconosVector>(q0());
  state_q_[1] = std::make_shared<siconos::algebra::SiconosVector>(velocity0());

  p_[1] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  p_[1]->setZero();

  x_size_ = 2 * ndof_;

  // Mass ...
  // Deal with 'plugged' mass later
  mass_storage_ = std::make_unique<siconos::algebra::SiconosSparseMatrix>(ndof_, ndof_);
  hasConstantMass_ = true;
  hasMass_ = true;
  computemass_ = nullptr;
  hasLUMass_ = false;
  siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
      mass_storage_, [&](auto& matrix) { FEModel_->computeMassMatrix(matrix, materials_); },
      "mass_storage_");

  stiffnessMatrix_storage =
      std::make_unique<siconos::algebra::SiconosSparseMatrix>(ndof_, ndof_);
  siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
      stiffnessMatrix_storage,
      [&](auto& matrix) { FEModel_->computeStiffnessMatrix(matrix, materials_); },
      "stiffnessMatrix_storage");

  // dampingMatrix_ = std::make_shared<siconos::algebra::SiconosSparseMatrix>(ndof_, ndof_);
}

void siconos::mechanics::fem::FiniteElementLinearTIDS::applyDirichletBoundaryConditions(
    int physical_entity_tag, const std::vector<int>& node_dof_index, double imposedVelocity) {
  if (!boundaryConditions_)
    boundaryConditions_ = std::make_shared<siconos::modeling::BoundaryCondition>(
        siconos::modeling::BoundaryCondition::Indices{});

  FEModel_->applyDirichletBoundaryConditions(physical_entity_tag, node_dof_index,
                                             boundaryConditions_, imposedVelocity);

  reactionToBoundaryConditions_ = std::make_shared<siconos::algebra::SiconosVector>(
      boundaryConditions_->velocityIndices().size());
};

void siconos::mechanics::fem::FiniteElementLinearTIDS::applyNodalForces(
    int physical_entity_tag, const siconos::algebra::SiconosVector& nodal_forces) {
  if (!std::holds_alternative<siconos::algebra::OwnedDenseVector>(fext_storage_)) {
    fext_storage_ = std::make_unique<siconos::algebra::SiconosVector>(ndof_);
  }
  hasFext_ = true;
  hasConstantFext_ = true;  // set later by applyNodalForces
  //  We should probably use something like setComputeFext(applyNodalForce)
  computefext_ = nullptr;

  use_fext([&](auto& fext) {
    fext.setZero();
    FEModel_->applyNodalForces(physical_entity_tag, nodal_forces, fext);
  });
};

double siconos::mechanics::fem::FiniteElementLinearTIDS::elasticPotentialEnergy() const {
  siconos::algebra::SiconosVector tmp{ndof_};
  useStiffness([&](const auto& K) { tmp = K * *state_q_[0]; });
  return 0.5 * state_q_[0]->dot(tmp);
}

void siconos::mechanics::fem::FiniteElementLinearTIDS::display(bool brief) const {
  std::cout << "===== FiniteElementLinearTIDS display ===== " << std::endl;
  LagrangianSparseLinearTIDS::display();
  FEModel_->display(brief);
}

std::vector<double> siconos::mechanics::fem::FiniteElementLinearTIDS::computeStrainTensor() const {
  std::vector<double> epsilon;
  auto femodel = FEModel_;
  if (!femodel) return epsilon;

  auto q_ptr = this->q();
  const auto& q0_vec = this->q0();
  if (!q_ptr) return epsilon;
  const auto& q_vec = *q_ptr;
  siconos::algebra::SiconosVector u_vec = q_vec - q0_vec;

  // Get material properties from the first material
  double E = 0.0;
  double nu = 0.0;
  if (!materials_.empty()) {
    E = materials_.begin()->second.elasticYoungModulus();
    nu = materials_.begin()->second.Poisson_s_ratio();
  }

  // D matrix for plane stress (used only for stress, not strain)
  // Strain is computed directly from B * u

  for (auto& elem : femodel->elements()) {
    auto nodes = elem->nodes();
    if (nodes.size() != 3) continue;  // T3 only for now

    // Get global DOF indices for this element
    std::vector<siconos::algebra::Index> dofs;
    for (auto& node : nodes) {
      auto node_dofs = node->global_dof_index();
      if (node_dofs.size() >= 2) {
        dofs.push_back(node_dofs[0]);
        dofs.push_back(node_dofs[1]);
      }
    }
    if (dofs.size() != 6) continue;

    // Assemble element displacement vector
    std::vector<double> u_element(6);
    for (int i = 0; i < 6; i++) {
      u_element[i] = u_vec(dofs[i]);
    }

    // Compute B matrix
    siconos::algebra::SiconosDenseMatrix Be(3, 6);
    FEModel_->computeElementaryBMatrix_direct(*elem, Be, 1.0);  // thickness = 1.0

    // Compute strain = B * u
    // Be is stored as Eigen matrix internally
    Eigen::Map<Eigen::MatrixXd> Be_map(Be.data(), Be.rows(), Be.cols());
    Eigen::Map<Eigen::VectorXd> u_vec_elem(u_element.data(), 6);
    Eigen::VectorXd eps_vec = Be_map * u_vec_elem;

    epsilon.push_back(eps_vec[0]);  // exx
    epsilon.push_back(eps_vec[1]);  // eyy
    epsilon.push_back(eps_vec[2]);  // exy
  }

  return epsilon;
}

std::vector<double> siconos::mechanics::fem::FiniteElementLinearTIDS::computeStressTensor() const {
  std::vector<double> sigma;
  auto femodel = FEModel_;
  if (!femodel) return sigma;

  auto q_ptr = this->q();
  const auto& q0_vec = this->q0();
  if (!q_ptr) return sigma;
  const auto& q_vec = *q_ptr;
  siconos::algebra::SiconosVector u_vec = q_vec - q0_vec;

  // Get material properties from the first material
  double E = 0.0;
  double nu = 0.0;
  if (!materials_.empty()) {
    E = materials_.begin()->second.elasticYoungModulus();
    nu = materials_.begin()->second.Poisson_s_ratio();
  }

  // D matrix for plane stress
  double D11 = E / (1.0 - nu * nu);
  double D12 = E * nu / (1.0 - nu * nu);
  double D33 = E / (2.0 * (1.0 + nu));

  for (auto& elem : femodel->elements()) {
    auto nodes = elem->nodes();
    if (nodes.size() != 3) continue;  // T3 only for now

    // Get global DOF indices for this element
    std::vector<siconos::algebra::Index> dofs;
    for (auto& node : nodes) {
      auto node_dofs = node->global_dof_index();
      if (node_dofs.size() >= 2) {
        dofs.push_back(node_dofs[0]);
        dofs.push_back(node_dofs[1]);
      }
    }
    if (dofs.size() != 6) continue;

    // Assemble element displacement vector
    std::vector<double> u_element(6);
    for (int i = 0; i < 6; i++) {
      u_element[i] = q_vec(dofs[i]);
    }

    // Compute B matrix
    siconos::algebra::SiconosDenseMatrix Be(3, 6);
    FEModel_->computeElementaryBMatrix_direct(*elem, Be, 1.0);  // thickness = 1.0

    // Compute strain = B * u
    Eigen::Map<Eigen::MatrixXd> Be_map(Be.data(), Be.rows(), Be.cols());
    Eigen::Map<Eigen::VectorXd> u_vec(u_element.data(), 6);
    Eigen::VectorXd eps_vec = Be_map * u_vec;

    // Compute stress = D * strain
    Eigen::Matrix3d D_mat;
    D_mat << D11, D12, 0,
             D12, D11, 0,
             0,   0,   D33;
    Eigen::Vector3d sigma_vec = D_mat * eps_vec;

    sigma.push_back(sigma_vec[0]);  // sxx
    sigma.push_back(sigma_vec[1]);  // syy
    sigma.push_back(sigma_vec[2]);  // sxy
  }

  return sigma;
}

std::vector<double> siconos::mechanics::fem::FiniteElementLinearTIDS::computeStrainTensor(
    const siconos::algebra::SiconosVector& displacement) const {
  std::vector<double> epsilon;
  auto femodel = FEModel_;
  if (!femodel) return epsilon;

  for (auto& elem : femodel->elements()) {
    auto nodes = elem->nodes();
    if (nodes.size() != 3) continue;  // T3 only for now

    // Get global DOF indices for this element
    std::vector<siconos::algebra::Index> dofs;
    for (auto& node : nodes) {
      auto node_dofs = node->global_dof_index();
      if (node_dofs.size() >= 2) {
        dofs.push_back(node_dofs[0]);
        dofs.push_back(node_dofs[1]);
      }
    }
    if (dofs.size() != 6) continue;

    // Assemble element displacement vector
    std::vector<double> u_element(6);
    for (int i = 0; i < 6; i++) {
      u_element[i] = displacement(dofs[i]);
    }

    // Compute B matrix
    siconos::algebra::SiconosDenseMatrix Be(3, 6);
    FEModel_->computeElementaryBMatrix_direct(*elem, Be, 1.0);  // thickness = 1.0

    // Compute strain = B * u
    Eigen::Map<Eigen::MatrixXd> Be_map(Be.data(), Be.rows(), Be.cols());
    Eigen::Map<Eigen::VectorXd> u_vec(u_element.data(), 6);
    Eigen::VectorXd eps_vec = Be_map * u_vec;

    epsilon.push_back(eps_vec[0]);  // exx
    epsilon.push_back(eps_vec[1]);  // eyy
    epsilon.push_back(eps_vec[2]);  // exy
  }

  return epsilon;
}

std::vector<double> siconos::mechanics::fem::FiniteElementLinearTIDS::computeStressTensor(
    const siconos::algebra::SiconosVector& displacement) const {
  std::vector<double> sigma;
  auto femodel = FEModel_;
  if (!femodel) return sigma;

  // Get material properties from the first material
  double E = 0.0;
  double nu = 0.0;
  if (!materials_.empty()) {
    E = materials_.begin()->second.elasticYoungModulus();
    nu = materials_.begin()->second.Poisson_s_ratio();
  }

  // D matrix for plane stress
  double D11 = E / (1.0 - nu * nu);
  double D12 = E * nu / (1.0 - nu * nu);
  double D33 = E / (2.0 * (1.0 + nu));

  for (auto& elem : femodel->elements()) {
    auto nodes = elem->nodes();
    if (nodes.size() != 3) continue;  // T3 only for now

    // Get global DOF indices for this element
    std::vector<siconos::algebra::Index> dofs;
    for (auto& node : nodes) {
      auto node_dofs = node->global_dof_index();
      if (node_dofs.size() >= 2) {
        dofs.push_back(node_dofs[0]);
        dofs.push_back(node_dofs[1]);
      }
    }
    if (dofs.size() != 6) continue;

    // Assemble element displacement vector
    std::vector<double> u_element(6);
    for (int i = 0; i < 6; i++) {
      u_element[i] = displacement(dofs[i]);
    }

    // Compute B matrix
    siconos::algebra::SiconosDenseMatrix Be(3, 6);
    FEModel_->computeElementaryBMatrix_direct(*elem, Be, 1.0);  // thickness = 1.0

    // Compute strain = B * u
    Eigen::Map<Eigen::MatrixXd> Be_map(Be.data(), Be.rows(), Be.cols());
    Eigen::Map<Eigen::VectorXd> u_vec(u_element.data(), 6);
    Eigen::VectorXd eps_vec = Be_map * u_vec;

    // Compute stress = D * strain
    Eigen::Matrix3d D_mat;
    D_mat << D11, D12, 0,
             D12, D11, 0,
             0,   0,   D33;
    Eigen::Vector3d sigma_vec = D_mat * eps_vec;

    sigma.push_back(sigma_vec[0]);  // sxx
    sigma.push_back(sigma_vec[1]);  // syy
    sigma.push_back(sigma_vec[2]);  // sxy
  }

  return sigma;
}
