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

#include "FiniteElement.hpp"

#include <iostream>

#include "FENode.hpp"
#include "FETypes.hpp"
#include "Mesh.hpp"  // MElement
#include "SiconosException.hpp"
#include "Tools.hpp"  // enum_to_string

siconos::mechanics::fem::FElement::FElement(std::shared_ptr<MElement> mesh_elem,
                                            const std::vector<std::shared_ptr<FENode>>& nodes)
    : num_{mesh_elem->num()}, type_{mesh_elem->type()}, mElement_{mesh_elem} {
  for (auto node : nodes) nodes_.push_back(node);
  ndof_ = nodes_.size() * number_of_dof_per_node(type_);

  switch (type_) {
    case FiniteElementType::T3:
    case FiniteElementType::B2:
      dimStress_ =
          3;  // T3: sigma_xx, sigma_yy, sigma_xy // B2: axial, curvature GP1, curvature GP2
      break;
    case FiniteElementType::TH4:
    case FiniteElementType::B3:
      dimStress_ = 6;  // TH4: sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz //
                       // B2: axial, curvature_x GP1, curvature_x GP2, twisting, ,
                       // curvature_y GP1, curvature_y GP2
      break;
    default:
      dimStress_ = 0;
  }

  if (is_valid_beam_element(type_)) {
    // Compute Te, length ...
    initialize_beam_element();
  }
}

std::span<const std::vector<double>> siconos::mechanics::fem::FElement::GaussPoints(
    int order) const {
  FiniteElementType type = mElement_->type();
  switch (type) {
    case FiniteElementType::T3:
      if (order == 1)
        return GaussPointsT3_1;
      else if (order == 2)
        return GaussPointsT3_2;
      break;
    case FiniteElementType::TH4:
      if (order == 1)
        return GaussPointsTH4_1;
      else if (order == 2)
        return GaussPointsTH4_2;
      break;
    case FiniteElementType::B2:
    case FiniteElementType::B3:
      assert(order == 3);
      return GaussPointsB2_3;
      break;

    default:
      throw("FElement::GaussPoints(). element type not recognized");
  }
  return GaussPointsEmpty;
}

// Note FP : todo, template this class with number of nodes in the element ?
void siconos::mechanics::fem::FElement::shapeFunctionIso2D(double ksi, double eta,
                                                           std::vector<double>& N,
                                                           std::vector<double>& Nksi,
                                                           std::vector<double>& Neta) const {
  FiniteElementType type = mElement_->type();
  switch (type) {
    case FiniteElementType::T3: {
      N[0] = 1.0 - ksi - eta;
      N[1] = ksi;
      N[2] = eta;
      Nksi[0] = -1.0;
      Nksi[1] = 1.0;
      Nksi[2] = 0.0;
      Neta[0] = -1.0;
      Neta[1] = 0.0;
      Neta[2] = 1.0;
      break;
    }
    default:
      THROW_EXCEPTION("FElement::shapeFunctionIso2D(). element type not recognized");
  }
}

void siconos::mechanics::fem::FElement::shapeFunctionIso3D(double ksi, double eta, double zeta,
                                                           std::vector<double>& N,
                                                           std::vector<double>& Nksi,
                                                           std::vector<double>& Neta,
                                                           std::vector<double>& Nzeta) const {
  FiniteElementType type = mElement_->type();
  switch (type) {
    case FiniteElementType::TH4: {
      N[0] = 1.0 - ksi - eta - zeta;
      N[1] = ksi;
      N[2] = eta;
      N[3] = zeta;

      Nksi[0] = -1.0;
      Nksi[1] = 1.0;
      Nksi[2] = 0.0;
      Nksi[3] = 0.0;

      Neta[0] = -1.0;
      Neta[1] = 0.0;
      Neta[2] = 1.0;
      Neta[3] = 0.0;

      Nzeta[0] = -1.0;
      Nzeta[1] = 0.0;
      Nzeta[2] = 0.0;
      Nzeta[3] = 1.0;

      break;
    }
    default:
      throw("FElement::shapeFunctionIso3D(). element type not recognized");
  }
}

void siconos::mechanics::fem::FElement::display() {
  std::cout << " - FElement - number: " << mElement_->num()
            << "            - type: " << siconos::tools::enum_to_string(mElement_->type())
            << "            - ndof: " << ndof_
            << "            - number of nodes: " << nodes_.size();
  std::cout << std::endl;
  for (std::shared_ptr<FENode> n : nodes_) {
    n->display();
  }
};

void siconos::mechanics::fem::FElement::initialize_beam_element() {
  assert(is_valid_beam_element(type_));

  length_ =
      sqrt(pow(nodes_[0]->x() - nodes_[1]->x(), 2) + pow(nodes_[0]->y() - nodes_[1]->y(), 2) +
           pow(nodes_[0]->z() - nodes_[1]->z(), 2));

  Te_.resize(ndof_, ndof_);
  Te_.setZero();
  if (type_ == FiniteElementType::B2) {
    double c = (nodes_[1]->x() - nodes_[0]->x()) / length();
    double s = (nodes_[1]->y() - nodes_[0]->y()) / length();
    Te_(0, 0) = c;
    Te_(0, 1) = s;
    Te_(1, 0) = -s;
    Te_(1, 1) = c;
    Te_(2, 2) = 1;
    Te_(3, 3) = c;
    Te_(3, 4) = s;
    Te_(4, 3) = -s;
    Te_(4, 4) = c;
    Te_(5, 5) = 1;
  } else if (type_ == FiniteElementType::B3) {
    siconos::algebra::SiconosMatrix33 vx;
    vx.setZero();
    double lx = (nodes_[1]->x() - nodes_[0]->x()) / length();
    double ly = (nodes_[1]->y() - nodes_[0]->y()) / length();
    double lz = (nodes_[1]->z() - nodes_[0]->z()) / length();
    double s = sqrt(lz * lz + ly * ly);
    vx(0, 1) = ly;
    vx(0, 2) = lz;
    vx(1, 0) = -ly;
    vx(2, 0) = -lz;

    double coeff = 1.0 / (1.0 + lx);
    siconos::algebra::SiconosMatrix33 R = siconos::algebra::identity33 + vx + coeff * vx * vx;

    for (int a = 0; a < 2; a++) {
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          Te_(6 * a + i, 6 * a + j) = R(i, j);
        }
      }
    }
    for (int a = 0; a < 2; a++) {
      for (int i = 0; i < 3; i++) {
        Te_(3 + 6 * a + i, 3 + 6 * a + i) = 1.;
      }
    }
  }

  Te_transpose_ = Te_.transpose();
}