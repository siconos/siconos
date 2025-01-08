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

#include "Cable2d3DR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::fem::cable::Cable2d3DR::initialize(siconos::modeling::Interaction& inter) {
  auto qSize = inter.getSizeOfDS();

  if (!jacobianhOver_q_internal_storage_) {
    jacobianhOver_q_internal_storage_ = std::make_unique<std::vector<double>>(2 * qSize);
  }
  jacobianhOver_q_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobianhOver_q_internal_storage_->data(), 2, qSize);
}

void siconos::fem::cable::Cable2d3DR::computeh(const siconos::algebra::BlockVector& q,
                                               Eigen::Ref<siconos::algebra::SiconosVector> y) {
  DEBUG_BEGIN("Cable2d3DR::computeh(...)\n");

  // LagrangianScleronomousR::computeh(q, z, y);
  auto& position = *((q.getAllVect())[0]);
  _Pc1->setValue(0, position(_node_dof_index));
  _Pc1->setValue(1, position(_node_dof_index + 1));
  _Pc1->setValue(2, position(_node_dof_index + 2));
  y.setValue(0, distance());
  DEBUG_EXPR(y.display(););
  DEBUG_EXPR(display(););
  DEBUG_END("Cable2d3DR::computeh(...)\n")
}

void siconos::fem::cable::Cable2d3DR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector& q) {
  DEBUG_BEGIN(
      "Cable2d3DR::computeJacobianhOver_q(const siconos::algebra::BlockVector& q, "
      "siconos::algebra::BlockVector& z \n");

  double Nx = _Normal->getValue(0);
  double Ny = _Normal->getValue(1);
  double Nz = _Normal->getValue(2);

  double Tx = _Tangent->getValue(0);
  double Ty = _Tangent->getValue(1);
  double Tz = _Tangent->getValue(2);

  DEBUG_PRINTF("N_x = %4.2e,\t N_y = %4.2e,\t N_z= %4.2e\n", Nx, Ny, Nz);
  DEBUG_PRINTF("T_x = %4.2e,\t T_y = %4.2e,\t T_z= %4.2e\n", Tx, Ty, Tz);

  jacobianhOver_q_view_->setValue(0, _node_dof_index, Nx);
  jacobianhOver_q_view_->setValue(0, _node_dof_index + 1, Ny);
  jacobianhOver_q_view_->setValue(0, _node_dof_index + 2, Nz);

  jacobianhOver_q_view_->setValue(1, _node_dof_index, Tx);
  jacobianhOver_q_view_->setValue(1, _node_dof_index + 1, Ty);
  jacobianhOver_q_view_->setValue(1, _node_dof_index + 2, Tz);

  if (q.size() == 6) {
    DEBUG_PRINT("take into account second ds\n");
    THROW_EXCEPTION("Cable2d3DR is not implemented for cable/cable contact");
  }

  DEBUG_END(
      "Cable2d3DR::computeJacobianhOver_q(const siconos::algebra::BlockVector& q, "
      "siconos::algebra::BlockVector& z) \n");
}

double siconos::fem::cable::Cable2d3DR::distance() const {
  DEBUG_BEGIN("Cable2d3DR::distance(...)\n")
  siconos::algebra::SiconosVector dpc(*_Pc2 - *_Pc1);
  DEBUG_EXPR(_Pc1->display(););
  DEBUG_EXPR(_Pc2->display(););
  DEBUG_EXPR(dpc.display(););
  DEBUG_END("Cable2d3DR::distance(...)\n")
  return dpc.norm2() * (_Normal->dot(dpc) >= 0 ? -1 : 1);
}

void siconos::fem::cable::Cable2d3DR::updateContactPoint(
    std::shared_ptr<siconos::algebra::SiconosVector> pc1,
    std::shared_ptr<siconos::algebra::SiconosVector> pc2,
    std::shared_ptr<siconos::algebra::SiconosVector> normal,
    std::shared_ptr<siconos::algebra::SiconosVector> tangent) {
  _Pc1 = pc1;
  _Pc2 = pc2;
  _Normal = normal;
  _Tangent = tangent;
};

/** update the contact points from references
 */
void siconos::fem::cable::Cable2d3DR::updateContactPoint(
    siconos::algebra::SiconosVector& pc1, siconos::algebra::SiconosVector& pc2,
    siconos::algebra::SiconosVector& normal, siconos::algebra::SiconosVector& tangent) {
  *_Pc1 = pc1;
  *_Pc2 = pc2;
  *_Normal = normal;
  *_Tangent = tangent;
};

void siconos::fem::cable::Cable2d3DR::updateContactPoint(double pc1[3], double pc2[3],
                                                         double normal[3], double tangent[3]) {
  _Pc1->setValue(0, pc1[0]);
  _Pc1->setValue(1, pc1[1]);
  _Pc1->setValue(2, pc1[2]);
  _Pc2->setValue(0, pc2[0]);
  _Pc2->setValue(1, pc2[1]);
  _Pc2->setValue(2, pc2[2]);
  _Normal->setValue(0, normal[0]);
  _Normal->setValue(1, normal[1]);
  _Normal->setValue(2, normal[2]);
  _Tangent->setValue(0, tangent[0]);
  _Tangent->setValue(1, tangent[1]);
  _Tangent->setValue(2, normal[2]);
};

void siconos::fem::cable::Cable2d3DR::display() const {
  LagrangianR::display();

  std::cout << " _node_dof_index :" << _node_dof_index << "\n";

  std::cout << " _Pc1: \n";
  if (_Pc1)
    _Pc1->display();
  else
    std::cout << " nullptr\n";

  std::cout << " _Pc2 :\n";
  if (_Pc2)
    _Pc2->display();
  else
    std::cout << " nullptr\n";

  std::cout << " _Normal:\n";
  if (_Normal)
    _Normal->display();
  else
    std::cout << " nullptr\n";

  std::cout << " _Tangent:\n";
  if (_Tangent)
    _Tangent->display();
  else
    std::cout << " nullptr\n";
}
