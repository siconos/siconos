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

#include "Lagrangian2d1DR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::modeling::Lagrangian2d1DR::initialize(Interaction& inter) {
  // proj_with_q  jacobianhOver_q_Proj =
  // std::make_shared<siconos::algebra::SiconosMatrix>(jacobianhOver_q_->rows(),jacobianhOver_q_->cols()));

  if ((inter.getSizeOfDS() != 3) and (inter.getSizeOfDS() != 6)) {
    THROW_EXCEPTION(
        "siconos::modeling::Lagrangian2d1DR::initialize(Interaction& inter). The size of ds "
        "must of size 3 or 6");
  }

  auto qSize = 3 * (inter.getSizeOfDS() / 3);
  if (!jacobianhOver_q_internal_storage_) {
    jacobianhOver_q_internal_storage_ = std::make_unique<std::vector<double>>(qSize);
  }
  jacobianhOver_q_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobianhOver_q_internal_storage_->data(), 1, qSize);
}

void siconos::modeling::Lagrangian2d1DR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector& q) {
  DEBUG_BEGIN(
      "siconos::modeling::Lagrangian2d1DR::computeJacobianhOver_q(Interaction& inter, "
      "siconos::algebra::BlockVector q0 \n");

  double Nx = (*_Nc)(0);
  double Ny = (*_Nc)(1);
  double Px = (*_Pc1)(0);
  double Py = (*_Pc1)(1);
  double G1x = q(0);
  double G1y = q(1);

  jacobianhOver_q_view_->setValue(0, 0, Nx);
  jacobianhOver_q_view_->setValue(0, 1, Ny);
  jacobianhOver_q_view_->setValue(0, 2, (G1y - Py) * Nx - (G1x - Px) * Ny);

  if (q.size() == 6) {
    DEBUG_PRINT("take into account second ds\n");
    double G2x = q(3);
    double G2y = q(4);

    jacobianhOver_q_view_->setValue(0, 3, -Nx);
    jacobianhOver_q_view_->setValue(0, 4, -Ny);
    jacobianhOver_q_view_->setValue(0, 5, -((G2y - Py) * Nx - (G2x - Px) * Ny));
  }
  DEBUG_EXPR(siconos::algebra::print(*jacobianhOver_q_););
  DEBUG_END(
      "siconos::modeling::Lagrangian2d1DR::computeJacobianhOver_q(Interaction& inter, "
      "siconos::algebra::BlockVector q0) \n");
}

double siconos::modeling::Lagrangian2d1DR::distance() const {
  DEBUG_BEGIN("siconos::modeling::Lagrangian2d1DR::distance(...)\n")
  siconos::algebra::SiconosVector dpc(*_Pc2 - *_Pc1);
  DEBUG_END("siconos::modeling::Lagrangian2d1DR::distance(...)\n")
  return dpc.norm() * (_Nc->dot(dpc) >= 0 ? -1 : 1);
}

void siconos::modeling::Lagrangian2d1DR::computeh(
    const siconos::algebra::BlockVector& q, Eigen::Ref<siconos::algebra::SiconosVector> y) {
  DEBUG_BEGIN("siconos::modeling::Lagrangian2d1DR::computeh(...)\n");
  DEBUG_EXPR(siconos::algebra::print(q));

  LagrangianScleronomousR::computeh(q, y);
  y(0) = distance();
  DEBUG_EXPR(siconos::algebra::print(y););
  DEBUG_EXPR(display(););
  DEBUG_END("siconos::modeling::Lagrangian2d1DR::computeh(...)\n")
}
void siconos::modeling::Lagrangian2d1DR::display() const {
  LagrangianR::display();

  std::cout << " _Pc1 :" << std::endl;
  if (_Pc1)
    siconos::algebra::print(*_Pc1);
  else
    std::cout << " nullptr :" << std::endl;

  std::cout << " _Pc2 :" << std::endl;
  if (_Pc2)
    siconos::algebra::print(*_Pc2);
  else
    std::cout << " nullptr :" << std::endl;

  std::cout << " _Nc :" << std::endl;
  if (_Nc)
    siconos::algebra::print(*_Nc);
  else
    std::cout << " nullptr :" << std::endl;
  std::cout << " _relNc :" << std::endl;
  if (_relNc)
    siconos::algebra::print(*_relNc);
  else
    std::cout << " nullptr :" << std::endl;
}
