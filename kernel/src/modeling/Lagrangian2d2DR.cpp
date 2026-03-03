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
#include "Lagrangian2d2DR.hpp"

#include <iostream>

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::modeling::Lagrangian2d2DR::initialize(Interaction& inter) {
  // proj_with_q  jacobianhOver_q_Proj =
  // std::make_shared<siconos::algebra::siconos::algebra::SiconosMatrix>(jacobianhOver_q_->rows(),jacobianhOver_q_->cols()));

  if ((inter.getSizeOfDS() != 3) and (inter.getSizeOfDS() != 6)) {
    THROW_EXCEPTION(
        "siconos::modeling::Lagrangian2d2DR::initialize(Interaction& inter). The size of ds "
        "must of size 3 or 6");
  }
  auto qSize = 3 * (inter.getSizeOfDS() / 3);
  if (!jacobianhOver_q_internal_storage_) {
    jacobianhOver_q_internal_storage_ = std::make_unique<std::vector<double>>(2 * qSize);
  }
  jacobianhOver_q_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobianhOver_q_internal_storage_->data(), 2, qSize);
}

double siconos::modeling::Lagrangian2d2DR::distance() const {
  DEBUG_BEGIN("siconos::modeling::Lagrangian2d2DR::distance(...)\n")
  siconos::algebra::SiconosVector dpc(*_Pc2 - *_Pc1);
  DEBUG_EXPR(siconos::algebra::print(*_Pc1););
  DEBUG_EXPR(siconos::algebra::print(*_Pc2););
  DEBUG_EXPR(siconos::algebra::print(dpc););
  DEBUG_END("siconos::modeling::Lagrangian2d2DR::distance(...)\n")
  return dpc.norm() * (_Nc->dot(dpc) >= 0 ? -1 : 1);
}

void siconos::modeling::Lagrangian2d2DR::computeh(
    const siconos::algebra::BlockVector& q, Eigen::Ref<siconos::algebra::SiconosVector> y) {
  DEBUG_BEGIN("siconos::modeling::Lagrangian2d2DR::computeh(...)\n");
  DEBUG_EXPR(siconos::algebra::print(q));

  DEBUG_EXPR(siconos::algebra::print(*_Pc1););
  DEBUG_EXPR(siconos::algebra::print(*_Pc2););
  DEBUG_EXPR(siconos::algebra::print(*_Nc););

  LagrangianScleronomousR::computeh(q, y);
  y(0) = distance();

  DEBUG_EXPR(siconos::algebra::print(y););
  DEBUG_EXPR(display(););
  DEBUG_END("siconos::modeling::Lagrangian2d2DR::computeh(...)\n");
  // getchar();
}

void siconos::modeling::Lagrangian2d2DR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector& q) {
  DEBUG_BEGIN(
      "siconos::modeling::Lagrangian2d2DR::computeJacobianhOver_q(Interaction& inter, "
      "Ssiconos::algebra::BlockVector q0 \n");

  double Nx = (*_Nc)(0);
  double Ny = (*_Nc)(1);
  double Px = (*_Pc1)(0);
  double Py = (*_Pc1)(1);
  double G1x = q(0);
  double G1y = q(1);

  /* construct tangent vector */
  double Tx = -Ny;
  double Ty = Nx;

  double lever_arm_x = Px - G1x;
  double lever_arm_y = Py - G1y;
  DEBUG_PRINTF("N_x = %4.2e,\t N_ y = %4.2e\n", Nx, Ny);
  DEBUG_PRINTF("lever_arm_x = %4.2e,\t lever_arm_ y = %4.2e\n", lever_arm_x, lever_arm_y);

  // jacobianhOver_q_->setValue(0,0,Nx);
  // jacobianhOver_q_->setValue(0,1,Ny);
  // jacobianhOver_q_->setValue(0,2,lever_arm_x*Ny - lever_arm_y*Nx );

  // jacobianhOver_q_->setValue(1,0,Tx);
  // jacobianhOver_q_->setValue(1,1,Ty);
  // jacobianhOver_q_->setValue(1,2,lever_arm_x*Ty - lever_arm_y*Tx );

  double* array = &*jacobianhOver_q_view_->data();
  array[0] = Nx;
  array[2] = Ny;
  array[4] = lever_arm_x * Ny - lever_arm_y * Nx;

  array[1] = Tx;
  array[3] = Ty;
  array[5] = lever_arm_x * Ty - lever_arm_y * Tx;

  array[0] = Nx;
  array[2] = Ny;
  array[4] = lever_arm_x * Ny - lever_arm_y * Nx;

  array[1] = Tx;
  array[3] = Ty;
  array[5] = lever_arm_x * Ty - lever_arm_y * Tx;

  if (q.size() == 6) {
    DEBUG_PRINT("take into account second ds\n");
    double G2x = q(3);
    double G2y = q(4);
    lever_arm_x = Px - G2x;
    lever_arm_y = Py - G2y;

    DEBUG_PRINTF("lever_arm_x = %4.2e,\t lever_arm_ y = %4.2e\n", lever_arm_x, lever_arm_y);

    // jacobianhOver_q_->setValue(0,3,-Nx);
    // jacobianhOver_q_->setValue(0,4,-Ny);
    // jacobianhOver_q_->setValue(0,5,lever_arm_y * Nx - lever_arm_x*Ny);

    // jacobianhOver_q_->setValue(1,3,-Tx);
    // jacobianhOver_q_->setValue(1,4,-Ty);
    // jacobianhOver_q_->setValue(1,5,lever_arm_y * Tx - lever_arm_x*Ty);
    array[6] = -Nx;
    array[8] = -Ny;
    array[10] = lever_arm_y * Nx - lever_arm_x * Ny;

    array[7] = -Tx;
    array[9] = -Ty;
    array[11] = lever_arm_y * Tx - lever_arm_x * Ty;
  }
  DEBUG_EXPR(siconos::algebra::print(*jacobianhOver_q_););
  DEBUG_END(
      "siconos::modeling::Lagrangian2d2DR::computeJacobianhOver_q(Interaction& inter, "
      "siconos::algebra::BlockVector q0) \n");
}

void siconos::modeling::Lagrangian2d2DR::display() const {
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

  // std::cout << " _relPc1 :" << std::endl;
  // if(_relPc1)
  //   siconos::algebra::print(*_relPc1);
  // else
  //   std::cout << " nullptr :" << std::endl;

  // std::cout << " _relPc2 :" << std::endl;
  // if(_relPc2)
  //   siconos::algebra::print(*_relPc2);
  // else
  //   std::cout << " nullptr :" << std::endl;

  std::cout << " _Nc :" << std::endl;
  if (_Nc)
    siconos::algebra::print(*_Nc);
  else
    std::cout << " nullptr :" << std::endl;
  // std::cout << " _relNc :" << std::endl;
  // if(_relNc)
  //   siconos::algebra::print(*_relNc);
  // else
  //   std::cout << " nullptr :" << std::endl;
}

// void siconos::modeling::Lagrangian2d2DR::computeOutput(double time, Interaction& inter,
// siconos::algebra::blocks::size_type derivativeNumber)
// {

//   DEBUG_PRINTF("siconos::modeling::Lagrangian2d2DR::computeOutput(double time, Interaction&
//   inter, InteractionProperties& interProp, siconos::algebra::blocks::size_type derivativeNumber) with time = %f and
//   derivativeNumber = %i\n", time, derivativeNumber);
//   std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& ds_vars =
//   inter.read_dynamical_systems_variables(); siconos::algebra::SiconosVector& y =
//   *inter.y(derivativeNumber); if(derivativeNumber == 0)
//   {
//     computeh(*ds_vars[tools::enum_to_index(LagrangianR::ds_var::q0)], *ds_vars[LagrangianR::z], y);
//   }
//   else
//   {
//     computeJacobianhOver_q(*ds_vars[tools::enum_to_index(LagrangianR::ds_var::q0)], *ds_vars[LagrangianR::z]);
//     if(derivativeNumber == 1)
//     {
//       assert(jacobianhOver_q_);

//       // direct prod to save time
//       //siconos::algebra::prod(*jacobianhOver_q_, *ds_vars[tools::enum_to_index(ds_var::q1)], y);

//       double * A = &*jacobianhOver_q_->data();
//       std::shared_ptr<siconos::algebra::BlockVector> v = ds_vars[tools::enum_to_index(ds_var::q1)];
//       double *  v_ds_1 = v->vector(0)->data();

//       y(0)= A[0]* v_ds_1[0] + A[2]* v_ds_1[1]  + A[4]* v_ds_1[2];
//       y(1)= A[1]* v_ds_1[0] + A[3]* v_ds_1[1]  + A[5]* v_ds_1[2];

//       if (v ->numberOfBlocks() >1 )
//       {
//         double *  v_ds_2 = v->vector(1)->data();
//         y(0) += A[6]* v_ds_2[0] + A[8]* v_ds_2[1]  + A[10]* v_ds_2[2];
//         y(1) += A[7]* v_ds_2[0] + A[9]* v_ds_2[1]  + A[11]* v_ds_2[2];

//       }

//     // }
//     //   else
//     //   {
//     //     y(0)= A[0]* (*v)(0) + A[2]* (*v)(1) + A[4]* (*v)(2)
//     //         + A[6]* (*v)(3) + A[8]* (*v)(4) + A[10]* (*v)(5);
//     //     y(1)= A[1]* (*v)(0) + A[3]* (*v)(1) + A[5]* (*v)(2)
//     //         + A[7]* (*v)(3) + A[9]* (*v)(4) + A[11]* (*v)(5);

//     //   }
//       // for (siconos::algebra::blocks::size_type i =0; i < 2; i++)
//       // {
//       //   y(i)= A[i]* (*v)(0);
//       //   for (siconos::algebra::blocks::size_type j =1; j<v->size(); j++)
//       //   {
//       //     y(i) += A[i+j*2] * (*v)(j);
//       //   }
//       // }
//     }
//     else
//       THROW_EXCEPTION("siconos::modeling::Lagrangian2d2DR::computeOutput(double time,
//       Interaction& inter, InteractionProperties& interProp, siconos::algebra::blocks::size_type derivativeNumber),
//       index out of range");
//   }
// }

// void LagrangianScleronomousR::computeInput(double time, Interaction& inter, siconos::algebra::blocks::size_type
// level)
// {
//   DEBUG_BEGIN("void LagrangianScleronomousR::computeInput(double time, Interaction& inter,
//   InteractionProperties& interProp, siconos::algebra::blocks::size_type level) \n");

//   DEBUG_PRINTF("level = %i\n", level);
//   const auto& ds_vars = inter.read_dynamical_systems_variables();
//   computeJacobianhOver_q(*ds_vars[tools::enum_to_index(LagrangianR::ds_var::q0)], *ds_vars[LagrangianR::z]);
//   // get lambda of the concerned interaction
//   siconos::algebra::SiconosVector& lambda = *inter.lambda(level);
//   DEBUG_EXPR(siconos::algebra::print(lambda););
//   DEBUG_EXPR(siconos::algebra::print(*jacobianhOver_q_););
//   // data[name] += trans(G) * lambda
//   //siconos::algebra::prod(lambda, *jacobianhOver_q_, *ds_vars[tools::enum_to_index(ds_var::p0) + level],
//   false);

//   double * A = &*jacobianhOver_q_->data();
//   std::shared_ptr<siconos::algebra::BlockVector> v = ds_vars[tools::enum_to_index(ds_var::q1)];
//   int v_size= v->size();
//   for (siconos::algebra::blocks::size_type i =0; i < 2; i++)
//   {
//     y(i)= A[i]* (*v)(0);
//     for (siconos::algebra::blocks::size_type j =1; j<v->size(); j++)
//     {
//       y(i) += A[i+j*2] * (*v)(j);
//     }
//   }

//   DEBUG_END("void LagrangianScleronomousR::computeInput(double time, Interaction& inter,
//   InteractionProperties& interProp, siconos::algebra::blocks::size_type level) \n");
// }
