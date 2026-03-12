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

#include "NewtonEuler5DR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "RotationQuaternion.hpp"  // siconos::geometry::computeRotationMatrix
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "op3x3.h"  // numerics: orthobasefromvector
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

/*
See devNotes.pdf for details. A detailed documentation is available in DevNotes.pdf: chapter
'NewtonEulerR: computation of \nabla q H'. Subsection 'Case RFC3D: using the local frame local
velocities'
*/
void siconos::modeling::NewtonEuler5DR::initialize(Interaction& inter) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEuler5DR::siconos::modeling::NewtonEuler5DR::initialize("
      "Interaction& inter)\n");
  auto qSize = 7 * (inter.getSizeOfDS() / 6);
  /*keep only the distance.*/

  H_NE_internal_storage_.resize(5, qSize);
  H_NE_view_ =
      std::make_shared<siconos::algebra::MapType>(H_NE_internal_storage_.data(), 5, qSize);
  H_NE_view_->setZero();
  NewtonEulerR::initialize(inter);
  _rotationAbsoluteToContactFrame = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  _rotationBodyToAbsoluteFrame = std::make_shared<siconos::algebra::SiconosMatrix33>();
  _AUX1 = std::make_shared<siconos::algebra::SiconosMatrix33>();
  _AUX2 = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  _NPG1 = std::make_shared<siconos::algebra::SiconosMatrix33>();
  _NPG2 = std::make_shared<siconos::algebra::SiconosMatrix33>();

  //  _isContact=1;
  DEBUG_END(
      "siconos::modeling::NewtonEuler5DR::siconos::modeling::NewtonEuler5DR::initialize("
      "Interaction& inter)\n");
}

void siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1) {
  DEBUG_BEGIN("siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts()\n");
  double Nx = (*_Nc)(0);
  double Ny = (*_Nc)(1);
  double Nz = (*_Nc)(2);
  double Px = (*_Pc1)(0);
  double Py = (*_Pc1)(1);
  double Pz = (*_Pc1)(2);
  double G1x = q1(0);
  double G1y = q1(1);
  double G1z = q1(2);

  DEBUG_PRINT("contact normal:\n");
  DEBUG_EXPR(siconos::algebra::print(*_Nc););
  DEBUG_PRINTF("_Nc->norm() -1.0 = %e\n", _Nc->norm() - 1.0);
  DEBUG_PRINT("contact point :\n");
  DEBUG_EXPR(siconos::algebra::print(*_Pc1););
  DEBUG_PRINT("center of mass :\n");
  DEBUG_EXPR(siconos::algebra::print(q1););

  assert(_Nc->norm() > 0.0 && std::abs(_Nc->norm() - 1.0) < 1e-6 &&
         "siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts. Normal vector "
         "not consistent ");

  double t[6];
  double* pt = t;

  // 1 - Construction of the local contact frame from the normal vector

  if (orthoBaseFromVector(&Nx, &Ny, &Nz, pt, pt + 1, pt + 2, pt + 3, pt + 4, pt + 5))
    THROW_EXCEPTION(
        "siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts. Problem in "
        "calling orthoBaseFromVector");

  // 2 - Construction of the rotation matrix from the absolute frame to the local contact frame
  pt = t;
  _rotationAbsoluteToContactFrame->setValue(0, 0, Nx);
  _rotationAbsoluteToContactFrame->setValue(1, 0, *pt);
  _rotationAbsoluteToContactFrame->setValue(2, 0, *(pt + 3));
  _rotationAbsoluteToContactFrame->setValue(0, 1, Ny);
  _rotationAbsoluteToContactFrame->setValue(1, 1, *(pt + 1));
  _rotationAbsoluteToContactFrame->setValue(2, 1, *(pt + 4));
  _rotationAbsoluteToContactFrame->setValue(0, 2, Nz);
  _rotationAbsoluteToContactFrame->setValue(1, 2, *(pt + 2));
  _rotationAbsoluteToContactFrame->setValue(2, 2, *(pt + 5));
  DEBUG_PRINT("_rotationAbsoluteToContactFrame:\n");
  DEBUG_EXPR(siconos::algebra::print(*_rotationAbsoluteToContactFrame););

  // 3 - Construction of the lever arm matrix in  the absolute frame

  _NPG1->setZero();
  (*_NPG1)(0, 0) = 0;
  (*_NPG1)(0, 1) = -(G1z - Pz);
  (*_NPG1)(0, 2) = (G1y - Py);
  (*_NPG1)(1, 0) = (G1z - Pz);
  (*_NPG1)(1, 1) = 0;
  (*_NPG1)(1, 2) = -(G1x - Px);
  (*_NPG1)(2, 0) = -(G1y - Py);
  (*_NPG1)(2, 1) = (G1x - Px);
  (*_NPG1)(2, 2) = 0;

  DEBUG_PRINT("lever arm skew matrix :\n");
  DEBUG_EXPR(siconos::algebra::print(*_NPG1););

  /* The Jacobian matrix (H) is given by the product
   * H = _rotationAbsoluteToContactFrame
   * for the translation part and
   * H = _rotationAbsoluteToContactFrame * leverArmMatrix * _rotationBodyToAbsoluteFrame
   * for the rotation part and
   */

  // 4 - Compute the rotation matrix from the body-fixed frame to the absolute frame
  siconos::geometry::computeRotationMatrix(q1, *_rotationBodyToAbsoluteFrame);
  DEBUG_EXPR(siconos::algebra::print(*_rotationBodyToAbsoluteFrame););

  // 5 - compose the body lever arm matrix with the rotation matrix
  _AUX1->noalias() = *_NPG1 * *_rotationBodyToAbsoluteFrame;

  DEBUG_EXPR(siconos::algebra::print(*_rotationBodyToAbsoluteFrame););
  DEBUG_EXPR(siconos::algebra::print(*_AUX1););

  // 6 -  Rotate the resulting matric in the contact frame
  _AUX2->noalias() = *_rotationAbsoluteToContactFrame * *_AUX1;
  DEBUG_EXPR(siconos::algebra::print(*_rotationAbsoluteToContactFrame););
  DEBUG_EXPR(siconos::algebra::print(*_AUX2););

  // 7 - fill the Jacobian
  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      H_NE_prod_T_->setValue(ii, jj, (*_rotationAbsoluteToContactFrame)(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      H_NE_prod_T_->setValue(ii, jj, (*_AUX2)(ii, jj - 3));

  _AUX2->noalias() = *_rotationAbsoluteToContactFrame * *_rotationBodyToAbsoluteFrame;
  DEBUG_EXPR(siconos::algebra::print(*_AUX2););

  for (unsigned int ii = 3; ii < 5; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      H_NE_prod_T_->setValue(ii, jj, (*_AUX2)(ii - 2, jj - 3));

  DEBUG_EXPR(siconos::algebra::print(*jacobianhOver_q_T););

  // DEBUG_EXPR_WE(
  //   std::shared_ptr<siconos::algebra::SiconosMatrix> jaux =
  //   std::make_shared<siconos::algebra::SiconosMatrix>(*jacobianhOver_q_T)); jaux->trans();
  //   std::shared_ptr<siconos::algebra::SiconosVector> v =
  //   std::make_shared<siconos::algebra::SiconosVector>(3));
  //   std::shared_ptr<siconos::algebra::SiconosVector> vRes =
  //   std::make_shared<siconos::algebra::SiconosVector>(6)); v->setZero(); (*v)(0) = 1;
  //   *vRes= *jaux **v;
  //   siconos::algebra::print(*vRes);
  //   v->setZero();
  //   (*v)(1) = 1;
  //   *vRes= *jaux **v;
  //   siconos::algebra::print(*vRes);
  //   v->setZero();
  //   (*v)(2) = 1;
  //   *vRes= *jaux **v;
  //   siconos::algebra::print(*vRes);
  //   );
  DEBUG_END(
      "siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts(std::shared_ptr<"
      "siconos::algebra::SiconosVector> q1)\n");
}

void siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q2) {
  DEBUG_BEGIN("siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts()\n");

  double Nx = (*_Nc)(0);
  double Ny = (*_Nc)(1);
  double Nz = (*_Nc)(2);
  double Px = (*_Pc1)(0);
  double Py = (*_Pc1)(1);
  double Pz = (*_Pc1)(2);
  double G1x = q1(0);
  double G1y = q1(1);
  double G1z = q1(2);
  double G2x = q2(0);
  double G2y = q2(1);
  double G2z = q2(2);

  DEBUG_PRINT("contact normal:\n");
  DEBUG_EXPR(siconos::algebra::print(*_Nc););
  DEBUG_PRINT("contact point :\n");
  DEBUG_EXPR(siconos::algebra::print(*_Pc1););
  DEBUG_PRINT("center of mass :\n");
  DEBUG_EXPR(siconos::algebra::print(q1););

  double t[6];
  double* pt = t;
  if (orthoBaseFromVector(&Nx, &Ny, &Nz, pt, pt + 1, pt + 2, pt + 3, pt + 4, pt + 5))
    THROW_EXCEPTION(
        "siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts. Problem in "
        "calling orthoBaseFromVector");
  pt = t;

  _rotationAbsoluteToContactFrame->setValue(0, 0, Nx);
  _rotationAbsoluteToContactFrame->setValue(1, 0, *pt);
  _rotationAbsoluteToContactFrame->setValue(2, 0, *(pt + 3));
  _rotationAbsoluteToContactFrame->setValue(0, 1, Ny);
  _rotationAbsoluteToContactFrame->setValue(1, 1, *(pt + 1));
  _rotationAbsoluteToContactFrame->setValue(2, 1, *(pt + 4));
  _rotationAbsoluteToContactFrame->setValue(0, 2, Nz);
  _rotationAbsoluteToContactFrame->setValue(1, 2, *(pt + 2));
  _rotationAbsoluteToContactFrame->setValue(2, 2, *(pt + 5));

  _NPG1->setZero();

  (*_NPG1)(0, 0) = 0;
  (*_NPG1)(0, 1) = -(G1z - Pz);
  (*_NPG1)(0, 2) = (G1y - Py);
  (*_NPG1)(1, 0) = (G1z - Pz);
  (*_NPG1)(1, 1) = 0;
  (*_NPG1)(1, 2) = -(G1x - Px);
  (*_NPG1)(2, 0) = -(G1y - Py);
  (*_NPG1)(2, 1) = (G1x - Px);
  (*_NPG1)(2, 2) = 0;

  _NPG2->setZero();

  (*_NPG2)(0, 0) = 0;
  (*_NPG2)(0, 1) = -(G2z - Pz);
  (*_NPG2)(0, 2) = (G2y - Py);
  (*_NPG2)(1, 0) = (G2z - Pz);
  (*_NPG2)(1, 1) = 0;
  (*_NPG2)(1, 2) = -(G2x - Px);
  (*_NPG2)(2, 0) = -(G2y - Py);
  (*_NPG2)(2, 1) = (G2x - Px);
  (*_NPG2)(2, 2) = 0;

  siconos::geometry::computeRotationMatrix(q1, *_rotationBodyToAbsoluteFrame);
  _AUX1->noalias() = *_NPG1 * *_rotationBodyToAbsoluteFrame;
  _AUX2->noalias() = *_rotationAbsoluteToContactFrame * *_AUX1;

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      H_NE_prod_T_->setValue(ii, jj, (*_rotationAbsoluteToContactFrame)(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      H_NE_prod_T_->setValue(ii, jj, (*_AUX2)(ii, jj - 3));

  _AUX2->noalias() = *_rotationAbsoluteToContactFrame * *_rotationBodyToAbsoluteFrame;
  DEBUG_EXPR(siconos::algebra::print(*_AUX2););

  for (unsigned int ii = 3; ii < 5; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      H_NE_prod_T_->setValue(ii, jj, (*_AUX2)(ii - 2, jj - 3));

  siconos::geometry::computeRotationMatrix(q2, *_rotationBodyToAbsoluteFrame);
  _AUX1->noalias() = *_NPG2 * *_rotationBodyToAbsoluteFrame;
  _AUX2->noalias() = *_rotationAbsoluteToContactFrame * *_AUX1;

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      H_NE_prod_T_->setValue(ii, jj + 6, -(*_rotationAbsoluteToContactFrame)(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      H_NE_prod_T_->setValue(ii, jj + 6, -(*_AUX2)(ii, jj - 3));

  _AUX2->noalias() = *_rotationAbsoluteToContactFrame * *_rotationBodyToAbsoluteFrame;
  DEBUG_EXPR(siconos::algebra::print(*_AUX2););

  for (unsigned int ii = 3; ii < 5; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      H_NE_prod_T_->setValue(ii, jj, -(*_AUX2)(ii - 2, jj - 3));

  DEBUG_EXPR(siconos::algebra::print(*jacobianhOver_q_T););

  DEBUG_END(
      "siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts(std::shared_ptr<"
      "siconos::algebra::SiconosVector> q1, std::shared_ptr<siconos::algebra::SiconosVector> "
      "q2)\n");
}

void siconos::modeling::NewtonEuler5DR::computeH_NE_prod_T(
    const Interaction& inter, const siconos::algebra::BlockVector& q0) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEuler5DR::computeH_NE_prod_T(Interaction& inter,  "
      "std::shared_ptr<siconos::algebra::BlockVector> q0)\n");
  if (q0.numberOfBlocks() > 1) {
    RFC3DcomputeJachqTFromContacts(*q0.vector(0), *q0.vector(1));
  } else {
    RFC3DcomputeJachqTFromContacts(*q0.vector(0));
  }
  DEBUG_END(
      "siconos::modeling::NewtonEuler5DR::computeH_NE_prod_T(Interaction& inter,  "
      "std::shared_ptr<siconos::algebra::BlockVector> q0)\n");
}
