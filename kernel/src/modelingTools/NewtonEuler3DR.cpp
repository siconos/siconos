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

#include "NewtonEuler3DR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "RotationQuaternion.hpp"  // siconos::geometry::computeRotationMatrix
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"  // For mat prod
#include "SiconosVector.hpp"
#include "op3x3.h"  // numerics: orthobasefromvector

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

/*
See devNotes.pdf for details. A detailed documentation is available in DevNotes.pdf: chapter
'NewtonEulerR: computation of \nabla q H'. Subsection 'Case FC3D: using the local frame local
velocities'
*/
void siconos::modeling::NewtonEuler3DR::initialize(Interaction& inter) {
  unsigned int qSize = 7 * (inter.getSizeOfDS() / 6);
  /*keep only the distance.*/

  //  H_NE_internal_storage_ = std::make_unique<std::vector<double>>(3 * qSize);
  H_NE_internal_storage_.resize(3, qSize);
  H_NE_view_ =
      std::make_shared<siconos::algebra::MapType>(H_NE_internal_storage_.data(), 3, qSize);
  H_NE_view_->setZero();
  NewtonEulerR::initialize(inter);

  /* VA 12/04/2016 All of what follows should be put in WorkM*/
  _rotationAbsoluteToContactFrame = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  _rotationBodyToAbsoluteFrame = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  _AUX1 = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  _AUX2 = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  _NPG1 = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  _NPG2 = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  //  _isContact=1;
}

void siconos::modeling::NewtonEuler3DR::FC3DcomputeJachqTFromContacts(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& q1) {
  DEBUG_BEGIN("siconos::modeling::NewtonEuler3DR::FC3DcomputeJachqTFromContacts()\n");
  double Nx = _Nc->getValue(0);
  double Ny = _Nc->getValue(1);
  double Nz = _Nc->getValue(2);
  double Px = _Pc1->getValue(0);
  double Py = _Pc1->getValue(1);
  double Pz = _Pc1->getValue(2);
  double G1x = q1.getValue(0);
  double G1y = q1.getValue(1);
  double G1z = q1.getValue(2);

  DEBUG_PRINT("contact normal:\n");
  DEBUG_EXPR(_Nc->display(););
  DEBUG_PRINTF("_Nc->norm2() -1.0 = %e\n", _Nc->norm2() - 1.0);
  DEBUG_PRINT("contact point :\n");
  DEBUG_EXPR(_Pc1->display(););
  DEBUG_PRINT("center of mass :\n");
  DEBUG_EXPR(q1.display(););

  assert(_Nc->norm2() > 0.0 && std::abs(_Nc->norm2() - 1.0) < 1e-6 &&
         "siconos::modeling::NewtonEuler3DR::FC3DcomputeJachqTFromContacts. Normal vector not "
         "consistent ");

  double t[6];
  double* pt = t;

  // 1 - Construction of the local contact frame from the normal vector

  if (orthoBaseFromVector(&Nx, &Ny, &Nz, pt, pt + 1, pt + 2, pt + 3, pt + 4, pt + 5))
    THROW_EXCEPTION(
        "siconos::modeling::NewtonEuler3DR::FC3DcomputeJachqTFromContacts. Problem in calling "
        "orthoBaseFromVector");

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
  DEBUG_EXPR(_rotationAbsoluteToContactFrame->display(););

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
  DEBUG_EXPR(_NPG1->display(););

  /* The Jacobian matrix (H) is given by the product
   * H = _rotationAbsoluteToContactFrame
   * for the translation part and
   * H = _rotationAbsoluteToContactFrame * leverArmMatrix * _rotationBodyToAbsoluteFrame
   * for the rotation part and
   */

  // 4 - Compute the rotation matrix from the body-fixed frame to the absolute frame
  siconos::geometry::computeRotationMatrix(q1, *_rotationBodyToAbsoluteFrame);
  DEBUG_EXPR(_rotationBodyToAbsoluteFrame->display(););

  // 5 - compose the body lever arm matrix with the rotation matrix
  siconos::algebra::prod(*_NPG1, *_rotationBodyToAbsoluteFrame, *_AUX1, true);
  DEBUG_EXPR(_rotationBodyToAbsoluteFrame->display(););
  DEBUG_EXPR(_AUX1->display(););

  // 6 -  Rotate the resulting matric in the contact frame
  siconos::algebra::prod(*_rotationAbsoluteToContactFrame, *_AUX1, *_AUX2, true);
  DEBUG_EXPR(_rotationAbsoluteToContactFrame->display(););
  DEBUG_EXPR(_AUX2->display(););

  // 7 - fill the Jacobian

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      H_NE_prod_T_->setValue(ii, jj, _rotationAbsoluteToContactFrame->getValue(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      H_NE_prod_T_->setValue(ii, jj, _AUX2->getValue(ii, jj - 3));

  DEBUG_EXPR(jacobianhOver_q_T->display(););
  // DEBUG_EXPR_WE(
  //   std::shared_ptr<siconos::algebra::SiconosMatrix> jaux =
  //   std::make_shared<siconos::algebra::SiconosMatrix>(*jacobianhOver_q_T)); jaux->trans();
  //   std::shared_ptr<siconos::algebra::SiconosVector> v =
  //   std::make_shared<siconos::algebra::SiconosVector>(3));
  //   std::shared_ptr<siconos::algebra::SiconosVector> vRes =
  //   std::make_shared<siconos::algebra::SiconosVector>(6)); v->setZero(); v->setValue(0, 1);
  //   siconos::algebra::prod(*jaux, *v, *vRes, true);
  //   vRes->display();
  //   v->setZero();
  //   v->setValue(1, 1);
  //   siconos::algebra::prod(*jaux, *v, *vRes, true);
  //   vRes->display();
  //   v->setZero();
  //   v->setValue(2, 1);
  //   siconos::algebra::prod(*jaux, *v, *vRes, true);
  //   vRes->display();
  //   );
  DEBUG_END(
      "siconos::modeling::NewtonEuler3DR::FC3DcomputeJachqTFromContacts(std::shared_ptr<"
      "siconos::algebra::SiconosVector> q1)\n");
}

void siconos::modeling::NewtonEuler3DR::FC3DcomputeJachqTFromContacts(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& q1,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& q2) {
  double Nx = _Nc->getValue(0);
  double Ny = _Nc->getValue(1);
  double Nz = _Nc->getValue(2);
  double Px = _Pc1->getValue(0);
  double Py = _Pc1->getValue(1);
  double Pz = _Pc1->getValue(2);
  double G1x = q1.getValue(0);
  double G1y = q1.getValue(1);
  double G1z = q1.getValue(2);
  double G2x = q2.getValue(0);
  double G2y = q2.getValue(1);
  double G2z = q2.getValue(2);

  DEBUG_PRINT("contact normal:\n");
  DEBUG_EXPR(_Nc->display(););
  DEBUG_PRINT("contact point :\n");
  DEBUG_EXPR(_Pc1->display(););
  DEBUG_PRINT("center of mass :\n");
  DEBUG_EXPR(q1.display(););

  double t[6];
  double* pt = t;
  if (orthoBaseFromVector(&Nx, &Ny, &Nz, pt, pt + 1, pt + 2, pt + 3, pt + 4, pt + 5))
    THROW_EXCEPTION(
        "siconos::modeling::NewtonEuler3DR::FC3DcomputeJachqTFromContacts. Problem in calling "
        "orthoBaseFromVector");
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
  siconos::algebra::prod(*_NPG1, *_rotationBodyToAbsoluteFrame, *_AUX1, true);
  siconos::algebra::prod(*_rotationAbsoluteToContactFrame, *_AUX1, *_AUX2, true);

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      H_NE_prod_T_->setValue(ii, jj, _rotationAbsoluteToContactFrame->getValue(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      H_NE_prod_T_->setValue(ii, jj, _AUX2->getValue(ii, jj - 3));

  siconos::geometry::computeRotationMatrix(q2, *_rotationBodyToAbsoluteFrame);
  siconos::algebra::prod(*_NPG2, *_rotationBodyToAbsoluteFrame, *_AUX1, true);
  siconos::algebra::prod(*_rotationAbsoluteToContactFrame, *_AUX1, *_AUX2, true);

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      H_NE_prod_T_->setValue(ii, jj + 6, -_rotationAbsoluteToContactFrame->getValue(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      H_NE_prod_T_->setValue(ii, jj + 6, -_AUX2->getValue(ii, jj - 3));
}

void siconos::modeling::NewtonEuler3DR::computeH_NE_prod_T(
    const Interaction& inter, const siconos::algebra::BlockVector& q0) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEuler3DR::computeH_NE_prod_T(Interaction& inter,  "
      "std::shared_ptr<siconos::algebra::BlockVector> q0)\n");
  if (q0.numberOfBlocks() > 1) {
    FC3DcomputeJachqTFromContacts(*q0.vector(0), *q0.vector(1));
  } else {
    FC3DcomputeJachqTFromContacts(*q0.vector(0));
  }
  DEBUG_END(
      "siconos::modeling::NewtonEuler3DR::computeH_NE_prod_T(Interaction& inter,  "
      "std::shared_ptr<siconos::algebra::BlockVector> q0)\n");
}
