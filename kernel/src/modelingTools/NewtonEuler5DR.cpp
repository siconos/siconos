/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include "Interaction.hpp"
#include "BlockVector.hpp"
#include "RotationQuaternion.hpp"
#include "op3x3.h"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

/*
See devNotes.pdf for details. A detailed documentation is available in DevNotes.pdf: chapter 'NewtonEulerR: computation of \nabla q H'. Subsection 'Case RFC3D: using the local frame local velocities'
*/
void NewtonEuler5DR::initialize(Interaction& inter)
{
  DEBUG_BEGIN("NewtonEuler5DR::NewtonEuler5DR::initialize(Interaction& inter)\n");
  NewtonEuler1DR::initialize(inter);
  unsigned int qSize = 7 * (inter.getSizeOfDS() / 6);
  /*keep only the distance.*/
  _jachq.reset(new SimpleMatrix(5, qSize));

  _rotationAbsoluteToContactFrame.reset(new SimpleMatrix(3, 3));
  _AUX2.reset(new SimpleMatrix(3, 3));
  //  _isContact=1;
  DEBUG_END("NewtonEuler5DR::NewtonEuler5DR::initialize(Interaction& inter)\n");
}
void NewtonEuler5DR::RFC3DcomputeJachqTFromContacts(SP::SiconosVector q1)
{
  DEBUG_BEGIN("NewtonEuler5DR::RFC3DcomputeJachqTFromContacts(SP::SiconosVector q1)\n");
  double Nx = _Nc->getValue(0);
  double Ny = _Nc->getValue(1);
  double Nz = _Nc->getValue(2);
  double Px = _Pc1->getValue(0);
  double Py = _Pc1->getValue(1);
  double Pz = _Pc1->getValue(2);
  double G1x = q1->getValue(0);
  double G1y = q1->getValue(1);
  double G1z = q1->getValue(2);

  DEBUG_PRINT("contact normal:\n");
  DEBUG_EXPR(_Nc->display(););
  DEBUG_PRINTF("_Nc->norm2() -1.0 = %e\n",_Nc->norm2()-1.0);
  DEBUG_PRINT("contact point :\n");
  DEBUG_EXPR(_Pc1->display(););
  DEBUG_PRINT("center of mass :\n");
  DEBUG_EXPR(q1->display(););

  assert(_Nc->norm2() >0.0
         && std::abs(_Nc->norm2()-1.0) < 1e-6
         && "NewtonEuler5DR::RFC3DcomputeJachqTFromContacts. Normal vector not consistent ") ;

  double t[6];
  double * pt = t;

  // 1 - Construction of the local contact frame from the normal vector
  
  if (orthoBaseFromVector(&Nx, &Ny, &Nz, pt, pt + 1, pt + 2, pt + 3, pt + 4, pt + 5))
    RuntimeException::selfThrow("NewtonEuler5DR::RFC3DcomputeJachqTFromContacts. Problem in calling orthoBaseFromVector");

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
  _NPG1->zero();
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
  computeRotationMatrix(q1,_rotationBodyToAbsoluteFrame);
  DEBUG_EXPR(_rotationBodyToAbsoluteFrame->display(););

  // 5 - compose the body lever arm matrix with the rotation matrix 
  prod(*_NPG1, *_rotationBodyToAbsoluteFrame, *_AUX1, true);
  DEBUG_EXPR(_rotationBodyToAbsoluteFrame->display(););
  DEBUG_EXPR(_AUX1->display(););

  // 6 -  Rotate the resulting matric in the contact frame
  prod(*_rotationAbsoluteToContactFrame, *_AUX1, *_AUX2, true);
  DEBUG_EXPR(_rotationAbsoluteToContactFrame->display(););
  DEBUG_EXPR(_AUX2->display(););


  // 7 - fill the Jacobian
  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      _jachqT->setValue(ii, jj, _rotationAbsoluteToContactFrame->getValue(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      _jachqT->setValue(ii, jj, _AUX2->getValue(ii, jj - 3));

  prod(*_rotationAbsoluteToContactFrame, *_rotationBodyToAbsoluteFrame, *_AUX2, true);
  DEBUG_EXPR(_AUX2->display(););

  for (unsigned int ii = 3; ii < 5; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      _jachqT->setValue(ii, jj, _AUX2->getValue(ii-2, jj-3));

  DEBUG_EXPR(_jachqT->display(););

  // DEBUG_EXPR_WE(
  //   SP::SimpleMatrix jaux(new SimpleMatrix(*_jachqT));
  //   jaux->trans();
  //   SP::SiconosVector v(new SiconosVector(3));
  //   SP::SiconosVector vRes(new SiconosVector(6));
  //   v->zero();
  //   v->setValue(0, 1);
  //   prod(*jaux, *v, *vRes, true);
  //   vRes->display();
  //   v->zero();
  //   v->setValue(1, 1);
  //   prod(*jaux, *v, *vRes, true);
  //   vRes->display();
  //   v->zero();
  //   v->setValue(2, 1);
  //   prod(*jaux, *v, *vRes, true);
  //   vRes->display();
  //   );
  DEBUG_END("NewtonEuler5DR::RFC3DcomputeJachqTFromContacts(SP::SiconosVector q1)\n");
  //getchar();
}

void NewtonEuler5DR::RFC3DcomputeJachqTFromContacts(SP::SiconosVector q1, SP::SiconosVector q2)
{
  DEBUG_BEGIN("NewtonEuler5DR::RFC3DcomputeJachqTFromContacts(SP::SiconosVector q1, SP::SiconosVector q2)\n");
  double Nx = _Nc->getValue(0);
  double Ny = _Nc->getValue(1);
  double Nz = _Nc->getValue(2);
  double Px = _Pc1->getValue(0);
  double Py = _Pc1->getValue(1);
  double Pz = _Pc1->getValue(2);
  double G1x = q1->getValue(0);
  double G1y = q1->getValue(1);
  double G1z = q1->getValue(2);
  double G2x = q2->getValue(0);
  double G2y = q2->getValue(1);
  double G2z = q2->getValue(2);


  DEBUG_PRINT("contact normal:\n");
  DEBUG_EXPR(_Nc->display(););
  DEBUG_PRINT("contact point :\n");
  DEBUG_EXPR(_Pc1->display(););
  DEBUG_PRINT("center of mass :\n");
  DEBUG_EXPR(q1->display(););


  double t[6];
  double * pt = t;
  if(orthoBaseFromVector(&Nx, &Ny, &Nz, pt, pt + 1, pt + 2, pt + 3, pt + 4, pt + 5))
    RuntimeException::selfThrow("NewtonEuler5DR::RFC3DcomputeJachqTFromContacts. Problem in calling orthoBaseFromVector");
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

  _NPG1->zero();

  (*_NPG1)(0, 0) = 0;
  (*_NPG1)(0, 1) = -(G1z - Pz);
  (*_NPG1)(0, 2) = (G1y - Py);
  (*_NPG1)(1, 0) = (G1z - Pz);
  (*_NPG1)(1, 1) = 0;
  (*_NPG1)(1, 2) = -(G1x - Px);
  (*_NPG1)(2, 0) = -(G1y - Py);
  (*_NPG1)(2, 1) = (G1x - Px);
  (*_NPG1)(2, 2) = 0;

  _NPG2->zero();

  (*_NPG2)(0, 0) = 0;
  (*_NPG2)(0, 1) = -(G2z - Pz);
  (*_NPG2)(0, 2) = (G2y - Py);
  (*_NPG2)(1, 0) = (G2z - Pz);
  (*_NPG2)(1, 1) = 0;
  (*_NPG2)(1, 2) = -(G2x - Px);
  (*_NPG2)(2, 0) = -(G2y - Py);
  (*_NPG2)(2, 1) = (G2x - Px);
  (*_NPG2)(2, 2) = 0;





  computeRotationMatrix(q1,_rotationBodyToAbsoluteFrame);
  prod(*_NPG1, *_rotationBodyToAbsoluteFrame, *_AUX1, true);
  prod(*_rotationAbsoluteToContactFrame, *_AUX1, *_AUX2, true);


  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      _jachqT->setValue(ii, jj, _rotationAbsoluteToContactFrame->getValue(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      _jachqT->setValue(ii, jj, _AUX2->getValue(ii, jj - 3));

  prod(*_rotationAbsoluteToContactFrame, *_rotationBodyToAbsoluteFrame, *_AUX2, true);
  DEBUG_EXPR(_AUX2->display(););

  for (unsigned int ii = 3; ii < 5; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      _jachqT->setValue(ii, jj, _AUX2->getValue(ii-2, jj-3));



  computeRotationMatrix(q2,_rotationBodyToAbsoluteFrame);
  prod(*_NPG2, *_rotationBodyToAbsoluteFrame, *_AUX1, true);
  prod(*_rotationAbsoluteToContactFrame, *_AUX1, *_AUX2, true);

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      _jachqT->setValue(ii, jj + 6, -_rotationAbsoluteToContactFrame->getValue(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      _jachqT->setValue(ii, jj + 6, -_AUX2->getValue(ii, jj - 3));

  prod(*_rotationAbsoluteToContactFrame, *_rotationBodyToAbsoluteFrame, *_AUX2, true);
  DEBUG_EXPR(_AUX2->display(););

  for (unsigned int ii = 3; ii < 5; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      _jachqT->setValue(ii, jj, _AUX2->getValue(ii-2, jj-3));

  DEBUG_EXPR(_jachqT->display(););

  DEBUG_END("NewtonEuler5DR::RFC3DcomputeJachqTFromContacts(SP::SiconosVector q1, SP::SiconosVector q2)\n");

}

void NewtonEuler5DR::computeJachqT(Interaction& inter, SP::BlockVector q0)
{
  DEBUG_BEGIN("NewtonEuler5DR::computeJachqT(Interaction& inter,  SP::BlockVector q0)\n");
  if (q0->numberOfBlocks()>1)
  {
    RFC3DcomputeJachqTFromContacts((q0->getAllVect())[0], (q0->getAllVect())[1]);
  }
  else
  {
    RFC3DcomputeJachqTFromContacts((q0->getAllVect())[0]);
  }
  DEBUG_END("NewtonEuler5DR::computeJachqT(Interaction& inter,  SP::BlockVector q0)\n");


}
