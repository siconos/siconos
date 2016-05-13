/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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


#include "NewtonEulerFrom3DLocalFrameR.hpp"
#include "NewtonEulerDS.hpp"
#include <boost/math/quaternion.hpp>
#include "Interaction.hpp"
#include "BlockVector.hpp"

#include "op3x3.h"

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

/*
See devNotes.pdf for details. A detailed documentation is available in DevNotes.pdf: chapter 'NewtonEulerR: computation of \nabla q H'. Subsection 'Case FC3D: using the local frame local velocities'
*/
void NewtonEulerFrom3DLocalFrameR::initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  NewtonEulerFrom1DLocalFrameR::initComponents(inter, DSlink, workV, workM);
  unsigned int qSize = 7 * (inter.getSizeOfDS() / 6);
  /*keep only the distance.*/
  _jachq.reset(new SimpleMatrix(3, qSize));


  _Mabs_C.reset(new SimpleMatrix(3, 3));
  _AUX2.reset(new SimpleMatrix(3, 3));
  //  _isContact=1;
}
void NewtonEulerFrom3DLocalFrameR::FC3DcomputeJachqTFromContacts(SP::SiconosVector q1)
{

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
  DEBUG_PRINT("contact point :\n");
  DEBUG_EXPR(_Pc1->display(););
  DEBUG_PRINT("center of mass :\n");
  DEBUG_EXPR(q1->display(););

  assert(_Nc->norm2() >0.0 && "NewtonEulerFrom3DLocalFrameR::FC3DcomputeJachqTFromContacts. Normal vector not consistent ") ;

  double t[6];
  double * pt = t;
  orthoBaseFromVector(&Nx, &Ny, &Nz, pt, pt + 1, pt + 2, pt + 3, pt + 4, pt + 5);
  pt = t;
  _Mabs_C->setValue(0, 0, Nx);
  _Mabs_C->setValue(1, 0, *pt);
  _Mabs_C->setValue(2, 0, *(pt + 3));
  _Mabs_C->setValue(0, 1, Ny);
  _Mabs_C->setValue(1, 1, *(pt + 1));
  _Mabs_C->setValue(2, 1, *(pt + 4));
  _Mabs_C->setValue(0, 2, Nz);
  _Mabs_C->setValue(1, 2, *(pt + 2));
  _Mabs_C->setValue(2, 2, *(pt + 5));

  DEBUG_PRINT("_Mabs_C:\n");
  DEBUG_EXPR(_Mabs_C->display(););

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

//  d1->computeMObjToAbs();
//  SimpleMatrix& Mobj1_abs = *d1->MObjToAbs();

  computeMObjToAbs(q1,_MObjToAbs);

  DEBUG_PRINT("NewtonEulerFrom3DLocalFrameR::FC3DcomputeJachqTFromContacts, _MObjToAbs:");
  DEBUG_EXPR(_MObjToAbs->display(););



  prod(*_NPG1, *_MObjToAbs, *_AUX1, true);

  DEBUG_PRINT("NewtonEulerFrom3DLocalFrameR::FC3DcomputeJachqTFromContacts, *_MObjToAbs:");
  DEBUG_EXPR_WE(_MObjToAbs->display(););
  DEBUG_PRINT("NewtonEulerFrom3DLocalFrameR::FC3DcomputeJachqTFromContacts, AUX1:");
  DEBUG_EXPR_WE(_AUX1->display(););
  DEBUG_PRINT("NewtonEulerFrom3DLocalFrameR::FC3DcomputeJachqTFromContacts, AUX2:");
  DEBUG_EXPR(_AUX2->display(););


  prod(*_Mabs_C, *_AUX1, *_AUX2, true);

  DEBUG_PRINT("NewtonEulerFrom3DLocalFrameR::FC3DcomputeJachqTFromContacts, Mabs_C:");
  DEBUG_EXPR_WE(_Mabs_C->display(););
  DEBUG_PRINT("NewtonEulerFrom3DLocalFrameR::FC3DcomputeJachqTFromContacts, AUX1:");
  DEBUG_EXPR_WE(_AUX1->display(););
  DEBUG_PRINT("NewtonEulerFrom3DLocalFrameR::FC3DcomputeJachqTFromContacts, AUX2:");
  DEBUG_EXPR_WE(_AUX2->display(););


  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      _jachqT->setValue(ii, jj, _Mabs_C->getValue(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      _jachqT->setValue(ii, jj, _AUX2->getValue(ii, jj - 3));



  DEBUG_PRINT("NewtonEulerFrom3DLocalFrameR::FC3DcomputeJachqTFromContacts, _jahcqT:\n");
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
}

void NewtonEulerFrom3DLocalFrameR::FC3DcomputeJachqTFromContacts(SP::SiconosVector q1, SP::SiconosVector q2)
{
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
  orthoBaseFromVector(&Nx, &Ny, &Nz, pt, pt + 1, pt + 2, pt + 3, pt + 4, pt + 5);
  pt = t;
  _Mabs_C->setValue(0, 0, Nx);
  _Mabs_C->setValue(1, 0, *pt);
  _Mabs_C->setValue(2, 0, *(pt + 3));
  _Mabs_C->setValue(0, 1, Ny);
  _Mabs_C->setValue(1, 1, *(pt + 1));
  _Mabs_C->setValue(2, 1, *(pt + 4));
  _Mabs_C->setValue(0, 2, Nz);
  _Mabs_C->setValue(1, 2, *(pt + 2));
  _Mabs_C->setValue(2, 2, *(pt + 5));

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




//  d1->computeMObjToAbs();
//  SimpleMatrix& Mobj1_abs = *d1->MObjToAbs();

  computeMObjToAbs(q1,_MObjToAbs);

  prod(*_NPG1, *_MObjToAbs, *_AUX1, true);
  prod(*_Mabs_C, *_AUX1, *_AUX2, true);


  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      _jachqT->setValue(ii, jj, _Mabs_C->getValue(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      _jachqT->setValue(ii, jj, _AUX2->getValue(ii, jj - 3));

//  d2->computeMObjToAbs();
//  SimpleMatrix& Mobj2_abs = *d2->MObjToAbs();


  computeMObjToAbs(q2,_MObjToAbs);

  prod(*_NPG2, *_MObjToAbs, *_AUX1, true);
  prod(*_Mabs_C, *_AUX1, *_AUX2, true);

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      _jachqT->setValue(ii, jj + 6, -_Mabs_C->getValue(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      _jachqT->setValue(ii, jj + 6, -_AUX2->getValue(ii, jj - 3));

}

void NewtonEulerFrom3DLocalFrameR::computeJachqT(Interaction& inter, SP::BlockVector q0)
{
  DEBUG_BEGIN("NewtonEulerFrom3DLocalFrameR::computeJachqT(Interaction& inter,  SP::BlockVector q0)\n");
  if (q0->getNumberOfBlocks()>1)
  {
    FC3DcomputeJachqTFromContacts((q0->getAllVect())[0], (q0->getAllVect())[1]);
  }
  else
  {
    FC3DcomputeJachqTFromContacts((q0->getAllVect())[0]);
  }
  DEBUG_END("NewtonEulerFrom3DLocalFrameR::computeJachqT(Interaction& inter,  SP::BlockVector q0)\n");


}
