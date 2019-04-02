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


#include "NewtonEuler1DR.hpp"
#include "RotationQuaternion.hpp"
#include "Interaction.hpp"
#include "BlockVector.hpp"
#include <boost/math/quaternion.hpp>

//#define NERI_DEBUG

//#define NEFC3D_DEBUG
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"
/*
See devNotes.pdf for details. A detailed documentation is available in DevNotes.pdf: chapter 'NewtonEulerR: computation of \nabla q H'. Subsection 'Case FC3D: using the local frame local velocities'
*/
void NewtonEuler1DR::NIcomputeJachqTFromContacts(SP::SiconosVector q1)
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
#ifdef NEFC3D_DEBUG
  printf("contact normal:\n");
  _Nc->display();
  printf("point de contact :\n");
  _Pc1->display();
  printf("center of masse :\n");
  q1->display();
#endif
  _rotationAbsoluteToContactFrame->setValue(0, 0, Nx);
  _rotationAbsoluteToContactFrame->setValue(0, 1, Ny);
  _rotationAbsoluteToContactFrame->setValue(0, 2, Nz);

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


  computeRotationMatrix(q1,_rotationBodyToAbsoluteFrame);
  prod(*_NPG1, *_rotationBodyToAbsoluteFrame, *_AUX1, true);

  prod(*_rotationAbsoluteToContactFrame, *_AUX1, *_AUX2, true);


  for (unsigned int jj = 0; jj < 3; jj++)
    _jachqT->setValue(0, jj, _rotationAbsoluteToContactFrame->getValue(0, jj));

  for (unsigned int jj = 3; jj < 6; jj++)
    _jachqT->setValue(0, jj, _AUX2->getValue(0, jj - 3));

#ifdef NEFC3D_DEBUG
  printf("NewtonEuler1DR jhqt\n");
  _jachqT->display();
#endif
}

void NewtonEuler1DR::NIcomputeJachqTFromContacts(SP::SiconosVector q1, SP::SiconosVector q2)
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

  _rotationAbsoluteToContactFrame->setValue(0, 0, Nx);
  _rotationAbsoluteToContactFrame->setValue(0, 1, Ny);
  _rotationAbsoluteToContactFrame->setValue(0, 2, Nz);

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

  computeRotationMatrix(q1,_rotationBodyToAbsoluteFrame);
  prod(*_NPG1, *_rotationBodyToAbsoluteFrame, *_AUX1, true);
  prod(*_rotationAbsoluteToContactFrame, *_AUX1, *_AUX2, true);

  for (unsigned int jj = 0; jj < 3; jj++)
    _jachqT->setValue(0, jj, _rotationAbsoluteToContactFrame->getValue(0, jj));


  for (unsigned int jj = 3; jj < 6; jj++)
    _jachqT->setValue(0, jj, _AUX2->getValue(0, jj - 3));

  double G2x = q2->getValue(0);
  double G2y = q2->getValue(1);
  double G2z = q2->getValue(2);

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



  computeRotationMatrix(q2,_rotationBodyToAbsoluteFrame);
  prod(*_NPG2, *_rotationBodyToAbsoluteFrame, *_AUX1, true);

  prod(*_rotationAbsoluteToContactFrame, *_AUX1, *_AUX2, true);

  for (unsigned int jj = 0; jj < 3; jj++)
    _jachqT->setValue(0, jj + 6, -_rotationAbsoluteToContactFrame->getValue(0, jj));

  for (unsigned int jj = 3; jj < 6; jj++)
    _jachqT->setValue(0, jj + 6, -_AUX2->getValue(0, jj - 3));
}

void NewtonEuler1DR::initialize(Interaction& inter)
{
  NewtonEulerR::initialize(inter);
  //proj_with_q  _jachqProj.reset(new SimpleMatrix(_jachq->size(0),_jachq->size(1)));
  unsigned int qSize = 7 * (inter.getSizeOfDS() / 6);
  _jachq.reset(new SimpleMatrix(1, qSize));

  /* VA 12/04/2016 All of what follows should be put in WorkM*/
  _rotationAbsoluteToContactFrame.reset(new SimpleMatrix(1, 3));
  _rotationBodyToAbsoluteFrame.reset(new SimpleMatrix(3, 3));
  _AUX1.reset(new SimpleMatrix(3, 3));
  _AUX2.reset(new SimpleMatrix(1, 3));
  _NPG1.reset(new SimpleMatrix(3, 3));
  _NPG2.reset(new SimpleMatrix(3, 3));
  //  _isContact=1;
}




void NewtonEuler1DR::computeJachq(double time, Interaction& inter, SP::BlockVector q0)
{

  DEBUG_BEGIN("NewtonEuler1DR::computeJachq(double time, Interaction& inter, SP::BlockVector q0 ) \n");
  DEBUG_PRINTF("with time =  %f\n",time);
  DEBUG_PRINTF("with inter =  %p\n",&inter);


  _jachq->setValue(0, 0, _Nc->getValue(0));
  _jachq->setValue(0, 1, _Nc->getValue(1));
  _jachq->setValue(0, 2, _Nc->getValue(2));
  if (inter.has2Bodies())
  {
    _jachq->setValue(0, 7, -_Nc->getValue(0));
    _jachq->setValue(0, 8, -_Nc->getValue(1));
    _jachq->setValue(0, 9, -_Nc->getValue(2));
  }

  for (unsigned int iDS =0 ; iDS < q0->numberOfBlocks()  ; iDS++)
  {
    SP::SiconosVector q = (q0->getAllVect())[iDS];
    double sign = 1.0;
    DEBUG_PRINTF("NewtonEuler1DR::computeJachq : ds%d->q :", iDS);
    DEBUG_EXPR_WE(q->display(););

    ::boost::math::quaternion<double>    quatGP;
    if (iDS == 0)
    {
      ::boost::math::quaternion<double>    quatAux(0, _Pc1->getValue(0) - q->getValue(0), _Pc1->getValue(1) - q->getValue(1),
                                                   _Pc1->getValue(2) - q->getValue(2));
      quatGP = quatAux;
    }
    else
    {
      sign = -1.0;
      //cout<<"NewtonEuler1DR::computeJachq sign is -1 \n";
      ::boost::math::quaternion<double>    quatAux(0, _Pc2->getValue(0) - q->getValue(0), _Pc2->getValue(1) - q->getValue(1),
                                                   _Pc2->getValue(2) - q->getValue(2));
      quatGP = quatAux;
    }
    DEBUG_PRINTF("NewtonEuler1DR::computeJachq :GP :%lf, %lf, %lf\n", quatGP.R_component_2(), quatGP.R_component_3(), quatGP.R_component_4());
    DEBUG_PRINTF("NewtonEuler1DR::computeJachq :Q :%e,%e, %e, %e\n", q->getValue(3), q->getValue(4), q->getValue(5), q->getValue(6));
    ::boost::math::quaternion<double>    quatQ(q->getValue(3), q->getValue(4), q->getValue(5), q->getValue(6));
    ::boost::math::quaternion<double>    quatcQ(q->getValue(3), -q->getValue(4), -q->getValue(5), -q->getValue(6));
    ::boost::math::quaternion<double>    quat0(1, 0, 0, 0);
    ::boost::math::quaternion<double>    quatBuff;
    ::boost::math::quaternion<double>    _2qiquatGP;
    _2qiquatGP = quatGP;
    _2qiquatGP *= 2 * (q->getValue(3));
    quatBuff = (quatGP * quatQ) + (quatcQ * quatGP) - _2qiquatGP;

    DEBUG_PRINTF("NewtonEuler1DR::computeJachq :quattBuuf : %e,%e,%e \n", quatBuff.R_component_2(), quatBuff.R_component_3(), quatBuff.R_component_4());

    _jachq->setValue(0, 7 * iDS + 3, sign * (quatBuff.R_component_2()*_Nc->getValue(0) +
                                             quatBuff.R_component_3()*_Nc->getValue(1) + quatBuff.R_component_4()*_Nc->getValue(2)));
    //cout<<"WARNING NewtonEuler1DR set jachq \n";
    //_jachq->setValue(0,7*iDS+3,0);
    for (unsigned int i = 1; i < 4; i++)
    {
      ::boost::math::quaternion<double>    quatei(0, (i == 1) ? 1 : 0, (i == 2) ? 1 : 0, (i == 3) ? 1 : 0);
      _2qiquatGP = quatGP;
      _2qiquatGP *= 2 * (q->getValue(3 + i));
      quatBuff = quatei * quatcQ * quatGP - quatGP * quatQ * quatei - _2qiquatGP;
      _jachq->setValue(0, 7 * iDS + 3 + i, sign * (quatBuff.R_component_2()*_Nc->getValue(0) +
                                                   quatBuff.R_component_3()*_Nc->getValue(1) + quatBuff.R_component_4()*_Nc->getValue(2)));
    }
  }

  DEBUG_EXPR(_jachq->display(););
  DEBUG_END("NewtonEuler1DR::computeJachq(double time, Interaction& inter, SP::BlockVector q0 \n");

}

void NewtonEuler1DR::computeJachqT(Interaction& inter, SP::BlockVector q0 )
{
  DEBUG_BEGIN("NewtonEuler1DR::computeJachqT(Interaction& inter, SP::BlockVector q0 \n")

  if (q0->numberOfBlocks()>1)
  {
    NIcomputeJachqTFromContacts((q0->getAllVect())[0], (q0->getAllVect())[1]);
  }
  else
  {
    NIcomputeJachqTFromContacts((q0->getAllVect())[0]);
  }

  DEBUG_END("NewtonEuler1DR::computeJachqT(Interaction& inter, SP::BlockVector q0) \n");

}

double NewtonEuler1DR::distance() const
{
  SiconosVector dpc(*_Pc2 - *_Pc1);
  return dpc.norm2() * (inner_prod(*_Nc, dpc) >= 0 ? -1 : 1);
}

void NewtonEuler1DR::computeh(double time, BlockVector& q0,
                                            SiconosVector &y)
{
  // Contact points and normal are stored as relative to q1 and q2, if
  // no q2 then pc2 and normal are absolute.

  // Update pc1 based on q0 and relPc1
  SP::SiconosVector q1 = (q0.getAllVect())[0];
  ::boost::math::quaternion<double> qq1((*q1)(3), (*q1)(4), (*q1)(5), (*q1)(6));
  ::boost::math::quaternion<double> qpc1(0,(*_relPc1)(0),(*_relPc1)(1),(*_relPc1)(2));

  // apply q1 rotation and add
  qpc1 = qq1 * qpc1 / qq1;
  (*_Pc1)(0) = qpc1.R_component_2() + (*q1)(0);
  (*_Pc1)(1) = qpc1.R_component_3() + (*q1)(1);
  (*_Pc1)(2) = qpc1.R_component_4() + (*q1)(2);

  if (q0.numberOfBlocks() > 1)
  {
    // Update pc2 based on q0 and relPc2
    SP::SiconosVector q2 = (q0.getAllVect())[1];
    ::boost::math::quaternion<double> qq2((*q2)(3), (*q2)(4), (*q2)(5), (*q2)(6));
    ::boost::math::quaternion<double> qpc2(0,(*_relPc2)(0),(*_relPc2)(1),(*_relPc2)(2));

    // apply q2 rotation and add
    qpc2 = qq2 * qpc2 / qq2;
    (*_Pc2)(0) = qpc2.R_component_2() + (*q2)(0);
    (*_Pc2)(1) = qpc2.R_component_3() + (*q2)(1);
    (*_Pc2)(2) = qpc2.R_component_4() + (*q2)(2);

    // same for normal
    ::boost::math::quaternion<double> qnc(0, (*_relNc)(0), (*_relNc)(1), (*_relNc)(2));
    qnc = qq2 * qnc / qq2;
    (*_Nc)(0) = qnc.R_component_2();
    (*_Nc)(1) = qnc.R_component_3();
    (*_Nc)(2) = qnc.R_component_4();
  }
  else
  {
    *_Pc2 = *_relPc2;
    *_Nc = *_relNc;
  }

  NewtonEulerR::computeh(time, q0, y);
}
