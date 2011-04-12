/* Siconos-Kernel, Copyright INRIA 2005-2010.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */


#include "NewtonEulerRFC3D.hpp"
#include "NewtonEulerDS.hpp"
#include <boost/math/quaternion.hpp>
using namespace std;
/*
  See devNotes.pdf for details.
 */
void NewtonEulerRFC3D::initComponents()
{
  NewtonEulerRImpact::initComponents();
  /*keep only the distance.*/
  /*Warning, in current version, user of FC3D has to set _y and _yProj in the computeh */
  _yProj.reset(new SimpleVector(1));
  _Mabs_C.reset(new SimpleMatrix(3, 3));
  _AUX2.reset(new SimpleMatrix(3, 3));
}
void NewtonEulerRFC3D::FC3DcomputeJachqTFromContacts(SP::SimpleVector Pc, SP::SimpleVector Nc, SP::SiconosVector G1, SP::SiconosMatrix jhqT)
{

  double Nx = Nc->getValue(0);
  double Ny = Nc->getValue(1);
  double Nz = Nc->getValue(2);
  double Px = Pc->getValue(0);
  double Py = Pc->getValue(1);
  double Pz = Pc->getValue(2);
  double G1x = G1->getValue(0);
  double G1y = G1->getValue(1);
  double G1z = G1->getValue(2);
#ifdef NEFC3D_DEBUG
  printf("contact normal:\n");
  Nc->display();
  printf("point de contact :\n");
  Pc->display();
  printf("center of masse :\n");
  G1->display();
#endif
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



  double q0 = G1->getValue(3);
  double q1 = G1->getValue(4);
  double q2 = G1->getValue(5);
  double q3 = G1->getValue(6);

  ::boost::math::quaternion<double>    quatQ(q0, q1, q2, q3);
  ::boost::math::quaternion<double>    quatcQ(q0, -q1, -q2, -q3);
  ::boost::math::quaternion<double>    quatx(0, 1, 0, 0);
  ::boost::math::quaternion<double>    quaty(0, 0, 1, 0);
  ::boost::math::quaternion<double>    quatz(0, 0, 0, 1);


  ::boost::math::quaternion<double>    quatBuff;
  quatBuff = quatcQ * quatx * quatQ;
  _Mobj1_abs->setValue(0, 0, quatBuff.R_component_2());
  _Mobj1_abs->setValue(0, 1, quatBuff.R_component_3());
  _Mobj1_abs->setValue(0, 2, quatBuff.R_component_4());
  quatBuff = quatcQ * quaty * quatQ;
  _Mobj1_abs->setValue(1, 0, quatBuff.R_component_2());
  _Mobj1_abs->setValue(1, 1, quatBuff.R_component_3());
  _Mobj1_abs->setValue(1, 2, quatBuff.R_component_4());
  quatBuff = quatcQ * quatz * quatQ;
  _Mobj1_abs->setValue(2, 0, quatBuff.R_component_2());
  _Mobj1_abs->setValue(2, 1, quatBuff.R_component_3());
  _Mobj1_abs->setValue(2, 2, quatBuff.R_component_4());

  prod(*_NPG1, *_Mobj1_abs, *_AUX1, true);
  prod(*_Mabs_C, *_AUX1, *_AUX2, true);


  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      jhqT->setValue(ii, jj, _Mabs_C->getValue(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      jhqT->setValue(ii, jj, _AUX2->getValue(ii, jj - 3));

}

void NewtonEulerRFC3D::FC3DcomputeJachqTFromContacts(SP::SimpleVector Pc, SP::SimpleVector Nc, SP::SiconosVector G1, SP::SiconosVector G2, SP::SiconosMatrix jhqT)
{
  double Nx = Nc->getValue(0);
  double Ny = Nc->getValue(1);
  double Nz = Nc->getValue(2);
  double Px = Pc->getValue(0);
  double Py = Pc->getValue(1);
  double Pz = Pc->getValue(2);
  double G1x = G1->getValue(0);
  double G1y = G1->getValue(1);
  double G1z = G1->getValue(2);

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



  double q0 = G1->getValue(3);
  double q1 = G1->getValue(4);
  double q2 = G1->getValue(5);
  double q3 = G1->getValue(6);

  ::boost::math::quaternion<double>    quatQ(q0, q1, q2, q3);
  ::boost::math::quaternion<double>    quatcQ(q0, -q1, -q2, -q3);
  ::boost::math::quaternion<double>    quatx(0, 1, 0, 0);
  ::boost::math::quaternion<double>    quaty(0, 0, 1, 0);
  ::boost::math::quaternion<double>    quatz(0, 0, 0, 1);
  ::boost::math::quaternion<double>    quatBuff;
  quatBuff = quatcQ * quatx * quatQ;
  _Mobj1_abs->setValue(0, 0, quatBuff.R_component_2());
  _Mobj1_abs->setValue(0, 1, quatBuff.R_component_3());
  _Mobj1_abs->setValue(0, 2, quatBuff.R_component_4());
  quatBuff = quatcQ * quaty * quatQ;
  _Mobj1_abs->setValue(1, 0, quatBuff.R_component_2());
  _Mobj1_abs->setValue(1, 1, quatBuff.R_component_3());
  _Mobj1_abs->setValue(1, 2, quatBuff.R_component_4());
  quatBuff = quatcQ * quatz * quatQ;
  _Mobj1_abs->setValue(2, 0, quatBuff.R_component_2());
  _Mobj1_abs->setValue(2, 1, quatBuff.R_component_3());
  _Mobj1_abs->setValue(2, 2, quatBuff.R_component_4());


  prod(*_NPG1, *_Mobj1_abs, *_AUX1, true);
  prod(*_Mabs_C, *_AUX1, *_AUX2, true);


  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      jhqT->setValue(ii, jj, _Mabs_C->getValue(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      jhqT->setValue(ii, jj, _AUX2->getValue(ii, jj - 3));

  double G2x = G2->getValue(0);
  double G2y = G2->getValue(1);
  double G2z = G2->getValue(2);

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

  q0 = G2->getValue(3);
  q1 = G2->getValue(4);
  q2 = G2->getValue(5);
  q3 = G2->getValue(6);

  ::boost::math::quaternion<double>    quatQ2(q0, q1, q2, q3);
  ::boost::math::quaternion<double>    quatcQ2(q0, -q1, -q2, -q3);
  quatBuff = quatcQ2 * quatx * quatQ2;

  _Mobj2_abs->setValue(0, 0, quatBuff.R_component_2());
  _Mobj2_abs->setValue(0, 1, quatBuff.R_component_3());
  _Mobj2_abs->setValue(0, 2, quatBuff.R_component_4());
  quatBuff = quatcQ2 * quaty * quatQ2;
  _Mobj2_abs->setValue(1, 0, quatBuff.R_component_2());
  _Mobj2_abs->setValue(1, 1, quatBuff.R_component_3());
  _Mobj2_abs->setValue(1, 2, quatBuff.R_component_4());
  quatBuff = quatcQ2 * quatz * quatQ2;
  _Mobj2_abs->setValue(2, 0, quatBuff.R_component_2());
  _Mobj2_abs->setValue(2, 1, quatBuff.R_component_3());
  _Mobj2_abs->setValue(2, 2, quatBuff.R_component_4());

  prod(*_NPG2, *_Mobj2_abs, *_AUX1, true);
  prod(*_Mabs_C, *_AUX1, *_AUX2, true);

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      jhqT->setValue(ii, jj + 6, -_Mabs_C->getValue(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      jhqT->setValue(ii, jj + 6, -_AUX2->getValue(ii, jj - 3));

}
void NewtonEulerRFC3D::computeJachqT()
{
  DSIterator itDS = interaction()->dynamicalSystemsBegin();
  SP::NewtonEulerDS d1 =  boost::static_pointer_cast<NewtonEulerDS> (*itDS);
  SP::SiconosVector Q1 = d1->q();
  itDS++;
  if (itDS != interaction()->dynamicalSystemsEnd())
  {
    SP::NewtonEulerDS d2 =  boost::static_pointer_cast<NewtonEulerDS> (*itDS);
    SP::SiconosVector Q2 = d2->q();
    FC3DcomputeJachqTFromContacts(_Pc1, _Nc, Q1, Q2, _jachqT);
  }
  else
  {
    FC3DcomputeJachqTFromContacts(_Pc1, _Nc, Q1, _jachqT);
  }
}
