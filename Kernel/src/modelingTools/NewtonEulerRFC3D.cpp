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
}
void FC3DcomputeJachqTFromContacts(SP::SimpleVector Pc, SP::SimpleVector Nc, SP::SiconosVector G1, SP::SiconosMatrix jhqT)
{
  /*Matrix converting the contact coordinate to the absolute coordinate*/
  SP::SimpleMatrix MC_abs(new SimpleMatrix(3, 3));
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
  MC_abs->setValue(0, 0, Nx);
  MC_abs->setValue(0, 1, *pt);
  MC_abs->setValue(0, 2, *(pt + 3));
  MC_abs->setValue(1, 0, Ny);
  MC_abs->setValue(1, 1, *(pt + 1));
  MC_abs->setValue(1, 2, *(pt + 4));
  MC_abs->setValue(2, 0, Nz);
  MC_abs->setValue(2, 1, *(pt + 2));
  MC_abs->setValue(2, 2, *(pt + 5));
  SP::SimpleMatrix Mabs_C(new SimpleMatrix(*MC_abs));
  Mabs_C->trans();
  SP::SimpleMatrix NPG(new SimpleMatrix(3, 3));
  NPG->zero();

  (*NPG)(0, 0) = 0;
  (*NPG)(0, 1) = -(G1z - Pz);
  (*NPG)(0, 2) = (G1y - Py);
  (*NPG)(1, 0) = (G1z - Pz);
  (*NPG)(1, 1) = 0;
  (*NPG)(1, 2) = -(G1x - Px);
  (*NPG)(2, 0) = -(G1y - Py);
  (*NPG)(2, 1) = (G1x - Px);
  (*NPG)(2, 2) = 0;



  double q0 = G1->getValue(3);
  double q1 = G1->getValue(4);
  double q2 = G1->getValue(5);
  double q3 = G1->getValue(6);

  ::boost::math::quaternion<double>    quatQ(q0, q1, q2, q3);
  ::boost::math::quaternion<double>    quatcQ(q0, -q1, -q2, -q3);
  ::boost::math::quaternion<double>    quatx(0, 1, 0, 0);
  ::boost::math::quaternion<double>    quaty(0, 0, 1, 0);
  ::boost::math::quaternion<double>    quatz(0, 0, 0, 1);
  SP::SimpleMatrix Mabs_obj(new SimpleMatrix(3, 3));
  ::boost::math::quaternion<double>    quatBuff;
  quatBuff = quatcQ * quatx * quatQ;
  Mabs_obj->setValue(0, 0, quatBuff.R_component_2());
  Mabs_obj->setValue(1, 0, quatBuff.R_component_3());
  Mabs_obj->setValue(2, 0, quatBuff.R_component_4());
  quatBuff = quatcQ * quaty * quatQ;
  Mabs_obj->setValue(0, 1, quatBuff.R_component_2());
  Mabs_obj->setValue(1, 1, quatBuff.R_component_3());
  Mabs_obj->setValue(2, 1, quatBuff.R_component_4());
  quatBuff = quatcQ * quatz * quatQ;
  Mabs_obj->setValue(0, 2, quatBuff.R_component_2());
  Mabs_obj->setValue(1, 2, quatBuff.R_component_3());
  Mabs_obj->setValue(2, 2, quatBuff.R_component_4());
  SP::SimpleMatrix Mobj_abs(new SimpleMatrix(*Mabs_obj));
  Mobj_abs->trans();

#ifdef NEFC3D_DEBUG
  printf("first line of res %e %e %e\n", Mabs_obj->getValue(0, 0), Mabs_obj->getValue(0, 1), Mabs_obj->getValue(0, 2));
  printf("coord de Rabs_obj:\n");
  Mabs_obj->display();
#endif

  SP::SimpleMatrix AUX1(new SimpleMatrix(3, 3));
  SP::SimpleMatrix AUX2(new SimpleMatrix(3, 3));
  prod(*NPG, *Mobj_abs, *AUX1, true);
  prod(*Mabs_C, *AUX1, *AUX2, true);


  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      jhqT->setValue(ii, jj, Mabs_C->getValue(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      jhqT->setValue(ii, jj, AUX2->getValue(ii, jj - 3));

}

void FC3DcomputeJachqTFromContacts(SP::SimpleVector Pc, SP::SimpleVector Nc, SP::SiconosVector G1, SP::SiconosVector G2, SP::SiconosMatrix jhqT)
{
  /*Matrix converting the contact coordinate to the absolute coordinate*/
  SP::SimpleMatrix MC_abs(new SimpleMatrix(3, 3));
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
  MC_abs->setValue(0, 0, Nx);
  MC_abs->setValue(0, 1, *pt);
  MC_abs->setValue(0, 2, *(pt + 3));
  MC_abs->setValue(1, 0, Ny);
  MC_abs->setValue(1, 1, *(pt + 1));
  MC_abs->setValue(1, 2, *(pt + 4));
  MC_abs->setValue(2, 0, Nz);
  MC_abs->setValue(2, 1, *(pt + 2));
  MC_abs->setValue(2, 2, *(pt + 5));
  SP::SimpleMatrix Mabs_C(new SimpleMatrix(*MC_abs));
  Mabs_C->trans();
  SP::SimpleMatrix NPG(new SimpleMatrix(3, 3));
  NPG->zero();

  (*NPG)(0, 0) = 0;
  (*NPG)(0, 1) = -(G1z - Pz);
  (*NPG)(0, 2) = (G1y - Py);
  (*NPG)(1, 0) = (G1z - Pz);
  (*NPG)(1, 1) = 0;
  (*NPG)(1, 2) = -(G1x - Px);
  (*NPG)(2, 0) = -(G1y - Py);
  (*NPG)(2, 1) = (G1x - Px);
  (*NPG)(2, 2) = 0;



  double q0 = G1->getValue(3);
  double q1 = G1->getValue(4);
  double q2 = G1->getValue(5);
  double q3 = G1->getValue(6);

  ::boost::math::quaternion<double>    quatQ(q0, q1, q2, q3);
  ::boost::math::quaternion<double>    quatcQ(q0, -q1, -q2, -q3);
  ::boost::math::quaternion<double>    quatx(0, 1, 0, 0);
  ::boost::math::quaternion<double>    quaty(0, 0, 1, 0);
  ::boost::math::quaternion<double>    quatz(0, 0, 0, 1);
  SP::SimpleMatrix Mabs_obj(new SimpleMatrix(3, 3));
  ::boost::math::quaternion<double>    quatBuff;
  quatBuff = quatcQ * quatx * quatQ;
  Mabs_obj->setValue(0, 0, quatBuff.R_component_2());
  Mabs_obj->setValue(1, 0, quatBuff.R_component_3());
  Mabs_obj->setValue(2, 0, quatBuff.R_component_4());
  quatBuff = quatcQ * quaty * quatQ;
  Mabs_obj->setValue(0, 1, quatBuff.R_component_2());
  Mabs_obj->setValue(1, 1, quatBuff.R_component_3());
  Mabs_obj->setValue(2, 1, quatBuff.R_component_4());
  quatBuff = quatcQ * quatz * quatQ;
  Mabs_obj->setValue(0, 2, quatBuff.R_component_2());
  Mabs_obj->setValue(1, 2, quatBuff.R_component_3());
  Mabs_obj->setValue(2, 2, quatBuff.R_component_4());
  SP::SimpleMatrix Mobj_abs(new SimpleMatrix(*Mabs_obj));
  Mobj_abs->trans();

#ifdef NEFC3D_DEBUG
  printf("first line of res %e %e %e\n", Mabs_obj->getValue(0, 0), Mabs_obj->getValue(0, 1), Mabs_obj->getValue(0, 2));
  printf("coord de Rabs_obj:\n");
  Mabs_obj->display();
#endif

  SP::SimpleMatrix AUX1(new SimpleMatrix(3, 3));
  SP::SimpleMatrix AUX2(new SimpleMatrix(3, 3));
  prod(*NPG, *Mobj_abs, *AUX1, true);
  prod(*Mabs_C, *AUX1, *AUX2, true);


  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      jhqT->setValue(ii, jj, Mabs_C->getValue(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      jhqT->setValue(ii, jj, AUX2->getValue(ii, jj - 3));

  double G2x = G2->getValue(0);
  double G2y = G2->getValue(1);
  double G2z = G2->getValue(2);

  SP::SimpleMatrix NPG2(new SimpleMatrix(3, 3));
  NPG2->zero();
  (*NPG2)(0, 0) = 0;
  (*NPG2)(0, 1) = -(G2z - Pz);
  (*NPG2)(0, 2) = (G2y - Py);
  (*NPG2)(1, 0) = (G2z - Pz);
  (*NPG2)(1, 1) = 0;
  (*NPG2)(1, 2) = -(G2x - Px);
  (*NPG2)(2, 0) = -(G2y - Py);
  (*NPG2)(2, 1) = (G2x - Px);
  (*NPG2)(2, 2) = 0;

  SP::SimpleMatrix Mabs_obj2(new SimpleMatrix(3, 3));
  q0 = G2->getValue(3);
  q1 = G2->getValue(4);
  q2 = G2->getValue(5);
  q3 = G2->getValue(6);

  ::boost::math::quaternion<double>    quatQ2(q0, q1, q2, q3);
  ::boost::math::quaternion<double>    quatcQ2(q0, -q1, -q2, -q3);
  quatBuff = quatcQ2 * quatx * quatQ2;
  Mabs_obj2->setValue(0, 0, quatBuff.R_component_2());
  Mabs_obj2->setValue(1, 0, quatBuff.R_component_3());
  Mabs_obj2->setValue(2, 0, quatBuff.R_component_4());
  quatBuff = quatcQ2 * quaty * quatQ2;
  Mabs_obj2->setValue(0, 1, quatBuff.R_component_2());
  Mabs_obj2->setValue(1, 1, quatBuff.R_component_3());
  Mabs_obj2->setValue(2, 1, quatBuff.R_component_4());
  quatBuff = quatcQ2 * quatz * quatQ2;
  Mabs_obj2->setValue(0, 2, quatBuff.R_component_2());
  Mabs_obj2->setValue(1, 2, quatBuff.R_component_3());
  Mabs_obj2->setValue(2, 2, quatBuff.R_component_4());
  SP::SimpleMatrix Mobj2_abs(new SimpleMatrix(*Mabs_obj2));
  Mobj2_abs->trans();
  prod(*NPG2, *Mobj2_abs, *AUX1, true);
  prod(*Mabs_C, *AUX1, *AUX2, true);

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 0; jj < 3; jj++)
      jhqT->setValue(ii, jj + 6, -Mabs_C->getValue(ii, jj));

  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int jj = 3; jj < 6; jj++)
      jhqT->setValue(ii, jj + 6, -AUX2->getValue(ii, jj - 3));

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
