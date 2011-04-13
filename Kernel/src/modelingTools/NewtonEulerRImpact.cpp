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


#include "NewtonEulerRImpact.hpp"
#include <boost/math/quaternion.hpp>
#include "NewtonEulerDS.hpp"
//#define NERI_DEBUG
using namespace std;

//#define NEFC3D_DEBUG

void NewtonEulerRImpact::NIcomputeJachqTFromContacts(SP::SimpleVector Pc, SP::SimpleVector Nc, SP::SiconosVector G1, SP::SiconosMatrix jhqT)
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
  _Mabs_C->setValue(0, 0, Nx);
  _Mabs_C->setValue(0, 1, Ny);
  _Mabs_C->setValue(0, 2, Nz);

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


  for (unsigned int jj = 0; jj < 3; jj++)
    jhqT->setValue(0, jj, _Mabs_C->getValue(0, jj));

  for (unsigned int jj = 3; jj < 6; jj++)
    jhqT->setValue(0, jj, _AUX2->getValue(0, jj - 3));

#ifdef NEFC3D_DEBUG
  printf("NewtonEulerRImpact jhqt\n");
  jhqT->display();
#endif
}

void NewtonEulerRImpact::NIcomputeJachqTFromContacts(SP::SimpleVector Pc, SP::SimpleVector Nc, SP::SiconosVector G1, SP::SiconosVector G2, SP::SiconosMatrix jhqT)
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

  _Mabs_C->setValue(0, 0, Nx);
  _Mabs_C->setValue(0, 1, Ny);
  _Mabs_C->setValue(0, 2, Nz);

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



  for (unsigned int jj = 0; jj < 3; jj++)
    jhqT->setValue(0, jj, _Mabs_C->getValue(0, jj));


  for (unsigned int jj = 3; jj < 6; jj++)
    jhqT->setValue(0, jj, _AUX2->getValue(0, jj - 3));

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

  for (unsigned int jj = 0; jj < 3; jj++)
    jhqT->setValue(0, jj + 6, -_Mabs_C->getValue(0, jj));

  for (unsigned int jj = 3; jj < 6; jj++)
    jhqT->setValue(0, jj + 6, -_AUX2->getValue(0, jj - 3));
}
void  NewtonEulerRImpact::initComponents()
{
  NewtonEulerR::initComponents();
  _jachqProj.reset(new SimpleMatrix(_jachq->size(0), _jachq->size(1)));
  _Mabs_C.reset(new SimpleMatrix(1, 3));
  _AUX1.reset(new SimpleMatrix(3, 3));
  _AUX2.reset(new SimpleMatrix(1, 3));
  _NPG1.reset(new SimpleMatrix(3, 3));
  _NPG2.reset(new SimpleMatrix(3, 3));
  _Mobj1_abs.reset(new SimpleMatrix(3, 3));
  _Mobj2_abs.reset(new SimpleMatrix(3, 3));

}

void NewtonEulerRImpact::computeJachq(double t)
{
  DSIterator itDS = interaction()->dynamicalSystemsBegin();
  SP::DynamicalSystem aux = *itDS;
  //assert (&(*aux)==&(*_ds1));
  itDS++;

  bool has2Bodies = false;
  if (itDS != interaction()->dynamicalSystemsEnd())
    has2Bodies = true;
  _jachq->setValue(0, 0, _Nc->getValue(0));
  _jachq->setValue(0, 1, _Nc->getValue(1));
  _jachq->setValue(0, 2, _Nc->getValue(2));
  if (has2Bodies)
  {
    _jachq->setValue(0, 7, -_Nc->getValue(0));
    _jachq->setValue(0, 8, -_Nc->getValue(1));
    _jachq->setValue(0, 9, -_Nc->getValue(2));
  }
  SP::BlockVector BlockX = boost::static_pointer_cast<BlockVector>((data[q0]));
  for (int iDS = 0; iDS < 2; iDS++)
  {
    if (!has2Bodies && iDS == 1)
      continue;
    double sign = 1.0;
    SP::SiconosVector q = (BlockX->getAllVect())[iDS];
#ifdef NERI_DEBUG
    printf("NewtonEulerRImpact::computeJachq : ds%d->q :", iDS);
    q->display();
#endif
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
      //cout<<"NewtonEulerRImpact::computeJachq sign is -1 \n";
      ::boost::math::quaternion<double>    quatAux(0, _Pc2->getValue(0) - q->getValue(0), _Pc2->getValue(1) - q->getValue(1),
          _Pc2->getValue(2) - q->getValue(2));
      quatGP = quatAux;
    }
#ifdef NERI_DEBUG
    printf("NewtonEulerRImpact::computeJachq :GP :%lf, %lf, %lf\n", quatGP.R_component_2(), quatGP.R_component_3(), quatGP.R_component_4());
    printf("NewtonEulerRImpact::computeJachq :Q :%e,%e, %e, %e\n", q->getValue(3), q->getValue(4), q->getValue(5), q->getValue(6));
#endif
    ::boost::math::quaternion<double>    quatQ(q->getValue(3), q->getValue(4), q->getValue(5), q->getValue(6));
    ::boost::math::quaternion<double>    quatcQ(q->getValue(3), -q->getValue(4), -q->getValue(5), -q->getValue(6));
    ::boost::math::quaternion<double>    quat0(1, 0, 0, 0);
    ::boost::math::quaternion<double>    quatBuff;
    ::boost::math::quaternion<double>    _2qiquatGP;
    _2qiquatGP = quatGP;
    _2qiquatGP *= 2 * (q->getValue(3));
    quatBuff = (quatGP * quatQ) + (quatcQ * quatGP) - _2qiquatGP;
#ifdef NERI_DEBUG
    printf("NewtonEulerRImpact::computeJachq :quattBuuf : %e,%e,%e \n", quatBuff.R_component_2(), quatBuff.R_component_3(), quatBuff.R_component_4());
#endif
    _jachq->setValue(0, 7 * iDS + 3, sign * (quatBuff.R_component_2()*_Nc->getValue(0) +
                     quatBuff.R_component_3()*_Nc->getValue(1) + quatBuff.R_component_4()*_Nc->getValue(2)));
    //cout<<"WARNING NewtonEulerRImpact set jachq \n";
    //_jachq->setValue(0,7*iDS+3,0);
    for (int i = 1; i < 4; i++)
    {
      ::boost::math::quaternion<double>    quatei(0, (i == 1) ? 1 : 0, (i == 2) ? 1 : 0, (i == 3) ? 1 : 0);
      _2qiquatGP = quatGP;
      _2qiquatGP *= 2 * (q->getValue(3 + i));
      quatBuff = quatei * quatcQ * quatGP - quatGP * quatQ * quatei - _2qiquatGP;
      _jachq->setValue(0, 7 * iDS + 3 + i, sign * (quatBuff.R_component_2()*_Nc->getValue(0) +
                       quatBuff.R_component_3()*_Nc->getValue(1) + quatBuff.R_component_4()*_Nc->getValue(2)));
    }
  }
  //cout<<"WARNING NewtonEulerRImpact set jachq to zedro \n";
  //_jachq->setValue(0,4,0);_jachq->setValue(0,5,0);_jachq->setValue(0,6,0);
  *_jachqProj = *_jachq;
  //_jachqProj->setValue(0,4,0);_jachqProj->setValue(0,5,0);_jachqProj->setValue(0,6,0);
#ifdef NERI_DEBUG
  printf("NewtonEulerRImpact::computeJachq :");
  _jachq->display();
#endif

}
void NewtonEulerRImpact::computeJachqT()
{
  DSIterator itDS = interaction()->dynamicalSystemsBegin();
  SP::NewtonEulerDS d1 =  boost::static_pointer_cast<NewtonEulerDS> (*itDS);
  SP::SiconosVector Q1 = d1->q();
  itDS++;
  if (itDS != interaction()->dynamicalSystemsEnd())
  {
    SP::NewtonEulerDS d2 =  boost::static_pointer_cast<NewtonEulerDS> (*itDS);
    SP::SiconosVector Q2 = d2->q();
    NIcomputeJachqTFromContacts(_Pc1, _Nc, Q1, Q2, _jachqT);
  }
  else
  {
    NIcomputeJachqTFromContacts(_Pc1, _Nc, Q1, _jachqT);
  }
}
