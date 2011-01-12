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



void NIcomputeJachqTFromContacts(SP::SimpleVector Pc, SP::SimpleVector Nc, SP::SiconosVector G1, SP::SiconosMatrix jhqT)
{
  SP::SimpleMatrix M(new SimpleMatrix(1, 3));
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
  M->setValue(0, 0, Nx);
  M->setValue(0, 1, Ny);
  M->setValue(0, 2, Nz);
  //  M->display();
  SP::SimpleMatrix N(new SimpleMatrix(3, 6));
  N->zero();
  (*N)(0, 0) = 1;
  (*N)(1, 1) = 1;
  (*N)(2, 2) = 1;
  (*N)(0, 3) = 0;
  (*N)(0, 4) = Pz - G1z;
  (*N)(0, 5) = -(Py - G1y);
  (*N)(1, 3) = -(Pz - G1z);
  (*N)(1, 4) = 0;
  (*N)(1, 5) = Px - G1x;
  (*N)(2, 3) = Py - G1y;
  (*N)(2, 4) = -(Px - G1x);
  (*N)(2, 5) = 0;

  prod(*M, *N, *jhqT, true);

}

void NIcomputeJachqTFromContacts(SP::SimpleVector Pc, SP::SimpleVector Nc, SP::SiconosVector G1, SP::SiconosVector G2, SP::SiconosMatrix jhqT)
{
  SP::SimpleMatrix M(new SimpleMatrix(1, 3));
  double Nx = Nc->getValue(0);
  double Ny = Nc->getValue(1);
  double Nz = Nc->getValue(2);
  double Px = Pc->getValue(0);
  double Py = Pc->getValue(1);
  double Pz = Pc->getValue(2);
  double G1x = G1->getValue(0);
  double G1y = G1->getValue(1);
  double G1z = G1->getValue(2);
  double G2x = G2->getValue(0);
  double G2y = G2->getValue(1);
  double G2z = G2->getValue(2);
  M->setValue(0, 0, Nx);
  M->setValue(0, 1, Ny);
  M->setValue(0, 2, Nz);
  //cout<<"M display\n";
  //M->display();
  SP::SimpleMatrix N(new SimpleMatrix(3, 12));
  N->zero();
  (*N)(0, 0) = 1;
  (*N)(1, 1) = 1;
  (*N)(2, 2) = 1;
  (*N)(0, 3) = 0;
  (*N)(0, 4) = Pz - G1z;
  (*N)(0, 5) = -(Py - G1y);
  (*N)(1, 3) = -(Pz - G1z);
  (*N)(1, 4) = 0;
  (*N)(1, 5) = Px - G1x;
  (*N)(2, 3) = Py - G1y;
  (*N)(2, 4) = -(Px - G1x);
  (*N)(2, 5) = 0;

  (*N)(0, 6) = -1;
  (*N)(1, 7) = -1;
  (*N)(2, 8) = -1;
  (*N)(0, 9) = 0;
  (*N)(0, 10) = -(Pz - G2z);
  (*N)(0, 11) = (Py - G2y);
  (*N)(1, 9) = (Pz - G2z);
  (*N)(1, 10) = 0;
  (*N)(1, 11) = -(Px - G2x);
  (*N)(2, 9) = -(Py - G2y);
  (*N)(2, 10) = (Px - G2x);
  (*N)(2, 11) = 0;
  prod(*M, *N, *jhqT, true);
  //cout<<"jhqt\n";
  //jhqT->display();

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
    ::boost::math::quaternion<float>    quatGP;
    if (iDS == 0)
    {
      ::boost::math::quaternion<float>    quatAux(0, _Pc1->getValue(0) - q->getValue(0), _Pc1->getValue(1) - q->getValue(1),
          _Pc1->getValue(2) - q->getValue(2));
      quatGP = quatAux;
    }
    else
    {
      sign = -1.0;
      //cout<<"NewtonEulerRImpact::computeJachq sign is -1 \n";
      ::boost::math::quaternion<float>    quatAux(0, _Pc2->getValue(0) - q->getValue(0), _Pc2->getValue(1) - q->getValue(1),
          _Pc2->getValue(2) - q->getValue(2));
      quatGP = quatAux;
    }
#ifdef NERI_DEBUG
    printf("NewtonEulerRImpact::computeJachq :GP :%lf, %lf, %lf\n", quatGP.R_component_2(), quatGP.R_component_3(), quatGP.R_component_4());
    printf("NewtonEulerRImpact::computeJachq :Q :%e,%e, %e, %e\n", q->getValue(3), q->getValue(4), q->getValue(5), q->getValue(6));
#endif
    ::boost::math::quaternion<float>    quatQ(q->getValue(3), q->getValue(4), q->getValue(5), q->getValue(6));
    ::boost::math::quaternion<float>    quatcQ(q->getValue(3), -q->getValue(4), -q->getValue(5), -q->getValue(6));
    ::boost::math::quaternion<float>    quat0(1, 0, 0, 0);
    ::boost::math::quaternion<float>    quatBuff;
    quatBuff = (quatGP * quatQ) + (quatcQ * quatGP);
#ifdef NERI_DEBUG
    printf("NewtonEulerRImpact::computeJachq :quattBuuf : %e,%e,%e \n", quatBuff.R_component_2(), quatBuff.R_component_3(), quatBuff.R_component_4());
#endif
    _jachq->setValue(0, 7 * iDS + 3, sign * (quatBuff.R_component_2()*_Nc->getValue(0) +
                     quatBuff.R_component_3()*_Nc->getValue(1) + quatBuff.R_component_4()*_Nc->getValue(2)));
    //cout<<"WARNING NewtonEulerRImpact set jachq \n";
    //_jachq->setValue(0,7*iDS+3,0);
    for (int i = 1; i < 4; i++)
    {
      ::boost::math::quaternion<float>    quatei(0, (i == 1) ? 1 : 0, (i == 2) ? 1 : 0, (i == 3) ? 1 : 0);
      quatBuff = quatei * quatcQ * quatGP - quatGP * quatQ * quatei;
      _jachq->setValue(0, 7 * iDS + 3 + i, sign * (quatBuff.R_component_2()*_Nc->getValue(0) +
                       quatBuff.R_component_3()*_Nc->getValue(1) + quatBuff.R_component_4()*_Nc->getValue(2)));
    }
  }
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
