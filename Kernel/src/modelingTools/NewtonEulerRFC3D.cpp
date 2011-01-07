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
using namespace std;
/*
  m \dot V = Fext + R
  I \dot \homega + \homega I \homega = Mext + R*PG ,
  with * vectoriel product, R reaction in the globla frame. P the point of contact.
  r is the reaction in the local frame.
  M^t r=R
            |nx t1x t2x|
  with M^t =|ny t1y t2y|
            |nz t1z t2z|


  |  R | |1 0 0     |
  |    |=|0 1 0     |
  |R*PG| |0 0 1     |*R=N^t *R = N^t * M^t r
         |0 -PGz PGy|
         |PGz 0 -PGx|
         |-PGy PGx 0|


so jachqT=M*N



 */
void NewtonEulerRFC3D::initComponents()
{
  NewtonEulerRImpact::initComponents();
  /*keep only the distance.*/
  /*Warning, in current version, user of FC3D has to set _y and _yProj in the computeh */
  _yProj.reset(new SimpleVector(1));
}
void computeJachqTFromContacts(SP::SimpleVector Pc, SP::SimpleVector Nc, SP::SiconosVector G1, SP::SiconosMatrix jhqT)
{
  SP::SimpleMatrix M(new SimpleMatrix(3, 3));
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
  M->setValue(0, 0, Nx);
  M->setValue(0, 1, Ny);
  M->setValue(0, 2, Nz);
  M->setValue(1, 0, *pt);
  M->setValue(1, 1, *(pt + 1));
  M->setValue(1, 2, *(pt + 2));
  M->setValue(2, 0, *(pt + 3));
  M->setValue(2, 1, *(pt + 4));
  M->setValue(2, 2, *(pt + 5));
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
  //  cout<<"jhqt\n";
  //  jhqT->display();

}

void computeJachqTFromContacts(SP::SimpleVector Pc, SP::SimpleVector Nc, SP::SiconosVector G1, SP::SiconosVector G2, SP::SiconosMatrix jhqT)
{
  SP::SimpleMatrix M(new SimpleMatrix(3, 3));
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
  double t[6];
  double * pt = t;
  //cout<<"computeJachqTFromContacts Nx Ny Nz"<<Nx<<" "<<Ny<<" "<<Nz<<endl;
  orthoBaseFromVector(&Nx, &Ny, &Nz, pt, pt + 1, pt + 2, pt + 3, pt + 4, pt + 5);
  pt = t;
  M->setValue(0, 0, Nx);
  M->setValue(0, 1, Ny);
  M->setValue(0, 2, Nz);
  M->setValue(1, 0, *pt);
  M->setValue(1, 1, *(pt + 1));
  M->setValue(1, 2, *(pt + 2));
  M->setValue(2, 0, *(pt + 3));
  M->setValue(2, 1, *(pt + 4));
  M->setValue(2, 2, *(pt + 5));
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
    computeJachqTFromContacts(_Pc1, _Nc, Q1, Q2, _jachqT);
  }
  else
  {
    computeJachqTFromContacts(_Pc1, _Nc, Q1, _jachqT);
  }
}
