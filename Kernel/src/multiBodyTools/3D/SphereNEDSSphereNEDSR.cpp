/* Siconos-Example version 3.0.0, Copyright INRIA 2005-2010.
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
 * Foundation, Inc., 51 Franklin St, Fifth FLOOR, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 *
 */

#include <math.h>
#include "SphereNEDSSphereNEDSR.hpp"

SphereNEDSSphereNEDSR::SphereNEDSSphereNEDSR(double r,
    double rr)
  : NewtonEulerRFC3D()
{
  r1 = r;
  r2 = rr;
  r1pr2 = r1 + r2;
}

double SphereNEDSSphereNEDSR::distance(double x1, double y1, double z1, double r1,
                                       double x2, double y2, double z2, double r2)
{
  double dx = x1 - x2;
  double dy = y1 - y2;
  double dz = z1 - z2;

  return (sqrt(dx * dx + dy * dy + dz * dz) - r1pr2);
}


void SphereNEDSSphereNEDSR::computeh(double)
{


  double q_0 = (*data[q0])(0);
  double q_1 = (*data[q0])(1);
  double q_2 = (*data[q0])(2);
  double q_7 = (*data[q0])(7);
  double q_8 = (*data[q0])(8);
  double q_9 = (*data[q0])(9);

  SP::SiconosVector y = interaction()->y(0);

  y->setValue(0, distance(q_0, q_1, q_2, r1, q_7, q_8, q_9, r2));
  _Pc->setValue(0, (r1 * q_0 + r2 * q_7) / (r1 + r2));
  _Pc->setValue(1, (r1 * q_1 + r2 * q_8) / (r1 + r2));
  _Pc->setValue(2, (r1 * q_2 + r2 * q_9) / (r1 + r2));
  _Nc->setValue(0, (q_0 - q_7) / (y->getValue(0) + r1pr2));
  _Nc->setValue(1, (q_1 - q_8) / (y->getValue(0) + r1pr2));
  _Nc->setValue(2, (q_2 - q_9) / (y->getValue(0) + r1pr2));
  //std::cout<<" SphereNEDSSphereNEDSR::computeh dist:"<<y->getValue(0)<<"\n";
  //std::cout<<"_Pc:\n";
  //_Pc->display();
  //std::cout<<"_Nc:\n";
  //_Nc->display();
};


