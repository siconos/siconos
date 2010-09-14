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
#include "SphereNEDSPlanR.hpp"

SphereNEDSPlanR::SphereNEDSPlanR(double r, double A, double B, double C, double D)
  : NewtonEulerRFC3D(), r(r), A(A), B(B), C(C), D(D)
{
  nN = sqrt(A * A + B * B + C * C);

  n1 = A / nN;
  n2 = B / nN;
  n3 = C / nN;



}

double SphereNEDSPlanR::distance(double x, double y, double z, double rad)
{

  return (fabs(A * x + B * y + C * z + D) / nN - rad);
}


void SphereNEDSPlanR::computeh(double)
{

  double q_0 = (*data[q0])(0);
  double q_1 = (*data[q0])(1);
  double q_2 = (*data[q0])(2);

  SP::SiconosVector y = interaction()->y(0);

  y->setValue(0, distance(q_0, q_1, q_2, r));
  _Pc->setValue(0, q_0 - r * n1);
  _Pc->setValue(1, q_1 - r * n2);
  _Pc->setValue(2, q_2 - r * n3);
  _Nc->setValue(0, n1);
  _Nc->setValue(1, n2);
  _Nc->setValue(2, n3);
};
