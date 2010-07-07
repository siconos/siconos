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
  : NewtonEulerR()
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

};

void SphereNEDSSphereNEDSR::computeJachq(double)
{

  double r, nN, nU;
  /* u ^ v  = n */
  double u1, u2, u3, v1, v2, v3, n1, n2, n3, r1u1, r1u2, r1u3, r1v1, r1v2, r1v3, r2u1, r2u2, r2u3, r2v1, r2v2, r2v3;

  SimpleMatrix *g = (SimpleMatrix *)_jachq.get();

  double q_0 = (*data[q0])(0);
  double q_1 = (*data[q0])(1);
  double q_2 = (*data[q0])(2);

  double q_7 = (*data[q0])(7);
  double q_8 = (*data[q0])(8);
  double q_9 = (*data[q0])(9);

  double a1 = (*data[q0])(3);
  double b1 = (*data[q0])(4);
  double c1 = (*data[q0])(5);
  double d1 = (*data[q0])(6);

  double a2 = (*data[q0])(10);
  double b2 = (*data[q0])(11);
  double c2 = (*data[q0])(12);
  double d2 = (*data[q0])(13);


  double A = q_0 - q_7;
  double B = q_1 - q_8;
  double C = q_2 - q_9;

  double x0 = pow(B, 2);
  double x1 = pow(A, 2);
  double x2 = pow(C, 2);
  double x3 = x0 + x1 + x2;
  double x4 = pow(x3, (-1.0 / 2.0));
  double x5 = B - C;
  double x6 = A + C;
  double x7 = A + B;
  double x8 = pow(x7, 2);
  double x9 = pow(x6, 2);
  double x10 = pow(x5, 2);
  double x11 = x10 + x8 + x9;
  double x12 = pow(x11, (-1.0 / 2.0));
  double x13 = -x6;
  double x14 = 2 * A * r1 * x12 * x4 * x7;
  double x15 = 2 * C * r1 * x12 * x13 * x4;
  double x16 = -2 * B * r1 * x12 * x4 * x7;
  double x17 = x15 + x16;
  double x18 = 2 * B * r1 * x12 * x4 * x5;
  double x19 = 2 * C * r1 * x12 * x4 * x5;
  double x20 = x14 - x19;
  double x21 = -2 * A * r1 * x12 * x13 * x4;
  double x22 = x18 + x21;
  double x23 = -2 * B * r2 * x12 * x4 * x7;
  double x24 = 2 * C * r2 * x12 * x13 * x4;
  double x25 = x23 + x24;
  double x26 = 2 * A * r2 * x12 * x4 * x7;
  double x27 = 2 * C * r2 * x12 * x4 * x5;
  double x28 = x26 - x27;
  double x29 = -2 * A * r2 * x12 * x13 * x4;
  double x30 = 2 * B * r2 * x12 * x4 * x5;
  double x31 = x29 + x30;
  double x32 = C * x12 * x4 * x5;
  double x33 = A * x12 * x13 * x4;
  double x34 = B * x12 * x4 * x5;
  double x35 = x33 - x34;
  double x36 = -A * x12 * x4 * x7;
  double x37 = x32 + x36;
  double x38 = B * x12 * x4 * x7;
  double x39 = -C * x12 * x13 * x4;
  double x40 = x38 + x39;
  double x41 = 2 * A * r1 * x35 * x4;
  double x42 = 2 * C * r1 * x4 * x40;
  double x43 = x41 - x42;
  double x44 = 2 * C * r1 * x37 * x4;
  double x45 = -2 * B * r1 * x35 * x4;
  double x46 = x44 + x45;
  double x47 = -2 * A * r1 * x37 * x4;
  double x48 = 2 * B * r1 * x4 * x40;
  double x49 = x47 + x48;
  double x50 = C * x12 * x13 * x4;
  double x51 = x38 - x50;
  double x52 = A * x12 * x4 * x7;
  double x53 = x32 - x52;
  double x54 = -x34;
  double x55 = x33 + x54;
  double x56 = -2 * C * r2 * x37 * x4;
  double x57 = 2 * A * r2 * x37 * x4;
  double x58 = 2 * C * r2 * x4 * x40;
  double x59 = 2 * A * r2 * x35 * x4;
  double x60 = x58 - x59;
  double x61 = 2 * B * r2 * x35 * x4;
  double x62 = x56 + x61;
  double x63 = -2 * B * r2 * x4 * x40;
  double x64 = x57 + x63;
  (*g)(0, 0) = A * x4;
  (*g)(0, 1) = B * x4;
  (*g)(0, 2) = C * x4;
  (*g)(0, 3) = 0;
  (*g)(0, 4) = 0;
  (*g)(0, 5) = 0;
  (*g)(0, 6) = 0;
  (*g)(0, 7) = -A * x4;
  (*g)(0, 8) = -B * x4;
  (*g)(0, 9) = -C * x4;
  (*g)(0, 10) = 0;
  (*g)(0, 11) = 0;
  (*g)(0, 12) = 0;
  (*g)(0, 13) = 0;
  (*g)(1, 0) = x12 * x5;
  (*g)(1, 1) = x12 * x13;
  (*g)(1, 2) = x12 * x7;
  (*g)(1, 3) = a1 * x17 - b1 * x20 - c1 * x22;
  (*g)(1, 4) = a1 * (x14 - 2 * r1 * x32) + b1 * x17 - d1 * (x18 - 2 * r1 * x33);
  (*g)(1, 5) = a1 * x22 + c1 * x17 + d1 * x20;
  (*g)(1, 6) = b1 * x22 + d1 * x17 - c1 * x20;
  (*g)(1, 7) = -x12 * x5;
  (*g)(1, 8) = x12 * x6;
  (*g)(1, 9) = x12 * (-A - B);
  (*g)(1, 10) = a2 * x25 - b2 * x28 - c2 * x31;
  (*g)(1, 11) = a2 * x28 + b2 * x25 - d2 * x31;
  (*g)(1, 12) = a2 * x31 + c2 * x25 + d2 * x28;
  (*g)(1, 13) = b2 * x31 + d2 * x25 - c2 * x28;
  (*g)(2, 0) = x51;
  (*g)(2, 1) = x53;
  (*g)(2, 2) = x35;
  (*g)(2, 3) = a1 * x46 - b1 * x43 - c1 * x49;
  (*g)(2, 4) = a1 * x43 + b1 * x46 - d1 * x49;
  (*g)(2, 5) = a1 * x49 + c1 * (x44 - 2 * B * r1 * x35 * x4) + d1 * (x41 - 2 * C * r1 * x4 * x40);
  (*g)(2, 6) = b1 * x49 + d1 * x46 - c1 * x43;
  (*g)(2, 7) = x51;
  (*g)(2, 8) = x53;
  (*g)(2, 9) = x35;
  (*g)(2, 10) = a2 * x62 - b2 * x60 - c2 * x64;
  (*g)(2, 11) = a2 * (-2 * A * r2 * x4 * x55 + 2 * C * r2 * x4 * x51) + b2 * (x56 + 2 * B * r2 * x4 * x55) - d2 * (x57 - 2 * B * r2 * x4 * x51);
  (*g)(2, 12) = a2 * x64 + c2 * x62 + d2 * x60;
  (*g)(2, 13) = b2 * x64 + d2 * (x61 - 2 * C * r2 * x37 * x4) - c2 * (x58 - 2 * A * r2 * x35 * x4);
};
