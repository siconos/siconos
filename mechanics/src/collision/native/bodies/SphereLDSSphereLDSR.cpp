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

#include <cmath>
#include "SphereLDSSphereLDSR.hpp"
#include <BlockVector.hpp>
#include "SimpleMatrix.hpp"

#include <op3x3.h>

SphereLDSSphereLDSR::SphereLDSSphereLDSR(double r, double rr) : LagrangianScleronomousR()
{
  r1 = r;
  r2 = rr;
  r1pr2 = r1 + r2;
}

double SphereLDSSphereLDSR::distance(double x1, double y1, double z1, double r1, double x2, double y2, double z2, double r2)
{
  double dx = x1 - x2;
  double dy = y1 - y2;
  double dz = z1 - z2;

  return (sqrt(dx * dx + dy * dy + dz * dz) - r1pr2);
}


void SphereLDSSphereLDSR::computeh(SiconosVector& q, SiconosVector& z, SiconosVector& y)
{

  y.setValue(0, distance(q(0), q(1), q(2), r1,
                         q(6), q(7), q(8), r2));
  y.setValue(1, 0.);

};

void SphereLDSSphereLDSR::computeJachq(SiconosVector& q, SiconosVector& z)
{

  double A, B, C;
  /* u ^ v  = n */
  double u1, u2, u3, v1, v2, v3, n1, n2, n3;

  SimpleMatrix *g = (SimpleMatrix *)_jachq.get();

  double q_0 = q(0);
  double q_1 = q(1);
  double q_2 = q(2);

  double theta1 = q(3);
  double phi1 = q(4);


  double q_6 = q(6);
  double q_7 = q(7);
  double q_8 = q(8);


  double theta2 = q(9);
  double phi2 = q(10);


  A = -(q_6 - q_0);
  B = -(q_7 - q_1);
  C = -(q_8 - q_2);

  n1 = A;
  n2 = B;
  n3 = C;

  if (orthoBaseFromVector(&n1, &n2, &n3, &u1, &u2, &u3, &v1, &v2, &v3))
        RuntimeException::selfThrow("SphereLDSSphereLDSR::computeJachq. Problem in calling orthoBaseFromVector");

  double x0 = 2 * n1 * r1 * u2;
  double x1 = 2 * n2 * r1 * u1;
  double x2 = x0 - x1;
  double x3 = sin(phi1);
  double x4 = 2 * n2 * r1 * u3;
  double x5 = 2 * n3 * r1 * u2;
  double x6 = x4 - x5;
  double x7 = cos(phi1);
  double x8 = sin(theta1);
  double x9 = 2 * n3 * r1 * u1;
  double x10 = 2 * n1 * r1 * u3;
  double x11 = x9 - x10;
  double x12 = 2 * n2 * r2 * u1;
  double x13 = 2 * n1 * r2 * u2;
  double x14 = x12 - x13;
  double x15 = sin(phi2);
  double x16 = 2 * n3 * r2 * u2;
  double x17 = 2 * n2 * r2 * u3;
  double x18 = x16 - x17;
  double x19 = cos(phi2);
  double x20 = sin(theta2);
  double x21 = 2 * n1 * r2 * u3;
  double x22 = 2 * n3 * r2 * u1;
  double x23 = x21 - x22;
  double x24 = cos(theta1);
  double x25 = 2 * n1 * r1 * v2;
  double x26 = 2 * n2 * r1 * v1;
  double x27 = x25 - x26;
  double x28 = 2 * n2 * r1 * v3;
  double x29 = 2 * n3 * r1 * v2;
  double x30 = x28 - x29;
  double x31 = 2 * n3 * r1 * v1;
  double x32 = 2 * n1 * r1 * v3;
  double x33 = x31 - x32;
  double x34 = cos(theta2);
  double x35 = 2 * n1 * r2 * v2;
  double x36 = 2 * n2 * r2 * v1;
  double x37 = x35 - x36;
  double x38 = 2 * n2 * r2 * v3;
  double x39 = 2 * n3 * r2 * v2;
  double x40 = x38 - x39;
  double x41 = 2 * n3 * r2 * v1;
  double x42 = 2 * n1 * r2 * v3;
  double x43 = x41 - x42;
  double x44 = x3 * x8;
  double x45 = x7 * x8;
  double x46 = x19 * x20;
  double x47 = x15 * x20;
  double x48 = x18 * x47;
  double x49 = x30 * x44;
  double x50 = x33 * x45;
  double x51 = x44 * x6;
  double x52 = x11 * x45;
  double x53 = x40 * x47;
  double x54 = x23 * x46;
  double x55 = x43 * x46;
  double x56 = x2 * x24;
  double x57 = -u2;
  double x58 = x6 * x7;
  double x59 = x14 * x34;
  double x60 = -n3;
  double x61 = -u1;
  double x62 = -n1;
  double x63 = x24 * x27;
  double x64 = x34 * x37;
  double x65 = x11 * x3;
  double x66 = x3 * x33;
  double x67 = -u3;
  double x68 = x19 * x40;
  double x69 = x30 * x7;
  double x70 = -n2;
  double x71 = x15 * x23;
  double x72 = x15 * x43;
  double x73 = x18 * x19;
  double x74 = x48 + x59;
  double x75 = x49 + x63;
  double x76 = x66 + x69;
  double x77 = x68 + x72;
  double x78 = x71 + x73;
  double x79 = x51 + x56;
  double x80 = x58 + x65;
  double x81 = x53 + x64;
  (*g)(0, 0) = n1;
  (*g)(0, 1) = n2;
  (*g)(0, 2) = n3;
  (*g)(0, 3) = 0;
  (*g)(0, 4) = 0;
  (*g)(0, 5) = 0;
  (*g)(0, 6) = x62;
  (*g)(0, 7) = x70;
  (*g)(0, 8) = x60;
  (*g)(0, 9) = 0;
  (*g)(0, 10) = 0;
  (*g)(0, 11) = 0;
  (*g)(1, 0) = u1;
  (*g)(1, 1) = u2;
  (*g)(1, 2) = u3;
  (*g)(1, 3) = x80;
  (*g)(1, 4) = x2;
  (*g)(1, 5) = x79 - x52;
  (*g)(1, 6) = x61;
  (*g)(1, 7) = x57;
  (*g)(1, 8) = x67;
  (*g)(1, 9) = x78;
  (*g)(1, 10) = x14;
  (*g)(1, 11) = x74 - x54;
  (*g)(2, 0) = v1;
  (*g)(2, 1) = v2;
  (*g)(2, 2) = v3;
  (*g)(2, 3) = x76;
  (*g)(2, 4) = x27;
  (*g)(2, 5) = x75 - x50;
  (*g)(2, 6) = v1;
  (*g)(2, 7) = v2;
  (*g)(2, 8) = v3;
  (*g)(2, 9) = x77;
  (*g)(2, 10) = x37;
  (*g)(2, 11) = x81 - x55;
};
