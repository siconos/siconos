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
#include "SphereLDSPlanR.hpp"
#include <BlockVector.hpp>
#include "SimpleMatrix.hpp"

#include <op3x3.h>

SphereLDSPlanR::SphereLDSPlanR(double r, double A, double B, double C, double D)
  : LagrangianScleronomousR(), r(r), A(A), B(B), C(C), D(D)
{
  n1 = A;
  n2 = B;
  n3 = C;

  nN = sqrt(A * A + B * B + C * C);

  if (orthoBaseFromVector(&n1, &n2, &n3, &u1, &u2, &u3, &v1, &v2, &v3))
      RuntimeException::selfThrow("SphereLDSPlanR::SphereLDSPlanR. Problem in calling orthoBaseFromVector");
  // r*u & r *v

  ru1 = r * u1;
  ru2 = r * u2;
  ru3 = r * u3;

  rv1 = r * v1;
  rv2 = r * v2;
  rv3 = r * v3;

}

double SphereLDSPlanR::distance(double x, double y, double z, double rad)
{

  return (fabs(A * x + B * y + C * z + D) / nN - rad);
}


void SphereLDSPlanR::computeh(SiconosVector& q, SiconosVector& z, SiconosVector& y)
{

  double q_0 = q(0);
  double q_1 = q(1);
  double q_2 = q(2);

  y.setValue(0, distance(q_0, q_1, q_2, r));

};

void normalize(SP::SiconosVector, unsigned int);

void SphereLDSPlanR::computeJachq(SiconosVector& q, SiconosVector& z)
{
  SimpleMatrix *g = (SimpleMatrix *)_jachq.get();

  double theta = q(3);
  double phi   = q(4);

  double cthe = cos(theta);
  double sthe = sin(theta);
  double cphi = cos(phi);
  double sphi = sin(phi);

  (*g)(0, 0) = n1;
  (*g)(1, 0) = u1;
  (*g)(2, 0) = v1;
  (*g)(0, 1) = n2;
  (*g)(1, 1) = u2;
  (*g)(2, 1) = v2;
  (*g)(0, 2) = n3;
  (*g)(1, 2) = u3;
  (*g)(2, 2) = v3;
  (*g)(0, 3) = 0;
  (*g)(1, 3) = -rv1 * cphi - rv2 * sphi;
  (*g)(2, 3) = ru1 * cphi + ru2 * sphi;
  (*g)(0, 4) = 0;
  (*g)(1, 4) = -rv3;
  (*g)(2, 4) = ru3;
  (*g)(0, 5) = 0;
  (*g)(1, 5) = -rv3 * cthe + rv2 * cphi * sthe - rv1 * sphi * sthe;
  (*g)(2, 5) = ru3 * cthe + ru1 * sphi * sthe - ru2 * cphi * sthe;

}

