/* Siconos-Example version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 *
 */

#include <math.h>
#include "SphereSphereR.h"

SphereSphereR::SphereSphereR(double r, double rr) : LagrangianScleronomousR()
{
  r1 = r;
  r2 = rr;
  r1pr2 = r1 + r2;
}

double SphereSphereR::distance(double x1, double y1, double z1, double r1, double x2, double y2, double z2, double r2)
{
  double dx = x1 - x2;
  double dy = y1 - y2;
  double dz = z1 - z2;

  return (sqrt(dx * dx + dy * dy + dz * dz) - r1pr2);
}


void SphereSphereR::computeH(double)
{

  // Warning: temporary method to have contiguous values in memory,
  // copy of block to simple.
  *workX = *data[q0];
  double *q = &(*workX)(0);

  SP::SiconosVector y = getInteractionPtr()->getYPtr(0);

  y->setValue(0, distance(q[0], q[1], q[2], r1, q[6], q[7], q[8], r2));
  y->setValue(1, 0.);

};

void SphereSphereR::computeJacH(double, unsigned int)
{

  SimpleMatrix *g = JacH[0].get();

  double q_0 = (*data[q0])(0);
  double q_1 = (*data[q0])(1);
  double q_2 = (*data[q0])(2);

  double q_6 = (*data[q0])(6);
  double q_7 = (*data[q0])(7);
  double q_8 = (*data[q0])(8);

  double r, A, B, C, D, nN, nU;
  /* u ^ v  = n */
  double u1, u2, u3, v1, v2, v3, n1, n2, n3, r1u1, r1u2, r1u3, r1v1, r1v2, r1v3, r2u1, r2u2, r2u3, r2v1, r2v2, r2v3;

  A = -(q_6 - q_0);
  B = -(q_7 - q_1);
  C = -(q_8 - q_2);

  nN = sqrt(A * A + B * B + C * C);

  n1 = A / nN;
  n2 = B / nN;
  n3 = C / nN;

  nU = sqrt((A + B) * (A + B) + (B - C) * (B - C) + (A + C) * (A + C));
  u1 = (B - C) / nU;
  u2 = -(A + C) / nU;
  u3 = (A + B) / nU;

  // v = n /\ u

  v1 = B * (A + B) / (nN * nU) - C * (-A - C) / (nN * nU);
  v2 = C * (B - C) / (nN * nU) - A * (A + B) / (nN * nU);
  v3 = A * (-A - C) / (nN * nU) - B * (B - C) / (nN * nU);

  r1u1 = r1 * u1;
  r1u2 = r1 * u2;
  r1u3 = r1 * u3;

  r1v1 = r1 * v1;
  r1v2 = r1 * v2;
  r1v3 = r1 * v3;

  r2u1 = r2 * u1;
  r2u2 = r2 * u2;
  r2u3 = r2 * u3;

  r2v1 = r2 * v1;
  r2v2 = r2 * v2;
  r2v3 = r2 * v3;

  (*g)(0, 0) = n1;
  (*g)(1, 0) = u1;
  (*g)(2, 0) = v1;
  (*g)(0, 1) = n2;
  (*g)(1, 1) = u2;
  (*g)(2, 1) = v2;
  (*g)(0, 2) = n3;
  (*g)(1, 2) = u3;
  (*g)(2, 2) = v3;
  (*g)(0, 3) = 0;          // -r n /\ n = 0
  (*g)(1, 3) = -r1v1;      // -r n /\ u = -r v
  (*g)(2, 3) = r1u1;       // -r n /\ v = r u
  (*g)(0, 4) = 0;
  (*g)(1, 4) = -r1v2;
  (*g)(2, 4) = r1u2;
  (*g)(0, 5) = 0;
  (*g)(1, 5) = -r1v3;
  (*g)(2, 5) = r1u3;

  (*g)(0, 6) = -n1;
  (*g)(1, 6) = -u1;
  (*g)(2, 6) = v1;
  (*g)(0, 7) = -n2;
  (*g)(1, 7) = -u2;
  (*g)(2, 7) = v2;
  (*g)(0, 8) = -n3;
  (*g)(1, 8) = -u3;
  (*g)(2, 8) = v3;
  (*g)(0, 9) = 0;          // -r n /\ n = 0
  (*g)(1, 9) = r2v1;       // -r n /\ u = -r v
  (*g)(2, 9) = -r2u1;      // -r n /\ v = r u
  (*g)(0, 10) = 0;
  (*g)(1, 10) = r2v2;
  (*g)(2, 10) = -r2u2;
  (*g)(0, 11) = 0;
  (*g)(1, 11) = r2v3;
  (*g)(2, 11) = -r2u3;

}

