/* Siconos-sample , Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

#include <math.h>

// parameters according to Table 1
// geometrical characteristics
double l1 = 0.1530;
double l2 = 0.3060;
double a = 0.05;
double b = 0.025;
double c = 0.001;
double d = 2. * (b + c);

// inertial properties
double m1 = 0.038;
double m2 = 0.038;
double m3 = 0.0760;
double J1 = 7.4e-5;
double J2 = 5.9e-4;
double J3 = 2.7e-6;

// force elements
double gravity = 9.81;

// plugins for smooth equations of motion
extern "C" void mass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  // columnwise definition of mass matrix
  mass[0] = J1 + (0.25 * m1 + m2 + m3) * l1 * l1;
  mass[1] = (0.5 * m2 + m3) * l1 * l2 * cos(q[1] - q[0]);
  mass[2] = 0.;

  mass[3] = (0.5 * m2 + m3) * l1 * l2 * cos(q[1] - q[0]);
  mass[4] = J2 + (0.25 * m2 + m3) * l2 * l2;
  mass[5] = 0.;

  mass[6] = 0.;
  mass[7] = 0.;
  mass[8] = J3;
}

extern "C" void NNL(unsigned int sizeOfq, const double *q, const double *velocity, double *NNL, unsigned int sizeZ, double* z)
{
  // nonlinear inertia terms (negative in h according to LagrangianDS)
  NNL[0] = (0.5 * m2 + m3) * l1 * l2 * sin(q[0] - q[1]) * velocity[1] * velocity[1];
  NNL[1] = -(0.5 * m2 + m3) * l1 * l2 * sin(q[0] - q[1]) * velocity[0] * velocity[0];
  NNL[2] = 0.;
}

extern "C" void jacobianNNLq(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  // Jacobian of nonlinear inertia terms with respect to q (columnwise)
  jacob[0] = (0.5 * m2 + m3) * l1 * l2 * cos(q[0] - q[1]) * velocity[1] * velocity[1];
  jacob[1] = -(0.5 * m2 + m3) * l1 * l2 * cos(q[0] - q[1]) * velocity[0] * velocity[0];
  jacob[2] = 0.;

  jacob[3] = -(0.5 * m2 + m3) * l1 * l2 * cos(q[0] - q[1]) * velocity[1] * velocity[1];
  jacob[4] = (0.5 * m2 + m3) * l1 * l2 * cos(q[0] - q[1]) * velocity[0] * velocity[0];
  jacob[5] = 0.;

  jacob[6] = 0.;
  jacob[7] = 0.;
  jacob[8] = 0.;
}

extern "C" void jacobianNNLqDot(unsigned int sizeOfq, const double *q, const  double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  // Jacobian of nonlinear inertia terms with respect to velocity (columnwise)
  jacob[0] =  0.;
  jacob[1] =  -2. * (0.5 * m2 + m3) * l1 * l2 * sin(q[0] - q[1]) * velocity[0];
  jacob[2] =  0.;

  jacob[3] = 2. * (0.5 * m2 + m3) * l1 * l2 * sin(q[0] - q[1]) * velocity[1];
  jacob[4] = 0.;
  jacob[5] = 0.;

  jacob[6] = 0.;
  jacob[7] = 0.;
  jacob[8] = 0.;
}

extern "C" void FInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *fInt, unsigned int sizeZ, double* z)
{
  // internal forces (negative in h according to LagrangianDS)
  fInt[0] = (0.5 * m1 + m2 + m3) * gravity * l1 * cos(q[0]);
  fInt[1] = (0.5 * m2 + m3) * gravity * l2 * cos(q[1]);
  fInt[2] = 0.;
}

extern "C" void jacobianFIntq(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  // Jacobian of internal forces with respect to q (columnwise)
  jacob[0] = -(0.5 * m1 + m2 + m3) * gravity * l1 * sin(q[0]);
  jacob[1] = 0.;
  jacob[2] = 0.;

  jacob[3] = 0.;
  jacob[4] = -(0.5 * m2 + m3) * gravity * l2 * sin(q[1]);
  jacob[5] = 0.;

  jacob[6] = 0.;
  jacob[7] = 0.;
  jacob[8] = 0.;
}

extern "C" void jacobianFIntqDot(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  // Jacobian of internal forces with respect to velocity (columnwise)
  jacob[0] = 0.;
  jacob[1] = 0.;
  jacob[2] = 0.;

  jacob[3] = 0.;
  jacob[4] = 0.;
  jacob[5] = 0.;

  jacob[6] = 0.;
  jacob[7] = 0.;
  jacob[8] = 0.;
}

extern "C" void g1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* g, unsigned int sizeZ, double* z)
{
  g[0] = 0.5 * d - (l1 * sin(q[0]) + l2 * sin(q[1]) - a * sin(q[2]) + b * cos(q[2]));
}

extern "C" void W1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* W, unsigned int sizeZ, double* z)
{
  W[0] = -l1 * cos(q[0]);
  W[1] = -l2 * cos(q[1]);
  W[2] = a * cos(q[2]) + b * sin(q[2]);
}

extern "C" void g2(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* g, unsigned int sizeZ, double* z)
{
  g[0] = 0.5 * d - (l1 * sin(q[0]) + l2 * sin(q[1]) + a * sin(q[2]) + b * cos(q[2]));
}

extern "C" void W2(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* W, unsigned int sizeZ, double* z)
{
  W[0] = -l1 * cos(q[0]);
  W[1] = -l2 * cos(q[1]);
  W[2] = -a * cos(q[2]) + b * sin(q[2]);
}

extern "C" void g3(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* g, unsigned int sizeZ, double* z)
{
  g[0] = 0.5 * d + l1 * sin(q[0]) + l2 * sin(q[1]) - a * sin(q[2]) - b * cos(q[2]);
}

extern "C" void W3(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* W, unsigned int sizeZ, double* z)
{
  W[0] = l1 * cos(q[0]);
  W[1] = l2 * cos(q[1]);
  W[2] = -a * cos(q[2]) + b * sin(q[2]);
}

extern "C" void g4(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* g, unsigned int sizeZ, double* z)
{
  g[0] = 0.5 * d + l1 * sin(q[0]) + l2 * sin(q[1]) + a * sin(q[2]) - b * cos(q[2]);
}

extern "C" void W4(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* W, unsigned int sizeZ, double* z)
{
  W[0] = l1 * cos(q[0]);
  W[1] = l2 * cos(q[1]);
  W[2] = a * cos(q[2]) + b * sin(q[2]);
}
