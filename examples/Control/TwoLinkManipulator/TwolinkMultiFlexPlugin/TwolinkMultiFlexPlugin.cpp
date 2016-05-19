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

#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>

#define PI 3.14159265

double l1 = 0.5;//length of the first link
double l2 = 0.5;//length of the second link
double m1 = 1; //mass of the first link
double m2 = 1; //mass of the second link
double I1 = 0.5;// the moment of inertia of the first link about the axis that passes through the center of mass (parallel to the Z axis)
double I2 = 0.5;// the moment of inertia of the second link about the axis that passes through the center of mass (parallel to the Z axis)
double g = 9.8;//gravitational acceleration
double gamma2 = 1;
double gamma1 = 10;
double Kf = 0.5;
double P = 10;
double ep = 0.5;
double Del = 0.5;
double delta = 0.2;
double alpha = 10;
double K1 = 2000;
double K2 = 2000;
double J1 = 0.1;
double J2 = 0.1;


SICONOS_EXPORT void mass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  //inertia matrix using (\theta)
  mass[0]  = m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  mass[1]  = I2 + m2 * l2 * l2 / 4 + m2 * l1 * l2 * cos(q[1]) / 2;
  mass[2]  = 0;
  mass[3]  = 0;
  mass[4]  = I2 + m2 * l2 * l2 / 4 + m2 * l1 * l2 * cos(q[1]) / 2;
  mass[5]  = I2 + m2 * l2 * l2 / 4;
  mass[6]  = 0;
  mass[7]  = 0;
  mass[8]  = 0;
  mass[9]  = 0;
  mass[10]  = J1;
  mass[11]  = 0;
  mass[12]  = 0;
  mass[13]  = 0;
  mass[14]  = 0;
  mass[15]  = J2;
}

SICONOS_EXPORT void FGyr(unsigned int sizeOfq, const double *q, const double *velocity, double *FGyr, unsigned int sizeZ, double* z)
{
  FGyr[0] = -m2 * l1 * l2 * sin(q[1]) * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2) + K1 * (q[0] - q[2]);
  FGyr[1] = m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[0] / 2 + g * m2 * l2 * cos(q[0] + q[1]) / 2 + K2 * (q[1] - q[3]);
  FGyr[2] = K1 * (q[2] - q[0]);
  FGyr[3] = K2 * (q[3] - q[1]);
}

SICONOS_EXPORT void jacobianFGyrq(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeOfZ, double* z)
{
  jacob[0] = -g * (l1 * sin(q[0]) * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) / 2) + K1;
  jacob[1] = -g * m2 * l2 * sin(q[0] + q[1]) / 2;
  jacob[2] = -K1;
  jacob[3] = 0;
  jacob[4] = -m2 * l1 * l2 * cos(q[1]) * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) - g * m2 * l2 * sin(q[0] + q[1]) / 2;
  jacob[5] = m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] / 2 - g * m2 * l2 * sin(q[0] + q[1]) / 2 + K2;
  jacob[6] = 0;
  jacob[7] = -K2;
  jacob[8] = -K1;
  jacob[9] = 0;
  jacob[10] = K1;
  jacob[11] = 0;
  jacob[12] = 0;
  jacob[13] = -K2;
  jacob[14] = 0;
  jacob[15] = K2;
}

SICONOS_EXPORT void jacobianVFGyr(unsigned int sizeOfq, const double *q, const  double *velocity, double *jacob, unsigned int sizeOfZ, double* z)
{
  jacob[0] = -m2 * l1 * l2 * sin(q[1]) * velocity[1];
  jacob[1] = m2 * l1 * l2 * sin(q[1]) * velocity[0];
  jacob[2] = 0;
  jacob[3] = 0;
  jacob[4] = -m2 * l1 * l2 * sin(q[1]) * (velocity[0] + velocity[1]);
  jacob[5] = 0;
  jacob[6] = 0;
  jacob[7] = 0;
  jacob[8] = 0;
  jacob[9] = 0;
  jacob[10] = 0;
  jacob[11] = 0;
  jacob[12] = 0;
  jacob[13] = 0;
  jacob[14] = 0;
  jacob[15] = 0;
}

SICONOS_EXPORT void U(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
{

  double m11 = m1 * l1 * l1 / 4 + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = m2 * l2 * l2 / 4 + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = m2 * l2 * l2 / 4 + I2;
  double detm = m11 * m22 - m12 * m12;
  double ddetm = (m12 - m22) * m2 * l1 * l2 * sin(q[1]) * velocity[1];

  double FGyr0 = -m2 * l1 * l2 * sin(q[1]) * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2) + K1 * (q[0] - q[2]);
  double FGyr1 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[0] / 2 + g * m2 * l2 * cos(q[0] + q[1]) / 2 + K2 * (q[1] - q[3]);

  // acceleration in q[0],q[1] coordinates
  double dv0 = -(m22 * FGyr0 - m12 * FGyr1) / detm;
  double dv1 = -(m11 * FGyr1 - m12 * FGyr0) / detm;

  double ddv0 = -((m22 * (-m2 * l1 * l2 * cos(q[1]) * velocity[1] * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) - m2 * l1 * l2 * sin(q[1]) * (dv0 * velocity[1] + dv1 * velocity[0] + dv1 * velocity[1]) - g * (l1 * sin(q[0]) * velocity[0] * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2) + K1 * (velocity[0] - velocity[2])) + m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr1 / 2 - m12 * (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] * velocity[1] / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * dv0 + K2 * (velocity[1] - velocity[3]))) * detm - (m22 * FGyr0 - m12 * FGyr1) * ddetm) / (detm * detm);
  double ddv1 = -((m11 * (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] * velocity[1] / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * dv0 + K2 * (velocity[1] - velocity[3])) + m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr0 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr1 + m12 * (m2 * l1 * l2 * cos(q[1]) * velocity[1] * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) + m2 * l1 * l2 * sin(q[1]) * (dv0 * velocity[1] + dv1 * velocity[0] + dv1 * velocity[1]) + g * m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2 - K1 * (velocity[0] - velocity[2]))) * detm - (m11 * FGyr1 - m12 * FGyr0) * ddetm) / (detm * detm);

  // generalized coordinates q=(x,y)^T
  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

  // time derivatives of x and y
  double x1 = -l1 * sin(q[0]) * velocity[0] - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]);
  double y1 = l1 * cos(q[0]) * velocity[0] + l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]);

  // acceleration in x,y coordinates
  double x2 =  -l1 * cos(q[0]) * velocity[0] * velocity[0] - l1 * sin(q[0]) * dv0 - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - l2 * sin(q[0] + q[1]) * (dv0 + dv1);
  double y2 =  l1 * cos(q[0]) * dv0 - l1 * sin(q[0]) * velocity[0] * velocity[0] + l2 * cos(q[0] + q[1]) * (dv0 + dv1) - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]);

  double x3 =  l1 * sin(q[0]) * velocity[0] * velocity[0] * velocity[0] - 2 * l1 * cos(q[0]) * velocity[0] * dv0 - l1 * cos(q[0]) * velocity[0] * dv0 - l1 * sin(q[0]) * ddv0 + l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - 2 * l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - l2 * sin(q[0] + q[1]) * (ddv0 + ddv1);
  double y3 = -l1 * sin(q[0]) * velocity[0] * dv0 + l1 * cos(q[0]) * ddv0 - l1 * cos(q[0]) * velocity[0] * velocity[0] * velocity[0] - 2 * l1 * sin(q[0]) * velocity[0] * dv0 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) + l2 * cos(q[0] + q[1]) * (ddv0 + ddv1) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - 2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1);

  // the gradient of the transformation \theta->q
  double grad11 = -y;
  double grad12 = -l2 * sin(q[0] + q[1]);
  double grad21 = x;
  double grad22 = l2 * cos(q[0] + q[1]);
  double det = grad11 * grad22 - grad12 * grad21;

  // d=(grad)^{-1}
  double d11 = grad22 / det;
  double d12 = -grad12 / det;
  double d21 = -grad21 / det;
  double d22 = grad11 / det;

  // time derivative of d

  double dd11 =  -(sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd12 = (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd21 =  -dd11 + (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));
  double dd22 =  -dd12 - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));

  // time derivative of dd

  double ddd11 = (-(sin(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] + cos(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd12 = ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + sin(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] - sin(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd21 =  -ddd11 + ((cos(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + sin(q[0]) * sin(q[1]) * dv0 + cos(q[0]) * cos(q[1]) * dv1) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd22 =  -ddd12 - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]));

  // time derivative of ddd

  double dddd11 = ((-(sin(q[0] + q[1]) * (ddv0 + ddv1) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * ddv1 + 3 * cos(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - 3 * cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * dv1 + cos(q[0] + q[1]) * cos(q[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * velocity[1] - velocity[1] * velocity[1] * velocity[1]) - sin(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - velocity[1] * velocity[1]) + sin(q[0] + q[1]) * cos(q[1]) * ((dv0 + dv1) * velocity[1] - (velocity[0] + velocity[1]) * dv1)) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - (-(sin(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] + cos(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd12 = (((cos(q[0] + q[1]) * (ddv0 + ddv1) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * ddv1 - 3 * sin(q[0] + q[1]) * sin(q[1]) * ((velocity[0] + velocity[1]) * (dv0 + dv1) - velocity[1] * dv1) + cos(q[0] + q[1]) * cos(q[1]) * ((dv0 + dv1) * velocity[1] - (velocity[0] + velocity[1]) * dv1) - (cos(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) + sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + sin(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] - sin(q[0] + q[1]) * cos(q[1]) * dv1) * (sin(q[1]) * sin(q[1])) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * sin(q[1]) * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd21 = -dddd11 + (((sin(q[0]) * sin(q[1]) * ddv0 + cos(q[0]) * cos(q[1]) * ddv1 + 3 * (velocity[0] * dv0 - velocity[1] * dv1) + sin(q[0]) * cos(q[1]) * (velocity[1] * dv0 - velocity[0] * dv1) + (cos(q[0]) * cos(q[1]) * velocity[1] - sin(q[0]) * sin(q[1]) * velocity[0]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + sin(q[0]) * sin(q[1]) * dv0 + cos(q[0]) * cos(q[1]) * dv1) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd22 = -dddd12 - (((-(sin(q[0]) * cos(q[1]) * velocity[1] + cos(q[0]) * sin(q[1]) * velocity[0]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + cos(q[0]) * sin(q[1]) * ddv0 - sin(q[0]) * cos(q[1]) * ddv1 + 3 * (velocity[1] * dv1 - velocity[0] * dv0) * sin(q[0]) * sin(q[1]) - cos(q[0]) * cos(q[1]) * (velocity[1] * dv0 - velocity[0] * dv1)) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1])) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));

  // M*d
  double mc11 = m11 * d11 + m12 * d21;
  double mc12 = m11 * d12 + m12 * d22;
  double mc21 = m12 * d11 + m22 * d21;
  double mc22 = m12 * d12 + m22 * d22;

  // M*dd
  double a11 = m11 * dd11 + m12 * dd21;
  double a12 = m11 * dd12 + m12 * dd22;
  double a21 = m12 * dd11 + m22 * dd21;
  double a22 = m12 * dd12 + m22 * dd22;

  double da11 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd11 + dd21 / 2) + m11 * ddd11 + m12 * ddd21;
  double da12 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2) + m11 * ddd12 + m12 * ddd22;
  double da21 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 / 2 + m12 * ddd11 + m22 * ddd21;
  double da22 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 / 2 + m12 * ddd12 + m22 * ddd22;

  double dda11 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (dd11 + dd21 / 2) + 2 * sin(q[1]) * velocity[1] * (ddd11 + ddd21 / 2)) + m11 * dddd11 + m12 * dddd21;
  double dda12 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (dd12 + dd22 / 2) + 2 * sin(q[1]) * velocity[1] * (ddd12 + ddd22 / 2)) + m11 * dddd12 + m12 * dddd22;
  double dda21 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd11 / 2 + sin(q[1]) * velocity[1] * ddd11) + m12 * dddd11 + m22 * dddd21;
  double dda22 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd12 / 2 + sin(q[1]) * velocity[1] * ddd12) + m12 * dddd12 + m22 * dddd22;

  double dmc11 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (d11 + d21 / 2) + m11 * dd11 + m12 * dd21;
  double dmc12 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (d12 + d22 / 2) + m11 * dd12 + m12 * dd22;
  double dmc21 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * d11 / 2 + m12 * dd11 + m22 * dd21;
  double dmc22 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * d12 / 2 + m12 * dd12 + m22 * dd22;

  double ddmc11 = -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (d11 + d21 / 2) - 2 * m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd11 + dd21 / 2) + m11 * ddd11 + m12 * ddd21;
  double ddmc12 = -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (d12 + d22 / 2) - 2 * m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2) + m11 * ddd12 + m12 * ddd22;
  double ddmc21 =  -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * d11 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 + m12 * ddd11 + m22 * ddd21;
  double ddmc22 =  -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * d12 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 + m12 * ddd12 + m22 * ddd22;

  double qd1 = 0.7 + 0.5 * cos(2 * PI * time / P + PI / 2);
  double qd2 = 0.5 * sin(2 * PI * time / P + PI / 2);
  double qd11 = -(2 * PI / P) * 0.5 * sin(2 * PI * time / P + PI / 2);
  double qd12 = (2 * PI / P) * 0.5 * cos(2 * PI * time / P + PI / 2);
  double qd21 = -(2 * PI / P) * (2 * PI / P) * 0.5 * cos(2 * PI * time / P + PI / 2);
  double qd22 = -(2 * PI / P) * (2 * PI / P) * 0.5 * sin(2 * PI * time / P + PI / 2);
  double qd31 = (2 * PI / P) * (2 * PI / P) * (2 * PI / P) * 0.5 * sin(2 * PI * time / P + PI / 2);
  double qd32 = -(2 * PI / P) * (2 * PI / P) * (2 * PI / P) * 0.5 * cos(2 * PI * time / P + PI / 2);
  double qd41 = (2 * PI / P) * (2 * PI / P) * (2 * PI / P) * (2 * PI / P) * 0.5 * cos(2 * PI * time / P + PI / 2);
  double qd42 = (2 * PI / P) * (2 * PI / P) * (2 * PI / P) * (2 * PI / P) * 0.5 * sin(2 * PI * time / P + PI / 2);

  double qr11 = qd11 - gamma2 * (x - qd1);
  double qr12 = qd12 - gamma2 * (y - qd2);

  double qr21 = qd21 - gamma2 * (x1 - qd11);
  double qr22 = qd22 - gamma2 * (y1 - qd12);

  double qr31 = qd31 - gamma2 * (x2 - qd21);
  double qr32 = qd32 - gamma2 * (y2 - qd22);

  double qr41 = qd41 - gamma2 * (x3 - qd31);
  double qr42 = qd42 - gamma2 * (y3 - qd32);

  double s1 = x1 - qr11;
  double s2 = y1 - qr12;

  double ki1 = sqrt((qd1 * qd1 + qd2 * qd2 + l1 * l1 + l2 * l2) * (qd1 * qd1 + qd2 * qd2 + l1 * l1 + l2 * l2) - 2 * ((qd1 * qd1 + qd2 * qd2) * (qd1 * qd1 + qd2 * qd2) + l1 * l1 * l1 * l1 + l2 * l2 * l2 * l2));
  double ki2 = qd1 * qd1 + qd2 * qd2 + l1 * l1 - l2 * l2;
  double ki3 = qd1 * qd1 + qd2 * qd2 - l1 * l1 - l2 * l2;

  double xd = atan2(qd2, qd1) + atan2(ki1, ki2);
  double yd = -atan2(ki1, ki3);
  double xd1 = (qd1 * qd12 - qd2 * qd11) / (qd1 * qd1 + qd2 * qd2) + ((qd1 * qd11 + qd2 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2)) / (ki1 * (qd1 * qd1 + qd2 * qd2));
  double yd1 = 2 * (qd1 * qd11 + qd2 * qd12) / ki1;
  double xd2 = (((qd1 * qd21 + qd11 * qd11 + qd2 * qd22 + qd12 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2) - 2 * (qd1 * qd11 + qd2 * qd12) * (qd1 * qd11 + qd2 * qd12)) * ki1 * (qd1 * qd1 + qd2 * qd2) - (qd1 * qd11 + qd2 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2) * (ki1 * 2 * (qd1 * qd11 + qd2 * qd12) - (qd1 * qd1 + qd2 * qd2) * 2 * (qd1 * qd11 + qd2 * qd12) * ki3 / ki1)) / (ki1 * (qd1 * qd1 + qd2 * qd2) * ki1 * (qd1 * qd1 + qd2 * qd2)) + ((qd1 * qd22 - qd2 * qd21) * (qd1 * qd1 + qd2 * qd2) - (qd1 * qd12 - qd2 * qd11) * 2 * (qd1 * qd11 + qd2 * qd12)) / ((qd1 * qd1 + qd2 * qd2) * (qd1 * qd1 + qd2 * qd2));
  double yd2 = 2 * ((qd1 * qd21 + qd11 * qd11 + qd2 * qd22 + qd12 * qd12) * ki1 * ki1 + (qd1 * qd11 + qd2 * qd12) * 2 * (qd1 * qd11 + qd2 * qd12) * ki3) / (ki1 * ki1 * ki1);

  double a1 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2);
  double a2 = a11 * qr11 + a12 * qr12;
  double a3 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * qr11 + d12 * qr12) / 2;
  double a4 = a21 * qr11 + a22 * qr12;

  double da1 = -m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * dv1 * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2);
  double da2 = a11 * qr21 + a12 * qr22 + da11 * qr11 + da12 * qr12;
  double da3 = m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[1] * (d11 * qr11 + d12 * qr12) / 2 + m2 * l1 * l2 * sin(q[1]) * dv0 * (d11 * qr11 + d12 * qr12) / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * (dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) / 2;
  double da4 = a21 * qr21 + a22 * qr22 + da21 * qr11 + da22 * qr12;

  double dda1 = (m2 * l1 * l2 * sin(q[1]) * velocity[1] * velocity[1] * velocity[1] - 2 * m2 * l1 * l2 * cos(q[1]) * velocity[1] * dv1) * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] * ((dd11 * qr11 + d11 * qr21 + dd12 * qr12 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - (m2 * l1 * l2 * cos(q[1]) * velocity[1] * dv1 + m2 * l1 * l2 * sin(q[1]) * ddv1) * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * dv1 * ((dd11 * qr11 + d11 * qr21 + dd12 * qr12 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - (m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] + m2 * l1 * l2 * sin(q[1]) * dv1) * ((dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((ddd11 * qr11 + dd11 * qr21 + ddd12 * qr12 + dd12 * qr22 + dd11 * qr21 + d11 * qr31 + dd12 * qr22 + d12 * qr32) + (ddd21 * qr11 + dd21 * qr21 + ddd22 * qr12 + dd22 * qr22 + dd21 * qr21 + d21 * qr31 + dd22 * qr22 + d22 * qr32) / 2);
  double dda2 = da11 * qr21 + a11 * qr31 + da12 * qr22 + a12 * qr32 + da11 * qr21 + dda11 * qr11 + da12 * qr22 + dda12 * qr12;
  double dda3 = (m2 * l1 * l2 * sin(q[1]) * ddv0 + 2 * m2 * l1 * l2 * cos(q[1]) * dv0 * velocity[1] + m2 * l1 * l2 * cos(q[1]) * velocity[0] * dv1 - m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[1] * velocity[1]) * (d11 * qr11 + d12 * qr12) / 2 + (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[1] + m2 * l1 * l2 * sin(q[1]) * dv0) * (dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + m2 * l1 * l2 * sin(q[1]) * velocity[0] * (ddd11 * qr11 + ddd12 * qr12 + dd11 * qr21 + dd12 * qr22 + dd11 * qr21 + dd12 * qr22 + d11 * qr31 + d12 * qr32) / 2;
  double dda4 = da21 * qr21 + a21 * qr31 + da22 * qr22 + a22 * qr32 + dda21 * qr11 + da21 * qr21 + dda22 * qr12 + da22 * qr22;

  double T01 = mc11 * qr21 + mc12 * qr22;
  double T02 = a1 + a2 + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2);
  double T03 = grad11 * gamma1 * s1 + grad21 * gamma1 * s2;

  double dT01 = dmc11 * qr21 + dmc12 * qr22 + mc11 * qr31 + mc12 * qr32;
  double dT02 = da1 + da2 - g * (l1 * sin(q[0]) * velocity[0] * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2);
  double dT03 = -y1 * gamma1 * s1 + x1 * gamma1 * s2 + grad11 * gamma1 * (x2 - qr21) + grad21 * gamma1 * (y2 - qr22);

  double ddT01 = ddmc11 * qr21 + ddmc12 * qr22 + 2 * dmc11 * qr31 + 2 * dmc12 * qr32 + mc11 * qr41 + mc12 * qr42;
  double ddT02 = dda1 + dda2 - g * (l1 * sin(q[0]) * dv0 * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (dv0 + dv1) / 2 + l1 * cos(q[0]) * velocity[0] * velocity[0] * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) / 2);
  double ddT03 = -y2 * gamma1 * s1 + x2 * gamma1 * s2 - 2 * y1 * gamma1 * (x2 - qr21) + 2 * x1 * gamma1 * (y2 - qr22) + grad11 * gamma1 * (x3 - qr31) + grad21 * gamma1 * (y3 - qr32);

  double T11 = mc21 * qr21 + mc22 * qr22;
  double T12 = a3 + a4 + g * m2 * l2 * cos(q[0] + q[1]) / 2;
  double T13 = grad12 * gamma1 * s1 + grad22 * gamma1 * s2;

  double dT11 = dmc21 * qr21 + dmc22 * qr22 + mc21 * qr31 + mc22 * qr32;
  double dT12 = da3 + da4 - g * m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2;
  double dT13 = -l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * s1 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * s2 + grad12 * gamma1 * (x2 - qr21) + grad22 * gamma1 * (y2 - qr22);

  double ddT11 = ddmc21 * qr21 + ddmc22 * qr22 + 2 * dmc21 * qr31 + 2 * dmc22 * qr32 + mc21 * qr41 + mc22 * qr42;
  double ddT12 = dda3 + dda4 - g * m2 * l2 * (sin(q[0] + q[1]) * (dv0 + dv1) / 2 + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) / 2);
  double ddT13 = l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * gamma1 * s1 - l2 * cos(q[0] + q[1]) * (dv0 + dv1) * gamma1 * s1 - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (x2 - qr21) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * gamma1 * s2 - l2 * sin(q[0] + q[1]) * (dv0 + dv1) * gamma1 * s2 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (y2 - qr22) + grad12 * gamma1 * (x3 - qr31) + grad22 * gamma1 * (y3 - qr32) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (x2 - qr21) - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (y2 - qr22);

  double thetad1 = xd + (T01 + T02 - T03) / K1;
  double thetad2 = yd + (T11 + T12 - T13) / K2;
  double thetad11 = xd1 + (dT01 + dT02 - dT03) / K1;
  double thetad12 = yd1 + (dT11 + dT12 - dT13) / K2;
  double thetad21 = xd2 + (ddT01 + ddT02 - ddT03) / K1;
  double thetad22 = yd2 + (ddT11 + ddT12 - ddT13) / K2;

  double s11 = (velocity[0] - xd1) + gamma2 * (q[0] - xd);
  double s12 = (velocity[1] - yd1) + gamma2 * (q[1] - yd);
  double s21 = (velocity[2] - thetad11) + gamma2 * (q[2] - thetad1);
  double s22 = (velocity[3] - thetad12) + gamma2 * (q[3] - thetad2);

  double thetar21 = thetad21 - gamma2 * (velocity[2] - thetad11);
  double thetar22 = thetad22 - gamma2 * (velocity[3] - thetad12);


  // control law
  U[0] = 0;//-(T01+T02-T03);//0;
  U[1] = 0;//-(T11+T12-T13);//0;
  U[2] = -(J1 * thetar21 + K1 * (thetad1 - xd) - gamma1 * s21);
  U[3] = -(J2 * thetar22 + K2 * (thetad2 - yd) - gamma1 * s22);

  // double V1 =(s1*s1*(d11*mc11+d21*mc21)+2*s1*s2*(d12*mc11+d22*mc21)+s2*s2*(d12*mc12+d22*mc22))/2+(J1*s21*s21+J2*s22*s22)/2;//
  double V1 = (s11 * s11 * m11 + 2 * s11 * s12 * m12 + s12 * s12 * m22) / 2 + (J1 * s21 * s21 + J2 * s22 * s22) / 2;
  z[6] = V1 + gamma1 * gamma2 * ((q[0] - xd) * (q[0] - xd) + (q[1] - yd) * (q[1] - yd) + (q[2] - thetad1) * (q[2] - thetad1) + (q[3] - thetad2) * (q[3] - thetad2)) + (q[0] - xd - q[2] + thetad1) * K1 * (q[0] - xd - q[2] + thetad1) / 2 + (q[1] - yd - q[3] + thetad2) * K2 * (q[1] - yd - q[3] + thetad2) / 2;

  z[9] = V1;
  z[11] = P;
  z[13] = -1;
  z[14] = qd1;
  z[15] = qd2;
  z[16] = U[2];
  z[17] = U[3];
  z[19] = s1;
  z[20] = s2;
  z[21] = qr11;
  z[22] = qr12;
  z[23] = qr21;
  z[24] = qr22;
}


SICONOS_EXPORT void U1(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
{
  double m11 = m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = (m2 * l2 * l2 / 4) + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = (m2 * l2 * l2 / 4) + I2;
  double detm = m11 * m22 - m12 * m12;
  double ddetm = (m12 - m22) * m2 * l1 * l2 * sin(q[1]) * velocity[1];

  double FGyr0 = -m2 * l1 * l2 * sin(q[1]) * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2) + K1 * (q[0] - q[2]);
  double FGyr1 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[0] / 2 + g * m2 * l2 * cos(q[0] + q[1]) / 2 + K2 * (q[1] - q[3]);

  // acceleration in q[0],q[1] coordinates
  double dv0 = -(m22 * FGyr0 - m12 * FGyr1) / detm;
  double dv1 = -(m11 * FGyr1 - m12 * FGyr0) / detm;

  double ddv0 = -((m22 * (-m2 * l1 * l2 * cos(q[1]) * velocity[1] * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) - m2 * l1 * l2 * sin(q[1]) * (dv0 * velocity[1] + dv1 * velocity[0] + dv1 * velocity[1]) - g * (l1 * sin(q[0]) * velocity[0] * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2) + K1 * (velocity[0] - velocity[2])) + m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr1 / 2 - m12 * (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] * velocity[1] / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * dv0 + K2 * (velocity[1] - velocity[3]))) * detm - (m22 * FGyr0 - m12 * FGyr1) * ddetm) / (detm * detm);
  double ddv1 = -((m11 * (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] * velocity[1] / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * dv0 + K2 * (velocity[1] - velocity[3])) + m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr0 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr1 + m12 * (m2 * l1 * l2 * cos(q[1]) * velocity[1] * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) + m2 * l1 * l2 * sin(q[1]) * (dv0 * velocity[1] + dv1 * velocity[0] + dv1 * velocity[1]) + g * m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2 - K1 * (velocity[0] - velocity[2]))) * detm - (m11 * FGyr1 - m12 * FGyr0) * ddetm) / (detm * detm);

  // generalized coordinates q=(x,y)^T
  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

  // time derivatives of x and y
  double x1 = -l1 * sin(q[0]) * velocity[0] - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]);
  double y1 = l1 * cos(q[0]) * velocity[0] + l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]);

  // acceleration in x,y coordinates
  double x2 = -l1 * cos(q[0]) * velocity[0] * velocity[0] - l1 * sin(q[0]) * dv0 - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - l2 * sin(q[0] + q[1]) * (dv0 + dv1);
  double y2 =  l1 * cos(q[0]) * dv0 - l1 * sin(q[0]) * velocity[0] * velocity[0] + l2 * cos(q[0] + q[1]) * (dv0 + dv1) - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]);

  double x3 = l1 * sin(q[0]) * velocity[0] * velocity[0] * velocity[0] - 2 * l1 * cos(q[0]) * velocity[0] * dv0 - l1 * cos(q[0]) * velocity[0] * dv0 - l1 * sin(q[0]) * ddv0 + l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - 2 * l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - l2 * sin(q[0] + q[1]) * (ddv0 + ddv1);
  double y3 = -l1 * sin(q[0]) * velocity[0] * dv0 + l1 * cos(q[0]) * ddv0 - l1 * cos(q[0]) * velocity[0] * velocity[0] * velocity[0] - 2 * l1 * sin(q[0]) * velocity[0] * dv0 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) + l2 * cos(q[0] + q[1]) * (ddv0 + ddv1) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - 2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1);

  // the gradient of the transformation \theta->q
  double grad11 = -y;
  double grad12 = -l2 * sin(q[0] + q[1]);
  double grad21 = x;
  double grad22 = l2 * cos(q[0] + q[1]);
  double det = grad11 * grad22 - grad12 * grad21;

  // d=(grad)^{-1}
  double d11 = grad22 / det;
  double d12 = -grad12 / det;
  double d21 = -grad21 / det;
  double d22 = grad11 / det;

  // time derivative of d

  double dd11 =  -(sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd12 = (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd21 =  -dd11 + (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));
  double dd22 =  -dd12 - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));

  // time derivative of dd

  double ddd11 = (-(sin(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] + cos(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd12 = ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + sin(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] - sin(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd21 =  -ddd11 + ((cos(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + sin(q[0]) * sin(q[1]) * dv0 + cos(q[0]) * cos(q[1]) * dv1) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd22 =  -ddd12 - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]));

  // time derivative of ddd

  double dddd11 = ((-(sin(q[0] + q[1]) * (ddv0 + ddv1) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * ddv1 + 3 * cos(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - 3 * cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * dv1 + cos(q[0] + q[1]) * cos(q[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * velocity[1] - velocity[1] * velocity[1] * velocity[1]) - sin(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - velocity[1] * velocity[1]) + sin(q[0] + q[1]) * cos(q[1]) * ((dv0 + dv1) * velocity[1] - (velocity[0] + velocity[1]) * dv1)) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - (-(sin(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] + cos(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd12 = (((cos(q[0] + q[1]) * (ddv0 + ddv1) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * ddv1 - 3 * sin(q[0] + q[1]) * sin(q[1]) * ((velocity[0] + velocity[1]) * (dv0 + dv1) - velocity[1] * dv1) + cos(q[0] + q[1]) * cos(q[1]) * ((dv0 + dv1) * velocity[1] - (velocity[0] + velocity[1]) * dv1) - (cos(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) + sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + sin(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] - sin(q[0] + q[1]) * cos(q[1]) * dv1) * (sin(q[1]) * sin(q[1])) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * sin(q[1]) * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd21 = -dddd11 + (((sin(q[0]) * sin(q[1]) * ddv0 + cos(q[0]) * cos(q[1]) * ddv1 + 3 * (velocity[0] * dv0 - velocity[1] * dv1) + sin(q[0]) * cos(q[1]) * (velocity[1] * dv0 - velocity[0] * dv1) + (cos(q[0]) * cos(q[1]) * velocity[1] - sin(q[0]) * sin(q[1]) * velocity[0]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + sin(q[0]) * sin(q[1]) * dv0 + cos(q[0]) * cos(q[1]) * dv1) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd22 = -dddd12 - (((-(sin(q[0]) * cos(q[1]) * velocity[1] + cos(q[0]) * sin(q[1]) * velocity[0]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + cos(q[0]) * sin(q[1]) * ddv0 - sin(q[0]) * cos(q[1]) * ddv1 + 3 * (velocity[1] * dv1 - velocity[0] * dv0) * sin(q[0]) * sin(q[1]) - cos(q[0]) * cos(q[1]) * (velocity[1] * dv0 - velocity[0] * dv1)) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1])) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));

  // M*d
  double mc11 = m11 * d11 + m12 * d21;
  double mc12 = m11 * d12 + m12 * d22;
  double mc21 = m12 * d11 + m22 * d21;
  double mc22 = m12 * d12 + m22 * d22;

  // M*dd
  double a11 = m11 * dd11 + m12 * dd21;
  double a12 = m11 * dd12 + m12 * dd22;
  double a21 = m12 * dd11 + m22 * dd21;
  double a22 = m12 * dd12 + m22 * dd22;

  double da11 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd21 / 2 + m11 * ddd11 + m12 * ddd21;
  double da12 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2) + m11 * ddd12 + m12 * ddd22;
  double da21 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 / 2 + m12 * ddd11 + m22 * ddd21;
  double da22 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 / 2 + m12 * ddd12 + m22 * ddd22;

  double dda11 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd11 + sin(q[1]) * velocity[1] * ddd11 + (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd21 / 2 + sin(q[1]) * velocity[1] * ddd21 / 2 + sin(q[1]) * velocity[1] * (ddd11 + ddd21 / 2)) + m11 * dddd11 + m12 * dddd21;
  double dda12 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (dd12 + dd22 / 2) + sin(q[1]) * velocity[1] * (ddd12 + ddd22 / 2) + sin(q[1]) * velocity[1] * (ddd12 + ddd22 / 2)) + m11 * dddd12 + m12 * dddd22;
  double dda21 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd11 / 2 + sin(q[1]) * velocity[1] * ddd11) + m12 * dddd11 + m22 * dddd21;
  double dda22 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd12 / 2 + sin(q[1]) * velocity[1] * ddd12) + m12 * dddd12 + m22 * dddd22;

  double dmc11 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (d11 + d21 / 2) + m11 * dd11 + m12 * dd21;
  double dmc12 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (d12 + d22 / 2) + m11 * dd12 + m12 * dd22;
  double dmc21 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * d11 / 2 + m12 * dd11 + m22 * dd21;
  double dmc22 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * d12 / 2 + m12 * dd12 + m22 * dd22;

  double ddmc11 = -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (d11 + d21 / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd11 + dd21 / 2) + m11 * ddd11 + m12 * ddd21 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd11 + dd21 / 2);
  double ddmc12 = -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (d12 + d22 / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2) + m11 * ddd12 + m12 * ddd22 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2);
  double ddmc21 =  -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * d11 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 / 2 + m12 * ddd11 + m22 * ddd21 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 / 2;
  double ddmc22 =  -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * d12 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 / 2 + m12 * ddd12 + m22 * ddd22 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 / 2;

  //double t0 = z[8];
  double t1 = z[8] + delta + Del;
  double t2 = z[8] + delta;


  double qd1 = 0;
  double qd11 = 0;
  double qd21 = 0;
  double qd31 = 0;
  double qd41 = 0;

  double a = -sqrt(z[7]) * alpha;
  double b0 = 0.5 * sin(PI / 2 + 2 * PI * z[8] / P);
  double b1 = 0.5 * (2 * PI / P) * cos(PI / 2 + 2 * PI * z[8] / P);
  double b2 = -0.5 * (2 * PI / P) * (2 * PI / P) * sin(PI / 2 + 2 * PI * z[8] / P);
  double b3 = -0.5 * (2 * PI / P) * (2 * PI / P) * (2 * PI / P) * cos(PI / 2 + 2 * PI * z[8] / P);

  double p0 = -((delta + Del) * (delta + Del) * (delta + Del) * b3 + 3 * (delta + Del) * (delta + Del) * b2 + 6 * (delta + Del) * b1 - 6 * (a - b0)) / (6 * (delta + Del) * (delta + Del) * (delta + Del) * (delta + Del));
  double p1 = ((delta + Del) * (delta + Del) * (delta + Del) * b3 + 6 * (delta + Del) * (delta + Del) * b2 + 18 * (delta + Del) * b1 - 24 * (a - b0)) / (6 * (delta + Del) * (delta + Del) * (delta + Del) * (delta + Del) * (delta + Del));
  double p2 = -((delta + Del) * (delta + Del) * (delta + Del) * b3 + 9 * (delta + Del) * (delta + Del) * b2 + 36 * (delta + Del) * b1 - 60 * (a - b0)) / (6 * (delta + Del) * (delta + Del) * (delta + Del) * (delta + Del) * (delta + Del) * (delta + Del));
  double p3 = ((delta + Del) * (delta + Del) * (delta + Del) * b3 + 12 * (delta + Del) * (delta + Del) * b2 + 60 * (delta + Del) * b1 - 120 * (a - b0)) / (6 * (delta + Del) * (delta + Del) * (delta + Del) * (delta + Del) * (delta + Del) * (delta + Del) * (delta + Del));

  double c0 = z[5];
  double c1 = -0.5 * (2 * PI / P) * sin(PI / 2 + 2 * PI * z[8] / P);
  double c2 = -0.5 * (2 * PI / P) * (2 * PI / P) * cos(PI / 2 + 2 * PI * z[8] / P);
  double c3 = 0.5 * (2 * PI / P) * (2 * PI / P) * (2 * PI / P) * sin(PI / 2 + 2 * PI * z[8] / P);
  double w0 = -(delta * delta * delta * c3 + 3 * delta * delta * c2 + 6 * delta * c1) / (6 * delta * delta * delta * delta);
  double w1 = (delta * delta * delta * c3 + 6 * delta * delta * c2 + 18 * delta * c1) / (6 * delta * delta * delta * delta * delta);
  double w2 = -(delta * delta * delta * c3 + 9 * delta * delta * c2 + 36 * delta * c1) / (6 * delta * delta * delta * delta * delta * delta);
  double w3 = (delta * delta * delta * c3 + 12 * delta * delta * c2 + 60 * delta * c1) / (6 * delta * delta * delta * delta * delta * delta * delta);

  double qd2 = 0;
  double qd12 = 0;
  double qd22 = 0;
  double qd32 = 0;
  double qd42 = 0;
  if (time <= t1)
  {
    if (time <= t2)
    {
      qd1 = c0 + c1 * (time - z[8]) + c2 * (time - z[8]) * (time - z[8]) / 2 + c3 * (time - z[8]) * (time - z[8]) * (time - z[8]) / 6 + w0 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) + w1 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - t2) + w2 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - t2) * (time - t2) + w3 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - t2) * (time - t2) * (time - t2);
      qd11 = c1 + c2 * (time - z[8]) + c3 * (time - z[8]) * (time - z[8]) / 2 + 4 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (w0 + w1 * (time - t2) + w2 * (time - t2) * (time - t2) + w3 * (time - t2) * (time - t2) * (time - t2)) + (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (w1 + 2 * w2 * (time - t2) + 3 * w3 * (time - t2) * (time - t2));
      qd21 = c2 + c3 * (time - z[8]) + 12 * (time - z[8]) * (time - z[8]) * (w0 + w1 * (time - t2) + w2 * (time - t2) * (time - t2) + w3 * (time - t2) * (time - t2) * (time - t2)) + 8 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (w1 + 2 * w2 * (time - t2) + 3 * w3 * (time - t2) * (time - t2)) + (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (2 * w2 + 6 * w3 * (time - t2));
      qd31 = c3 + 24 * (time - z[8]) * (w0 + w1 * (time - t2) + w2 * (time - t2) * (time - t2) + w3 * (time - t2) * (time - t2) * (time - t2)) + 36 * (time - z[8]) * (time - z[8]) * (w1 + 2 * w2 * (time - t2) + 3 * w3 * (time - t2) * (time - t2)) + 12 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (2 * w2 + 6 * w3 * (time - t2)) + (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * 6 * w3;
    }
    else
    {
      qd1 = z[5];
      qd11 = 0;
      qd21 = 0;
      qd31 = 0;
    }
    qd2 = b0 + b1 * (time - z[8]) + b2 * (time - z[8]) * (time - z[8]) / 2 + b3 * (time - z[8]) * (time - z[8]) * (time - z[8]) / 6 + p0 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) + p1 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - t1) + p2 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - t1) * (time - t1) + p3 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - t1) * (time - t1) * (time - t1);
    qd12 = b1 + b2 * (time - z[8]) + b3 * (time - z[8]) * (time - z[8]) / 2 + 4 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (p0 + p1 * (time - t1) + p2 * (time - t1) * (time - t1) + p3 * (time - t1) * (time - t1) * (time - t1)) + (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (p1 + 2 * p2 * (time - t1) + 3 * p3 * (time - t1) * (time - t1));
    qd22 = b2 + b3 * (time - z[8]) + 12 * (time - z[8]) * (time - z[8]) * (p0 + p1 * (time - t1) + p2 * (time - t1) * (time - t1) + p3 * (time - t1) * (time - t1) * (time - t1)) + 8 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (p1 + 2 * p2 * (time - t1) + 3 * p3 * (time - t1) * (time - t1)) + (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (2 * p2 + 6 * p3 * (time - t1));
    qd32 = b3 + 24 * (time - z[8]) * (p0 + p1 * (time - t1) + p2 * (time - t1) * (time - t1) + p3 * (time - t1) * (time - t1) * (time - t1)) + 36 * (time - z[8]) * (time - z[8]) * (p1 + 2 * p2 * (time - t1) + 3 * p3 * (time - t1) * (time - t1)) + 12 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (2 * p2 + 6 * p3 * (time - t1)) + (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * 6 * p3;
  }
  else
  {
    qd1 = z[5];
    qd11 = 0;
    qd21 = 0;
    qd31 = 0;
    qd2 = -sqrt(z[7]) * alpha;
    qd12 = 0;
    qd22 = 0;
    qd32 = 0;
  }


  double qr11 = qd11 - gamma2 * (x - qd1);
  double qr12 = qd12 - gamma2 * (y - qd2);

  double qr21 = qd21 - gamma2 * (x1 - qd11);
  double qr22 = qd22 - gamma2 * (y1 - qd12);

  double qr31 = qd31 - gamma2 * (x2 - qd21);
  double qr32 = qd32 - gamma2 * (y2 - qd22);

  double qr41 = qd41 - gamma2 * (x3 - qd31);
  double qr42 = qd42 - gamma2 * (y3 - qd32);

  double s1 = x1 - qr11;
  double s2 = y1 - qr12;

  double a1 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2);
  double a2 = a11 * qr11 + a12 * qr12;
  double a3 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * qr11 + d12 * qr12) / 2;
  double a4 = a21 * qr11 + a22 * qr12;

  double da1 = -m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * dv1 * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2);
  double da2 = a11 * qr21 + a12 * qr22 + da11 * qr11 + da12 * qr12;
  double da3 = m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[1] * (d11 * qr11 + d12 * qr12) / 2 + m2 * l1 * l2 * sin(q[1]) * dv0 * (d11 * qr11 + d12 * qr12) / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * (dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) / 2;
  double da4 = a21 * qr21 + a22 * qr22 + da21 * qr11 + da22 * qr12;

  double dda1 = (m2 * l1 * l2 * sin(q[1]) * velocity[1] * velocity[1] * velocity[1] - 2 * m2 * l1 * l2 * cos(q[1]) * velocity[1] * dv1) * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] * ((dd11 * qr11 + d11 * qr21 + dd12 * qr12 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - (m2 * l1 * l2 * cos(q[1]) * velocity[1] * dv1 + m2 * l1 * l2 * sin(q[1]) * ddv1) * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * dv1 * ((dd11 * qr11 + d11 * qr21 + dd12 * qr12 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - (m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] + m2 * l1 * l2 * sin(q[1]) * dv1) * ((dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((ddd11 * qr11 + dd11 * qr21 + ddd12 * qr12 + dd12 * qr22 + dd11 * qr21 + d11 * qr31 + dd12 * qr22 + d12 * qr32) + (ddd21 * qr11 + dd21 * qr21 + ddd22 * qr12 + dd22 * qr22 + dd21 * qr21 + d21 * qr31 + dd22 * qr22 + d22 * qr32) / 2);
  double dda2 = da11 * qr21 + a11 * qr31 + da12 * qr22 + a12 * qr32 + da11 * qr21 + dda11 * qr11 + da12 * qr22 + dda12 * qr12;
  double dda3 = (m2 * l1 * l2 * sin(q[1]) * ddv0 + 2 * m2 * l1 * l2 * cos(q[1]) * dv0 * velocity[1] + m2 * l1 * l2 * cos(q[1]) * velocity[0] * dv1 - m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[1] * velocity[1]) * (d11 * qr11 + d12 * qr12) / 2 + (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[1] + m2 * l1 * l2 * sin(q[1]) * dv0) * (dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + m2 * l1 * l2 * sin(q[1]) * velocity[0] * (ddd11 * qr11 + ddd12 * qr12 + dd11 * qr21 + dd12 * qr22 + dd11 * qr21 + dd12 * qr22 + d11 * qr31 + d12 * qr32) / 2;
  double dda4 = da21 * qr21 + a21 * qr31 + da22 * qr22 + a22 * qr32 + dda21 * qr11 + da21 * qr21 + dda22 * qr12 + da22 * qr22;

  double T01 = mc11 * qr21 + mc12 * qr22;
  double T02 = a1 + a2 + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2);
  double T03 = grad11 * gamma1 * s1 + grad21 * gamma1 * s2;

  double dT01 = dmc11 * qr21 + dmc12 * qr22 + mc11 * qr31 + mc12 * qr32;
  double dT02 = da1 + da2 - g * (l1 * sin(q[0]) * velocity[0] * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2);
  double dT03 = -y1 * gamma1 * s1 + x1 * gamma1 * s2 + grad11 * gamma1 * (x2 - qr21) + grad21 * gamma1 * (y2 - qr22);

  double ddT01 = ddmc11 * qr21 + ddmc12 * qr22 + 2 * dmc11 * qr31 + 2 * dmc12 * qr32 + mc11 * qr41 + mc12 * qr42;
  double ddT02 = dda1 + dda2 - g * (l1 * sin(q[0]) * dv0 * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (dv0 + dv1) / 2 + l1 * cos(q[0]) * velocity[0] * velocity[0] * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) / 2);
  double ddT03 = -y2 * gamma1 * s1 + x2 * gamma1 * s2 - 2 * y1 * gamma1 * (x2 - qr21) + 2 * x1 * gamma1 * (y2 - qr22) + grad11 * gamma1 * (x3 - qr31) + grad21 * gamma1 * (y3 - qr32);

  double T11 = mc21 * qr21 + mc22 * qr22;
  double T12 = a3 + a4 + g * m2 * l2 * cos(q[0] + q[1]) / 2;
  double T13 = grad12 * gamma1 * s1 + grad22 * gamma1 * s2;

  double dT11 = dmc21 * qr21 + dmc22 * qr22 + mc21 * qr31 + mc22 * qr32;
  double dT12 = da3 + da4 - g * m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2;
  double dT13 = -l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * s1 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * s2 + grad12 * gamma1 * (x2 - qr21) + grad22 * gamma1 * (y2 - qr22);

  double ddT11 = ddmc21 * qr21 + ddmc22 * qr22 + 2 * dmc21 * qr31 + 2 * dmc22 * qr32 + mc21 * qr41 + mc22 * qr42;
  double ddT12 = dda3 + dda4 - g * m2 * l2 * (sin(q[0] + q[1]) * (dv0 + dv1) / 2 + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) / 2);
  double ddT13 = l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * gamma1 * s1 - l2 * cos(q[0] + q[1]) * (dv0 + dv1) * gamma1 * s1 - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (x2 - qr21) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * gamma1 * s2 - l2 * sin(q[0] + q[1]) * (dv0 + dv1) * gamma1 * s2 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (y2 - qr22) + grad12 * gamma1 * (x3 - qr31) + grad22 * gamma1 * (y3 - qr32) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (x2 - qr21) - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (y2 - qr22);

  double ki1 = sqrt((qd1 * qd1 + qd2 * qd2 + l1 * l1 + l2 * l2) * (qd1 * qd1 + qd2 * qd2 + l1 * l1 + l2 * l2) - 2 * ((qd1 * qd1 + qd2 * qd2) * (qd1 * qd1 + qd2 * qd2) + l1 * l1 * l1 * l1 + l2 * l2 * l2 * l2));
  double ki2 = qd1 * qd1 + qd2 * qd2 + l1 * l1 - l2 * l2;
  double ki3 = qd1 * qd1 + qd2 * qd2 - l1 * l1 - l2 * l2;

  double xd = atan2(qd2, qd1) + atan2(ki1, ki2);
  double yd = -atan2(ki1, ki3);
  double xd1 = (qd1 * qd12 - qd2 * qd11) / (qd1 * qd1 + qd2 * qd2) + ((qd1 * qd11 + qd2 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2)) / (ki1 * (qd1 * qd1 + qd2 * qd2));
  double yd1 = 2 * (qd1 * qd11 + qd2 * qd12) / ki1;
  double xd2 = (((qd1 * qd21 + qd11 * qd11 + qd2 * qd22 + qd12 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2) - 2 * (qd1 * qd11 + qd2 * qd12) * (qd1 * qd11 + qd2 * qd12)) * ki1 * (qd1 * qd1 + qd2 * qd2) - (qd1 * qd11 + qd2 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2) * (ki1 * 2 * (qd1 * qd11 + qd2 * qd12) - (qd1 * qd1 + qd2 * qd2) * 2 * (qd1 * qd11 + qd2 * qd12) * ki3 / ki1)) / (ki1 * (qd1 * qd1 + qd2 * qd2) * ki1 * (qd1 * qd1 + qd2 * qd2)) + ((qd1 * qd22 - qd2 * qd21) * (qd1 * qd1 + qd2 * qd2) - (qd1 * qd12 - qd2 * qd11) * 2 * (qd1 * qd11 + qd2 * qd12)) / ((qd1 * qd1 + qd2 * qd2) * (qd1 * qd1 + qd2 * qd2));
  double yd2 = 2 * ((qd1 * qd21 + qd11 * qd11 + qd2 * qd22 + qd12 * qd12) * ki1 * ki1 + (qd1 * qd11 + qd2 * qd12) * 2 * (qd1 * qd11 + qd2 * qd12) * ki3) / (ki1 * ki1 * ki1);

  double thetad1 = xd + (T01 + T02 - T03) / K1;
  double thetad2 = yd + (T11 + T12 - T13) / K2;
  double thetad11 = xd1 + (dT01 + dT02 - dT03) / K1;
  double thetad12 = yd1 + (dT11 + dT12 - dT13) / K2;
  double thetad21 = xd2 + (ddT01 + ddT02 - ddT03) / K1;
  double thetad22 = yd2 + (ddT11 + ddT12 - ddT13) / K2;

  double s21 = (velocity[2] - thetad11) + gamma2 * (q[2] - thetad1);
  double s22 = (velocity[3] - thetad12) + gamma2 * (q[3] - thetad2);

  double thetar21 = thetad21 - gamma2 * (velocity[2] - thetad11);
  double thetar22 = thetad22 - gamma2 * (velocity[3] - thetad12);


  // control law
  U[0] = 0;//-(T01+T02-T03);//0;
  U[1] = 0;//-(T11+T12-T13);//0;
  U[2] = -(J1 * thetar21 + K1 * (thetad1 - xd) - gamma1 * s21);
  U[3] = -(J2 * thetar22 + K2 * (thetad2 - yd) - gamma1 * s22);

  double V1 = (s1 * s1 * (d11 * mc11 + d21 * mc21) + 2 * s1 * s2 * (d12 * mc11 + d22 * mc21) + s2 * s2 * (d12 * mc12 + d22 * mc22)) / 2 + (J1 * s21 * s21 + J2 * s22 * s22) / 2;
  z[6] = V1 + gamma1 * gamma2 * ((q[0] - xd) * (q[0] - xd) + (q[1] - yd) * (q[1] - yd) + (q[2] - thetad1) * (q[2] - thetad1) + (q[3] - thetad2) * (q[3] - thetad2)) + (q[0] - xd - q[2] + thetad1) * K1 * (q[0] - xd - q[2] + thetad1) / 2 + (q[1] - yd - q[3] + thetad2) * K2 * (q[1] - yd - q[3] + thetad2) / 2;


  // z[9] = V1;
  z[11] = P;
  z[14] = qd1;
  z[15] = qd2;
  z[16] = U[2];
  z[17] = U[3];
  z[19] = s1;
  z[20] = s2;
  z[21] = qr11;
  z[22] = qr12;
  z[23] = qr21;
  z[24] = qr22;
}

SICONOS_EXPORT void U2(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
{
  double m11 = m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = (m2 * l2 * l2 / 4) + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = (m2 * l2 * l2 / 4) + I2;
  double detm = m11 * m22 - m12 * m12;
  double ddetm = (m12 - m22) * m2 * l1 * l2 * sin(q[1]) * velocity[1];

  double FGyr0 = -m2 * l1 * l2 * sin(q[1]) * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2) + K1 * (q[0] - q[2]);
  double FGyr1 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[0] / 2 + g * m2 * l2 * cos(q[0] + q[1]) / 2 + K2 * (q[1] - q[3]);

  // acceleration in q[0],q[1] coordinates
  double dv0 = -(m22 * FGyr0 - m12 * FGyr1) / detm;
  double dv1 = -(m11 * FGyr1 - m12 * FGyr0) / detm;

  double ddv0 = -((m22 * (-m2 * l1 * l2 * cos(q[1]) * velocity[1] * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) - m2 * l1 * l2 * sin(q[1]) * (dv0 * velocity[1] + dv1 * velocity[0] + dv1 * velocity[1]) - g * (l1 * sin(q[0]) * velocity[0] * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2) + K1 * (velocity[0] - velocity[2])) + m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr1 / 2 - m12 * (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] * velocity[1] / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * dv0 + K2 * (velocity[1] - velocity[3]))) * detm - (m22 * FGyr0 - m12 * FGyr1) * ddetm) / (detm * detm);
  double ddv1 = -((m11 * (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] * velocity[1] / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * dv0 + K2 * (velocity[1] - velocity[3])) + m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr0 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr1 + m12 * (m2 * l1 * l2 * cos(q[1]) * velocity[1] * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) + m2 * l1 * l2 * sin(q[1]) * (dv0 * velocity[1] + dv1 * velocity[0] + dv1 * velocity[1]) + g * m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2 - K1 * (velocity[0] - velocity[2]))) * detm - (m11 * FGyr1 - m12 * FGyr0) * ddetm) / (detm * detm);

  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

  double x1 = -l1 * sin(q[0]) * velocity[0] - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]);
  double y1 = l1 * cos(q[0]) * velocity[1] + l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]);

  double x2 = -l1 * cos(q[0]) * velocity[0] * velocity[0] - l1 * sin(q[0]) * dv0 - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - l2 * sin(q[0] + q[1]) * (dv0 + dv1);
  double y2 = l1 * cos(q[0]) * dv0 - l1 * sin(q[0]) * velocity[0] * velocity[0] + l2 * cos(q[0] + q[1]) * (dv0 + dv1) - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]);

  double x3 = l1 * sin(q[0]) * velocity[0] * velocity[0] * velocity[0] - 2 * l1 * cos(q[0]) * velocity[0] * dv0 - l1 * cos(q[0]) * velocity[0] * dv0 - l1 * sin(q[0]) * ddv0 + l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - 2 * l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - l2 * sin(q[0] + q[1]) * (ddv0 + ddv1);
  double y3 = -l1 * sin(q[0]) * velocity[0] * dv0 + l1 * cos(q[0]) * ddv0 - l1 * cos(q[0]) * velocity[0] * velocity[0] * velocity[0] - 2 * l1 * sin(q[0]) * velocity[0] * dv0 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) + l2 * cos(q[0] + q[1]) * (ddv0 + ddv1) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - 2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1);

  double grad11 = -y;
  double grad12 = -l2 * sin(q[0] + q[1]);
  double grad21 = x;
  double grad22 = l2 * cos(q[0] + q[1]);
  double det = grad11 * grad22 - grad12 * grad21;

  double d11 = grad22 / det;
  double d12 = -grad12 / det;
  double d21 = -grad21 / det;
  double d22 = grad11 / det;

  double dd11 =  -(sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd12 = (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd21 =  -dd11 + (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));
  double dd22 =  -dd12 - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));

  double ddd11 = (-(sin(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] + cos(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd12 = ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + sin(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] - sin(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd21 =  -ddd11 + ((cos(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + sin(q[0]) * sin(q[1]) * dv0 + cos(q[0]) * cos(q[1]) * dv1) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd22 =  -ddd12 - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]));

  double dddd11 = ((-(sin(q[0] + q[1]) * (ddv0 + ddv1) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * ddv1 + 3 * cos(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - 3 * cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * dv1 + cos(q[0] + q[1]) * cos(q[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * velocity[1] - velocity[1] * velocity[1] * velocity[1]) - sin(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - velocity[1] * velocity[1]) + sin(q[0] + q[1]) * cos(q[1]) * ((dv0 + dv1) * velocity[1] - (velocity[0] + velocity[1]) * dv1)) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - (-(sin(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] + cos(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd12 = (((cos(q[0] + q[1]) * (ddv0 + ddv1) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * ddv1 - 3 * sin(q[0] + q[1]) * sin(q[1]) * ((velocity[0] + velocity[1]) * (dv0 + dv1) - velocity[1] * dv1) + cos(q[0] + q[1]) * cos(q[1]) * ((dv0 + dv1) * velocity[1] - (velocity[0] + velocity[1]) * dv1) - (cos(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) + sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + sin(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] - sin(q[0] + q[1]) * cos(q[1]) * dv1) * (sin(q[1]) * sin(q[1])) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * sin(q[1]) * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd21 = -dddd11 + (((sin(q[0]) * sin(q[1]) * ddv0 + cos(q[0]) * cos(q[1]) * ddv1 + 3 * (velocity[0] * dv0 - velocity[1] * dv1) + sin(q[0]) * cos(q[1]) * (velocity[1] * dv0 - velocity[0] * dv1) + (cos(q[0]) * cos(q[1]) * velocity[1] - sin(q[0]) * sin(q[1]) * velocity[0]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + sin(q[0]) * sin(q[1]) * dv0 + cos(q[0]) * cos(q[1]) * dv1) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd22 = -dddd12 - (((-(sin(q[0]) * cos(q[1]) * velocity[1] + cos(q[0]) * sin(q[1]) * velocity[0]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + cos(q[0]) * sin(q[1]) * ddv0 - sin(q[0]) * cos(q[1]) * ddv1 + 3 * (velocity[1] * dv1 - velocity[0] * dv0) * sin(q[0]) * sin(q[1]) - cos(q[0]) * cos(q[1]) * (velocity[1] * dv0 - velocity[0] * dv1)) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1])) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));

  double a11 = m11 * dd11 + m12 * dd21;
  double a12 = m11 * dd12 + m12 * dd22;
  double a21 = m12 * dd11 + m22 * dd21;
  double a22 = m12 * dd12 + m22 * dd22;

  double da11 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd21 / 2 + m11 * ddd11 + m12 * ddd21;
  double da12 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2) + m11 * ddd12 + m12 * ddd22;
  double da21 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 / 2 + m12 * ddd11 + m22 * ddd21;
  double da22 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 / 2 + m12 * ddd12 + m22 * ddd22;

  double dda11 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd11 + sin(q[1]) * velocity[1] * ddd11 + (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd21 / 2 + sin(q[1]) * velocity[1] * ddd21 / 2 + sin(q[1]) * velocity[1] * (ddd11 + ddd21 / 2)) + m11 * dddd11 + m12 * dddd21;
  double dda12 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (dd12 + dd22 / 2) + sin(q[1]) * velocity[1] * (ddd12 + ddd22 / 2) + sin(q[1]) * velocity[1] * (ddd12 + ddd22 / 2)) + m11 * dddd12 + m12 * dddd22;
  double dda21 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd11 / 2 + sin(q[1]) * velocity[1] * ddd11) + m12 * dddd11 + m22 * dddd21;
  double dda22 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd12 / 2 + sin(q[1]) * velocity[1] * ddd12) + m12 * dddd12 + m22 * dddd22;

  double mc11 = m11 * d11 + m12 * d21;
  double mc12 = m11 * d12 + m12 * d22;
  double mc21 = m12 * d11 + m22 * d21;
  double mc22 = m12 * d12 + m22 * d22;

  double dmc11 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (d11 + d21 / 2) + m11 * dd11 + m12 * dd21;
  double dmc12 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (d12 + d22 / 2) + m11 * dd12 + m12 * dd22;
  double dmc21 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * d11 / 2 + m12 * dd11 + m22 * dd21;
  double dmc22 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * d12 / 2 + m12 * dd12 + m22 * dd22;

  double ddmc11 = -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (d11 + d21 / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd11 + dd21 / 2) + m11 * ddd11 + m12 * ddd21 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd11 + dd21 / 2);
  double ddmc12 = -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (d12 + d22 / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2) + m11 * ddd12 + m12 * ddd22 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2);
  double ddmc21 =  -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * d11 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 / 2 + m12 * ddd11 + m22 * ddd21 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 / 2;
  double ddmc22 =  -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * d12 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 / 2 + m12 * ddd12 + m22 * ddd22 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 / 2;

  double qr11 = -gamma2 * (x - z[5]);
  double qr12 = -gamma2 * y;

  double qr21 = -gamma2 * x1;
  double qr22 = -gamma2 * y1;

  double qr31 = -gamma2 * x2;
  double qr32 = -gamma2 * y2;

  double qr41 = -gamma2 * x3;
  double qr42 = -gamma2 * y3;

  double s1 = x1 - qr11;
  double s2 = y1 - qr12;
  double s2bar = y1 + gamma2 * (y + sqrt(z[7]) * alpha);

  double a1 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2);
  double a2 = a11 * qr11 + a12 * qr12;
  double a3 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * qr11 + d12 * qr12) / 2;
  double a4 = a21 * qr11 + a22 * qr12;

  double da1 = -m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * dv1 * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2);
  double da2 = a11 * qr21 + a12 * qr22 + da11 * qr11 + da12 * qr12;
  double da3 = m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[1] * (d11 * qr11 + d12 * qr12) / 2 + m2 * l1 * l2 * sin(q[1]) * dv0 * (d11 * qr11 + d12 * qr12) / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * (dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) / 2;
  double da4 = a21 * qr21 + a22 * qr22 + da21 * qr11 + da22 * qr12;

  double dda1 = (m2 * l1 * l2 * sin(q[1]) * velocity[1] * velocity[1] * velocity[1] - 2 * m2 * l1 * l2 * cos(q[1]) * velocity[1] * dv1) * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] * ((dd11 * qr11 + d11 * qr21 + dd12 * qr12 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - (m2 * l1 * l2 * cos(q[1]) * velocity[1] * dv1 + m2 * l1 * l2 * sin(q[1]) * ddv1) * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * dv1 * ((dd11 * qr11 + d11 * qr21 + dd12 * qr12 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - (m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] + m2 * l1 * l2 * sin(q[1]) * dv1) * ((dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((ddd11 * qr11 + dd11 * qr21 + ddd12 * qr12 + dd12 * qr22 + dd11 * qr21 + d11 * qr31 + dd12 * qr22 + d12 * qr32) + (ddd21 * qr11 + dd21 * qr21 + ddd22 * qr12 + dd22 * qr22 + dd21 * qr21 + d21 * qr31 + dd22 * qr22 + d22 * qr32) / 2);
  double dda2 = da11 * qr21 + a11 * qr31 + da12 * qr22 + a12 * qr32 + da11 * qr21 + dda11 * qr11 + da12 * qr22 + dda12 * qr12;
  double dda3 = (m2 * l1 * l2 * sin(q[1]) * ddv0 + 2 * m2 * l1 * l2 * cos(q[1]) * dv0 * velocity[1] + m2 * l1 * l2 * cos(q[1]) * velocity[0] * dv1 - m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[1] * velocity[1]) * (d11 * qr11 + d12 * qr12) / 2 + (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[1] + m2 * l1 * l2 * sin(q[1]) * dv0) * (dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + m2 * l1 * l2 * sin(q[1]) * velocity[0] * (ddd11 * qr11 + ddd12 * qr12 + dd11 * qr21 + dd12 * qr22 + dd11 * qr21 + dd12 * qr22 + d11 * qr31 + d12 * qr32) / 2;
  double dda4 = da21 * qr21 + a21 * qr31 + da22 * qr22 + a22 * qr32 + dda21 * qr11 + da21 * qr21 + dda22 * qr12 + da22 * qr22;


  double T01 = mc11 * qr21 + mc12 * qr22;
  double T02 = a1 + a2 + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2);
  double T03 = grad11 * gamma1 * s1 + grad21 * gamma1 * s2bar;

  double dT01 = dmc11 * qr21 + dmc12 * qr22 + mc11 * qr31 + mc12 * qr32;
  double dT02 = da1 + da2 - g * (l1 * sin(q[0]) * velocity[0] * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2);
  double dT03 = -y1 * gamma1 * s1 + x1 * gamma1 * s2bar + grad11 * gamma1 * (x2 - qr21) + grad21 * gamma1 * (y2 - qr22);

  double ddT01 = ddmc11 * qr21 + ddmc12 * qr22 + 2 * dmc11 * qr31 + 2 * dmc12 * qr32 + mc11 * qr41 + mc12 * qr42;
  double ddT02 = dda1 + dda2 - g * (l1 * sin(q[0]) * dv0 * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (dv0 + dv1) / 2 + l1 * cos(q[0]) * velocity[0] * velocity[0] * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) / 2);
  double ddT03 = -y2 * gamma1 * s1 + x2 * gamma1 * s2bar - 2 * y1 * gamma1 * (x2 - qr21) + 2 * x1 * gamma1 * (y2 - qr22) + grad11 * gamma1 * (x3 - qr31) + grad21 * gamma1 * (y3 - qr32);

  double T11 = mc21 * qr21 + mc22 * qr22;
  double T12 = a3 + a4 + g * m2 * l2 * cos(q[0] + q[1]) / 2;
  double T13 = grad12 * gamma1 * s1 + grad22 * gamma1 * s2bar;

  double dT11 = dmc21 * qr21 + dmc22 * qr22 + mc21 * qr31 + mc22 * qr32;
  double dT12 = da3 + da4 - g * m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2;
  double dT13 = -l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * s1 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * s2bar + grad12 * gamma1 * (x2 - qr21) + grad22 * gamma1 * (y2 + gamma2 * y1);

  double ddT11 = ddmc21 * qr21 + ddmc22 * qr22 + 2 * dmc21 * qr31 + 2 * dmc22 * qr32 + mc21 * qr41 + mc22 * qr42;
  double ddT12 = dda3 + dda4 - g * m2 * l2 * (sin(q[0] + q[1]) * (dv0 + dv1) / 2 + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) / 2);
  double ddT13 = l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * gamma1 * s1 - l2 * cos(q[0] + q[1]) * (dv0 + dv1) * gamma1 * s1 - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (x2 - qr21) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * gamma1 * s2 - l2 * sin(q[0] + q[1]) * (dv0 + dv1) * gamma1 * s2 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (y2 - qr22) + grad12 * gamma1 * (x3 - qr31) + grad22 * gamma1 * (y3 - qr32) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (x2 - qr21) - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (y2 - qr22);

  double qd1 = z[5];
  double qd2 = 0;
  double qd11 = 0;
  double qd12 = 0;
  double qd21 = 0;
  double qd22 = 0;

  double ki1 = sqrt((qd1 * qd1 + qd2 * qd2 + l1 * l1 + l2 * l2) * (qd1 * qd1 + qd2 * qd2 + l1 * l1 + l2 * l2) - 2 * ((qd1 * qd1 + qd2 * qd2) * (qd1 * qd1 + qd2 * qd2) + l1 * l1 * l1 * l1 + l2 * l2 * l2 * l2));
  double ki2 = qd1 * qd1 + qd2 * qd2 + l1 * l1 - l2 * l2;
  double ki3 = qd1 * qd1 + qd2 * qd2 - l1 * l1 - l2 * l2;

  double xd = atan2(qd2, qd1) + atan2(ki1, ki2);
  double yd = -atan2(ki1, ki3);
  double xd1 = (qd1 * qd12 - qd2 * qd11) / (qd1 * qd1 + qd2 * qd2) + ((qd1 * qd11 + qd2 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2)) / (ki1 * (qd1 * qd1 + qd2 * qd2));
  double yd1 = 2 * (qd1 * qd11 + qd2 * qd12) / ki1;
  double xd2 = (((qd1 * qd21 + qd11 * qd11 + qd2 * qd22 + qd12 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2) - 2 * (qd1 * qd11 + qd2 * qd12) * (qd1 * qd11 + qd2 * qd12)) * ki1 * (qd1 * qd1 + qd2 * qd2) - (qd1 * qd11 + qd2 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2) * (ki1 * 2 * (qd1 * qd11 + qd2 * qd12) - (qd1 * qd1 + qd2 * qd2) * 2 * (qd1 * qd11 + qd2 * qd12) * ki3 / ki1)) / (ki1 * (qd1 * qd1 + qd2 * qd2) * ki1 * (qd1 * qd1 + qd2 * qd2)) + ((qd1 * qd22 - qd2 * qd21) * (qd1 * qd1 + qd2 * qd2) - (qd1 * qd12 - qd2 * qd11) * 2 * (qd1 * qd11 + qd2 * qd12)) / ((qd1 * qd1 + qd2 * qd2) * (qd1 * qd1 + qd2 * qd2));
  double yd2 = 2 * ((qd1 * qd21 + qd11 * qd11 + qd2 * qd22 + qd12 * qd12) * ki1 * ki1 + (qd1 * qd11 + qd2 * qd12) * 2 * (qd1 * qd11 + qd2 * qd12) * ki3) / (ki1 * ki1 * ki1);

  double thetad1 = xd + (T01 + T02 - T03) / K1;
  double thetad2 = yd + (T11 + T12 - T13) / K2;
  double thetad11 = xd1 + (dT01 + dT02 - dT03) / K1;
  double thetad12 = yd1 + (dT11 + dT12 - dT13) / K2;
  double thetad21 = xd2 + (ddT01 + ddT02 - ddT03) / K1;
  double thetad22 = yd2 + (ddT11 + ddT12 - ddT13) / K2;

  double s21 = (velocity[2] - thetad11) + gamma2 * (q[2] - thetad1);
  double s22 = (velocity[3] - thetad12) + gamma2 * (q[3] - thetad2);

  double thetar21 = thetad21 - gamma2 * (velocity[2] - thetad11);
  double thetar22 = thetad22 - gamma2 * (velocity[3] - thetad12);

  U[0] = 0;//-(T01+T02-T03);//0;
  U[1] = 0;//-(T11+T12-T13);//0;
  U[2] = -(thetar21 + K1 * (thetad1 - xd) - gamma1 * s21);
  U[3] = -(thetar22 + K2 * (thetad2 - yd) - gamma1 * s22);

  double V1 = (s1 * s1 * (d11 * mc11 + d21 * mc21) + 2 * s1 * s2 * (d12 * mc11 + d22 * mc21) + s2 * s2 * (d12 * mc12 + d22 * mc22)) / 2 + J1 * s21 * s21 + J2 * s22 * s22;
  z[6] = V1 + gamma1 * gamma2 * ((q[0] - xd) * (q[0] - xd) + (q[1] - yd) * (q[1] - yd) + (q[2] - thetad1) * (q[2] - thetad1) + (q[3] - thetad2) * (q[3] - thetad2)) + (q[0] - xd - q[2] + thetad1) * K1 * (q[0] - xd - q[2] + thetad1) / 2 + (q[1] - yd - q[3] + thetad2) * K2 * (q[1] - yd - q[3] + thetad2) / 2;

  z[11] = P;
  z[14] = qd1;
  z[15] = -sqrt(z[7]) * alpha;
  z[16] = U[2];
  z[17] = U[3];
  z[19] = s1;
  z[20] = s2;
  z[21] = qr11;
  z[22] = qr12;
  z[23] = qr21;
  z[24] = qr22;
}

SICONOS_EXPORT void U3(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
{
  double m11 = m1 * l1 * l1 / 4 + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = m2 * l2 * l2 / 4 + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = m2 * l2 * l2 / 4 + I2;
  double detm = m11 * m22 - m12 * m12;
  double ddetm = (m12 - m22) * m2 * l1 * l2 * sin(q[1]) * velocity[1];

  double FGyr0 = -m2 * l1 * l2 * sin(q[1]) * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2) + K1 * (q[0] - q[2]);
  double FGyr1 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[0] / 2 + g * m2 * l2 * cos(q[0] + q[1]) / 2 + K2 * (q[1] - q[3]);

  // acceleration in q[0],q[1] coordinates
  double dv0 = -(m22 * FGyr0 - m12 * FGyr1) / detm;
  double dv1 = -(m11 * FGyr1 - m12 * FGyr0) / detm;

  double ddv0 = -((m22 * (-m2 * l1 * l2 * cos(q[1]) * velocity[1] * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) - m2 * l1 * l2 * sin(q[1]) * (dv0 * velocity[1] + dv1 * velocity[0] + dv1 * velocity[1]) - g * (l1 * sin(q[0]) * velocity[0] * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2) + K1 * (velocity[0] - velocity[2])) + m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr1 / 2 - m12 * (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] * velocity[1] / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * dv0 + K2 * (velocity[1] - velocity[3]))) * detm - (m22 * FGyr0 - m12 * FGyr1) * ddetm) / (detm * detm);
  double ddv1 = -((m11 * (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] * velocity[1] / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * dv0 + K2 * (velocity[1] - velocity[3])) + m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr0 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr1 + m12 * (m2 * l1 * l2 * cos(q[1]) * velocity[1] * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) + m2 * l1 * l2 * sin(q[1]) * (dv0 * velocity[1] + dv1 * velocity[0] + dv1 * velocity[1]) + g * m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2 - K1 * (velocity[0] - velocity[2]))) * detm - (m11 * FGyr1 - m12 * FGyr0) * ddetm) / (detm * detm);

  // generalized coordinates q=(x,y)^T
  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

  // time derivatives of x and y
  double x1 = -l1 * sin(q[0]) * velocity[0] - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]);
  double y1 = l1 * cos(q[0]) * velocity[0] + l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]);

  // acceleration in x,y coordinates
  double x2 =  -l1 * cos(q[0]) * velocity[0] * velocity[0] - l1 * sin(q[0]) * dv0 - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - l2 * sin(q[0] + q[1]) * (dv0 + dv1);
  double y2 =  l1 * cos(q[0]) * dv0 - l1 * sin(q[0]) * velocity[0] * velocity[0] + l2 * cos(q[0] + q[1]) * (dv0 + dv1) - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]);

  double x3 =  l1 * sin(q[0]) * velocity[0] * velocity[0] * velocity[0] - 2 * l1 * cos(q[0]) * velocity[0] * dv0 - l1 * cos(q[0]) * velocity[0] * dv0 - l1 * sin(q[0]) * ddv0 + l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - 2 * l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - l2 * sin(q[0] + q[1]) * (ddv0 + ddv1);
  double y3 = -l1 * sin(q[0]) * velocity[0] * dv0 + l1 * cos(q[0]) * ddv0 - l1 * cos(q[0]) * velocity[0] * velocity[0] * velocity[0] - 2 * l1 * sin(q[0]) * velocity[0] * dv0 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) + l2 * cos(q[0] + q[1]) * (ddv0 + ddv1) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - 2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1);

  // the gradient of the transformation \theta->q
  double grad11 = -y;
  double grad12 = -l2 * sin(q[0] + q[1]);
  double grad21 = x;
  double grad22 = l2 * cos(q[0] + q[1]);
  double det = grad11 * grad22 - grad12 * grad21;

  // d=(grad)^{-1}
  double d11 = grad22 / det;
  double d12 = -grad12 / det;
  double d21 = -grad21 / det;
  double d22 = grad11 / det;

  // time derivative of d

  double dd11 =  -(sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd12 = (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd21 =  -dd11 + (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));
  double dd22 =  -dd12 - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));

  // time derivative of dd

  double ddd11 = (-(sin(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] + cos(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd12 = ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + sin(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] - sin(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd21 =  -ddd11 + ((cos(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + sin(q[0]) * sin(q[1]) * dv0 + cos(q[0]) * cos(q[1]) * dv1) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd22 =  -ddd12 - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]));

  // time derivative of ddd

  double dddd11 = ((-(sin(q[0] + q[1]) * (ddv0 + ddv1) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * ddv1 + 3 * cos(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - 3 * cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * dv1 + cos(q[0] + q[1]) * cos(q[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * velocity[1] - velocity[1] * velocity[1] * velocity[1]) - sin(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - velocity[1] * velocity[1]) + sin(q[0] + q[1]) * cos(q[1]) * ((dv0 + dv1) * velocity[1] - (velocity[0] + velocity[1]) * dv1)) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - (-(sin(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] + cos(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd12 = (((cos(q[0] + q[1]) * (ddv0 + ddv1) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * ddv1 - 3 * sin(q[0] + q[1]) * sin(q[1]) * ((velocity[0] + velocity[1]) * (dv0 + dv1) - velocity[1] * dv1) + cos(q[0] + q[1]) * cos(q[1]) * ((dv0 + dv1) * velocity[1] - (velocity[0] + velocity[1]) * dv1) - (cos(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) + sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + sin(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] - sin(q[0] + q[1]) * cos(q[1]) * dv1) * (sin(q[1]) * sin(q[1])) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * sin(q[1]) * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd21 = -dddd11 + (((sin(q[0]) * sin(q[1]) * ddv0 + cos(q[0]) * cos(q[1]) * ddv1 + 3 * (velocity[0] * dv0 - velocity[1] * dv1) + sin(q[0]) * cos(q[1]) * (velocity[1] * dv0 - velocity[0] * dv1) + (cos(q[0]) * cos(q[1]) * velocity[1] - sin(q[0]) * sin(q[1]) * velocity[0]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + sin(q[0]) * sin(q[1]) * dv0 + cos(q[0]) * cos(q[1]) * dv1) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd22 = -dddd12 - (((-(sin(q[0]) * cos(q[1]) * velocity[1] + cos(q[0]) * sin(q[1]) * velocity[0]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + cos(q[0]) * sin(q[1]) * ddv0 - sin(q[0]) * cos(q[1]) * ddv1 + 3 * (velocity[1] * dv1 - velocity[0] * dv0) * sin(q[0]) * sin(q[1]) - cos(q[0]) * cos(q[1]) * (velocity[1] * dv0 - velocity[0] * dv1)) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1])) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));

  // M*d
  double mc11 = m11 * d11 + m12 * d21;
  double mc12 = m11 * d12 + m12 * d22;
  double mc21 = m12 * d11 + m22 * d21;
  double mc22 = m12 * d12 + m22 * d22;

  // M*dd
  double a11 = m11 * dd11 + m12 * dd21;
  double a12 = m11 * dd12 + m12 * dd22;
  double a21 = m12 * dd11 + m22 * dd21;
  double a22 = m12 * dd12 + m22 * dd22;

  double da11 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd11 + dd21 / 2) + m11 * ddd11 + m12 * ddd21;
  double da12 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2) + m11 * ddd12 + m12 * ddd22;
  double da21 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 / 2 + m12 * ddd11 + m22 * ddd21;
  double da22 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 / 2 + m12 * ddd12 + m22 * ddd22;

  double dda11 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (dd11 + dd21 / 2) + 2 * sin(q[1]) * velocity[1] * (ddd11 + ddd21 / 2)) + m11 * dddd11 + m12 * dddd21;
  double dda12 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (dd12 + dd22 / 2) + 2 * sin(q[1]) * velocity[1] * (ddd12 + ddd22 / 2)) + m11 * dddd12 + m12 * dddd22;
  double dda21 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd11 / 2 + sin(q[1]) * velocity[1] * ddd11) + m12 * dddd11 + m22 * dddd21;
  double dda22 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd12 / 2 + sin(q[1]) * velocity[1] * ddd12) + m12 * dddd12 + m22 * dddd22;

  double dmc11 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (d11 + d21 / 2) + m11 * dd11 + m12 * dd21;
  double dmc12 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (d12 + d22 / 2) + m11 * dd12 + m12 * dd22;
  double dmc21 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * d11 / 2 + m12 * dd11 + m22 * dd21;
  double dmc22 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * d12 / 2 + m12 * dd12 + m22 * dd22;

  double ddmc11 = -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (d11 + d21 / 2) - 2 * m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd11 + dd21 / 2) + m11 * ddd11 + m12 * ddd21;
  double ddmc12 = -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (d12 + d22 / 2) - 2 * m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2) + m11 * ddd12 + m12 * ddd22;
  double ddmc21 =  -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * d11 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 + m12 * ddd11 + m22 * ddd21;
  double ddmc22 =  -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * d12 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 + m12 * ddd12 + m22 * ddd22;

  double td = z[8] + P / 5;

  double c0 = z[5];
  double c1 = 0.7 + (sqrt(z[7]) / 10) * alpha;
  double w0 = (c1 - c0) / ((P / 5) * (P / 5) * (P / 5) * (P / 5));
  double w1 = -4 * (c1 - c0) / ((P / 5) * (P / 5) * (P / 5) * (P / 5) * (P / 5));
  double w2 = 10 * (c1 - c0) / ((P / 5) * (P / 5) * (P / 5) * (P / 5) * (P / 5) * (P / 5));
  double w3 = -20 * (c1 - c0) / ((P / 5) * (P / 5) * (P / 5) * (P / 5) * (P / 5) * (P / 5) * (P / 5));

  double qd1 = 0;
  double qd11 = 0;
  double qd21 = 0;
  double qd31 = 0;
  double qd41 = 0;

  if (time <= td)
  {
    qd1 = c0 + w0 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) + w1 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - td) + w2 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - td) * (time - td) + w3 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - td) * (time - td) * (time - td);
    qd11 = (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (w1 + 2 * w2 * (time - td) + 3 * w3 * (time - td) * (time - td)) + 4 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (w0 + w1 * (time - td) + w2 * (time - td) * (time - td) + w3 * (time - td) * (time - td) * (time - td));
    qd21 = (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (2 * w2 + 6 * w3 * (time - td)) + 8 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (w1 + 2 * w2 * (time - td) + 3 * w3 * (time - td) * (time - td)) + 12 * (time - z[8]) * (time - z[8]) * (w0 + w1 * (time - td) + w2 * (time - td) * (time - td) + w3 * (time - td) * (time - td) * (time - td));
    qd31 = (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * 6 * w3 + 12 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (2 * w2 + 6 * w3 * (time - td)) + 36 * (time - z[8]) * (time - z[8]) * (w1 + 2 * w2 * (time - td) + 3 * w3 * (time - td) * (time - td)) + 24 * (time - z[8]) * (w0 + w1 * (time - td) + w2 * (time - td) * (time - td) + w3 * (time - td) * (time - td) * (time - td));
  }
  else
  {
    qd1 = 0.7 + (sqrt(z[7]) / 10) * alpha;
    qd11 = 0;
    qd21 = 0;
    qd31 = 0;
  }
  double qd2 = 0;
  double qd12 = 0;
  double qd22 = 0;
  double qd32 = 0;
  double qd42 = 0;

  double qr11 = qd11 - gamma2 * (x - qd1);
  double qr12 = -gamma2 * y;

  double qr21 = qd21 - gamma2 * (x1 - qd11);
  double qr22 = -gamma2 * y1;

  double qr31 = qd31 - gamma2 * (x2 - qd21);
  double qr32 = -gamma2 * y2;

  double qr41 = qd41 - gamma2 * (x3 - qd31);
  double qr42 = -gamma2 * y3;

  double s1 = x1 - qr11;
  double s2 = y1 - qr12;

  double ki1 = sqrt((qd1 * qd1 + qd2 * qd2 + l1 * l1 + l2 * l2) * (qd1 * qd1 + qd2 * qd2 + l1 * l1 + l2 * l2) - 2 * ((qd1 * qd1 + qd2 * qd2) * (qd1 * qd1 + qd2 * qd2) + l1 * l1 * l1 * l1 + l2 * l2 * l2 * l2));
  double ki2 = qd1 * qd1 + qd2 * qd2 + l1 * l1 - l2 * l2;
  double ki3 = qd1 * qd1 + qd2 * qd2 - l1 * l1 - l2 * l2;

  double xd = atan2(qd2, qd1) + atan2(ki1, ki2);
  double yd = -atan2(ki1, ki3);
  double xd1 = (qd1 * qd12 - qd2 * qd11) / (qd1 * qd1 + qd2 * qd2) + ((qd1 * qd11 + qd2 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2)) / (ki1 * (qd1 * qd1 + qd2 * qd2));
  double yd1 = 2 * (qd1 * qd11 + qd2 * qd12) / ki1;
  double xd2 = (((qd1 * qd21 + qd11 * qd11 + qd2 * qd22 + qd12 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2) - 2 * (qd1 * qd11 + qd2 * qd12) * (qd1 * qd11 + qd2 * qd12)) * ki1 * (qd1 * qd1 + qd2 * qd2) - (qd1 * qd11 + qd2 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2) * (ki1 * 2 * (qd1 * qd11 + qd2 * qd12) - (qd1 * qd1 + qd2 * qd2) * 2 * (qd1 * qd11 + qd2 * qd12) * ki3 / ki1)) / (ki1 * (qd1 * qd1 + qd2 * qd2) * ki1 * (qd1 * qd1 + qd2 * qd2)) + ((qd1 * qd22 - qd2 * qd21) * (qd1 * qd1 + qd2 * qd2) - (qd1 * qd12 - qd2 * qd11) * 2 * (qd1 * qd11 + qd2 * qd12)) / ((qd1 * qd1 + qd2 * qd2) * (qd1 * qd1 + qd2 * qd2));
  double yd2 = 2 * ((qd1 * qd21 + qd11 * qd11 + qd2 * qd22 + qd12 * qd12) * ki1 * ki1 + (qd1 * qd11 + qd2 * qd12) * 2 * (qd1 * qd11 + qd2 * qd12) * ki3) / (ki1 * ki1 * ki1);

  double a1 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2);
  double a2 = a11 * qr11 + a12 * qr12;
  double a3 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * qr11 + d12 * qr12) / 2;
  double a4 = a21 * qr11 + a22 * qr12;

  double da1 = -m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * dv1 * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2);
  double da2 = a11 * qr21 + a12 * qr22 + da11 * qr11 + da12 * qr12;
  double da3 = m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[1] * (d11 * qr11 + d12 * qr12) / 2 + m2 * l1 * l2 * sin(q[1]) * dv0 * (d11 * qr11 + d12 * qr12) / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * (dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) / 2;
  double da4 = a21 * qr21 + a22 * qr22 + da21 * qr11 + da22 * qr12;

  double dda1 = (m2 * l1 * l2 * sin(q[1]) * velocity[1] * velocity[1] * velocity[1] - 2 * m2 * l1 * l2 * cos(q[1]) * velocity[1] * dv1) * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] * ((dd11 * qr11 + d11 * qr21 + dd12 * qr12 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - (m2 * l1 * l2 * cos(q[1]) * velocity[1] * dv1 + m2 * l1 * l2 * sin(q[1]) * ddv1) * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * dv1 * ((dd11 * qr11 + d11 * qr21 + dd12 * qr12 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - (m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] + m2 * l1 * l2 * sin(q[1]) * dv1) * ((dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((ddd11 * qr11 + dd11 * qr21 + ddd12 * qr12 + dd12 * qr22 + dd11 * qr21 + d11 * qr31 + dd12 * qr22 + d12 * qr32) + (ddd21 * qr11 + dd21 * qr21 + ddd22 * qr12 + dd22 * qr22 + dd21 * qr21 + d21 * qr31 + dd22 * qr22 + d22 * qr32) / 2);
  double dda2 = da11 * qr21 + a11 * qr31 + da12 * qr22 + a12 * qr32 + da11 * qr21 + dda11 * qr11 + da12 * qr22 + dda12 * qr12;
  double dda3 = (m2 * l1 * l2 * sin(q[1]) * ddv0 + 2 * m2 * l1 * l2 * cos(q[1]) * dv0 * velocity[1] + m2 * l1 * l2 * cos(q[1]) * velocity[0] * dv1 - m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[1] * velocity[1]) * (d11 * qr11 + d12 * qr12) / 2 + (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[1] + m2 * l1 * l2 * sin(q[1]) * dv0) * (dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + m2 * l1 * l2 * sin(q[1]) * velocity[0] * (ddd11 * qr11 + ddd12 * qr12 + dd11 * qr21 + dd12 * qr22 + dd11 * qr21 + dd12 * qr22 + d11 * qr31 + d12 * qr32) / 2;
  double dda4 = da21 * qr21 + a21 * qr31 + da22 * qr22 + a22 * qr32 + dda21 * qr11 + da21 * qr21 + dda22 * qr12 + da22 * qr22;

  double T01 = mc11 * qr21 + mc12 * qr22;
  double T02 = a1 + a2 + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2);
  double T03 = grad11 * gamma1 * s1 + grad21 * gamma1 * s2;

  double dT01 = dmc11 * qr21 + dmc12 * qr22 + mc11 * qr31 + mc12 * qr32;
  double dT02 = da1 + da2 - g * (l1 * sin(q[0]) * velocity[0] * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2);
  double dT03 = -y1 * gamma1 * s1 + x1 * gamma1 * s2 + grad11 * gamma1 * (x2 - qr21) + grad21 * gamma1 * (y2 - qr22);

  double ddT01 = ddmc11 * qr21 + ddmc12 * qr22 + 2 * dmc11 * qr31 + 2 * dmc12 * qr32 + mc11 * qr41 + mc12 * qr42;
  double ddT02 = dda1 + dda2 - g * (l1 * sin(q[0]) * dv0 * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (dv0 + dv1) / 2 + l1 * cos(q[0]) * velocity[0] * velocity[0] * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) / 2);
  double ddT03 = -y2 * gamma1 * s1 + x2 * gamma1 * s2 - 2 * y1 * gamma1 * (x2 - qr21) + 2 * x1 * gamma1 * (y2 - qr22) + grad11 * gamma1 * (x3 - qr31) + grad21 * gamma1 * (y3 - qr32);

  double T11 = mc21 * qr21 + mc22 * qr22;
  double T12 = a3 + a4 + g * m2 * l2 * cos(q[0] + q[1]) / 2;
  double T13 = grad12 * gamma1 * s1 + grad22 * gamma1 * s2;

  double dT11 = dmc21 * qr21 + dmc22 * qr22 + mc21 * qr31 + mc22 * qr32;
  double dT12 = da3 + da4 - g * m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2;
  double dT13 = -l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * s1 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * s2 + grad12 * gamma1 * (x2 - qr21) + grad22 * gamma1 * (y2 - qr22);

  double ddT11 = ddmc21 * qr21 + ddmc22 * qr22 + 2 * dmc21 * qr31 + 2 * dmc22 * qr32 + mc21 * qr41 + mc22 * qr42;
  double ddT12 = dda3 + dda4 - g * m2 * l2 * (sin(q[0] + q[1]) * (dv0 + dv1) / 2 + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) / 2);
  double ddT13 = l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * gamma1 * s1 - l2 * cos(q[0] + q[1]) * (dv0 + dv1) * gamma1 * s1 - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (x2 - qr21) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * gamma1 * s2 - l2 * sin(q[0] + q[1]) * (dv0 + dv1) * gamma1 * s2 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (y2 - qr22) + grad12 * gamma1 * (x3 - qr31) + grad22 * gamma1 * (y3 - qr32) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (x2 - qr21) - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (y2 - qr22);

  double wd = (td - time + fabs(td - time)) / 2;

  double ld = (m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d12 - d22 / 2) - a21 - (mc21 / mc11) * (m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 - d21 / 2) - a11)) * s1 - gamma1 * mc21 * s1 / mc11 + (1 + Kf) * wd;

  double thetad1 = xd + (T01 + T02 - T03 + grad21 * (Kf * z[4] - ld)) / K1;
  double thetad2 = yd + (T11 + T12 - T13 + grad22 * (Kf * z[4] - ld)) / K2;
  double thetad11 = xd1 + (dT01 + dT02 - dT03) / K1;
  double thetad12 = yd1 + (dT11 + dT12 - dT13) / K2;
  double thetad21 = xd2 + (ddT01 + ddT02 - ddT03) / K1;
  double thetad22 = yd2 + (ddT11 + ddT12 - ddT13) / K2;

  double s11 = (velocity[0] - xd1) + gamma2 * (q[0] - xd);
  double s12 = (velocity[1] - yd1) + gamma2 * (q[1] - yd);
  double s21 = (velocity[2] - thetad11) + gamma2 * (q[2] - thetad1);
  double s22 = (velocity[3] - thetad12) + gamma2 * (q[3] - thetad2);

  double thetar21 = thetad21 - gamma2 * (velocity[2] - thetad11);
  double thetar22 = thetad22 - gamma2 * (velocity[3] - thetad12);

  // control law
  U[0] = 0;//-(T01+T02-T03+grad21*(Kf*1000*z[4]-ld));//0;
  U[1] = 0;//-(T11+T12-T13+grad22*(Kf*1000*z[4]-ld));// 0;
  U[2] = -(J1 * thetar21 + K1 * (thetad1 - xd) - gamma1 * s21);
  U[3] = -(J2 * thetar22 + K2 * (thetad2 - yd) - gamma1 * s22);

  double V1 = (s11 * s11 * m11 + 2 * s11 * s12 * m12 + s12 * s12 * m22) / 2 + (J1 * s21 * s21 + J2 * s22 * s22) / 2;
  z[6] = V1 + gamma1 * gamma2 * ((q[0] - xd) * (q[0] - xd) + (q[1] - yd) * (q[1] - yd) + (q[2] - thetad1) * (q[2] - thetad1) + (q[3] - thetad2) * (q[3] - thetad2)) + (q[0] - xd - q[2] + thetad1) * K1 * (q[0] - xd - q[2] + thetad1) / 2 + (q[1] - yd - q[3] + thetad2) * K2 * (q[1] - yd - q[3] + thetad2) / 2;

  z[11] = P;
  z[14] = qd1;
  z[15] = qd2;
  z[16] = U[2];
  z[17] = U[3];
  z[19] = qd21;
  z[20] = qd22;
  z[21] = qr11;
  z[22] = qr12;
  z[23] = qr21;
  z[24] = qr22;
}

SICONOS_EXPORT void U4(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
{
  double m11 = m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = (m2 * l2 * l2 / 4) + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = (m2 * l2 * l2 / 4) + I2;
  double detm = m11 * m22 - m12 * m12;
  double ddetm = (m12 - m22) * m2 * l1 * l2 * sin(q[1]) * velocity[1];

  double FGyr0 = -m2 * l1 * l2 * sin(q[1]) * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2) + K1 * (q[0] - q[2]);
  double FGyr1 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[0] / 2 + g * m2 * l2 * cos(q[0] + q[1]) / 2 + K2 * (q[1] - q[3]);

  // acceleration in q[0],q[1] coordinates
  double dv0 = -(m22 * FGyr0 - m12 * FGyr1) / detm;
  double dv1 = -(m11 * FGyr1 - m12 * FGyr0) / detm;

  double ddv0 = -((m22 * (-m2 * l1 * l2 * cos(q[1]) * velocity[1] * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) - m2 * l1 * l2 * sin(q[1]) * (dv0 * velocity[1] + dv1 * velocity[0] + dv1 * velocity[1]) - g * (l1 * sin(q[0]) * velocity[0] * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2) + K1 * (velocity[0] - velocity[2])) + m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr1 / 2 - m12 * (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] * velocity[1] / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * dv0 + K2 * (velocity[1] - velocity[3]))) * detm - (m22 * FGyr0 - m12 * FGyr1) * ddetm) / (detm * detm);
  double ddv1 = -((m11 * (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] * velocity[1] / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * dv0 + K2 * (velocity[1] - velocity[3])) + m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr0 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr1 + m12 * (m2 * l1 * l2 * cos(q[1]) * velocity[1] * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) + m2 * l1 * l2 * sin(q[1]) * (dv0 * velocity[1] + dv1 * velocity[0] + dv1 * velocity[1]) + g * m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2 - K1 * (velocity[0] - velocity[2]))) * detm - (m11 * FGyr1 - m12 * FGyr0) * ddetm) / (detm * detm);

  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

  double x1 = -l1 * sin(q[0]) * velocity[0] - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]);
  double y1 = l1 * cos(q[0]) * velocity[1] + l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]);

  double x2 = -l1 * cos(q[0]) * velocity[0] * velocity[0] - l1 * sin(q[0]) * dv0 - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - l2 * sin(q[0] + q[1]) * (dv0 + dv1);
  double y2 = l1 * cos(q[0]) * dv0 - l1 * sin(q[0]) * velocity[0] * velocity[0] + l2 * cos(q[0] + q[1]) * (dv0 + dv1) - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]);

  double x3 = l1 * sin(q[0]) * velocity[0] * velocity[0] * velocity[0] - 2 * l1 * cos(q[0]) * velocity[0] * dv0 - l1 * cos(q[0]) * velocity[0] * dv0 - l1 * sin(q[0]) * ddv0 + l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - 2 * l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - l2 * sin(q[0] + q[1]) * (ddv0 + ddv1);
  double y3 = -l1 * sin(q[0]) * velocity[0] * dv0 + l1 * cos(q[0]) * ddv0 - l1 * cos(q[0]) * velocity[0] * velocity[0] * velocity[0] - 2 * l1 * sin(q[0]) * velocity[0] * dv0 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) + l2 * cos(q[0] + q[1]) * (ddv0 + ddv1) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - 2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1);

  double grad11 = -y;
  double grad12 = -l2 * sin(q[0] + q[1]);
  double grad21 = x;
  double grad22 = l2 * cos(q[0] + q[1]);
  double det = grad11 * grad22 - grad12 * grad21;

  double d11 = grad22 / det;
  double d12 = -grad12 / det;
  double d21 = -grad21 / det;
  double d22 = grad11 / det;

  double dd11 =  -(sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd12 = (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd21 =  -dd11 + (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));
  double dd22 =  -dd12 - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));

  double ddd11 = (-(sin(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] + cos(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd12 = ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + sin(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] - sin(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd21 =  -ddd11 + ((cos(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + sin(q[0]) * sin(q[1]) * dv0 + cos(q[0]) * cos(q[1]) * dv1) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd22 =  -ddd12 - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]));

  // time derivative of ddd

  double dddd11 = ((-(sin(q[0] + q[1]) * (ddv0 + ddv1) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * ddv1 + 3 * cos(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - 3 * cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * dv1 + cos(q[0] + q[1]) * cos(q[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * velocity[1] - velocity[1] * velocity[1] * velocity[1]) - sin(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - velocity[1] * velocity[1]) + sin(q[0] + q[1]) * cos(q[1]) * ((dv0 + dv1) * velocity[1] - (velocity[0] + velocity[1]) * dv1)) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - (-(sin(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] + cos(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd12 = (((cos(q[0] + q[1]) * (ddv0 + ddv1) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * ddv1 - 3 * sin(q[0] + q[1]) * sin(q[1]) * ((velocity[0] + velocity[1]) * (dv0 + dv1) - velocity[1] * dv1) + cos(q[0] + q[1]) * cos(q[1]) * ((dv0 + dv1) * velocity[1] - (velocity[0] + velocity[1]) * dv1) - (cos(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) + sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + sin(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] - sin(q[0] + q[1]) * cos(q[1]) * dv1) * (sin(q[1]) * sin(q[1])) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * sin(q[1]) * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd21 = -dddd11 + (((sin(q[0]) * sin(q[1]) * ddv0 + cos(q[0]) * cos(q[1]) * ddv1 + 3 * (velocity[0] * dv0 - velocity[1] * dv1) + sin(q[0]) * cos(q[1]) * (velocity[1] * dv0 - velocity[0] * dv1) + (cos(q[0]) * cos(q[1]) * velocity[1] - sin(q[0]) * sin(q[1]) * velocity[0]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + sin(q[0]) * sin(q[1]) * dv0 + cos(q[0]) * cos(q[1]) * dv1) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd22 = -dddd12 - (((-(sin(q[0]) * cos(q[1]) * velocity[1] + cos(q[0]) * sin(q[1]) * velocity[0]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + cos(q[0]) * sin(q[1]) * ddv0 - sin(q[0]) * cos(q[1]) * ddv1 + 3 * (velocity[1] * dv1 - velocity[0] * dv0) * sin(q[0]) * sin(q[1]) - cos(q[0]) * cos(q[1]) * (velocity[1] * dv0 - velocity[0] * dv1)) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1])) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));

  double mc11 = m11 * d11 + m12 * d21;
  double mc12 = m11 * d12 + m12 * d22;
  double mc21 = m12 * d11 + m22 * d21;
  double mc22 = m12 * d12 + m22 * d22;

  double a11 = m11 * dd11 + m12 * dd21;
  double a12 = m11 * dd12 + m12 * dd22;
  double a21 = m12 * dd11 + m22 * dd21;
  double a22 = m12 * dd12 + m22 * dd22;

  double da11 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd21 / 2 + m11 * ddd11 + m12 * ddd21;
  double da12 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2) + m11 * ddd12 + m12 * ddd22;
  double da21 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 / 2 + m12 * ddd11 + m22 * ddd21;
  double da22 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 / 2 + m12 * ddd12 + m22 * ddd22;

  double dda11 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd11 + sin(q[1]) * velocity[1] * ddd11 + (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd21 / 2 + sin(q[1]) * velocity[1] * ddd21 / 2 + sin(q[1]) * velocity[1] * (ddd11 + ddd21 / 2)) + m11 * dddd11 + m12 * dddd21;
  double dda12 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (dd12 + dd22 / 2) + sin(q[1]) * velocity[1] * (ddd12 + ddd22 / 2) + sin(q[1]) * velocity[1] * (ddd12 + ddd22 / 2)) + m11 * dddd12 + m12 * dddd22;
  double dda21 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd11 / 2 + sin(q[1]) * velocity[1] * ddd11) + m12 * dddd11 + m22 * dddd21;
  double dda22 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd12 / 2 + sin(q[1]) * velocity[1] * ddd12) + m12 * dddd12 + m22 * dddd22;

  double dmc11 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (d11 + d21 / 2) + m11 * dd11 + m12 * dd21;
  double dmc12 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (d12 + d22 / 2) + m11 * dd12 + m12 * dd22;
  double dmc21 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * d11 / 2 + m12 * dd11 + m22 * dd21;
  double dmc22 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * d12 / 2 + m12 * dd12 + m22 * dd22;

  double ddmc11 = -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (d11 + d21 / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd11 + dd21 / 2) + m11 * ddd11 + m12 * ddd21 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd11 + dd21 / 2);
  double ddmc12 = -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (d12 + d22 / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2) + m11 * ddd12 + m12 * ddd22 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2);
  double ddmc21 =  -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * d11 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 / 2 + m12 * ddd11 + m22 * ddd21 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 / 2;
  double ddmc22 =  -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * d12 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 / 2 + m12 * ddd12 + m22 * ddd22 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 / 2;

  double qd1 = 0.7;
  double qd11 = 0;
  double qd21 = 0;
  double qd31 = 0;
  double qd41 = 0;
  double qd2 = 0;
  double qd12 = 0;
  double qd22 = 0;

  double qr11 = qd11 - gamma2 * (x - qd1);
  double qr12 = -gamma2 * y;

  double qr21 = qd21 - gamma2 * (x1 - qd11);
  double qr22 = -gamma2 * y1;

  double qr31 = qd31 - gamma2 * (x2 - qd21);
  double qr32 = -gamma2 * y2;

  double qr41 = qd41 - gamma2 * (x3 - qd31);
  double qr42 = -gamma2 * y3;

  double s1 = x1 - qr11;
  double s2 = y1 - qr12;
  double s1bar = x1 + gamma2 * (x - 0.7 - (sqrt(z[7]) / 10) * alpha);

  double a1 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2);
  double a2 = a11 * qr11 + a12 * qr12;
  double a3 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * qr11 + d12 * qr12) / 2;
  double a4 = a21 * qr11 + a22 * qr12;

  double da1 = -m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * dv1 * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2);
  double da2 = a11 * qr21 + a12 * qr22 + da11 * qr11 + da12 * qr12;
  double da3 = m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[1] * (d11 * qr11 + d12 * qr12) / 2 + m2 * l1 * l2 * sin(q[1]) * dv0 * (d11 * qr11 + d12 * qr12) / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * (dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) / 2;
  double da4 = a21 * qr21 + a22 * qr22 + da21 * qr11 + da22 * qr12;

  double dda1 = (m2 * l1 * l2 * sin(q[1]) * velocity[1] * velocity[1] * velocity[1] - 2 * m2 * l1 * l2 * cos(q[1]) * velocity[1] * dv1) * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] * ((dd11 * qr11 + d11 * qr21 + dd12 * qr12 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - (m2 * l1 * l2 * cos(q[1]) * velocity[1] * dv1 + m2 * l1 * l2 * sin(q[1]) * ddv1) * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * dv1 * ((dd11 * qr11 + d11 * qr21 + dd12 * qr12 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - (m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] + m2 * l1 * l2 * sin(q[1]) * dv1) * ((dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((ddd11 * qr11 + dd11 * qr21 + ddd12 * qr12 + dd12 * qr22 + dd11 * qr21 + d11 * qr31 + dd12 * qr22 + d12 * qr32) + (ddd21 * qr11 + dd21 * qr21 + ddd22 * qr12 + dd22 * qr22 + dd21 * qr21 + d21 * qr31 + dd22 * qr22 + d22 * qr32) / 2);
  double dda2 = da11 * qr21 + a11 * qr31 + da12 * qr22 + a12 * qr32 + da11 * qr21 + dda11 * qr11 + da12 * qr22 + dda12 * qr12;
  double dda3 = (m2 * l1 * l2 * sin(q[1]) * ddv0 + 2 * m2 * l1 * l2 * cos(q[1]) * dv0 * velocity[1] + m2 * l1 * l2 * cos(q[1]) * velocity[0] * dv1 - m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[1] * velocity[1]) * (d11 * qr11 + d12 * qr12) / 2 + (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[1] + m2 * l1 * l2 * sin(q[1]) * dv0) * (dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + m2 * l1 * l2 * sin(q[1]) * velocity[0] * (ddd11 * qr11 + ddd12 * qr12 + dd11 * qr21 + dd12 * qr22 + dd11 * qr21 + dd12 * qr22 + d11 * qr31 + d12 * qr32) / 2;
  double dda4 = da21 * qr21 + a21 * qr31 + da22 * qr22 + a22 * qr32 + dda21 * qr11 + da21 * qr21 + dda22 * qr12 + da22 * qr22;

  double T01 = mc11 * qr21 + mc12 * qr22;
  double T02 = a1 + a2;
  double T03 = grad11 * gamma1 * s1bar + grad21 * gamma1 * s2;

  double dT01 = dmc11 * qr21 + dmc12 * qr22 + mc11 * qr31 + mc12 * qr32;
  double dT02 = da1 + da2;
  double dT03 = -y1 * gamma1 * s1bar + x1 * gamma1 * s2 + grad11 * gamma1 * (x2 - qr21) + grad21 * gamma1 * (y2 - qr22);

  double ddT01 = ddmc11 * qr21 + ddmc12 * qr22 + 2 * dmc11 * qr31 + 2 * dmc12 * qr32 + mc11 * qr41 + mc12 * qr42;
  double ddT02 = dda1 + dda2;
  double ddT03 = -y2 * gamma1 * s1bar + x2 * gamma1 * s2 - 2 * y1 * gamma1 * (x2 - qr21) + 2 * x1 * gamma1 * (y2 - qr22) + grad11 * gamma1 * (x3 - qr31) + grad21 * gamma1 * (y3 - qr32);

  double T11 = mc21 * qr21 + mc22 * qr22;
  double T12 = a3 + a4;
  double T13 = grad12 * gamma1 * s1bar + grad22 * gamma1 * s2;

  double dT11 = dmc21 * qr21 + dmc22 * qr22 + mc21 * qr31 + mc22 * qr32;
  double dT12 = da3 + da4;
  double dT13 = -l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * s1bar - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * s2 + grad12 * gamma1 * (x2 - qr21) + grad22 * gamma1 * (y2 - qr22);

  double ddT11 = ddmc21 * qr21 + ddmc22 * qr22 + 2 * dmc21 * qr31 + 2 * dmc22 * qr32 + mc21 * qr41 + mc22 * qr42;
  double ddT12 = dda3 + dda4;
  double ddT13 = l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * gamma1 * s1 - l2 * cos(q[0] + q[1]) * (dv0 + dv1) * gamma1 * s1 - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (x2 - qr21) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * gamma1 * s2 - l2 * sin(q[0] + q[1]) * (dv0 + dv1) * gamma1 * s2 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (y2 - qr22) + grad12 * gamma1 * (x3 - qr31) + grad22 * gamma1 * (y3 - qr32) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (x2 - qr21) - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (y2 - qr22);

  double ki1 = sqrt((qd1 * qd1 + qd2 * qd2 + l1 * l1 + l2 * l2) * (qd1 * qd1 + qd2 * qd2 + l1 * l1 + l2 * l2) - 2 * ((qd1 * qd1 + qd2 * qd2) * (qd1 * qd1 + qd2 * qd2) + l1 * l1 * l1 * l1 + l2 * l2 * l2 * l2));
  double ki2 = qd1 * qd1 + qd2 * qd2 + l1 * l1 - l2 * l2;
  double ki3 = qd1 * qd1 + qd2 * qd2 - l1 * l1 - l2 * l2;

  double xd = atan2(qd2, qd1) + atan2(ki1, ki2);
  double yd = -atan2(ki1, ki3);
  double xd1 = (qd1 * qd12 - qd2 * qd11) / (qd1 * qd1 + qd2 * qd2) + ((qd1 * qd11 + qd2 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2)) / (ki1 * (qd1 * qd1 + qd2 * qd2));
  double yd1 = 2 * (qd1 * qd11 + qd2 * qd12) / ki1;
  double xd2 = (((qd1 * qd21 + qd11 * qd11 + qd2 * qd22 + qd12 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2) - 2 * (qd1 * qd11 + qd2 * qd12) * (qd1 * qd11 + qd2 * qd12)) * ki1 * (qd1 * qd1 + qd2 * qd2) - (qd1 * qd11 + qd2 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2) * (ki1 * 2 * (qd1 * qd11 + qd2 * qd12) - (qd1 * qd1 + qd2 * qd2) * 2 * (qd1 * qd11 + qd2 * qd12) * ki3 / ki1)) / (ki1 * (qd1 * qd1 + qd2 * qd2) * ki1 * (qd1 * qd1 + qd2 * qd2)) + ((qd1 * qd22 - qd2 * qd21) * (qd1 * qd1 + qd2 * qd2) - (qd1 * qd12 - qd2 * qd11) * 2 * (qd1 * qd11 + qd2 * qd12)) / ((qd1 * qd1 + qd2 * qd2) * (qd1 * qd1 + qd2 * qd2));
  double yd2 = 2 * ((qd1 * qd21 + qd11 * qd11 + qd2 * qd22 + qd12 * qd12) * ki1 * ki1 + (qd1 * qd11 + qd2 * qd12) * 2 * (qd1 * qd11 + qd2 * qd12) * ki3) / (ki1 * ki1 * ki1);

  double thetad1 = xd + (T01 + T02 - T03 + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2)) / K1;
  double thetad2 = yd + (T11 + T12 - T13 + g * m2 * l2 * cos(q[0] + q[1]) / 2) / K2;
  double thetad11 = xd1 + (dT01 + dT02 - dT03 - g * (l1 * sin(q[0]) * velocity[0] * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2)) / K1;
  double thetad12 = yd1 + (dT11 + dT12 - dT13 - g * m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2) / K2;
  double thetad21 = xd2 + (ddT01 + ddT02 - ddT03 - g * (l1 * sin(q[0]) * dv0 * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (dv0 + dv1) / 2 + l1 * cos(q[0]) * velocity[0] * velocity[0] * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) / 2)) / K1;
  double thetad22 = yd2 + (ddT11 + ddT12 - ddT13 - g * m2 * l2 * (sin(q[0] + q[1]) * (dv0 + dv1) / 2 + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) / 2)) / K2;

  double s21 = (velocity[2] - thetad11) + gamma2 * (q[2] - thetad1);
  double s22 = (velocity[3] - thetad12) + gamma2 * (q[3] - thetad2);

  double thetar21 = thetad21 - gamma2 * (velocity[2] - thetad11);
  double thetar22 = thetad22 - gamma2 * (velocity[3] - thetad12);

  // control law
  U[0] = 0;//-(T01+T02-T03);// 0;
  U[1] = 0;//-(T11+T12-T13);//0;
  U[2] = -(J1 * thetar21 + K1 * (thetad1 - xd) - gamma1 * s21);
  U[3] = -(J2 * thetar22 + K2 * (thetad2 - yd) - gamma1 * s22);

  double V1 = (s1 * s1 * (d11 * mc11 + d21 * mc21) + 2 * s1 * s2 * (d12 * mc11 + d22 * mc21) + s2 * s2 * (d12 * mc12 + d22 * mc22)) / 2 + (J1 * s21 * s21 + J2 * s22 * s22) / 2;
  z[6] = V1 + gamma1 * gamma2 * ((q[0] - xd) * (q[0] - xd) + (q[1] - yd) * (q[1] - yd) + (q[2] - thetad1) * (q[2] - thetad1) + (q[3] - thetad2) * (q[3] - thetad2)) + (q[0] - xd - q[2] + thetad1) * K1 * (q[0] - xd - q[2] + thetad1) / 2 + (q[1] - yd - q[3] + thetad2) * K2 * (q[1] - yd - q[3] + thetad2) / 2;

  z[11] = P;
  z[14] = qd1;
  z[15] = qd2;
  z[16] = U[2];
  z[17] = U[3];
  z[19] = s1;
  z[20] = s2;
  z[21] = qr11;
  z[22] = qr12;
  z[23] = qr21;
  z[24] = qr22;
}

SICONOS_EXPORT void U5(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
{
  double m11 = m1 * l1 * l1 / 4 + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = m2 * l2 * l2 / 4 + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = m2 * l2 * l2 / 4 + I2;
  double detm = m11 * m22 - m12 * m12;
  double ddetm = (m12 - m22) * m2 * l1 * l2 * sin(q[1]) * velocity[1];

  double FGyr0 = -m2 * l1 * l2 * sin(q[1]) * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2) + K1 * (q[0] - q[2]);
  double FGyr1 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[0] / 2 + g * m2 * l2 * cos(q[0] + q[1]) / 2 + K2 * (q[1] - q[3]);

  // acceleration in q[0],q[1] coordinates
  double dv0 = -(m22 * FGyr0 - m12 * FGyr1) / detm;
  double dv1 = -(m11 * FGyr1 - m12 * FGyr0) / detm;

  double ddv0 = -((m22 * (-m2 * l1 * l2 * cos(q[1]) * velocity[1] * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) - m2 * l1 * l2 * sin(q[1]) * (dv0 * velocity[1] + dv1 * velocity[0] + dv1 * velocity[1]) - g * (l1 * sin(q[0]) * velocity[0] * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2) + K1 * (velocity[0] - velocity[2])) + m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr1 / 2 - m12 * (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] * velocity[1] / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * dv0 + K2 * (velocity[1] - velocity[3]))) * detm - (m22 * FGyr0 - m12 * FGyr1) * ddetm) / (detm * detm);
  double ddv1 = -((m11 * (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] * velocity[1] / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * dv0 + K2 * (velocity[1] - velocity[3])) + m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr0 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr1 + m12 * (m2 * l1 * l2 * cos(q[1]) * velocity[1] * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) + m2 * l1 * l2 * sin(q[1]) * (dv0 * velocity[1] + dv1 * velocity[0] + dv1 * velocity[1]) + g * m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2 - K1 * (velocity[0] - velocity[2]))) * detm - (m11 * FGyr1 - m12 * FGyr0) * ddetm) / (detm * detm);

  // generalized coordinates q=(x,y)^T
  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

  // time derivatives of x and y
  double x1 = -l1 * sin(q[0]) * velocity[0] - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]);
  double y1 = l1 * cos(q[0]) * velocity[0] + l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]);

  // acceleration in x,y coordinates
  double x2 =  -l1 * cos(q[0]) * velocity[0] * velocity[0] - l1 * sin(q[0]) * dv0 - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - l2 * sin(q[0] + q[1]) * (dv0 + dv1);
  double y2 =  l1 * cos(q[0]) * dv0 - l1 * sin(q[0]) * velocity[0] * velocity[0] + l2 * cos(q[0] + q[1]) * (dv0 + dv1) - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]);

  double x3 =  l1 * sin(q[0]) * velocity[0] * velocity[0] * velocity[0] - 2 * l1 * cos(q[0]) * velocity[0] * dv0 - l1 * cos(q[0]) * velocity[0] * dv0 - l1 * sin(q[0]) * ddv0 + l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - 2 * l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - l2 * sin(q[0] + q[1]) * (ddv0 + ddv1);
  double y3 = -l1 * sin(q[0]) * velocity[0] * dv0 + l1 * cos(q[0]) * ddv0 - l1 * cos(q[0]) * velocity[0] * velocity[0] * velocity[0] - 2 * l1 * sin(q[0]) * velocity[0] * dv0 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) + l2 * cos(q[0] + q[1]) * (ddv0 + ddv1) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - 2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1);

  // the gradient of the transformation \theta->q
  double grad11 = -y;
  double grad12 = -l2 * sin(q[0] + q[1]);
  double grad21 = x;
  double grad22 = l2 * cos(q[0] + q[1]);
  double det = grad11 * grad22 - grad12 * grad21;

  // d=(grad)^{-1}
  double d11 = grad22 / det;
  double d12 = -grad12 / det;
  double d21 = -grad21 / det;
  double d22 = grad11 / det;

  // time derivative of d

  double dd11 =  -(sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd12 = (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd21 =  -dd11 + (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));
  double dd22 =  -dd12 - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));

  // time derivative of dd

  double ddd11 = (-(sin(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] + cos(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd12 = ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + sin(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] - sin(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd21 =  -ddd11 + ((cos(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + sin(q[0]) * sin(q[1]) * dv0 + cos(q[0]) * cos(q[1]) * dv1) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd22 =  -ddd12 - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]));

  // time derivative of ddd

  double dddd11 = ((-(sin(q[0] + q[1]) * (ddv0 + ddv1) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * ddv1 + 3 * cos(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - 3 * cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * dv1 + cos(q[0] + q[1]) * cos(q[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * velocity[1] - velocity[1] * velocity[1] * velocity[1]) - sin(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - velocity[1] * velocity[1]) + sin(q[0] + q[1]) * cos(q[1]) * ((dv0 + dv1) * velocity[1] - (velocity[0] + velocity[1]) * dv1)) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - (-(sin(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] + cos(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd12 = (((cos(q[0] + q[1]) * (ddv0 + ddv1) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * ddv1 - 3 * sin(q[0] + q[1]) * sin(q[1]) * ((velocity[0] + velocity[1]) * (dv0 + dv1) - velocity[1] * dv1) + cos(q[0] + q[1]) * cos(q[1]) * ((dv0 + dv1) * velocity[1] - (velocity[0] + velocity[1]) * dv1) - (cos(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) + sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + sin(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] - sin(q[0] + q[1]) * cos(q[1]) * dv1) * (sin(q[1]) * sin(q[1])) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * sin(q[1]) * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd21 = -dddd11 + (((sin(q[0]) * sin(q[1]) * ddv0 + cos(q[0]) * cos(q[1]) * ddv1 + 3 * (velocity[0] * dv0 - velocity[1] * dv1) + sin(q[0]) * cos(q[1]) * (velocity[1] * dv0 - velocity[0] * dv1) + (cos(q[0]) * cos(q[1]) * velocity[1] - sin(q[0]) * sin(q[1]) * velocity[0]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + sin(q[0]) * sin(q[1]) * dv0 + cos(q[0]) * cos(q[1]) * dv1) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd22 = -dddd12 - (((-(sin(q[0]) * cos(q[1]) * velocity[1] + cos(q[0]) * sin(q[1]) * velocity[0]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + cos(q[0]) * sin(q[1]) * ddv0 - sin(q[0]) * cos(q[1]) * ddv1 + 3 * (velocity[1] * dv1 - velocity[0] * dv0) * sin(q[0]) * sin(q[1]) - cos(q[0]) * cos(q[1]) * (velocity[1] * dv0 - velocity[0] * dv1)) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1])) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));

  // M*d
  double mc11 = m11 * d11 + m12 * d21;
  double mc12 = m11 * d12 + m12 * d22;
  double mc21 = m12 * d11 + m22 * d21;
  double mc22 = m12 * d12 + m22 * d22;

  // M*dd
  double a11 = m11 * dd11 + m12 * dd21;
  double a12 = m11 * dd12 + m12 * dd22;
  double a21 = m12 * dd11 + m22 * dd21;
  double a22 = m12 * dd12 + m22 * dd22;

  double da11 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd11 + dd21 / 2) + m11 * ddd11 + m12 * ddd21;
  double da12 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2) + m11 * ddd12 + m12 * ddd22;
  double da21 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 / 2 + m12 * ddd11 + m22 * ddd21;
  double da22 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 / 2 + m12 * ddd12 + m22 * ddd22;

  double dda11 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (dd11 + dd21 / 2) + 2 * sin(q[1]) * velocity[1] * (ddd11 + ddd21 / 2)) + m11 * dddd11 + m12 * dddd21;
  double dda12 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (dd12 + dd22 / 2) + 2 * sin(q[1]) * velocity[1] * (ddd12 + ddd22 / 2)) + m11 * dddd12 + m12 * dddd22;
  double dda21 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd11 / 2 + sin(q[1]) * velocity[1] * ddd11) + m12 * dddd11 + m22 * dddd21;
  double dda22 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd12 / 2 + sin(q[1]) * velocity[1] * ddd12) + m12 * dddd12 + m22 * dddd22;

  double dmc11 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (d11 + d21 / 2) + m11 * dd11 + m12 * dd21;
  double dmc12 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (d12 + d22 / 2) + m11 * dd12 + m12 * dd22;
  double dmc21 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * d11 / 2 + m12 * dd11 + m22 * dd21;
  double dmc22 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * d12 / 2 + m12 * dd12 + m22 * dd22;

  double ddmc11 = -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (d11 + d21 / 2) - 2 * m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd11 + dd21 / 2) + m11 * ddd11 + m12 * ddd21;
  double ddmc12 = -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (d12 + d22 / 2) - 2 * m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2) + m11 * ddd12 + m12 * ddd22;
  double ddmc21 =  -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * d11 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 + m12 * ddd11 + m22 * ddd21;
  double ddmc22 =  -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * d12 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 + m12 * ddd12 + m22 * ddd22;

  double k = trunc(time / P) + 1;
  double td = k * P - z[8];

  double c0 = 0.5 / (td * td * td * td);
  double c1 = -2 / (td * td * td * td * td);
  double c2 = 5 / (td * td * td * td * td * td);
  double c3 = -10 / (td * td * td * td * td * td * td);

  double qd1 = 0.7;
  double qd11 = 0;
  double qd21 = 0;
  double qd31 = 0;
  double qd41 = 0;

  double qd2 = c0 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) + c1 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - k * P) + c2 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - k * P) * (time - k * P) + c3 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - k * P) * (time - k * P) * (time - k * P);
  double qd12 = 4 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (c0 + c1 * (time - k * P) + c2 * (time - k * P) * (time - k * P) + c3 * (time - k * P) * (time - k * P) * (time - k * P)) + (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (c1 + 2 * c2 * (time - k * P) + 3 * c3 * (time - k * P) * (time - k * P));
  double qd22 = 12 * (time - z[8]) * (time - z[8]) * (c0 + c1 * (time - k * P) + c2 * (time - k * P) * (time - k * P) + c3 * (time - k * P) * (time - k * P) * (time - k * P)) + 8 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (c1 + 2 * c2 * (time - k * P) + 3 * c3 * (time - k * P) * (time - k * P)) + (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * (2 * c2 + 6 * c3 * (time - k * P));
  double qd32 = 24 * (time - z[8]) * (c0 + c1 * (time - k * P) + c2 * (time - k * P) * (time - k * P) + c3 * (time - k * P) * (time - k * P) * (time - k * P)) + 36 * (time - z[8]) * (time - z[8]) * (c1 + 2 * c2 * (time - k * P) + 3 * c3 * (time - k * P) * (time - k * P)) + 12 * (time - z[8]) * (time - z[8]) * (time - z[8]) * (2 * c2 + 6 * c3 * (time - k * P)) + (time - z[8]) * (time - z[8]) * (time - z[8]) * (time - z[8]) * 6 * c3;
  double qd42 = 0;

  double qr11 = qd11 - gamma2 * (x - qd1);
  double qr12 = qd12 - gamma2 * (y - qd2);

  double qr21 = qd21 - gamma2 * (x1 - qd11);
  double qr22 = qd22 - gamma2 * (y1 - qd12);

  double qr31 = qd31 - gamma2 * (x2 - qd21);
  double qr32 = qd32 - gamma2 * (y2 - qd22);

  double qr41 = qd41 - gamma2 * (x3 - qd31);
  double qr42 = qd42 - gamma2 * (y3 - qd32);

  double s1 = x1 - qr11;
  double s2 = y1 - qr12;

  double ki1 = sqrt((qd1 * qd1 + qd2 * qd2 + l1 * l1 + l2 * l2) * (qd1 * qd1 + qd2 * qd2 + l1 * l1 + l2 * l2) - 2 * ((qd1 * qd1 + qd2 * qd2) * (qd1 * qd1 + qd2 * qd2) + l1 * l1 * l1 * l1 + l2 * l2 * l2 * l2));
  double ki2 = qd1 * qd1 + qd2 * qd2 + l1 * l1 - l2 * l2;
  double ki3 = qd1 * qd1 + qd2 * qd2 - l1 * l1 - l2 * l2;

  double xd = atan2(qd2, qd1) + atan2(ki1, ki2);
  double yd = -atan2(ki1, ki3);
  double xd1 = (qd1 * qd12 - qd2 * qd11) / (qd1 * qd1 + qd2 * qd2) + ((qd1 * qd11 + qd2 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2)) / (ki1 * (qd1 * qd1 + qd2 * qd2));
  double yd1 = 2 * (qd1 * qd11 + qd2 * qd12) / ki1;
  double xd2 = (((qd1 * qd21 + qd11 * qd11 + qd2 * qd22 + qd12 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2) - 2 * (qd1 * qd11 + qd2 * qd12) * (qd1 * qd11 + qd2 * qd12)) * ki1 * (qd1 * qd1 + qd2 * qd2) - (qd1 * qd11 + qd2 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2) * (ki1 * 2 * (qd1 * qd11 + qd2 * qd12) - (qd1 * qd1 + qd2 * qd2) * 2 * (qd1 * qd11 + qd2 * qd12) * ki3 / ki1)) / (ki1 * (qd1 * qd1 + qd2 * qd2) * ki1 * (qd1 * qd1 + qd2 * qd2)) + ((qd1 * qd22 - qd2 * qd21) * (qd1 * qd1 + qd2 * qd2) - (qd1 * qd12 - qd2 * qd11) * 2 * (qd1 * qd11 + qd2 * qd12)) / ((qd1 * qd1 + qd2 * qd2) * (qd1 * qd1 + qd2 * qd2));
  double yd2 = 2 * ((qd1 * qd21 + qd11 * qd11 + qd2 * qd22 + qd12 * qd12) * ki1 * ki1 + (qd1 * qd11 + qd2 * qd12) * 2 * (qd1 * qd11 + qd2 * qd12) * ki3) / (ki1 * ki1 * ki1);

  double a1 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2);
  double a2 = a11 * qr11 + a12 * qr12;
  double a3 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * qr11 + d12 * qr12) / 2;
  double a4 = a21 * qr11 + a22 * qr12;

  double da1 = -m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * dv1 * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2);
  double da2 = a11 * qr21 + a12 * qr22 + da11 * qr11 + da12 * qr12;
  double da3 = m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[1] * (d11 * qr11 + d12 * qr12) / 2 + m2 * l1 * l2 * sin(q[1]) * dv0 * (d11 * qr11 + d12 * qr12) / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * (dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) / 2;
  double da4 = a21 * qr21 + a22 * qr22 + da21 * qr11 + da22 * qr12;

  double dda1 = (m2 * l1 * l2 * sin(q[1]) * velocity[1] * velocity[1] * velocity[1] - 2 * m2 * l1 * l2 * cos(q[1]) * velocity[1] * dv1) * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] * ((dd11 * qr11 + d11 * qr21 + dd12 * qr12 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - (m2 * l1 * l2 * cos(q[1]) * velocity[1] * dv1 + m2 * l1 * l2 * sin(q[1]) * ddv1) * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * dv1 * ((dd11 * qr11 + d11 * qr21 + dd12 * qr12 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - (m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] + m2 * l1 * l2 * sin(q[1]) * dv1) * ((dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((ddd11 * qr11 + dd11 * qr21 + ddd12 * qr12 + dd12 * qr22 + dd11 * qr21 + d11 * qr31 + dd12 * qr22 + d12 * qr32) + (ddd21 * qr11 + dd21 * qr21 + ddd22 * qr12 + dd22 * qr22 + dd21 * qr21 + d21 * qr31 + dd22 * qr22 + d22 * qr32) / 2);
  double dda2 = da11 * qr21 + a11 * qr31 + da12 * qr22 + a12 * qr32 + da11 * qr21 + dda11 * qr11 + da12 * qr22 + dda12 * qr12;
  double dda3 = (m2 * l1 * l2 * sin(q[1]) * ddv0 + 2 * m2 * l1 * l2 * cos(q[1]) * dv0 * velocity[1] + m2 * l1 * l2 * cos(q[1]) * velocity[0] * dv1 - m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[1] * velocity[1]) * (d11 * qr11 + d12 * qr12) / 2 + (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[1] + m2 * l1 * l2 * sin(q[1]) * dv0) * (dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + m2 * l1 * l2 * sin(q[1]) * velocity[0] * (ddd11 * qr11 + ddd12 * qr12 + dd11 * qr21 + dd12 * qr22 + dd11 * qr21 + dd12 * qr22 + d11 * qr31 + d12 * qr32) / 2;
  double dda4 = da21 * qr21 + a21 * qr31 + da22 * qr22 + a22 * qr32 + dda21 * qr11 + da21 * qr21 + dda22 * qr12 + da22 * qr22;

  double T01 = mc11 * qr21 + mc12 * qr22;
  double T02 = a1 + a2 + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2);
  double T03 = grad11 * gamma1 * s1 + grad21 * gamma1 * s2;

  double dT01 = dmc11 * qr21 + dmc12 * qr22 + mc11 * qr31 + mc12 * qr32;
  double dT02 = da1 + da2 - g * (l1 * sin(q[0]) * velocity[0] * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2);
  double dT03 = -y1 * gamma1 * s1 + x1 * gamma1 * s2 + grad11 * gamma1 * (x2 - qr21) + grad21 * gamma1 * (y2 - qr22);

  double ddT01 = ddmc11 * qr21 + ddmc12 * qr22 + 2 * dmc11 * qr31 + 2 * dmc12 * qr32 + mc11 * qr41 + mc12 * qr42;
  double ddT02 = dda1 + dda2 - g * (l1 * sin(q[0]) * dv0 * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (dv0 + dv1) / 2 + l1 * cos(q[0]) * velocity[0] * velocity[0] * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) / 2);
  double ddT03 = -y2 * gamma1 * s1 + x2 * gamma1 * s2 - 2 * y1 * gamma1 * (x2 - qr21) + 2 * x1 * gamma1 * (y2 - qr22) + grad11 * gamma1 * (x3 - qr31) + grad21 * gamma1 * (y3 - qr32);

  double T11 = mc21 * qr21 + mc22 * qr22;
  double T12 = a3 + a4 + g * m2 * l2 * cos(q[0] + q[1]) / 2;
  double T13 = grad12 * gamma1 * s1 + grad22 * gamma1 * s2;

  double dT11 = dmc21 * qr21 + dmc22 * qr22 + mc21 * qr31 + mc22 * qr32;
  double dT12 = da3 + da4 - g * m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2;
  double dT13 = -l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * s1 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * s2 + grad12 * gamma1 * (x2 - qr21) + grad22 * gamma1 * (y2 - qr22);

  double ddT11 = ddmc21 * qr21 + ddmc22 * qr22 + 2 * dmc21 * qr31 + 2 * dmc22 * qr32 + mc21 * qr41 + mc22 * qr42;
  double ddT12 = dda3 + dda4 - g * m2 * l2 * (sin(q[0] + q[1]) * (dv0 + dv1) / 2 + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) / 2);
  double ddT13 = l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * gamma1 * s1 - l2 * cos(q[0] + q[1]) * (dv0 + dv1) * gamma1 * s1 - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (x2 - qr21) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * gamma1 * s2 - l2 * sin(q[0] + q[1]) * (dv0 + dv1) * gamma1 * s2 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (y2 - qr22) + grad12 * gamma1 * (x3 - qr31) + grad22 * gamma1 * (y3 - qr32) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (x2 - qr21) - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (y2 - qr22);

  double ld = (m2 * l1 * l2 * sin(q[1]) * velocity[1] * d11 / 2 - a12 - (mc12 / mc22) * (m2 * l1 * l2 * sin(q[1]) * velocity[1] * d21 / 2 - a22)) * s2 - gamma1 * mc12 * s2 / mc22 - (k * P - time) * (1 + Kf);

  z[12] = (-mc11 * ld + fabs(mc11 * ld)) / 2;

  double thetad1 = xd + (T01 + T02 - T03 + grad11 * (Kf * z[18] - ld)) / K1;
  double thetad2 = yd + (T11 + T12 - T13 + grad12 * (Kf * z[18] - ld)) / K2;
  double thetad11 = xd1 + (dT01 + dT02 - dT03) / K1;
  double thetad12 = yd1 + (dT11 + dT12 - dT13) / K2;
  double thetad21 = xd2 + (ddT01 + ddT02 - ddT03) / K1;
  double thetad22 = yd2 + (ddT11 + ddT12 - ddT13) / K2;

  double s11 = (velocity[0] - xd1) + gamma2 * (q[0] - xd);
  double s12 = (velocity[1] - yd1) + gamma2 * (q[1] - yd);
  double s21 = (velocity[2] - thetad11) + gamma2 * (q[2] - thetad1);
  double s22 = (velocity[3] - thetad12) + gamma2 * (q[3] - thetad2);

  double thetar21 = thetad21 - gamma2 * (velocity[2] - thetad11);
  double thetar22 = thetad22 - gamma2 * (velocity[3] - thetad12);

  // control law
  U[0] = 0;//-(T01+T02-T03+grad11*(Kf*1000*z[18]-ld));//0;
  U[1] = 0;//-(T11+T12-T13+grad12*(Kf*1000*z[18]-ld));//0;
  U[2] = -(J1 * thetar21 + K1 * (thetad1 - xd) - gamma1 * s21);
  U[3] = -(J2 * thetar22 + K2 * (thetad2 - yd) - gamma1 * s22);

  double V1 = (s11 * s11 * m11 + 2 * s11 * s12 * m12 + s12 * s12 * m22) / 2 + (J1 * s21 * s21 + J2 * s22 * s22) / 2;
  z[6] = V1 + gamma1 * gamma2 * ((q[0] - xd) * (q[0] - xd) + (q[1] - yd) * (q[1] - yd) + (q[2] - thetad1) * (q[2] - thetad1) + (q[3] - thetad2) * (q[3] - thetad2)) + (q[0] - xd - q[2] + thetad1) * K1 * (q[0] - xd - q[2] + thetad1) / 2 + (q[1] - yd - q[3] + thetad2) * K2 * (q[1] - yd - q[3] + thetad2) / 2;

  z[11] = P;
  z[14] = qd1;
  z[15] = qd2;
  z[16] = U[2];
  z[17] = U[3];
  z[7] = 0;
  z[19] = s21;
  z[20] = s22;
  z[21] = qr11;
  z[22] = qr12;
  z[23] = qr21;
  z[24] = qr22;
}


SICONOS_EXPORT void U6(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
{
  double m11 = m1 * l1 * l1 / 4 + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = m2 * l2 * l2 / 4 + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = m2 * l2 * l2 / 4 + I2;
  double detm = m11 * m22 - m12 * m12;
  double ddetm = (m12 - m22) * m2 * l1 * l2 * sin(q[1]) * velocity[1];

  double FGyr0 = -m2 * l1 * l2 * sin(q[1]) * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2) + K1 * (q[0] - q[2]);
  double FGyr1 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[0] / 2 + g * m2 * l2 * cos(q[0] + q[1]) / 2 + K2 * (q[1] - q[3]);

  // acceleration in q[0],q[1] coordinates
  double dv0 = -(m22 * FGyr0 - m12 * FGyr1) / detm;
  double dv1 = -(m11 * FGyr1 - m12 * FGyr0) / detm;

  double ddv0 = -((m22 * (-m2 * l1 * l2 * cos(q[1]) * velocity[1] * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) - m2 * l1 * l2 * sin(q[1]) * (dv0 * velocity[1] + dv1 * velocity[0] + dv1 * velocity[1]) - g * (l1 * sin(q[0]) * velocity[0] * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2) + K1 * (velocity[0] - velocity[2])) + m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr1 / 2 - m12 * (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] * velocity[1] / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * dv0 + K2 * (velocity[1] - velocity[3]))) * detm - (m22 * FGyr0 - m12 * FGyr1) * ddetm) / (detm * detm);
  double ddv1 = -((m11 * (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] * velocity[1] / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * dv0 + K2 * (velocity[1] - velocity[3])) + m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr0 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * FGyr1 + m12 * (m2 * l1 * l2 * cos(q[1]) * velocity[1] * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2) + m2 * l1 * l2 * sin(q[1]) * (dv0 * velocity[1] + dv1 * velocity[0] + dv1 * velocity[1]) + g * m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2 - K1 * (velocity[0] - velocity[2]))) * detm - (m11 * FGyr1 - m12 * FGyr0) * ddetm) / (detm * detm);

  // generalized coordinates q=(x,y)^T
  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

  // time derivatives of x and y
  double x1 = -l1 * sin(q[0]) * velocity[0] - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]);
  double y1 = l1 * cos(q[0]) * velocity[0] + l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]);

  // acceleration in x,y coordinates
  double x2 =  -l1 * cos(q[0]) * velocity[0] * velocity[0] - l1 * sin(q[0]) * dv0 - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - l2 * sin(q[0] + q[1]) * (dv0 + dv1);
  double y2 =  l1 * cos(q[0]) * dv0 - l1 * sin(q[0]) * velocity[0] * velocity[0] + l2 * cos(q[0] + q[1]) * (dv0 + dv1) - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]);

  double x3 =  l1 * sin(q[0]) * velocity[0] * velocity[0] * velocity[0] - 2 * l1 * cos(q[0]) * velocity[0] * dv0 - l1 * cos(q[0]) * velocity[0] * dv0 - l1 * sin(q[0]) * ddv0 + l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - 2 * l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - l2 * sin(q[0] + q[1]) * (ddv0 + ddv1);
  double y3 = -l1 * sin(q[0]) * velocity[0] * dv0 + l1 * cos(q[0]) * ddv0 - l1 * cos(q[0]) * velocity[0] * velocity[0] * velocity[0] - 2 * l1 * sin(q[0]) * velocity[0] * dv0 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) + l2 * cos(q[0] + q[1]) * (ddv0 + ddv1) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - 2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1);

  // the gradient of the transformation \theta->q
  double grad11 = -y;
  double grad12 = -l2 * sin(q[0] + q[1]);
  double grad21 = x;
  double grad22 = l2 * cos(q[0] + q[1]);
  double det = grad11 * grad22 - grad12 * grad21;

  // d=(grad)^{-1}
  double d11 = grad22 / det;
  double d12 = -grad12 / det;
  double d21 = -grad21 / det;
  double d22 = grad11 / det;

  // time derivative of d

  double dd11 =  -(sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd12 = (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd21 =  -dd11 + (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));
  double dd22 =  -dd12 - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));

  // time derivative of dd

  double ddd11 = (-(sin(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] + cos(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd12 = ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + sin(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] - sin(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd21 =  -ddd11 + ((cos(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + sin(q[0]) * sin(q[1]) * dv0 + cos(q[0]) * cos(q[1]) * dv1) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double ddd22 =  -ddd12 - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]));

  // time derivative of ddd

  double dddd11 = ((-(sin(q[0] + q[1]) * (ddv0 + ddv1) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * ddv1 + 3 * cos(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) * (dv0 + dv1) - 3 * cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * dv1 + cos(q[0] + q[1]) * cos(q[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * velocity[1] - velocity[1] * velocity[1] * velocity[1]) - sin(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - velocity[1] * velocity[1]) + sin(q[0] + q[1]) * cos(q[1]) * ((dv0 + dv1) * velocity[1] - (velocity[0] + velocity[1]) * dv1)) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - (-(sin(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] + cos(q[0] + q[1]) * cos(q[1]) * dv1) * sin(q[1]) + (sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd12 = (((cos(q[0] + q[1]) * (ddv0 + ddv1) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * ddv1 - 3 * sin(q[0] + q[1]) * sin(q[1]) * ((velocity[0] + velocity[1]) * (dv0 + dv1) - velocity[1] * dv1) + cos(q[0] + q[1]) * cos(q[1]) * ((dv0 + dv1) * velocity[1] - (velocity[0] + velocity[1]) * dv1) - (cos(q[0] + q[1]) * sin(q[1]) * (velocity[0] + velocity[1]) + sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * ((velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * (dv0 + dv1) * sin(q[1]) + sin(q[0] + q[1]) * sin(q[1]) * velocity[1] * velocity[1] - sin(q[0] + q[1]) * cos(q[1]) * dv1) * (sin(q[1]) * sin(q[1])) - (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) * 2 * sin(q[1]) * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd21 = -dddd11 + (((sin(q[0]) * sin(q[1]) * ddv0 + cos(q[0]) * cos(q[1]) * ddv1 + 3 * (velocity[0] * dv0 - velocity[1] * dv1) + sin(q[0]) * cos(q[1]) * (velocity[1] * dv0 - velocity[0] * dv1) + (cos(q[0]) * cos(q[1]) * velocity[1] - sin(q[0]) * sin(q[1]) * velocity[0]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + sin(q[0]) * sin(q[1]) * dv0 + cos(q[0]) * cos(q[1]) * dv1) * sin(q[1]) - (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));
  double dddd22 = -dddd12 - (((-(sin(q[0]) * cos(q[1]) * velocity[1] + cos(q[0]) * sin(q[1]) * velocity[0]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1]) + cos(q[0]) * sin(q[1]) * ddv0 - sin(q[0]) * cos(q[1]) * ddv1 + 3 * (velocity[1] * dv1 - velocity[0] * dv0) * sin(q[0]) * sin(q[1]) - cos(q[0]) * cos(q[1]) * (velocity[1] * dv0 - velocity[0] * dv1)) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) * 2 * (cos(q[1]) * dv1 - sin(q[1]) * velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - ((cos(q[0]) * sin(q[1]) * dv0 - sin(q[0]) * cos(q[1]) * dv1 - sin(q[0]) * sin(q[1]) * (velocity[0] * velocity[0] - velocity[1] * velocity[1])) * sin(q[1]) - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1])) * 2 * cos(q[1]) * velocity[1]) * 3 * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]) * sin(q[1]));

  // M*d
  double mc11 = m11 * d11 + m12 * d21;
  double mc12 = m11 * d12 + m12 * d22;
  double mc21 = m12 * d11 + m22 * d21;
  double mc22 = m12 * d12 + m22 * d22;

  // M*dd
  double a11 = m11 * dd11 + m12 * dd21;
  double a12 = m11 * dd12 + m12 * dd22;
  double a21 = m12 * dd11 + m22 * dd21;
  double a22 = m12 * dd12 + m22 * dd22;

  double da11 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd11 + dd21 / 2) + m11 * ddd11 + m12 * ddd21;
  double da12 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2) + m11 * ddd12 + m12 * ddd22;
  double da21 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 / 2 + m12 * ddd11 + m22 * ddd21;
  double da22 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 / 2 + m12 * ddd12 + m22 * ddd22;

  double dda11 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (dd11 + dd21 / 2) + 2 * sin(q[1]) * velocity[1] * (ddd11 + ddd21 / 2)) + m11 * dddd11 + m12 * dddd21;
  double dda12 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (dd12 + dd22 / 2) + 2 * sin(q[1]) * velocity[1] * (ddd12 + ddd22 / 2)) + m11 * dddd12 + m12 * dddd22;
  double dda21 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd11 / 2 + sin(q[1]) * velocity[1] * ddd11) + m12 * dddd11 + m22 * dddd21;
  double dda22 = -m2 * l1 * l2 * ((cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * dd12 / 2 + sin(q[1]) * velocity[1] * ddd12) + m12 * dddd12 + m22 * dddd22;

  double dmc11 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (d11 + d21 / 2) + m11 * dd11 + m12 * dd21;
  double dmc12 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * (d12 + d22 / 2) + m11 * dd12 + m12 * dd22;
  double dmc21 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * d11 / 2 + m12 * dd11 + m22 * dd21;
  double dmc22 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * d12 / 2 + m12 * dd12 + m22 * dd22;

  double ddmc11 = -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (d11 + d21 / 2) - 2 * m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd11 + dd21 / 2) + m11 * ddd11 + m12 * ddd21;
  double ddmc12 = -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * (d12 + d22 / 2) - 2 * m2 * l1 * l2 * sin(q[1]) * velocity[1] * (dd12 + dd22 / 2) + m11 * ddd12 + m12 * ddd22;
  double ddmc21 =  -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * d11 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd11 + m12 * ddd11 + m22 * ddd21;
  double ddmc22 =  -m2 * l1 * l2 * (cos(q[1]) * velocity[1] * velocity[1] + sin(q[1]) * dv1) * d12 / 2 - m2 * l1 * l2 * sin(q[1]) * velocity[1] * dd12 + m12 * ddd12 + m22 * ddd22;
  double k = trunc(z[8] / P);
  double t0 = k * P;
  double t1 = k * P + ep;
  //double b = -z[10]*ep*P/(2*PI);

  double b0 = 0.5 * cos(PI / 2 + 2 * PI * (k * P + ep) / P) / (ep * ep * ep * ep);
  double b1 = (-0.5 * (2 * PI / P) * sin(PI / 2 + 2 * PI * (k * P + ep) / P) * ep - 4 * 0.5 * cos(PI / 2 + 2 * PI * (k * P + ep) / P)) / (ep * ep * ep * ep * ep);
  double b2 = (-0.5 * ep * ep * (2 * PI / P) * (2 * PI / P) * cos(PI / 2 + 2 * PI * (k * P + ep) / P) + 8 * 0.5 * (2 * PI / P) * sin(PI / 2 + 2 * PI * (k * P + ep) / P) * ep + 20 * 0.5 * cos(PI / 2 + 2 * PI * (k * P + ep) / P)) / (2 * ep * ep * ep * ep * ep * ep);
  double b3 = (0.5 * ep * ep * ep * (2 * PI / P) * (2 * PI / P) * (2 * PI / P) * sin(PI / 2 + 2 * PI * (k * P + ep) / P) + 12 * 0.5 * ep * ep * (2 * PI / P) * (2 * PI / P) * cos(PI / 2 + 2 * PI * (k * P + ep) / P) - 60 * 0.5 * (2 * PI / P) * sin(PI / 2 + 2 * PI * (k * P + ep) / P) * ep - 120 * 0.5 * cos(PI / 2 + 2 * PI * (k * P + ep) / P)) / (6 * ep * ep * ep * ep * ep * ep * ep);

  double c0 = (0.5 * sin(PI / 2 + 2 * PI * (k * P + ep) / P) - 0.5) / (ep * ep * ep * ep);
  double c1 = (0.5 * ep * (2 * PI / P) * cos(PI / 2 + 2 * PI * (k * P + ep) / P) - 4 * (0.5 * sin(PI / 2 + 2 * PI * (k * P + ep) / P) - 0.5)) / (ep * ep * ep * ep * ep);
  double c2 = (-0.5 * ep * ep * (2 * PI / P) * (2 * PI / P) * sin(PI / 2 + 2 * PI * (k * P + ep) / P) - 8 * 0.5 * ep * (2 * PI / P) * cos(PI / 2 + 2 * PI * (k * P + ep) / P) + 20 * (0.5 * sin(PI / 2 + 2 * PI * (k * P + ep) / P) - 0.5)) / (2 * ep * ep * ep * ep * ep * ep);
  double c3 = (-0.5 * ep * ep * ep * (2 * PI / P) * (2 * PI / P) * (2 * PI / P) * cos(PI / 2 + 2 * PI * (k * P + ep) / P) + 12 * 0.5 * ep * ep * (2 * PI / P) * (2 * PI / P) * sin(PI / 2 + 2 * PI * (k * P + ep) / P) + 60 * 0.5 * ep * (2 * PI / P) * cos(PI / 2 + 2 * PI * (k * P + ep) / P) - 120 * (0.5 * sin(PI / 2 + 2 * PI * (k * P + ep) / P) - 0.5)) / (6 * ep * ep * ep * ep * ep * ep * ep);

  double qd1 = 0.7 + b0 * (time - t0) * (time - t0) * (time - t0) * (time - t0) + b1 * (time - t0) * (time - t0) * (time - t0) * (time - t0) * (time - t1) + b2 * (time - t0) * (time - t0) * (time - t0) * (time - t0) * (time - t1) * (time - t1) + b3 * (time - t0) * (time - t0) * (time - t0) * (time - t0) * (time - t1) * (time - t1) * (time - t1);
  double qd11 = 4 * (time - t0) * (time - t0) * (time - t0) * (b0 + b1 * (time - t1) + b2 * (time - t1) * (time - t1) + b3 * (time - t1) * (time - t1) * (time - t1)) + (time - t0) * (time - t0) * (time - t0) * (time - t0) * (b1 + 2 * b2 * (time - t1) + 3 * b3 * (time - t1) * (time - t1));
  double qd21 = 12 * (time - t0) * (time - t0) * (b0 + b1 * (time - t1) + b2 * (time - t1) * (time - t1) + b3 * (time - t1) * (time - t1) * (time - t1)) + 8 * (time - t0) * (time - t0) * (time - t0) * (b1 + 2 * b2 * (time - t1) + 3 * b3 * (time - t1) * (time - t1)) + (time - t0) * (time - t0) * (time - t0) * (time - t0) * (2 * b2 + 6 * b3 * (time - t1)); //-5*(time-t1)/ep;
  double qd31 = 24 * (time - t0) * (b0 + b1 * (time - t1) + b2 * (time - t1) * (time - t1) + b3 * (time - t1) * (time - t1) * (time - t1)) + 36 * (time - t0) * (time - t0) * (b1 + 2 * b2 * (time - t1) + 3 * b3 * (time - t1) * (time - t1)) + 12 * (time - t0) * (time - t0) * (time - t0) * (2 * b2 + 6 * b3 * (time - t1)) + (time - t0) * (time - t0) * (time - t0) * (time - t0) * 6 * b3;
  double qd41 = 0.5 * (2 * PI / P) * (2 * PI / P) * (2 * PI / P) * (2 * PI / P) * cos(PI / 2 + 2 * PI * (k * P + ep) / P);

  double qd2 = 0.5 + c0 * (time - t0) * (time - t0) * (time - t0) * (time - t0) + c1 * (time - t0) * (time - t0) * (time - t0) * (time - t0) * (time - t1) + c2 * (time - t0) * (time - t0) * (time - t0) * (time - t0) * (time - t1) * (time - t1) + c3 * (time - t0) * (time - t0) * (time - t0) * (time - t0) * (time - t1) * (time - t1) * (time - t1);
  double qd12 = 4 * (time - t0) * (time - t0) * (time - t0) * (c0 + c1 * (time - t1) + c2 * (time - t1) * (time - t1) + c3 * (time - t1) * (time - t1) * (time - t1)) + (time - t0) * (time - t0) * (time - t0) * (time - t0) * (c1 + 2 * c2 * (time - t1) + 3 * c3 * (time - t1) * (time - t1));
  double qd22 = 12 * (time - t0) * (time - t0) * (c0 + c1 * (time - t1) + c2 * (time - t1) * (time - t1) + c3 * (time - t1) * (time - t1) * (time - t1)) + 8 * (time - t0) * (time - t0) * (time - t0) * (c1 + 2 * c2 * (time - t1) + 3 * c3 * (time - t1) * (time - t1)) + (time - t0) * (time - t0) * (time - t0) * (time - t0) * (2 * c2 + 6 * c3 * (time - t1));
  double qd32 = 24 * (time - t0) * (c0 + c1 * (time - t1) + c2 * (time - t1) * (time - t1) + c3 * (time - t1) * (time - t1) * (time - t1)) + 36 * (time - t0) * (time - t0) * (c1 + 2 * c2 * (time - t1) + 3 * c3 * (time - t1) * (time - t1)) + 12 * (time - t0) * (time - t0) * (time - t0) * (2 * c2 + 6 * c3 * (time - t1)) + (time - t0) * (time - t0) * (time - t0) * (time - t0) * 6 * c3;
  double qd42 = 0.5 * (2 * PI / P) * (2 * PI / P) * (2 * PI / P) * (2 * PI / P) * sin(PI / 2 + 2 * PI * time / P);

  z[13] = (time - t1); //td*ep;

  double qr11 = qd11 - gamma2 * (x - qd1);
  double qr12 = qd12 - gamma2 * (y - qd2);

  double qr21 = qd21 - gamma2 * (x1 - qd11);
  double qr22 = qd22 - gamma2 * (y1 - qd12);

  double qr31 = qd31 - gamma2 * (x2 - qd21);
  double qr32 = qd32 - gamma2 * (y2 - qd22);

  double qr41 = qd41 - gamma2 * (x3 - qd31);
  double qr42 = qd42 - gamma2 * (y3 - qd32);

  double s1 = x1 - qr11;
  double s2 = y1 - qr12;

  double ki1 = sqrt((qd1 * qd1 + qd2 * qd2 + l1 * l1 + l2 * l2) * (qd1 * qd1 + qd2 * qd2 + l1 * l1 + l2 * l2) - 2 * ((qd1 * qd1 + qd2 * qd2) * (qd1 * qd1 + qd2 * qd2) + l1 * l1 * l1 * l1 + l2 * l2 * l2 * l2));
  double ki2 = qd1 * qd1 + qd2 * qd2 + l1 * l1 - l2 * l2;
  double ki3 = qd1 * qd1 + qd2 * qd2 - l1 * l1 - l2 * l2;

  double xd = atan2(qd2, qd1) + atan2(ki1, ki2);
  double yd = -atan2(ki1, ki3);
  double xd1 = (qd1 * qd12 - qd2 * qd11) / (qd1 * qd1 + qd2 * qd2) + ((qd1 * qd11 + qd2 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2)) / (ki1 * (qd1 * qd1 + qd2 * qd2));
  double yd1 = 2 * (qd1 * qd11 + qd2 * qd12) / ki1;
  double xd2 = (((qd1 * qd21 + qd11 * qd11 + qd2 * qd22 + qd12 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2) - 2 * (qd1 * qd11 + qd2 * qd12) * (qd1 * qd11 + qd2 * qd12)) * ki1 * (qd1 * qd1 + qd2 * qd2) - (qd1 * qd11 + qd2 * qd12) * (l1 * l1 - qd1 * qd1 - qd2 * qd2 - l2 * l2) * (ki1 * 2 * (qd1 * qd11 + qd2 * qd12) - (qd1 * qd1 + qd2 * qd2) * 2 * (qd1 * qd11 + qd2 * qd12) * ki3 / ki1)) / (ki1 * (qd1 * qd1 + qd2 * qd2) * ki1 * (qd1 * qd1 + qd2 * qd2)) + ((qd1 * qd22 - qd2 * qd21) * (qd1 * qd1 + qd2 * qd2) - (qd1 * qd12 - qd2 * qd11) * 2 * (qd1 * qd11 + qd2 * qd12)) / ((qd1 * qd1 + qd2 * qd2) * (qd1 * qd1 + qd2 * qd2));
  double yd2 = 2 * ((qd1 * qd21 + qd11 * qd11 + qd2 * qd22 + qd12 * qd12) * ki1 * ki1 + (qd1 * qd11 + qd2 * qd12) * 2 * (qd1 * qd11 + qd2 * qd12) * ki3) / (ki1 * ki1 * ki1);

  double a1 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2);
  double a2 = a11 * qr11 + a12 * qr12;
  double a3 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * qr11 + d12 * qr12) / 2;
  double a4 = a21 * qr11 + a22 * qr12;

  double da1 = -m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * dv1 * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2);
  double da2 = a11 * qr21 + a12 * qr22 + da11 * qr11 + da12 * qr12;
  double da3 = m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[1] * (d11 * qr11 + d12 * qr12) / 2 + m2 * l1 * l2 * sin(q[1]) * dv0 * (d11 * qr11 + d12 * qr12) / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * (dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) / 2;
  double da4 = a21 * qr21 + a22 * qr22 + da21 * qr11 + da22 * qr12;

  double dda1 = (m2 * l1 * l2 * sin(q[1]) * velocity[1] * velocity[1] * velocity[1] - 2 * m2 * l1 * l2 * cos(q[1]) * velocity[1] * dv1) * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] * ((dd11 * qr11 + d11 * qr21 + dd12 * qr12 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - (m2 * l1 * l2 * cos(q[1]) * velocity[1] * dv1 + m2 * l1 * l2 * sin(q[1]) * ddv1) * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * dv1 * ((dd11 * qr11 + d11 * qr21 + dd12 * qr12 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - (m2 * l1 * l2 * cos(q[1]) * velocity[1] * velocity[1] + m2 * l1 * l2 * sin(q[1]) * dv1) * ((dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + (dd21 * qr11 + dd22 * qr12 + d21 * qr21 + d22 * qr22) / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((ddd11 * qr11 + dd11 * qr21 + ddd12 * qr12 + dd12 * qr22 + dd11 * qr21 + d11 * qr31 + dd12 * qr22 + d12 * qr32) + (ddd21 * qr11 + dd21 * qr21 + ddd22 * qr12 + dd22 * qr22 + dd21 * qr21 + d21 * qr31 + dd22 * qr22 + d22 * qr32) / 2);
  double dda2 = da11 * qr21 + a11 * qr31 + da12 * qr22 + a12 * qr32 + da11 * qr21 + dda11 * qr11 + da12 * qr22 + dda12 * qr12;
  double dda3 = (m2 * l1 * l2 * sin(q[1]) * ddv0 + 2 * m2 * l1 * l2 * cos(q[1]) * dv0 * velocity[1] + m2 * l1 * l2 * cos(q[1]) * velocity[0] * dv1 - m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[1] * velocity[1]) * (d11 * qr11 + d12 * qr12) / 2 + (m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[1] + m2 * l1 * l2 * sin(q[1]) * dv0) * (dd11 * qr11 + dd12 * qr12 + d11 * qr21 + d12 * qr22) + m2 * l1 * l2 * sin(q[1]) * velocity[0] * (ddd11 * qr11 + ddd12 * qr12 + dd11 * qr21 + dd12 * qr22 + dd11 * qr21 + dd12 * qr22 + d11 * qr31 + d12 * qr32) / 2;
  double dda4 = da21 * qr21 + a21 * qr31 + da22 * qr22 + a22 * qr32 + dda21 * qr11 + da21 * qr21 + dda22 * qr12 + da22 * qr22;

  double T01 = mc11 * qr21 + mc12 * qr22;
  double T02 = a1 + a2 + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2);
  double T03 = grad11 * gamma1 * s1 + grad21 * gamma1 * s2;

  double dT01 = dmc11 * qr21 + dmc12 * qr22 + mc11 * qr31 + mc12 * qr32;
  double dT02 = da1 + da2 - g * (l1 * sin(q[0]) * velocity[0] * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2);
  double dT03 = -y1 * gamma1 * s1 + x1 * gamma1 * s2 + grad11 * gamma1 * (x2 - qr21) + grad21 * gamma1 * (y2 - qr22);

  double ddT01 = ddmc11 * qr21 + ddmc12 * qr22 + 2 * dmc11 * qr31 + 2 * dmc12 * qr32 + mc11 * qr41 + mc12 * qr42;
  double ddT02 = dda1 + dda2 - g * (l1 * sin(q[0]) * dv0 * (m2 + m1 / 2) + m2 * l2 * sin(q[0] + q[1]) * (dv0 + dv1) / 2 + l1 * cos(q[0]) * velocity[0] * velocity[0] * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) / 2);
  double ddT03 = -y2 * gamma1 * s1 + x2 * gamma1 * s2 - 2 * y1 * gamma1 * (x2 - qr21) + 2 * x1 * gamma1 * (y2 - qr22) + grad11 * gamma1 * (x3 - qr31) + grad21 * gamma1 * (y3 - qr32);

  double T11 = mc21 * qr21 + mc22 * qr22;
  double T12 = a3 + a4 + g * m2 * l2 * cos(q[0] + q[1]) / 2;
  double T13 = grad12 * gamma1 * s1 + grad22 * gamma1 * s2;

  double dT11 = dmc21 * qr21 + dmc22 * qr22 + mc21 * qr31 + mc22 * qr32;
  double dT12 = da3 + da4 - g * m2 * l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) / 2;
  double dT13 = -l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * s1 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * s2 + grad12 * gamma1 * (x2 - qr21) + grad22 * gamma1 * (y2 - qr22);

  double ddT11 = ddmc21 * qr21 + ddmc22 * qr22 + 2 * dmc21 * qr31 + 2 * dmc22 * qr32 + mc21 * qr41 + mc22 * qr42;
  double ddT12 = dda3 + dda4 - g * m2 * l2 * (sin(q[0] + q[1]) * (dv0 + dv1) / 2 + cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) / 2);
  double ddT13 = l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * gamma1 * s1 - l2 * cos(q[0] + q[1]) * (dv0 + dv1) * gamma1 * s1 - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (x2 - qr21) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * (velocity[0] + velocity[1]) * gamma1 * s2 - l2 * sin(q[0] + q[1]) * (dv0 + dv1) * gamma1 * s2 - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (y2 - qr22) + grad12 * gamma1 * (x3 - qr31) + grad22 * gamma1 * (y3 - qr32) - l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (x2 - qr21) - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * gamma1 * (y2 - qr22);


  double thetad1 = xd + (T01 + T02 - T03) / K1;
  double thetad2 = yd + (T11 + T12 - T13) / K2;
  double thetad11 = xd1 + (dT01 + dT02 - dT03) / K1;
  double thetad12 = yd1 + (dT11 + dT12 - dT13) / K2;
  double thetad21 = xd2 + (ddT01 + ddT02 - ddT03) / K1;
  double thetad22 = yd2 + (ddT11 + ddT12 - ddT13) / K2;

  double s11 = (velocity[0] - xd1) + gamma2 * (q[0] - xd);
  double s12 = (velocity[1] - yd1) + gamma2 * (q[1] - yd);
  double s21 = (velocity[2] - thetad11) + gamma2 * (q[2] - thetad1);
  double s22 = (velocity[3] - thetad12) + gamma2 * (q[3] - thetad2);

  double thetar21 = thetad21 - gamma2 * (velocity[2] - thetad11);
  double thetar22 = thetad22 - gamma2 * (velocity[3] - thetad12);


  // control law
  U[0] = 0;//-(T01+T02-T03);//0;
  U[1] = 0;//-(T11+T12-T13);//0;
  U[2] = -(J1 * thetar21 + K1 * (thetad1 - xd) - gamma1 * s21);
  U[3] = -(J2 * thetar22 + K2 * (thetad2 - yd) - gamma1 * s22);

  double V1 = (s1 * s1 * (d11 * mc11 + d21 * mc21) + 2 * s1 * s2 * (d12 * mc11 + d22 * mc21) + s2 * s2 * (d12 * mc12 + d22 * mc22)) / 2 + (J1 * s21 * s21 + J2 * s22 * s22) / 2; //double V1 =(s11*s11*m11+2*s11*s12*m12+s12*s12*m22)/2+(J1*s21*s21+J2*s22*s22)/2;
  z[6] = V1 + gamma1 * gamma2 * ((q[0] - xd) * (q[0] - xd) + (q[1] - yd) * (q[1] - yd) + (q[2] - thetad1) * (q[2] - thetad1) + (q[3] - thetad2) * (q[3] - thetad2)) + (q[0] - xd - q[2] + thetad1) * K1 * (q[0] - xd - q[2] + thetad1) / 2 + (q[1] - yd - q[3] + thetad2) * K2 * (q[1] - yd - q[3] + thetad2) / 2;

  z[9] = V1;
  z[11] = P;
  z[14] = qd1;
  z[15] = qd2;
  z[16] = U[2];
  z[17] = U[3];
  z[19] = s1;
  z[20] = s12;
  z[21] = qr11;
  z[22] = qr12;
  z[23] = qr21;
  z[24] = qr22;

}


SICONOS_EXPORT void jacobFintQ(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *jacobFintQ, unsigned int sizeOfZ, double* z)
{
  // double m11 = m1*(l1*l1/4)+I1+I2+m2*(l1*l1+(l2*l2/4)+l1*l2*cos(q[1]));
  //   double m12 = (m2*l2*l2/4)+I2+m2*l1*l2*cos(q[1])/2;
  //   double m22 = (m2*l2*l2/4)+I2;

  //   double x = l1*cos(q[0])+l2*cos(q[0]+q[1]);
  //   double y = l1*sin(q[0])+l2*sin(q[0]+q[1]);

  //   double grad11 = -y;
  //   double grad12 = -l2*sin(q[0]+q[1]);
  //   double grad21 = x;
  //   double grad22 = l2*cos(q[0]+q[1]);
  //   double det = grad11*grad22-grad12*grad21;

  //   double d11 = grad22/det;
  //   double d12 = -grad12/det;
  //   double d21 = -grad21/det;
  //   double d22 = grad11/det;

  //   double dd11 =  -(sin(q[0]+q[1])*(velocity[0]+velocity[1])*sin(q[1])+cos(q[0]+q[1])*cos(q[1])*velocity[1])/(l1*sin(q[1])*sin(q[1]));
  //   double dd12 =  (cos(q[0]+q[1])*(velocity[0]+velocity[1])*sin(q[1])-sin(q[0]+q[1])*cos(q[1])*velocity[1])/(l1*sin(q[1])*sin(q[1]));
  //   double dd21 =  -dd11+(sin(q[0])*sin(q[1])*velocity[0]+cos(q[0])*cos(q[1])*velocity[1])/(l2*sin(q[1])*sin(q[1]));
  //   double dd22 =  -dd12-(cos(q[0])*sin(q[1])*velocity[0]-sin(q[0])*cos(q[1])*velocity[1])/(l2*sin(q[1])*sin(q[1]));

  //   double qr11 = z[21];
  //   double qr12 = z[22];
  //   double qr21 = z[23];
  //   double qr22 = z[24];

  //   double a11 = m11*dd11+m12*dd21;
  //   double a12 = m11*dd12+m12*dd22;
  //   double a21 = m12*dd11+m22*dd21;
  //   double a22 = m12*dd12+m22*dd22;

  //   double A1 = -m2*l1*l2*sin(q[1])*velocity[1]*((-sin(q[0]+q[1])*qr11+cos(q[0]+q[1])*qr12)/(2*l1*sin(q[1]))+(sin(q[0])*qr11-cos(q[0])*qr12)/(2*l2*sin(q[1])));
  //   double A2 = m2*l1*l2*sin(q[1])*velocity[0]*((-sin(q[0]+q[1])*qr11+cos(q[0]+q[1])*qr12)/(2*l1*sin(q[1])));
  //   double A3 = -m2*l1*l2*sin(q[1])*velocity[1]*((-cos(q[0])*qr11-sin(q[0])*qr12)/(2*l1*sin(q[1])*sin(q[1]))+(cos(q[0])*cos(q[1])*qr11+sin(q[0])*sin(q[1])*qr12)/(2*l2*sin(q[1])*sin(q[1])));
  //   double A4 = m2*l1*l2*sin(q[1])*velocity[0]*((-cos(q[0])*qr11-sin(q[0])*qr12)/(2*l1*sin(q[1])*sin(q[1])));

  //   double B1 = ((m11-m12)*sin(q[0]+q[1])/(l1*sin(q[1]))-m12*sin(q[0])/(l2*sin(q[1])))*qr21+((-m11+m12)*cos(q[0]+q[1])/(l1*sin(q[1]))+m12*cos(q[0])/(l2*sin(q[1])))*qr22;
  //   double B2 = ((m12-m22)*sin(q[0]+q[1])/(l1*sin(q[1]))+m22*sin(q[0])/(l2*sin(q[1])))*qr21+((-m12+m22)*cos(q[0]+q[1])/(l1*sin(q[1]))+m22*cos(q[0])/(l2*sin(q[1])))*qr22;
  //   double B3 = (m2*l1*l2*sin(q[1])*(d11+d21/2)+(m11-m12)*cos(q[0])/(l1*sin(q[1])*sin(q[1]))-m12*cos(q[0]-q[1])/(l2*sin(q[1])*sin(q[1])))*qr21+(m2*l1*l2*sin(q[1])*(d12+d22/2)+(m11-m12)*sin(q[0])/(l1*sin(q[1])*sin(q[1]))+m12*sin(q[0]-q[1])/(l2*sin(q[1])*sin(q[1])));
  //   double B4 = (m2*l1*l2*sin(q[1])*d11/2+(m12-m22)*cos(q[0])/(l1*sin(q[1])*sin(q[1]))-m22*cos(q[0]-q[1])/(l2*sin(q[1])*sin(q[1])))*qr21+(m2*l1*l2*sin(q[1])*(d12/2)+(m12-m22)*sin(q[0])/(l1*sin(q[1])*sin(q[1]))+m22*sin(q[0]-q[1])/(l2*sin(q[1])*sin(q[1])));

  //   double C1 = ((m11-m12)*(-(cos(q[0]+q[1])*(velocity[0]+velocity[1])*sin(q[1])-sin(q[0]+q[1])*cos(q[1])*velocity[1])/(l1*sin(q[1])*sin(q[1])*sin(q[1])))+m12*(cos(q[0])*sin(q[1])*velocity[0]-sin(q[0])*cos(q[1])*velocity[1])/(l2*sin(q[1])*sin(q[1])))*qr11+((m11-m12)*((-sin(q[0]+q[1])*(velocity[0]+velocity[1])*sin(q[1])-cos(q[0]+q[1])*cos(q[1])*velocity[1])/(l1*sin(q[1])*sin(q[1])))-m12*((-sin(q[0])*sin(q[1])*velocity[0]-cos(q[0])*cos(q[1])*velocity[1])/(l2*sin(q[1])*sin(q[1]))))*qr12;
  //   double C2 = ((m12-m22)*(-(cos(q[0]+q[1])*(velocity[0]+velocity[1])*sin(q[1])-sin(q[0]+q[1])*cos(q[1])*velocity[1])/(l1*sin(q[1])*sin(q[1])))+m22*(cos(q[0])*sin(q[1])*velocity[0]-sin(q[0])*cos(q[1])*velocity[1])/(l2*sin(q[1])*sin(q[1])))*qr11+((m12-m22)*((-sin(q[0]+q[1])*(velocity[0]+velocity[1])*sin(q[1])-cos(q[0]+q[1])*cos(q[1])*velocity[1])/(l1*sin(q[1])*sin(q[1])))-m22*((-sin(q[0])*sin(q[1])*velocity[0]-cos(q[0])*cos(q[1])*velocity[1])/(l2*sin(q[1])*sin(q[1]))))*qr12;
  //   double C3 = ((m11-m12)*(2*cos(q[0])*cos(q[1])*velocity[1]+sin(q[0])*sin(q[1])*velocity[0])/(l1*sin(q[1])*sin(q[1])*sin(q[1]))-m12*(sin(q[0])*cos(q[1])*sin(q[1])*velocity[0]+cos(q[0])*sin(q[1])*sin(q[1])*velocity[1]+2*cos(q[0])*cos(q[1])*cos(q[0])*velocity[1])/(l2*sin(q[1])*sin(q[1])*sin(q[1])))*qr11+((m11-m12)*(2*sin(q[0])*cos(q[1])*velocity[1]-cos(q[0])*sin(q[1])*velocity[0])/(l1*sin(q[1])*sin(q[1])*sin(q[1]))+m12*(cos(q[0])*cos(q[1])*sin(q[1])*velocity[0]-sin(q[0])*sin(q[1])*sin(q[1])*velocity[1]+2*sin(q[0])*cos(q[1])*cos(q[1])*velocity[1])/(l2*sin(q[1])*sin(q[1])*sin(q[1])))*qr12-m2*l1*l2*sin(q[1])*((dd11+dd21/2)*qr11+(dd12+dd22/2)*qr12);
  //   double C4 = ((m12-m22)*(2*cos(q[0])*cos(q[1])*velocity[1]+sin(q[0])*sin(q[1])*velocity[0])/(l1*sin(q[1])*sin(q[1])*sin(q[1]))-m22*(sin(q[0])*cos(q[1])*sin(q[1])*velocity[0]+cos(q[0])*sin(q[1])*sin(q[1])*velocity[1]+2*cos(q[0])*cos(q[1])*cos(q[0])*velocity[1])/(l2*sin(q[1])*sin(q[1])*sin(q[1])))*qr11+((m12-m22)*(2*sin(q[0])*cos(q[1])*velocity[1]-cos(q[0])*sin(q[1])*velocity[0])/(l1*sin(q[1])*sin(q[1])*sin(q[1]))+m22*(cos(q[0])*cos(q[1])*sin(q[1])*velocity[0]-sin(q[0])*sin(q[1])*sin(q[1])*velocity[1]+2*sin(q[0])*cos(q[1])*cos(q[1])*velocity[1])/(l2*sin(q[1])*sin(q[1])*sin(q[1])))*qr12-m2*l1*l2*sin(q[1])*(dd11*qr11+dd12*qr12)/2;

  jacobFintQ[0] = 0;
  jacobFintQ[1] = 0;
  jacobFintQ[2] = 0; //(1+gamma1*gamma2/K1)*(B1-A1-C1+gamma2*(m2*l1*l2*sin(q[1])*velocity[1]*((d11*y-d12*x)+(d21*y-d22*x)/2)-a11*y+a12*x-grad11*gamma1*y+grad22*gamma1*x));
  jacobFintQ[3] = 0; //(1+gamma1*gamma2/K2)*( B3-A3-C3+gamma2*(-m2*l1*l2*cos(q[1])*velocity[1]*((d11*qr11+d12*qr12)+(d21*qr11+d22*qr12)/2)-m2*l1*l2*sin(q[1])*velocity[1]*((d11*grad12+d12*grad22)+(d21*grad12+d22*grad22)/2)+a11*grad12+a12*grad22+grad11*gamma1*grad12+grad22*gamma1*grad22));
  jacobFintQ[4] = 0;
  jacobFintQ[5] = 0;
  jacobFintQ[6] = 0; //(1+gamma1*gamma2/K1)*(B2-A2-C2+gamma2*(-m2*l1*l2*sin(q[1])*velocity[0]*(d11*y-d12*x)/2-a21*y+a22*x-grad12*gamma1*y+grad22*gamma1*x));
  jacobFintQ[7] = 0; //(1+gamma1*gamma2/K2)*(B4-A4-C4+gamma2*(m2*l1*l2*cos(q[1])*velocity[0]*(d11*qr11+d12*qr12)/2+m2*l1*l2*sin(q[1])*velocity[0]*(d11*grad12+d12*grad22)/2+a21*grad12+a22*grad22+grad12*gamma1*grad12+grad22*gamma1*grad22));
  jacobFintQ[8] = 0;
  jacobFintQ[9] = 0;
  jacobFintQ[10] = 0; //gamma1*gamma2;
  jacobFintQ[11] = 0;
  jacobFintQ[12] = 0;
  jacobFintQ[13] = 0;
  jacobFintQ[14] = 0;
  jacobFintQ[15] = 0; //gamma1*gamma2;
}

SICONOS_EXPORT void jacobFintV(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacobFintV, unsigned int sizeOfZ, double* z)
{
  // double m11 = m1*(l1*l1/4)+I1+I2+m2*(l1*l1+(l2*l2/4)+l1*l2*cos(q[1]));
  //   double m12 = (m2*l2*l2/4)+I2+m2*l1*l2*cos(q[1])/2;
  //   double m22 = (m2*l2*l2/4)+I2;

  //   double x = l1*cos(q[0])+l2*cos(q[0]+q[1]);
  //   double y = l1*sin(q[0])+l2*sin(q[0]+q[1]);

  //   double grad11 = -y;
  //   double grad12 = -l2*sin(q[0]+q[1]);
  //   double grad21 = x;
  //   double grad22 = l2*cos(q[0]+q[1]);
  //   double det = grad11*grad22-grad12*grad21;

  //   double d11 = grad22/det;
  //   double d12 = -grad12/det;
  //   double d21 = -grad21/det;
  //   double d22 = grad11/det;

  //   double mc11 = m11*d11+m12*d21;
  //   double mc12 = m11*d12+m12*d22;
  //   double mc21 = m12*d11+m22*d21;
  //   double mc22 = m12*d12+m22*d22;

  //   double qr11 = z[21];
  //   double qr12 = z[22];

  jacobFintV[0] = 0;
  jacobFintV[1] = 0;
  jacobFintV[2] = 0; //(1+gamma1*gamma2/K1)*(gamma2*(-mc11*y+mc12*x)-grad11*gamma1*y+grad22*gamma1*x);
  jacobFintV[3] = 0; //(1+gamma1*gamma2/K2)*(gamma2*(mc11*grad12+mc12*grad22)+grad11*gamma1*grad12+grad21*gamma1*grad22+m2*l1*l2*sin(q[1])*((d11*qr11+d12*qr12)+(d21*qr11+d22*qr12)/2));;
  jacobFintV[4] = 0;
  jacobFintV[5] = 0;
  jacobFintV[6] = 0; //(1+gamma1*gamma2/K1)*(gamma2*(-mc21*y+mc22*x)-grad12*gamma1*y+grad22*gamma1*x-m2*l1*l2*sin(q[1])*(d11*qr11+d12*qr12)/2);
  jacobFintV[7] = 0; //(1+gamma1*gamma2/K2)*(gamma2*(mc21*grad12+mc22*grad22)+grad12*gamma1*grad12+grad22*gamma1*grad22);
  jacobFintV[8] = 0;
  jacobFintV[9] = 0;
  jacobFintV[10] = 0; //gamma1+gamma2*J1;
  jacobFintV[11] = 0;
  jacobFintV[12] = 0;
  jacobFintV[13] = 0;
  jacobFintV[14] = 0;
  jacobFintV[15] = 0; //gamma1+gamma2*J2;
}

SICONOS_EXPORT void h0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeOfZ, double* z)
{
  // printf("%e\n",q[0]);
  y[0] = 2 + l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  y[1] = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

}

SICONOS_EXPORT void G0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeOfZ, double* z)
{
  G[0] = 0;
  G[1] = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  G[2] = 0;
  G[3] = l2 * cos(q[0] + q[1]);
}


SICONOS_EXPORT void h01(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeOfZ, double* z)
{
  y[0] = 2 + l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
}

SICONOS_EXPORT void G01(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeOfZ, double* z)
{
  // Curious entry for the Jacobian
  G[0] = 0;
  G[1] = 0;
}


SICONOS_EXPORT void h02(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeOfZ, double* z)
{
  y[0] = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);
}

SICONOS_EXPORT void G02(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeOfZ, double* z)
{
  G[0] = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  G[1] = l2 * cos(q[0] + q[1]);
}



SICONOS_EXPORT void h3(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeOfZ, double* z)
{
  y[0] = 0.7 - l1 * cos(q[0]) - l2 * cos(q[0] + q[1]);
  y[1] = 2 + l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

}

SICONOS_EXPORT void G3(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeOfZ, double* z)
{
  G[0] = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);
  G[1] = 0;
  G[2] = l2 * sin(q[0] + q[1]);
  G[3] = 0;
}




SICONOS_EXPORT void h31(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeOfZ, double* z)
{
  y[0] = 0.7 - l1 * cos(q[0]) - l2 * cos(q[0] + q[1]);
}

SICONOS_EXPORT void G31(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeOfZ, double* z)
{
  G[0] = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);
  G[1] = l2 * sin(q[0] + q[1]);
}


SICONOS_EXPORT void h32(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeOfZ, double* z)
{
  y[0] = 2 + l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

}

SICONOS_EXPORT void G32(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeOfZ, double* z)
{
  G[0] = 0;
  G[1] = 0;
}
