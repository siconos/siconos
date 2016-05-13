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
double I1 = 1;// the moment of inertia of the first link about the axis that passes through the center of mass (parallel to the Z axis)
double I2 = 1;// the moment of inertia of the second link about the axis that passes through the center of mass (parallel to the Z axis)
double g = 9.8;//gravitational acceleration
double gamma2 = 15;
double gamma1 = 15;
double Kf = 0.5;
double P = 10;
double ep = 0.1;
double delta = 0.4;
double Del = 1;
double eps = 0.1;
double alpha = 100;


SICONOS_EXPORT void mass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  //inertia matrix using (\theta)
  mass[0]  = m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  mass[1]  = I2 + m2 * l2 * l2 / 4 + m2 * l1 * l2 * cos(q[1]) / 2;
  mass[2]  = I2 + m2 * l2 * l2 / 4 + m2 * l1 * l2 * cos(q[1]) / 2;
  mass[3]  = I2 + m2 * l2 * l2 / 4;
}

SICONOS_EXPORT void FGyr(unsigned int sizeOfq, const double *q, const double *velocity, double *FGyr, unsigned int sizeZ, double* z)
{
  FGyr[0] = -m2 * l1 * l2 * sin(q[1]) * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2);
  FGyr[1] = m2 * l1 * l2 * sin(q[1]) * velocity[0] * velocity[0] / 2;
}

SICONOS_EXPORT void jacobianFGyrq(unsigned int sizeOfq, const double *q, const double *velocity, double *jacobq, unsigned int sizeOfZ, double* z)
{
  jacobq[0] = 0;
  jacobq[1] = 0;
  jacobq[2] = -m2 * l1 * l2 * cos(q[1]) * (velocity[0] * velocity[1] + velocity[1] * velocity[1] / 2);
  jacobq[3] = m2 * l1 * l2 * cos(q[1]) * velocity[0] * velocity[0] / 2;
}

SICONOS_EXPORT void jacobianVFGyr(unsigned int sizeOfq, const double *q, const  double *velocity, double *jacobv, unsigned int sizeOfZ, double* z)
{
  jacobv[0] =   -m2 * l1 * l2 * sin(q[1]) * velocity[1];
  jacobv[1] =   m2 * l1 * l2 * sin(q[1]) * velocity[0];
  jacobv[2] =   -m2 * l1 * l2 * sin(q[1]) * (velocity[0] + velocity[1]);
  jacobv[3] =  0;
}

SICONOS_EXPORT void U(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
{

  double m11 = m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = (m2 * l2 * l2 / 4) + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = (m2 * l2 * l2 / 4) + I2;

  // generalized coordinates q=(x,y)^T
  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

  // time derivatives of x and y
  double x1 = -l1 * sin(q[0]) * velocity[0] - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]);
  double y1 = l1 * cos(q[0]) * velocity[0] + l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]);

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

  // M*d
  double mc11 = m11 * d11 + m12 * d21;
  double mc12 = m11 * d12 + m12 * d22;
  double mc21 = m12 * d11 + m22 * d21;
  double mc22 = m12 * d12 + m22 * d22;

  double qr11 = -(2 * PI / P) * 0.5 * sin(2 * PI * time / P + PI / 2) - gamma2 * (x - (0.7 + 0.5 * cos(2 * PI * time / P + PI / 2)));
  double qr12 = (2 * PI / P) * 0.5 * cos(2 * PI * time / P + PI / 2) - gamma2 * (y - 0.5 * sin(2 * PI * time / P + PI / 2));

  double qr21 = -(2 * PI / P) * (2 * PI / P) * 0.5 * cos(2 * PI * time / P + PI / 2) - gamma2 * (x1 + (2 * PI / P) * 0.5 * sin(2 * PI * time / P + PI / 2));
  double qr22 = -(2 * PI / P) * (2 * PI / P) * 0.5 * sin(2 * PI * time / P + PI / 2) - gamma2 * (y1 - (2 * PI / P) * 0.5 * cos(2 * PI * time / P + PI / 2));

  double s1 = x1 - qr11;
  double s2 = y1 - qr12;

  // M*dd
  double a11 = m11 * dd11 + m12 * dd21;
  double a12 = m11 * dd12 + m12 * dd22;
  double a21 = m12 * dd11 + m22 * dd21;
  double a22 = m12 * dd12 + m22 * dd22;

  double a1 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2);
  double a2 = a11 * qr11 + a12 * qr12;
  double a3 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * qr11 + d12 * qr12) / 2;
  double a4 = a21 * qr11 + a22 * qr12;

  double T01 = mc11 * qr21 + mc12 * qr22;
  double T02 = a1 + a2;
  double T03 = grad11 * gamma1 * s1 + grad21 * gamma1 * s2;

  double T11 = mc21 * qr21 + mc22 * qr22;
  double T12 = a3 + a4;
  double T13 = grad12 * gamma1 * s1 + grad22 * gamma1 * s2;

  // control law
  U[0] = -(T01 + T02 - T03);
  U[1] = -(T11 + T12 - T13);

  double V1 = (s1 * s1 * (d11 * mc11 + d21 * mc21) + 2 * s1 * s2 * (d12 * mc11 + d22 * mc21) + s2 * s2 * (d12 * mc12 + d22 * mc22)) / 2;
  z[6] = V1 + gamma1 * gamma2 * ((x - (0.7 + 0.5 * cos(2 * PI * time / P + PI / 2))) * (x - (0.7 + 0.5 * cos(2 * PI * time / P + PI / 2))) + (y - 0.5 * sin(2 * PI * time / P + PI / 2)) * (y - 0.5 * sin(2 * PI * time / P + PI / 2)));
  z[9] = V1;
  z[11] = P;
  z[14] = qr11;
  z[15] = qr12;
  z[18] = qr21;
  z[19] = qr22;
  z[22] = -U[0] + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2);
  z[23] = -U[1] + g * m2 * l2 * cos(q[0] + q[1]) / 2;
}


SICONOS_EXPORT void U10(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
{
  double m11 = m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = (m2 * l2 * l2 / 4) + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = (m2 * l2 * l2 / 4) + I2;

  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

  double x1 = -l1 * sin(q[0]) * velocity[0] - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]);
  double y1 = l1 * cos(q[0]) * velocity[0] + l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]);

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

  double mc11 = m11 * d11 + m12 * d21;
  double mc12 = m11 * d12 + m12 * d22;
  double mc21 = m12 * d11 + m22 * d21;
  double mc22 = m12 * d12 + m22 * d22;

  double qd1 = 0;
  double qd11 = 0;
  double qd21 = 0;
  double t2 = time - z[8];

  double qd2 = 0;
  double qd12 = 0;
  double qd22 = 0;
  double t3 = (time - z[8] - delta) / Del;

  double b0 = z[10];
  double b2 = -3 * b0 - 3 * sqrt(z[7]) * alpha;
  double b3 = 2 * b0 + 2 * sqrt(z[7]) * alpha;

  if (t2 < delta)
  {
    qd1 = 0.7 + 0.5 * cos(PI / 2 + 2 * PI * (z[8] + (t2 - delta) * (t2 - delta) * t2 / (delta * delta)) / P);
    qd11 = -(2 * PI / P) * 0.5 * sin(PI / 2 + 2 * PI * (z[8] + (t2 - delta) * (t2 - delta) * t2 / (delta * delta)) / P) * (2 * (t2 - delta) * t2 + (t2 - delta) * (t2 - delta)) / (delta * delta);
    qd21 = -(2 * PI / P) * (2 * PI / P) * 0.5 * cos(PI / 2 + 2 * PI * (z[8] + (t2 - delta) * (t2 - delta) * t2 / (delta * delta)) / P) * ((2 * (t2 - delta) * t2 + (t2 - delta) * (t2 - delta)) * (2 * (t2 - delta) * t2 + (t2 - delta) * (t2 - delta))) / (delta * delta * delta * delta) - (2 * PI / P) * 0.5 * sin(PI / 2 + 2 * PI * (z[8] + (t2 - delta) * (t2 - delta) * t2 / (delta * delta)) / P) * (2 * t2 + 4 * (t2 - delta)) / (delta * delta);
    qd2 = 0.5 * sin(PI / 2 + 2 * PI * (z[8] + (t2 - delta) * (t2 - delta) * t2 / (delta * delta)) / P);
    qd12 = (2 * PI / P) * 0.5 * cos(PI / 2 + 2 * PI * (z[8] + (t2 - delta) * (t2 - delta) * t2 / (delta * delta)) / P) * (2 * (t2 - delta) * t2 + (t2 - delta) * (t2 - delta)) / (delta * delta);
    qd22 = -(2 * PI / P) * (2 * PI / P) * 0.5 * sin(PI / 2 + 2 * PI * (z[8] + (t2 - delta) * (t2 - delta) * t2 / (delta * delta)) / P) * ((2 * (t2 - delta) * t2 + (t2 - delta) * (t2 - delta)) * (2 * (t2 - delta) * t2 + (t2 - delta) * (t2 - delta))) / (delta * delta * delta * delta) + (2 * PI / P) * 0.5 * cos(PI / 2 + 2 * PI * (z[8] + (t2 - delta) * (t2 - delta) * t2 / (delta * delta)) / P) * (2 * t2 + 4 * (t2 - delta)) / (delta * delta);
  }
  else
  {
    qd1 = z[5];
    qd11 = 0;
    qd21 = 0;
    qd2 = b3 * t3 * t3 * t3 + b2 * t3 * t3 + b0;
    qd12 = (3 * b3 * t3 * t3 + 2 * b2 * t3) / Del;
    qd22 = (6 * b3 * t3 + 2 * b2) / (Del * Del);
  }

  double qd2c = (fabs(qd2) + qd2) / 2;

  double qr11 = qd11 - gamma2 * (x - qd1);
  double qr12 = qd12 - gamma2 * (y - qd2);

  double qr21 = qd21 - gamma2 * (x1 - qd11);
  double qr22 = qd22 - gamma2 * (y1 - qd12);

  double s1 = x1 - qr11;
  double s2 = y1 - qr12;

  double s2ly = 0;
  if (qd2 >= 0)
    s2ly = y1 - qr12;
  else
    s2ly = y1 + gamma2 * y;

  double a11 = m11 * dd11 + m12 * dd21;
  double a12 = m11 * dd12 + m12 * dd22;
  double a21 = m12 * dd11 + m22 * dd21;
  double a22 = m12 * dd12 + m22 * dd22;

  double a1 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2);
  double a2 =  a11 * qr11 + a12 * qr12;
  double a3 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * qr11 + d12 * qr12) / 2;
  double a4 =  a21 * qr11 + a22 * qr12;

  double T01 = mc11 * qr21 + mc12 * qr22;
  double T02 = a1 + a2;
  double T03 = grad11 * gamma1 * s1 + grad21 * gamma1 * s2;

  double T11 = mc21 * qr21 + mc22 * qr22;
  double T12 = a3 + a4;
  double T13 = grad12 * gamma1 * s1 + grad22 * gamma1 * s2;

  U[0] = -(T01 + T02 - T03);
  U[1] = -(T11 + T12 - T13);

  double V1 = (s1 * s1 * (d11 * mc11 + d21 * mc21) + 2 * s1 * s2ly * (d12 * mc11 + d22 * mc21) + s2ly * s2ly * (d12 * mc12 + d22 * mc22)) / 2;
  z[6] = V1 + gamma1 * gamma2 * ((x - qd1) * (x - qd1) + (y - qd2c) * (y - qd2c));

  z[11] = P;
  z[14] = qr11;
  z[15] = qr12;
  z[18] = qr21;
  z[19] = qr22;
  z[22] = -U[0] + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2);
  z[23] = -U[1] + g * m2 * l2 * cos(q[0] + q[1]) / 2;
}

SICONOS_EXPORT void U11(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
{
  double m11 = m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = (m2 * l2 * l2 / 4) + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = (m2 * l2 * l2 / 4) + I2;

  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

  double x1 = -l1 * sin(q[0]) * velocity[0] - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]);
  double y1 = l1 * cos(q[0]) * velocity[0] + l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]);

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

  double mc11 = m11 * d11 + m12 * d21;
  double mc12 = m11 * d12 + m12 * d22;
  double mc21 = m12 * d11 + m22 * d21;
  double mc22 = m12 * d12 + m22 * d22;

  double qr11 = -gamma2 * (x - z[5]);
  double qr12 = -gamma2 * y;

  double qr21 = -gamma2 * x1;
  double qr22 = -gamma2 * y1;

  double s1 = x1 - qr11;
  double s2 = y1 - qr12;
  double s2bar = y1 + gamma2 * (y + sqrt(z[7]) * alpha);

  double a11 = m11 * dd11 + m12 * dd21;
  double a12 = m11 * dd12 + m12 * dd22;
  double a21 = m12 * dd11 + m22 * dd21;
  double a22 = m12 * dd12 + m22 * dd22;

  double a1 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2);
  double a2 =  a11 * qr11 + a12 * qr12;
  double a3 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * qr11 + d12 * qr12) / 2;
  double a4 =  a21 * qr11 + a22 * qr12;

  double T01 = mc11 * qr21 + mc12 * qr22;
  double T02 = a1 + a2;
  double T03 = grad11 * gamma1 * s1 + grad21 * gamma1 * s2bar;

  double T11 = mc21 * qr21 + mc22 * qr22;
  double T12 = a3 + a4;
  double T13 = grad12 * gamma1 * s1 + grad22 * gamma1 * s2bar;

  U[0] = -(T01 + T02 - T03);
  U[1] = -(T11 + T12 - T13);

  double V1 = (s1 * s1 * (d11 * mc11 + d21 * mc21) + 2 * s1 * s2 * (d12 * mc11 + d22 * mc21) + s2 * s2 * (d12 * mc12 + d22 * mc22)) / 2;
  z[6] = V1 + gamma1 * gamma2 * ((x - z[5]) * (x - z[5]) + y * y);
  z[11] = P;
  z[14] = qr11;
  z[15] = qr12;
  z[18] = qr21;
  z[19] = qr22;
  z[22] = -U[0] + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2);
  z[23] = -U[1] + g * m2 * l2 * cos(q[0] + q[1]) / 2;
}

SICONOS_EXPORT void U20(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
{
  double m11 = m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = (m2 * l2 * l2 / 4) + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = (m2 * l2 * l2 / 4) + I2;

  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

  double x1 = -l1 * sin(q[0]) * velocity[0] - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]);
  double y1 = l1 * cos(q[0]) * velocity[0] + l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]);

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

  double mc11 = m11 * d11 + m12 * d21;
  double mc12 = m11 * d12 + m12 * d22;
  double mc21 = m12 * d11 + m22 * d21;
  double mc22 = m12 * d12 + m22 * d22;

  double td = (time - z[8]) / (P / 5);

  double c0 = z[5];
  double c2 = -3 * c0 + 3 * (0.7 + sqrt(z[7]) * alpha);
  double c3 = 2 * c0 - 2 * (0.7 + sqrt(z[7]) * alpha);

  double qd1 = c3 * td * td * td + c2 * td * td + c0;
  double qd11 = (3 * c3 * td * td + 2 * c2 * td) / (P / 5);
  double qd21 = (6 * c3 * td + 2 * c2) / ((P / 5) * (P / 5));

  double qd1c = (qd1 + 0.7 - fabs(qd1 - 0.7)) / 2;

  double qr11 = qd11 - gamma2 * (x - qd1);
  double qr12 = -gamma2 * y;

  double qr21 = qd21 - gamma2 * (x1 - qd11);
  double qr22 = -gamma2 * y1;

  double s1 = x1 - qr11;
  double s2 = y1 - qr12;
  double s1ly = 0;
  if (qd1 <= 0.7)
    s1ly = x1 - qr11;
  else
    s1ly = x1 + gamma2 * (x - 0.7);

  double a11 = m11 * dd11 + m12 * dd21;
  double a12 = m11 * dd12 + m12 * dd22;
  double a21 = m12 * dd11 + m22 * dd21;
  double a22 = m12 * dd12 + m22 * dd22;

  double a1 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2);
  double a2 =  a11 * qr11 + a12 * qr12;
  double a3 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * qr11 + d12 * qr12) / 2;
  double a4 =  a21 * qr11 + a22 * qr12;

  double ld = (m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d12 - d22 / 2) - a21 - (mc21 / mc11) * (m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 - d21 / 2) - a11)) * s1 - gamma1 * mc21 * s1 / mc11 + 10;

  double T01 = mc11 * qr21 + mc12 * qr22;
  double T02 = a1 + a2;
  double T03 = grad11 * gamma1 * s1 + grad21 * gamma1 * s2;

  double T11 = mc21 * qr21 + mc22 * qr22;
  double T12 = a3 + a4;
  double T13 = grad12 * gamma1 * s1 + grad22 * gamma1 * s2;

  U[0] = -(T01 + T02 - T03 + grad21 * (Kf * z[4] - ld));
  U[1] = -(T11 + T12 - T13 + grad22 * (Kf * z[4] - ld));

  double V1 = (s1ly * s1ly * (d11 * mc11 + d21 * mc21) + 2 * s1ly * s2 * (d12 * mc11 + d22 * mc21) + s2 * s2 * (d12 * mc12 + d22 * mc22)) / 2;
  z[6] = V1 + gamma1 * gamma2 * ((x - qd1c) * (x - qd1c) + y * y);
  z[11] = P;
  z[14] = qr11;
  z[15] = qr12;
  z[18] = qr21;
  z[19] = qr22;
  z[22] = -U[0] + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2);
  z[23] = -U[1] + g * m2 * l2 * cos(q[0] + q[1]) / 2; //(k*P-time)*(1+Kf);
}

SICONOS_EXPORT void U21(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
{
  double m11 = m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = (m2 * l2 * l2 / 4) + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = (m2 * l2 * l2 / 4) + I2;

  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

  double x1 = -l1 * sin(q[0]) * velocity[0] - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]);
  double y1 = l1 * cos(q[0]) * velocity[0] + l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]);

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

  double mc11 = m11 * d11 + m12 * d21;
  double mc12 = m11 * d12 + m12 * d22;
  double mc21 = m12 * d11 + m22 * d21;
  double mc22 = m12 * d12 + m22 * d22;

  double qr11 = -gamma2 * (x - 0.7);
  double qr12 = -gamma2 * y;

  double qr21 = -gamma2 * x1;
  double qr22 = -gamma2 * y1;

  double s1 = x1 - qr11;
  double s2 = y1 - qr12;
  double s1bar = x1 + gamma2 * (x - 0.7 - sqrt(z[7]) * alpha);

  double a11 = m11 * dd11 + m12 * dd21;
  double a12 = m11 * dd12 + m12 * dd22;
  double a21 = m12 * dd11 + m22 * dd21;
  double a22 = m12 * dd12 + m22 * dd22;

  double a1 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2);
  double a2 =  a11 * qr11 + a12 * qr12;
  double a3 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * qr11 + d12 * qr12) / 2;
  double a4 =  a21 * qr11 + a22 * qr12;

  double ld = (m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d12 - d22 / 2) - a21 - (mc21 / mc11) * (m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 - d21 / 2) - a11)) * s1 - gamma1 * mc21 * s1 / mc11 + 10;

  double T01 = mc11 * qr21 + mc12 * qr22;
  double T02 = a1 + a2;
  double T03 = grad11 * gamma1 * s1bar + grad21 * gamma1 * s2;

  double T11 = mc21 * qr21 + mc22 * qr22;
  double T12 = a3 + a4;
  double T13 = grad12 * gamma1 * s1bar + grad22 * gamma1 * s2;

  U[0] = -(T01 + T02 - T03 + grad21 * (Kf * z[4] - ld));
  U[1] = -(T11 + T12 - T13 + grad22 * (Kf * z[4] - ld));

  double V1 = (s1 * s1 * (d11 * mc11 + d21 * mc21) + 2 * s1 * s2 * (d12 * mc11 + d22 * mc21) + s2 * s2 * (d12 * mc12 + d22 * mc22)) / 2;
  z[6] = V1 + gamma1 * gamma2 * ((x - 0.7) * (x - 0.7) + y * y);
  z[11] = P;
  z[14] = qr11;
  z[15] = qr12;
  z[18] = qr21;
  z[19] = qr22;
  z[22] = -U[0] + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2);
  z[23] = -U[1] + g * m2 * l2 * cos(q[0] + q[1]) / 2; //(k*P-time)*(1+Kf);
}

SICONOS_EXPORT void U22(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
{
  double m11 = m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = (m2 * l2 * l2 / 4) + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = (m2 * l2 * l2 / 4) + I2;

  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

  double x1 = -l1 * sin(q[0]) * velocity[0] - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]);
  double y1 = l1 * cos(q[0]) * velocity[0] + l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]);

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

  double mc11 = m11 * d11 + m12 * d21;
  double mc12 = m11 * d12 + m12 * d22;
  double mc21 = m12 * d11 + m22 * d21;
  double mc22 = m12 * d12 + m22 * d22;

  double k = trunc(time / P) + 1;
  double td = (time - z[8]) / (k * P - z[8]);

  double c2 = 1.5;
  double c3 = -1;

  double qd2 = c3 * td * td * td + c2 * td * td;
  double qd12 = (3 * c3 * td * td + 2 * c2 * td) / (k * P - z[8]);
  double qd22 = (6 * c3 * td + 2 * c2) / ((k * P - z[8]) * (k * P - z[8]));

  double qr11 = -gamma2 * (x - 0.7);
  double qr12 = qd12 - gamma2 * (y - qd2);

  double qr21 = -gamma2 * x1;
  double qr22 = qd22 - gamma2 * (y1 - qd12);

  double s1 = x1 - qr11;
  double s2 = y1 - qr12;

  double a11 = m11 * dd11 + m12 * dd21;
  double a12 = m11 * dd12 + m12 * dd22;
  double a21 = m12 * dd11 + m22 * dd21;
  double a22 = m12 * dd12 + m22 * dd22;

  double a1 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2);
  double a2 =  a11 * qr11 + a12 * qr12;
  double a3 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * qr11 + d12 * qr12) / 2;
  double a4 =  a21 * qr11 + a22 * qr12;

  double ld = (m2 * l1 * l2 * sin(q[1]) * velocity[1] * d11 / 2 - a12 - (mc12 / mc22) * (m2 * l1 * l2 * sin(q[1]) * velocity[1] * d21 / 2 - a22)) * s2 - gamma1 * mc12 * s2 / mc22 - (k * P - time) * (1 + Kf);
  z[12] = 1 + (-mc11 * ld + fabs(mc11 * ld)) / 2;

  double T01 = mc11 * qr21 + mc12 * qr22;
  double T02 = a1 + a2;
  double T03 = grad11 * gamma1 * s1 + grad21 * gamma1 * s2;

  double T11 = mc21 * qr21 + mc22 * qr22;
  double T12 = a3 + a4;
  double T13 = grad12 * gamma1 * s1 + grad22 * gamma1 * s2;

  U[0] = -(T01 + T02 - T03 + grad11 * (Kf * z[24] - ld));
  U[1] = -(T11 + T12 - T13 + grad12 * (Kf * z[24] - ld));

  double V1 = (s1 * s1 * (d11 * mc11 + d21 * mc21) + 2 * s1 * s2 * (d12 * mc11 + d22 * mc21) + s2 * s2 * (d12 * mc12 + d22 * mc22)) / 2;
  z[6] = V1 + gamma1 * gamma2 * ((x - 0.7) * (x - 0.7) + (y - qd2) * (y - qd2));
  z[11] = P;
  z[14] = qr11;
  z[15] = qr12;
  z[18] = qr21;
  z[19] = qr22;
  z[22] = -U[0] + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2);
  z[23] = -U[1] + g * m2 * l2 * cos(q[0] + q[1]) / 2; //(k*P-time)*(1+Kf);
}


SICONOS_EXPORT void U3(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
{
  double m11 = m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = (m2 * l2 * l2 / 4) + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = (m2 * l2 * l2 / 4) + I2;

  // generalized coordinates q=(x,y)^T
  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

  // time derivatives of x and y
  double x1 = -l1 * sin(q[0]) * velocity[0] - l2 * sin(q[0] + q[1]) * (velocity[0] + velocity[1]);
  double y1 = l1 * cos(q[0]) * velocity[0] + l2 * cos(q[0] + q[1]) * (velocity[0] + velocity[1]);

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

  double dd11 =  -(sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) + cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd12 = (cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]));
  double dd21 =  -dd11 + (sin(q[0]) * sin(q[1]) * velocity[0] + cos(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));
  double dd22 =  -dd12 - (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]));

  double mc11 = m11 * d11 + m12 * d21;
  double mc12 = m11 * d12 + m12 * d22;
  double mc21 = m12 * d11 + m22 * d21;
  double mc22 = m12 * d12 + m22 * d22;

  double k = trunc(z[8] / P);
  double td = (time - k * P) / ep;
  double b = -z[10] * ep * P / (4 * 0.5 * PI);

  double qd1 = 0.7 + 0.5 * cos(PI / 2 + 2 * PI * (k * P + (time - k * P) * sin(b * td + (PI / 2 - b) * td * td)) / P);
  double qd11 = -0.5 * (2 * PI / P) * sin(PI / 2 + 2 * PI * (k * P + (time - k * P) * sin(b * td + (PI / 2 - b) * td * td)) / P) * (sin(b * td + (PI / 2 - b) * td * td) + (time - k * P) * cos(b * td + (PI / 2 - b) * td * td) * (b / ep + 2 * (PI / 2 - b) * td / ep));
  double qd21 = -0.5 * (2 * PI / P) * (2 * PI / P) * cos(PI / 2 + 2 * PI * (k * P + (time - k * P) * sin(b * td + (PI / 2 - b) * td * td)) / P) * (sin(b * td + (PI / 2 - b) * td * td) + (time - k * P) * cos(b * td + (PI / 2 - b) * td * td) * (b / ep + 2 * (PI / 2 - b) * td / ep)) - 0.5 * (2 * PI / P) * sin(PI / 2 + 2 * PI * (k * P + (time - k * P) * sin(b * td + (PI / 2 - b) * td * td)) / P) * (2 * cos(b * td + (PI / 2 - b) * td * td) * (b / ep + 2 * (PI / 2 - b) * td / ep) + (time - k * P) * (2 * cos(b * td + (PI / 2 - b) * td * td) * (PI / 2 - b) / (ep * ep) - sin(b * td + (PI / 2 - b) * td * td) * (b / ep + 2 * (PI / 2 - b) * td / ep) * (b / ep + 2 * (PI / 2 - b) * td / ep)));

  double qd2 = 0.5 * sin(PI / 2 + 2 * PI * time / P);
  double qd12 = 0.5 * (2 * PI / P) * cos(PI / 2 + 2 * PI * time / P);
  double qd22 = -0.5 * (2 * PI / P) * (2 * PI / P) * sin(PI / 2 + 2 * PI * time / P);

  z[13] = (time - k * P); //td*ep;

  double qr11 = qd11 - gamma2 * (x - qd1);
  double qr12 = qd12 - gamma2 * (y - qd2);

  double qr21 = qd21 - gamma2 * (x1 - qd11);
  double qr22 = qd22 - gamma2 * (y1 - qd12);

  double s1 = x1 - qr11;
  double s2 = y1 - qr12;

  double a11 = m11 * dd11 + m12 * dd21;
  double a12 = m11 * dd12 + m12 * dd22;
  double a21 = m12 * dd11 + m22 * dd21;
  double a22 = m12 * dd12 + m22 * dd22;

  double a1 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2);
  double a2 = a11 * qr11 + a12 * qr12;
  double a3 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * qr11 + d12 * qr12) / 2;
  double a4 = a21 * qr11 + a22 * qr12;

  double T01 = mc11 * qr21 + mc12 * qr22;
  double T02 = a1 + a2;
  double T03 = grad11 * gamma1 * s1 + grad21 * gamma1 * s2;

  double T11 = mc21 * qr21 + mc22 * qr22;
  double T12 = a3 + a4;
  double T13 = grad12 * gamma1 * s1 + grad22 * gamma1 * s2;

  U[0] = -(T01 + T02 - T03);
  U[1] = -(T11 + T12 - T13);

  double V1 = (s1 * s1 * (d11 * mc11 + d21 * mc21) + 2 * s1 * s2 * (d12 * mc11 + d22 * mc21) + s2 * s2 * (d12 * mc12 + d22 * mc22)) / 2;
  z[6] = V1 + gamma1 * gamma2 * ((x - qd1) * (x - qd1) + (y - qd2) * (y - qd2));
  z[11] = P;
  z[7] = 0;
  z[14] = qr11;
  z[15] = qr12;
  z[18] = qr21;
  z[19] = qr22;
  z[22] = -U[0] + g * (l1 * cos(q[0]) * (m2 + m1 / 2) + m2 * l2 * cos(q[0] + q[1]) / 2);
  z[23] = -U[1] + g * m2 * l2 * cos(q[0] + q[1]) / 2;
}


SICONOS_EXPORT void jacobFintQ(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *jacobFintQ, unsigned int sizeOfZ, double* z)
{
  double m11 = m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = (m2 * l2 * l2 / 4) + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = (m2 * l2 * l2 / 4) + I2;

  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

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

  double qr11 = z[16];
  double qr12 = z[17];
  double qr21 = z[20];
  double qr22 = z[21];

  double a11 = m11 * dd11 + m12 * dd21;
  double a12 = m11 * dd12 + m12 * dd22;
  double a21 = m12 * dd11 + m22 * dd21;
  double a22 = m12 * dd12 + m22 * dd22;

  double A1 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((-sin(q[0] + q[1]) * qr11 + cos(q[0] + q[1]) * qr12) / (2 * l1 * sin(q[1])) + (sin(q[0]) * qr11 - cos(q[0]) * qr12) / (2 * l2 * sin(q[1])));
  double A2 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * ((-sin(q[0] + q[1]) * qr11 + cos(q[0] + q[1]) * qr12) / (2 * l1 * sin(q[1])));
  double A3 = -m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((-cos(q[0]) * qr11 - sin(q[0]) * qr12) / (2 * l1 * sin(q[1]) * sin(q[1])) + (cos(q[0]) * cos(q[1]) * qr11 + sin(q[0]) * sin(q[1]) * qr12) / (2 * l2 * sin(q[1]) * sin(q[1])));
  double A4 = m2 * l1 * l2 * sin(q[1]) * velocity[0] * ((-cos(q[0]) * qr11 - sin(q[0]) * qr12) / (2 * l1 * sin(q[1]) * sin(q[1])));

  double B1 = ((m11 - m12) * sin(q[0] + q[1]) / (l1 * sin(q[1])) - m12 * sin(q[0]) / (l2 * sin(q[1]))) * qr21 + ((-m11 + m12) * cos(q[0] + q[1]) / (l1 * sin(q[1])) + m12 * cos(q[0]) / (l2 * sin(q[1]))) * qr22;
  double B2 = ((m12 - m22) * sin(q[0] + q[1]) / (l1 * sin(q[1])) + m22 * sin(q[0]) / (l2 * sin(q[1]))) * qr21 + ((-m12 + m22) * cos(q[0] + q[1]) / (l1 * sin(q[1])) + m22 * cos(q[0]) / (l2 * sin(q[1]))) * qr22;
  double B3 = (m2 * l1 * l2 * sin(q[1]) * (d11 + d21 / 2) + (m11 - m12) * cos(q[0]) / (l1 * sin(q[1]) * sin(q[1])) - m12 * cos(q[0] - q[1]) / (l2 * sin(q[1]) * sin(q[1]))) * qr21 + (m2 * l1 * l2 * sin(q[1]) * (d12 + d22 / 2) + (m11 - m12) * sin(q[0]) / (l1 * sin(q[1]) * sin(q[1])) + m12 * sin(q[0] - q[1]) / (l2 * sin(q[1]) * sin(q[1])));
  double B4 = (m2 * l1 * l2 * sin(q[1]) * d11 / 2 + (m12 - m22) * cos(q[0]) / (l1 * sin(q[1]) * sin(q[1])) - m22 * cos(q[0] - q[1]) / (l2 * sin(q[1]) * sin(q[1]))) * qr21 + (m2 * l1 * l2 * sin(q[1]) * (d12 / 2) + (m12 - m22) * sin(q[0]) / (l1 * sin(q[1]) * sin(q[1])) + m22 * sin(q[0] - q[1]) / (l2 * sin(q[1]) * sin(q[1])));

  double C1 = ((m11 - m12) * (-(cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1]))) + m12 * (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]))) * qr11 + ((m11 - m12) * ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]))) - m12 * ((-sin(q[0]) * sin(q[1]) * velocity[0] - cos(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1])))) * qr12;
  double C2 = ((m12 - m22) * (-(cos(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - sin(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]))) + m22 * (cos(q[0]) * sin(q[1]) * velocity[0] - sin(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]))) * qr11 + ((m12 - m22) * ((-sin(q[0] + q[1]) * (velocity[0] + velocity[1]) * sin(q[1]) - cos(q[0] + q[1]) * cos(q[1]) * velocity[1]) / (l1 * sin(q[1]) * sin(q[1]))) - m22 * ((-sin(q[0]) * sin(q[1]) * velocity[0] - cos(q[0]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1])))) * qr12;
  double C3 = ((m11 - m12) * (2 * cos(q[0]) * cos(q[1]) * velocity[1] + sin(q[0]) * sin(q[1]) * velocity[0]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1])) - m12 * (sin(q[0]) * cos(q[1]) * sin(q[1]) * velocity[0] + cos(q[0]) * sin(q[1]) * sin(q[1]) * velocity[1] + 2 * cos(q[0]) * cos(q[1]) * cos(q[0]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]))) * qr11 + ((m11 - m12) * (2 * sin(q[0]) * cos(q[1]) * velocity[1] - cos(q[0]) * sin(q[1]) * velocity[0]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1])) + m12 * (cos(q[0]) * cos(q[1]) * sin(q[1]) * velocity[0] - sin(q[0]) * sin(q[1]) * sin(q[1]) * velocity[1] + 2 * sin(q[0]) * cos(q[1]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]))) * qr12 - m2 * l1 * l2 * sin(q[1]) * ((dd11 + dd21 / 2) * qr11 + (dd12 + dd22 / 2) * qr12);
  double C4 = ((m12 - m22) * (2 * cos(q[0]) * cos(q[1]) * velocity[1] + sin(q[0]) * sin(q[1]) * velocity[0]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1])) - m22 * (sin(q[0]) * cos(q[1]) * sin(q[1]) * velocity[0] + cos(q[0]) * sin(q[1]) * sin(q[1]) * velocity[1] + 2 * cos(q[0]) * cos(q[1]) * cos(q[0]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]))) * qr11 + ((m12 - m22) * (2 * sin(q[0]) * cos(q[1]) * velocity[1] - cos(q[0]) * sin(q[1]) * velocity[0]) / (l1 * sin(q[1]) * sin(q[1]) * sin(q[1])) + m22 * (cos(q[0]) * cos(q[1]) * sin(q[1]) * velocity[0] - sin(q[0]) * sin(q[1]) * sin(q[1]) * velocity[1] + 2 * sin(q[0]) * cos(q[1]) * cos(q[1]) * velocity[1]) / (l2 * sin(q[1]) * sin(q[1]) * sin(q[1]))) * qr12 - m2 * l1 * l2 * sin(q[1]) * (dd11 * qr11 + dd12 * qr12) / 2;

  jacobFintQ[0] = B1 - A1 - C1 + gamma2 * (m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * y - d12 * x) + (d21 * y - d22 * x) / 2) - a11 * y + a12 * x - grad11 * gamma1 * y + grad22 * gamma1 * x);
  jacobFintQ[1] = B2 - A2 - C2 + gamma2 * (-m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * y - d12 * x) / 2 - a21 * y + a22 * x - grad12 * gamma1 * y + grad22 * gamma1 * x);
  jacobFintQ[2] = B3 - A3 - C3 + gamma2 * (-m2 * l1 * l2 * cos(q[1]) * velocity[1] * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2) - m2 * l1 * l2 * sin(q[1]) * velocity[1] * ((d11 * grad12 + d12 * grad22) + (d21 * grad12 + d22 * grad22) / 2) + a11 * grad12 + a12 * grad22 + grad11 * gamma1 * grad12 + grad22 * gamma1 * grad22);
  jacobFintQ[3] = B4 - A4 - C4 + gamma2 * (m2 * l1 * l2 * cos(q[1]) * velocity[0] * (d11 * qr11 + d12 * qr12) / 2 + m2 * l1 * l2 * sin(q[1]) * velocity[0] * (d11 * grad12 + d12 * grad22) / 2 + a21 * grad12 + a22 * grad22 + grad12 * gamma1 * grad12 + grad22 * gamma1 * grad22);

}

SICONOS_EXPORT void jacobFintV(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacobFintV, unsigned int sizeOfZ, double* z)
{
  double m11 = m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q[1]));
  double m12 = (m2 * l2 * l2 / 4) + I2 + m2 * l1 * l2 * cos(q[1]) / 2;
  double m22 = (m2 * l2 * l2 / 4) + I2;

  double x = l1 * cos(q[0]) + l2 * cos(q[0] + q[1]);
  double y = l1 * sin(q[0]) + l2 * sin(q[0] + q[1]);

  double grad11 = -y;
  double grad12 = -l2 * sin(q[0] + q[1]);
  double grad21 = x;
  double grad22 = l2 * cos(q[0] + q[1]);
  double det = grad11 * grad22 - grad12 * grad21;

  double d11 = grad22 / det;
  double d12 = -grad12 / det;
  double d21 = -grad21 / det;
  double d22 = grad11 / det;

  double mc11 = m11 * d11 + m12 * d21;
  double mc12 = m11 * d12 + m12 * d22;
  double mc21 = m12 * d11 + m22 * d21;
  double mc22 = m12 * d12 + m22 * d22;

  double qr11 = z[16];
  double qr12 = z[17];

  jacobFintV[0] = gamma2 * (-mc11 * y + mc12 * x) - grad11 * gamma1 * y + grad22 * gamma1 * x;
  jacobFintV[1] = gamma2 * (-mc21 * y + mc22 * x) - grad12 * gamma1 * y + grad22 * gamma1 * x - m2 * l1 * l2 * sin(q[1]) * (d11 * qr11 + d12 * qr12) / 2;
  jacobFintV[2] = gamma2 * (mc11 * grad12 + mc12 * grad22) + grad11 * gamma1 * grad12 + grad21 * gamma1 * grad22 + m2 * l1 * l2 * sin(q[1]) * ((d11 * qr11 + d12 * qr12) + (d21 * qr11 + d22 * qr12) / 2);
  jacobFintV[3] = gamma2 * (mc21 * grad12 + mc22 * grad22) + grad12 * gamma1 * grad12 + grad22 * gamma1 * grad22;
}

SICONOS_EXPORT void h0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeOfZ, double* z)
{
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
