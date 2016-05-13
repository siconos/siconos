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
#include "Robot.h"
#include <stdio.h>
#include <math.h>

const unsigned int n0 = 7; // Real problem dimension

SICONOS_EXPORT void mass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  unsigned int n = sizeOfq;
  unsigned int n1 = n0 * n0;

  unsigned int i, j;

  double massTmp[n1];
  // mass set to zero
  for (i = 0; i < n1; i++)
    massTmp[i] = 0.0;

  double qTmp[n0];
  for (i = 1; i < n0; i++)
    qTmp[i] = 0.0;

  // compute mass matrix
  Inertia(massTmp, qTmp);

  // Motor Inertia effect
  massTmp[0]  += 0.75;
  massTmp[8]  += 0.75;
  massTmp[16] += 0.2125;
  massTmp[24] += 0.2125;
  massTmp[32] += 0.00575;
  massTmp[40] += 0.00575;
  massTmp[48] += 0.00575;

  // compute reduced mass matrix
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      mass[j + n * i] = massTmp[n0 * 2 * i + n0 + 2 * j + 1];
    }
  }
}

SICONOS_EXPORT void FGyr(unsigned int sizeOfq, const double *q, const double *velocity, double *FGyr, unsigned int sizeZ, double* z)
{
  unsigned int n = sizeOfq;

  unsigned int i;

  double FGyrTmp[n0];
  // FGyr set to zero
  for (i = 0; i < n0; i++)
    FGyrTmp[i] = 0.0;

  double qTmp[n0], vTmp[n0];
  for (i = 1; i < n0; i++)
  {
    qTmp[i] = 0.0;
    vTmp[i] = 0.0;
  }
  for (i = 0; i < n; i++)
  {
    qTmp[2 * i + 1] = q[i];
    vTmp[2 * i + 1] = velocity[i];
  }

  // compute mass matrix
  NLEffects(FGyrTmp, qTmp, vTmp);

  FGyrTmp[0] += 10 * vTmp[0];
  FGyrTmp[1] += 10 * vTmp[1];
  FGyrTmp[2] += 5 * vTmp[2];
  FGyrTmp[3] += 5 * vTmp[3];
  FGyrTmp[4] += 2 * vTmp[4];
  FGyrTmp[5] += 2 * vTmp[5];
  FGyrTmp[6] += 2 * vTmp[6];

  // compute reduced FGyr vector
  for (i = 0; i < n; i++)
    FGyr[i] = FGyrTmp[2 * i + 1];

}

SICONOS_EXPORT void jacobianFGyrq(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  unsigned int n = sizeOfq;
  unsigned int n1 = n0 * n0;

  unsigned int i, j;

  double jacoTmp[n1];
  // set to zero
  for (i = 0; i < n1; i++)
    jacoTmp[i] = 0.0;

  double qTmp[n0], vTmp[n0];

  for (i = 0; i < n0; i++)
  {
    qTmp[i] = 0.0;
    vTmp[i] = 0.0;
  }
  for (i = 0; i < n; i++)
  {
    qTmp[2 * i + 1] = q[i];
    vTmp[2 * i + 1] = velocity[i];
  }

  // compute jacobian matrix
  JacobianQNLEffects(jacoTmp, qTmp, vTmp);

  // compute reduced jacobian matrix
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      jacob[j + n * i] = jacoTmp[n0 * 2 * i + n0 + 2 * j + 1];
    }
  }
}

SICONOS_EXPORT void jacobianVFGyr(unsigned int sizeOfq, const double *q, const  double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  unsigned int n = sizeOfq;
  unsigned int n1 = n0 * n0;

  unsigned int i, j;

  double jacoTmp[n1];
  // set to zero
  for (i = 0; i < n1; i++)
    jacoTmp[i] = 0.0;

  double qTmp[n0], vTmp[n0];
  for (i = 0; i < n0; i++)
  {
    qTmp[i] = 0.0;
    vTmp[i] = 0.0;
  }
  for (i = 0; i < n; i++)
  {
    qTmp[2 * i + 1] = q[i];
    vTmp[2 * i + 1] = velocity[i];
  }

  // compute jacobian matrix
  JacobianVNLEffects(jacoTmp, qTmp, vTmp);

  // compute reduced jacobian matrix
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      jacob[j + n * i] = jacoTmp[n0 * 2 * i + n0 + 2 * j + 1];
    }
  }
}

SICONOS_EXPORT void h0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* z)
{
  y[0] = 0.45 * cos(q[0]);
}

SICONOS_EXPORT void G0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* z)
{
  G[0] = -0.45 * sin(q[0]);
  G[1] = 0.0;
}


SICONOS_EXPORT void h1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* z)
{
  y[0] = 0.45 * cos(q[0]) + 0.48 * cos(q[1] + q[0]);
}

SICONOS_EXPORT void G1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* z)
{
  G[0] = -0.45 * sin(q[0]) - 0.48 * sin(q[1] + q[0]);
  G[1] = -0.48 * sin(q[1] + q[0]);
}

SICONOS_EXPORT void h2(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* z)
{
  //  printf("BJBJBJBJ%f\n",y[0]);
  //printf("BJBJBJBsddssdJ%f\n",y[1]);
  y[0] = 0.45 * cos(q[0]);
  y[1] = 0.45 * cos(q[0]) + 0.48 * cos(q[1] + q[0]);
  //printf("AAAAAAAAABJBJBJBJ%f\n",y[0]);
  //  printf("AAAAAABJBJBJBsddssdJ%f\n",y[1]);
}

SICONOS_EXPORT void G2(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* z)
{
  G[0] = -0.45 * sin(q[0]);
  G[1] = -0.45 * sin(q[0]) - 0.48 * sin(q[1] + q[0]);
  G[2] = 0.0;
  G[3] = -0.48 * sin(q[1] + q[0]);
  G[4] = 0.0;
  G[5] = 0.0;
}

