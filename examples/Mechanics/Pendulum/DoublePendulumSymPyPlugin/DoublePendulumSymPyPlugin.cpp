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
#include <stdio.h>
#include <math.h>

#include "pendulum.h"

double l1 = 1.0 ;
double l2 = 1.0 ;

SICONOS_EXPORT void mass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{

  int n = sizeOfq;

  int n1 = n*n;
 
  // mass set to zero
  for (int i = 0; i < n1; i++)
  {

    mass[i] = 0.0;	

  }

  Inertia(mass, q);

}

SICONOS_EXPORT void FGyr(unsigned int sizeOfq, const double *q, const double *velocity, double *FGyr, unsigned int sizeZ, double* z)
{

  int n = sizeOfq;
 
  // FGyr set to zero
  for (int i = 0; i < n; i++)
  {

    FGyr[i] = 0.0;	

  }

  FGyrEffects(FGyr, q, velocity);

}

SICONOS_EXPORT void jacobianFGyrq(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{

  int n = sizeOfq;

  int n1 = n*n;
 
  // mass set to zero
  for (int i = 0; i < n1; i++)
  {

    jacob[i] = 0.0;	

  }

  JacobianQFGyr(jacob, q, velocity);

}

SICONOS_EXPORT void jacobianVFGyr(unsigned int sizeOfq, const double *q, const  double *velocity, double *jacob, unsigned int sizeZ, double* z)
{

  int n = sizeOfq;

  int n1 = n*n;
 
  // mass set to zero
  for (int i = 0; i < n1; i++)
  {

    jacob[i] = 0.0;	

  }

  JacobianVFGyr(jacob, q, velocity);

}

SICONOS_EXPORT void h0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* z)
{
  y[0] = l1 * sin(q[0]);
}

SICONOS_EXPORT void G0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* z)
{
  G[0] = l1 * cos(q[0]);
  G[1] = 0.0;
}

SICONOS_EXPORT void h1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* z)
{
  y[0] = l1 * sin(q[0]) + l2 * sin(q[1]);
}

SICONOS_EXPORT void G1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* z)
{
  G[0] = l1 * cos(q[0]);
  G[1] = l2 * cos(q[1]);
}

