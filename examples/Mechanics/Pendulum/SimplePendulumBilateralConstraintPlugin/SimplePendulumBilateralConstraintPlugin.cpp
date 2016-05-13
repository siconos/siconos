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
using namespace std;

#undef restrict
#define restrict __restrict

extern "C" double L;

SICONOS_EXPORT void h0(unsigned int sizeOfq,  double* restrict  q, unsigned int sizeOfY, double* restrict  y, unsigned int sizeZ, double* restrict  z)
{
  y[0] = pow(L, 2) - (pow(q[0], 2) + pow(q[1], 2));
}

SICONOS_EXPORT void G0(unsigned int sizeOfq,  double* restrict  q, unsigned int sizeOfY, double* restrict  G, unsigned int sizeZ, double* restrict  z)
{
  G[0] = -2.0 * q[0];
  G[1] = -2.0 * q[1];
}

SICONOS_EXPORT void G0dot(unsigned int sizeOfq,  double* restrict  q, unsigned int sizeOfqdot,  double* restrict  qdot, double* restrict  S, unsigned int sizeOfZ, double* restrict  z)
{
  S[0] = -2.0 * qdot[0];
  S[1] = -2.0 * qdot[1];
}
