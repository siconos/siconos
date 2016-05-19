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
#include <cstdio>
#include <cmath>

#pragma GCC diagnostic ignored "-Wmissing-declarations"

// BOUNCING BALL

#if defined(_MSC_VER)
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
#endif

const double m = 1; // ball mass
const double g = 9.8; // gravity
extern "C" DLLEXPORT void ballFExt(double time, unsigned int sizeOfq, double *fExt, unsigned int sizeZ, double* z)
{
  for (unsigned int i = 0; i < sizeOfq; i++)
    fExt[i] = 0.0;

  fExt[0] = -m * g;
}

// BALLBOWL

const double R = 0.5; // ball radius

extern "C" DLLEXPORT void groundFExt(double time, unsigned int sizeOfq, double *fExt, unsigned int sizeZ, double* z)
{
  for (unsigned int i = 0; i < sizeOfq; i++)
    fExt[i] = 0.0;
}

extern "C" DLLEXPORT void h0(unsigned int sizeOfq, double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* z)
{
  double R0 = 0.0;
  y[0] = q[0] + sqrt(R * R - q[1] * q[1]) - R0;
}

extern "C" DLLEXPORT void G0(unsigned int sizeOfq, double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* z)
{
  G[0] = 1.0;
  G[1] = -q[1] / (sqrt(R * R - q[1] * q[1]));
}


