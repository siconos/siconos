/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include <math.h>
const double omega = 1.4;
const double g = 10;
const double m = 1;

SICONOS_EXPORT void FExt(double time, unsigned int sizeOfq, double *fExt, unsigned int sizeZ, double* z)
{
  fExt[0] = sin(omega * time);
  fExt[1] = -m * g;
}
