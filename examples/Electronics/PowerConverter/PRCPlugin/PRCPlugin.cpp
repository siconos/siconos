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

#if defined(_MSC_VER)
#define _USE_MATH_DEFINES
#endif
#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <math.h>
#include <iostream>
#include <stdio.h>

//const double PI = 3.14159;


// ===== Dynamical System =====

// function to compute u
SICONOS_EXPORT void computeU(double time, unsigned int sizeU, double *U, unsigned int sizeZ, double* z)
{
  double f = 55000.0;
  if (time == 0.0)
    U[0] = 0;
  else
  {
    U[0] = z[0] * sin(2.0 * M_PI * f * time) / fabs(sin(2.0 * M_PI * f * time));
  }
}

