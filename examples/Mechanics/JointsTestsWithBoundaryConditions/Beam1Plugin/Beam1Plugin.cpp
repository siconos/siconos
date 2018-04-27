/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

// for M_PI
#if defined(_MSC_VER)
#define _USE_MATH_DEFINES
#endif
#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <math.h>
#include <stdio.h>

SICONOS_EXPORT void prescribedvelocity(double time, unsigned int sizeofprescribedvelocity, double *pv)
{
  /* the plugin implements v(t) = C + A cos(omega *t) */

  double C = 0.0 ;
  double omega = M_PI / 2.0;
  double A = 1.0;

  pv[0] =  C +  A * cos(omega * time);
  //printf("prescribed velocity = %e\n", pv[0]);
}
