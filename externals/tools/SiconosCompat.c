/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#ifdef _MSC_VER

// before MSVC 2013
#if _MSC_VER < 1800
#include <float.h>

float _sqrtf(float x)
{
  return sqrtf(x);
}
float _logf(float x)
{
  return logf(x);
}

extern "C" long int lroundf(float x)
{
  return (long int)floorl(x + .5);
}

// Classify floating point number - usually defined in math.h
extern "C" int __fpclassify(double x)
{
  return _fpclass(x);
}/* This is really bad --xhub */

#ifdef __cplusplus
namespace std
{
int isfinite(double x)
{
  return _finite(x);
}
}

#endif

#endif /* _MSC_VER < 1800 */



extern "C" double __cdecl __powidf2(double a, int b)
{
  const int recip = b < 0;
  double r = 1;
  while(1)
  {
    if(b & 1)
      r *= a;
    b /= 2;
    if(b == 0)
      break;
    a *= a;
  }
  return recip ? 1 / r : r;
}

#endif /* _MSC_VER */
