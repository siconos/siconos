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
#ifdef _WIN32
#define SICONOS_EXPORT extern "C" __declspec(dllexport)
#else
#define SICONOS_EXPORT extern "C"
#endif

const double m = 1; // ball mass
const double g = 9.8; // gravity

extern "C" double FextFunction(double time)
{
  double res = -0.0;
  return res;
}


SICONOS_EXPORT void ballFExt(double time, unsigned int sizeOfq, double *fExt, unsigned int sizeZ, double* z)
{
  for(unsigned int i = 0; i < sizeOfq; i++)
    fExt[i] = 0.0;
  fExt[0] = -m * g + FextFunction(time);
}

