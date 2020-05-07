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

#include <cmath>
#include "SphereNEDS.hpp"

SphereNEDS::SphereNEDS(double r, double m, SP::SiconosMatrix I,
                       SP::SiconosVector qinit,
                       SP::SiconosVector vinit)
  : NewtonEulerDS(qinit, vinit, m, I), radius(r)
{

  // note : _ndof = 3 in NewtonEuleurDS ? (=> _ndof = 6 ?)

  assert(qinit->size() == _qDim);
  assert(vinit->size() == 6);  // == _ndof

}

SphereNEDS::~SphereNEDS()
{}
