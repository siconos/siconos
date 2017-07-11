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
/*! \file NewtonEulerR.cpp

*/

#include "SiconosVector.hpp"
#include "BlockVector.hpp"
#include "NewtonEulerJointR.hpp"

void NewtonEulerJointR::projectVectorDoF(const SiconosVector& v,
                                         const BlockVector& q0,
                                         SiconosVector& ans, int axis,
                                         bool absoluteRef)
{
  SiconosVector ax(3);
  normalDoF(ax, q0, axis, absoluteRef);

  double L = v(0)*ax(0) + v(1)*ax(1) + v(2)*ax(2);
  ans(0) = ax(0) * L;
  ans(1) = ax(1) * L;
  ans(2) = ax(2) * L;
}

SP::SiconosVector NewtonEulerJointR::projectVectorDoF(const SiconosVector& v,
                                                      const BlockVector& q0,
                                                      int axis,
                                                      bool absoluteRef)
{
  SP::SiconosVector ans(std11::make_shared<SiconosVector>(3));
  projectVectorDoF(v, q0, *ans, axis, absoluteRef);
  return ans;
}

void NewtonEulerJointR::normalDoF(SiconosVector& ans,
                                  const BlockVector& q0, int axis,
                                  bool absoluteRef)
{
  _normalDoF(ans, q0, axis, absoluteRef);
}

SP::SiconosVector NewtonEulerJointR::normalDoF(const BlockVector& q0, int axis,
                                               bool absoluteRef)
{
  SP::SiconosVector ans(std11::make_shared<SiconosVector>(3));
  _normalDoF(*ans, q0, axis, absoluteRef);
  return ans;
}
