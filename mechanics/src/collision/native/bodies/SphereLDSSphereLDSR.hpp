/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

/*! \file SphereLDSSphereLDSR.hpp
  \brief Two spheres relation - Inherits from LagrangianScleronomousR
*/

#ifndef SphereLDSSphereLDSR_h
#define SphereLDSSphereLDSR_h

#include "MechanicsFwd.hpp"
#include "LagrangianScleronomousR.hpp"

class SphereLDSSphereLDSR : public LagrangianScleronomousR, public std::enable_shared_from_this<SphereLDSSphereLDSR>
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SphereLDSSphereLDSR);

  double r1, r2, r1pr2;

  SphereLDSSphereLDSR() {};

public:

  /** Constructor

  \param r1 disk1 radius
  \param r2 disk2 radius
  */
  SphereLDSSphereLDSR(double r1, double r2);

  double distance(double, double, double, double, double, double, double, double);

  using LagrangianScleronomousR::computeh;
  /** to compute the output y = h(t,q,z) of the Relation
      \param q coordinates of the dynamical systems involved in the relation
      \param z user defined parameters (optional)
      \param y the resulting vector
  */
  void computeh(const BlockVector& q, BlockVector& z, SiconosVector& y);

  /** to compute the jacobian of h(...). Set attribute _jachq (access: jacqhq())
      \param q coordinates of the dynamical systems involved in the relation
      \param z user defined parameters (optional)
  */
  void computeJachq(const BlockVector& q, BlockVector& z);

  /** visitors hook
   */
  ACCEPT_VISITORS();

};
#endif /* SphereLDSSphereLDSR_h */
