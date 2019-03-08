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

/*! \file SphereNEDSSphereNEDSR.hpp
  \brief Two spheres relation with the Newton Euler formalism and quaternions.
*/

#ifndef SphereNEDSSphereNEDSR_h
#define SphereNEDSSphereNEDSR_h

#include "MechanicsFwd.hpp"
#include "NewtonEuler3DR.hpp"
class SphereNEDSSphereNEDSR : public NewtonEuler3DR,
  public std11::enable_shared_from_this<SphereNEDSSphereNEDSR>
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SphereNEDSSphereNEDSR);

  double r1, r2, r1pr2;

  SphereNEDSSphereNEDSR() {};

public:

  /** Constructor

  \param r1 disk1 radius
  \param r2 disk2 radius
  */
  SphereNEDSSphereNEDSR(double r1, double r2);

  double distance(double, double, double, double, double, double, double, double);

  void computeh(double time, BlockVector& q0, SiconosVector& y);

  //void computeJachq(double);

  /** visitors hook
   */
  ACCEPT_VISITORS();

};
#endif /* SphereNEDSSphereNEDSR_h */
