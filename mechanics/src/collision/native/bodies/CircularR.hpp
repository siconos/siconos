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

/** \file CircleCircleR.hpp
 */
#ifndef CircularR_h
#define CircularR_h

#include "MechanicsFwd.hpp"
#include "Interaction.hpp"
#include "LagrangianScleronomousR.hpp"

/**  \class CircularR
 *   \brief Two circle relation - Inherits from LagrangianScleronomousR
 */

class CircularR : public LagrangianScleronomousR
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(CircularR);

  double _r1, _r2;

  CircularR() : LagrangianScleronomousR() {};

public:

  /** Constructor

  \param disk1 radius
  \param disk2 radius
  */
  CircularR(double r1, double r2): _r1(r1), _r2(r2) {};

  double getRadius1()
  {
    return _r1;
  };

  double getRadius2()
  {
    return _r2;
  };

  virtual double distance(double, double, double, double, double, double)
  {
    assert(0);
    return(0);
  };

  ~CircularR() {};

};
#endif /* CircularR_h */
