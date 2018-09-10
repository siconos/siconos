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

/** \file CircleCircleR.hpp
 *
 */

#ifndef CircleCircleR_h
#define CircleCircleR_h

#include "MechanicsFwd.hpp"
#include "CircularR.hpp"

/** \class CircleCircleR
 *  \brief Two disks relation - Inherits from LagrangianScleronomousR
 */
class CircleCircleR : public CircularR, public std11::enable_shared_from_this<CircleCircleR>
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(CircleCircleR);

  CircleCircleR() {};

public:

  /** Constructor
  \param rdisk1 radius
  \param rdisk2 radius
  */
  CircleCircleR(double rdisk1, double rdisk2);

  /** compute distance between 2 disks
      \param x1 x position of first disk
      \param y1 y position of first disk
      \param r1 radius of first disk
      \param x2 x position of second disk
      \param y2 y position of second disk
      \param r2 radius of second disk
      \return distance 
  */
  double distance(double x1, double y1, double r1,
                  double x2, double y2, double r2);

  using LagrangianScleronomousR::computeh;

  /** computeh implementation, see LagrangianScleronomousR
      \param q
      \param z
      \param y
  */
  void computeh(SiconosVector& q, SiconosVector& z, SiconosVector& y);

  /** computeh implementation, see LagrangianScleronomousR
      \param q
      \param z
  */
  void computeJachq(SiconosVector& q, SiconosVector& z);

  /** visitors hook
   */
  ACCEPT_VISITORS();

  ~CircleCircleR() {};

};
#endif /* CircleCircleR_h */
