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

/*! \file DiskDiskR.hpp
 */
#ifndef DiskDiskR_h
#define DiskDiskR_h

#include "MechanicsFwd.hpp"
#include <SiconosVisitor.hpp>
#include "CircularR.hpp"


/** \class DiskDiskR
 *  \brief Two disks relation - Inherits from LagrangianScleronomousR
 */
class DiskDiskR : public CircularR, public std11::enable_shared_from_this<DiskDiskR>
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(DiskDiskR);

  double r1pr2;

  DiskDiskR() : CircularR() {};

public:

  /** Constructor
  \param disk1 radius
  \param disk2 radius
  */
  DiskDiskR(double disk1, double disk2);

  double distance(double, double, double, double, double, double);

  void computeh(SiconosVector& q, SiconosVector& z, SiconosVector& y);

  void computeJachq(SiconosVector& q, SiconosVector& z);

  /** visitors hook
   */
  ACCEPT_VISITORS();

  ~DiskDiskR() {};

};
#endif /* DiskDiskR_h */
