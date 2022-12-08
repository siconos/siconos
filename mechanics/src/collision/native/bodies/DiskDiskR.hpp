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
class DiskDiskR : public CircularR, public std::enable_shared_from_this<DiskDiskR>
{
private:

  ACCEPT_SERIALIZATION(DiskDiskR);

  double r1pr2;

public:

  /** Constructor
   *
   *  \param disk1 radius
   *  \param disk2 radius
   */
  DiskDiskR(double disk1, double disk2);

  ~DiskDiskR() noexcept = default;

  double distance(double, double, double, double, double, double);

  /**
     to compute the output y = h(q,z) of the Relation
     
     \param q coordinates of the dynamical systems involved in the relation
     \param z user defined parameters (optional)
     \param y the resulting vector
  */
  void computeh(const BlockVector& q, BlockVector& z, SiconosVector& y);

  /**
     to compute the jacobian of h(...). Set attribute _jachq (access: jacqhq())
     
     \param q coordinates of the dynamical systems involved in the relation
     \param z user defined parameters (optional)
  */
  void computeJachq(const BlockVector& q, BlockVector& z);

  ACCEPT_VISITORS();

};
#endif /* DiskDiskR_h */
