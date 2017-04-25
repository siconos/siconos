/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2017 INRIA.
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
/*! \file NewtonEulerJointR.hpp

*/
#ifndef NewtonEulerJointRELATION_H
#define NewtonEulerJointRELATION_H

#include <MechanicsFwd.hpp>
#include <SiconosFwd.hpp>
#include <NewtonEulerR.hpp>

/** \class NewtonEulerJointR
 *  \brief This class implements an abstract Joint relation (articulation) between one or two Newton/Euler dynamical systems.
 */
class NewtonEulerJointR : public NewtonEulerR
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(NewtonEulerJointR);
  NewtonEulerJointR(): NewtonEulerR(), _allowSelfCollide(false) {};

  /** A flag determining whether this joint should block
   * "self-collision", i.e., if true, bodies connected by this joint
   * will not enter into unilateral contact. */
  bool _allowSelfCollide;

public:

  /** Return the value of the _allowSelfCollide flag. */
  bool allowSelfCollide() { return _allowSelfCollide; }

  /** Set the value of the _allowSelfCollide flag. */
  void setAllowSelfCollide(bool x) { _allowSelfCollide = x; }

  /** destructor
   */
  virtual ~NewtonEulerJointR() {};
};
#endif  //NewtonEulerJointRELATION_H
