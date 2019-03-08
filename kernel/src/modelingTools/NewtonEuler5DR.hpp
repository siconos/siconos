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
/*! \file NewtonEuler5DR.hpp
 */

#ifndef NEWTONEULERRELATIONRFC3D_H
#define NEWTONEULERRELATIONRFC3D_H

#include "NewtonEuler1DR.hpp"
/** NewtonEuler3DR
 *
 * This class is an interface for relation with impact and RFC3D.
 *
 */

class NewtonEuler5DR : public NewtonEuler1DR
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NewtonEuler3DR);

  void RFC3DcomputeJachqTFromContacts(SP::SiconosVector q1);
  void RFC3DcomputeJachqTFromContacts(SP::SiconosVector q1, SP::SiconosVector q2);

protected:

public:
  NewtonEuler5DR(): NewtonEuler1DR() {}

  /** destructor
  */
  virtual ~NewtonEuler5DR() {};

  /** initialize components specific to derived classes.
   * \param inter the interaction using this relation
   */
  virtual void initialize(Interaction& inter);

  /* Default implementation consists in multiplying jachq and T (see NewtonEulerR::computeJachqT)
   * but here we compute the operator from the the contact point locations
   * and the local frame at contact
   *  \param inter interaction that owns the relation
   *  \param q0  the block vector to the dynamical system position
   */
  virtual void computeJachqT(Interaction& inter, SP::BlockVector q0);

  ACCEPT_STD_VISITORS();
};
#endif // NEWTONEULERRELATIONRFC3D_H
