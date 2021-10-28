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
/*! \file NewtonEuler3DR.hpp

 */
#ifndef NEWTONEULERRELATIONFC3D_H
#define NEWTONEULERRELATIONFC3D_H

#include "NewtonEuler1DR.hpp"
/** NewtonEuler3DR
 *
 * This class is an interface for relation with impact and FC3D.
 * From NewtonEuler1DR, it inherits to the computation of the jacobian, this operator is use for the predictor of activation and deactivation of the Interaction.
 * The OSNSP is build using the matrix jachqT, that is computed from the point if contact pc1, pc2 and Nc.
 * Use this class consists in overload the method computeh, and children class has to set the menber pc1, pc2 and nc.
 *
 *
 */

class NewtonEuler3DR : public NewtonEuler1DR
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NewtonEuler3DR);

  void FC3DcomputeJachqTFromContacts(SP::SiconosVector q1);
  void FC3DcomputeJachqTFromContacts(SP::SiconosVector q1, SP::SiconosVector q2);

protected:

public:
  NewtonEuler3DR(): NewtonEuler1DR() {}

  /** destructor
  */
  virtual ~NewtonEuler3DR() {};

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
#endif // NEWTONEULERRELATIONFC3D_H
