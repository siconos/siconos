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
/*! \file CouplerJointR.hpp
*/
#ifndef CouplerJointRELATION_H
#define CouplerJointRELATION_H

#include <MechanicsFwd.hpp>
#include <SiconosFwd.hpp>
#include <NewtonEulerR.hpp>
#include <NewtonEulerJointR.hpp>
#include <Tools.hpp>

/** \class CouplerJointR
 *  \brief This class implements a coupling (equality) between two
 *  DoFs of any NewtonEulerJointR.  Can be used e.g. to implement a
 *  screw relation (coupled rotation and translation) based on
 *  CylindricalJointR.
 */
class CouplerJointR : public NewtonEulerJointR
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(CouplerJointR);

  SP::NewtonEulerJointR _joint;

  unsigned int _dof1, _dof2;
  double _ratio, _offset;

  /** Return the normal of the linear DoF axis.  \param axis must be 0 */
  virtual void _normalDoF(SiconosVector& ans, const BlockVector& q0, int axis,
                          bool absoluteRef=true);

public:

  /** Empty constructor. The relation may be initialized later by
   * setPoint, setAbsolute, and setBasePositions. */
  CouplerJointR();

  /** Initialize a joint friction for a common case: a single axis with a
   * single friction, either positive or negative. For use with
   * NewtonImpactNSL. */
  CouplerJointR(SP::NewtonEulerJointR joint, unsigned int dof1,
                unsigned int dof2, double ratio);

  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);

  virtual void computeJachq(double time, Interaction& inter, SP::BlockVector q0);

  /* Return the joint DoF index assigned to the first DoF. */
  unsigned int dof1() { return _dof1; }

  /* Return the joint DoF index assigned to the second DoF. */
  unsigned int dof2() { return _dof2; }

  /* Return the joint assigned to this friction relation. */
  SP::NewtonEulerJointR joint() { return _joint; }

  /* From provided position vectors, calculate initial offsets to be
   * maintained in the relation. */
  void setBasePositions(SP::SiconosVector q1,
                        SP::SiconosVector q2=SP::SiconosVector());

  /** Get the number of constraints defined in the joint
      \return the number of constraints
   */
  virtual unsigned int numberOfConstraints() { return 1; }

  /** Get the number of degrees of freedom defined in the joint
      \return the number of degrees of freedom (DoF)
   */
  virtual unsigned int numberOfDoF() { return 1; }

  /** Return the type of a degree of freedom of this joint.
      \return the type of the degree of freedom (DoF)
  */
  virtual DoF_Type typeOfDoF(unsigned int axis) {
    if (axis<1) return DOF_TYPE_LINEAR;
    else return DOF_TYPE_INVALID;
  };

  /** Compute the vector of linear and angular positions of the free axes */
  virtual void computehDoF(double time, BlockVector& q0, SiconosVector& y,
                           unsigned int axis);

  /** Compute the jacobian of linear and angular DoF with respect to some q */
  virtual void computeJachqDoF(double time, Interaction& inter,
                               SP::BlockVector q0, SimpleMatrix& jachq,
                               unsigned int axis);

  /** destructor
   */
  virtual ~CouplerJointR() {};
};
#endif  //CouplerJointRELATION_H
