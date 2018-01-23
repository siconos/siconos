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

  SP::NewtonEulerJointR _joint1, _joint2;
  SP::SiconosVector _ref1, _ref2;

  unsigned int _dof1, _dof2;
  unsigned int _ref1_index, _ref2_index;
  double _ratio, _offset;

  /** Return the normal of the linear DoF axis.  \param axis must be 0 */
  virtual void _normalDoF(SiconosVector& ans, const BlockVector& q0, int axis,
                          bool absoluteRef=true);

  /** An internal helper function to assign reference vectors during
   * computeh and computeJachq. */
  void makeBlockVectors(SP::SiconosVector q1, SP::SiconosVector q2,
                        BlockVector& q01, BlockVector& q02);

public:

  /** Empty constructor. The relation may be initialized later by
   * setPoint, setAbsolute, and setBasePositions. */
  CouplerJointR();

  /** Initialize a coupler. See setReferences() for an explanation of
   * the parameters. */
  CouplerJointR(SP::NewtonEulerJointR joint1, unsigned int dof1,
                SP::NewtonEulerJointR joint2, unsigned int dof2,
                double ratio,
                SP::SiconosVector ref1=SP::SiconosVector(), unsigned int ref1_index=0,
                SP::SiconosVector ref2=SP::SiconosVector(), unsigned int ref2_index=0);

  /** Initialize a coupler. See setReferences() for an explanation of
   * the parameters. */
  CouplerJointR(SP::NewtonEulerJointR joint1, unsigned int dof1,
                SP::NewtonEulerJointR joint2, unsigned int dof2,
                double ratio,
                SP::NewtonEulerDS refds1, unsigned int ref1_index,
                SP::NewtonEulerDS refds2, unsigned int ref2_index);

  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);

  virtual void computeJachq(double time, Interaction& inter, SP::BlockVector q0);

  /* Return the joint DoF index assigned to the first DoF. */
  unsigned int dof1() { return _dof1; }

  /* Return the joint DoF index assigned to the second DoF. */
  unsigned int dof2() { return _dof2; }

  /* Return the first reference joint assigned to this friction relation. */
  SP::NewtonEulerJointR joint1() { return _joint1; }

  /* Return the second reference joint assigned to this friction relation. */
  SP::NewtonEulerJointR joint2() { return _joint2; }

  /** Set reference joints and vectors.  This constraint maintains the
   *  equality theta2=theta1*ratio; theta1 is measured by joint1 with
   *  reference to some vector ref1 which must replace either the
   *  first or second DS of the current relation being defined.  If
   *  ref1 and ref2 are null, then no reference is used.  Typically
   *  ref1 and ref2 will be equal in order to implement gear ratios
   *  for example.  joint1 must be between ref1 and ds1 specified in
   *  setBasePositions(), while joint2 must be between ref2 and ds2.
   *
   *  \param joint1 The joint for the first reference measurement theta1.
   *  \param dof1 The degree of freedom index of joint1 to use for
   *              the first reference measurement theta1.
   *  \param ref1 The optional reference position for the first
   *              reference measurement theta1.
   *  \param ref1_index Must be 0 or 1, depending on where ref1
   *                    appears in joint1.
   *  \param joint2 The joint for the second reference measurement theta2.
   *  \param dof2 The degree of freedom index of joint2 to use for
   *              the second reference measurement theta2.
   *  \param ref2 The optional reference position for the second
   *              reference measurement theta2.
   *  \param ref2_index Must be 0 or 1, depending on where ref2
   *                    appears in joint2.
   */
  void setReferences(SP::NewtonEulerJointR joint1, unsigned int dof1,
                     SP::NewtonEulerJointR joint2, unsigned int dof2,
                     SP::SiconosVector ref1, unsigned int ref1_index,
                     SP::SiconosVector ref2, unsigned int ref2_index);

  /** Set reference joints and vectors.  This constraint maintains the
   *  equality theta2=theta1*ratio; theta1 is measured by joint1 with
   *  reference to some vector ref1 which must replace either the
   *  first or second DS of the current relation being defined.  If
   *  ref1 and ref2 are null, then no reference is used.  Typically
   *  ref1 and ref2 will be equal in order to implement gear ratios
   *  for example.  joint1 must be between ref1 and ds1 specified in
   *  setBasePositions(), while joint2 must be between ref2 and ds2.
   *  This alternative setReferences() takes NewtonEulerDS as
   *  parameters, but the reference vectors will be taken as
   *  refds1->q() and refds2->q() respectively.
   *
   *  \param joint1 The joint for the first reference measurement theta1.
   *  \param dof1 The degree of freedom index of joint1 to use for
   *              the first reference measurement theta1.
   *  \param ref1 The optional reference vector for the first
   *                reference measurement theta1.
   *  \param ref1_index Must be 0 or 1, depending on where ref1
   *                    appears in joint1.
   *  \param joint2 The joint for the second reference measurement theta2.
   *  \param dof2 The degree of freedom index of joint2 to use for
   *              the second reference measurement theta2.
   *  \param ref2 The optional reference vector for the second
   *                reference measurement theta2.
   *  \param ref2_index Must be 0 or 1, depending on where ref2
   *                    appears in joint2.
   */
  void setReferences(SP::NewtonEulerJointR joint1, unsigned int dof1,
                     SP::NewtonEulerJointR joint2, unsigned int dof2,
                     SP::NewtonEulerDS refds1, unsigned int ref1_index,
                     SP::NewtonEulerDS refds2, unsigned int ref2_index);

  /* Set the gear ratio. */
  void setRatio(double ratio);

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
