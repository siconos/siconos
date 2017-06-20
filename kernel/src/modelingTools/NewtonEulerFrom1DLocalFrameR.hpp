/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
/*! \file NewtonEulerFrom1DLocalFrameR.hpp

 */
#ifndef NEWTONEULERIMPACT_H
#define NEWTONEULERIMPACT_H

#include "NewtonEulerR.hpp"
#include "NewtonEulerDS.hpp"


/** NewtonEulerFrom1DLocalFrameR
 *
 * \author O. Bonnefon
 *  \version 3.0.0.
 *  \date Dec, 2010
 *
 * This class is an interface for a relation with impact.  It
 * implements the computation of the jacoboian of h from the points of
 * contacts and the normal.  Use this class consists in overloading
 * the method computeh, by setting the member pc1, pc2, nc and y.  The
 * matrix jachq is used both for the building of the OSNSP (with T)
 * and for the predictor of activation of deactivation of the Interaction.
 *
 */


class NewtonEulerFrom1DLocalFrameR : public NewtonEulerR
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NewtonEulerFrom1DLocalFrameR);

  /* Current Contact Points, may be updated within Newton loop based
   * on _relPc1, _relPc2. */
  SP::SiconosVector _Pc1;
  SP::SiconosVector _Pc2;

  /* Contact Points in coordinates relative to attached DS->q.  Set
   * these if _Pc1/_Pc2 are not calculated within the Newton loop. */
  SP::SiconosVector _relPc1;
  SP::SiconosVector _relPc2;

  /* Inward Normal at the contact.
   * \todo The meaning of "Inward" has to be explained carefully.
   */
  SP::SiconosVector _Nc;


  virtual void initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);

  /* Rotation matrix converting the absolute coordinate to the contact frame coordinate.
   * This matrix contains the unit vector(s)of the contact frame in row.
   */
  SP::SimpleMatrix _RotationAbsToContactFrame;

  /* Matrix converting */
  SP::SimpleMatrix _rotationMatrixAbsToBody;

  /* Cross product matrices that correspond the lever arm from
   * contact point to center of mass*/
  SP::SimpleMatrix _NPG1;
  SP::SimpleMatrix _NPG2;


  /*buffer matrices*/
  SP::SimpleMatrix _AUX1;
  SP::SimpleMatrix _AUX2;
private:
  void NIcomputeJachqTFromContacts(SP::SiconosVector q1);
  void NIcomputeJachqTFromContacts(SP::SiconosVector q1, SP::SiconosVector q2);
public:

  /** V.A. boolean _isOnCOntact ?? Why is it public members ?
  *  seems parametrize the projection algorithm
  *  the projection is done on the surface \f$y=0\f$ or on \f$y \geq 0\f$
  */
  bool _isOnContact;

  /** constructorx
  */
  NewtonEulerFrom1DLocalFrameR():
    NewtonEulerR(), _Pc1(new SiconosVector(3)), _Pc2(new SiconosVector(3)),
    _relPc1(new SiconosVector(3)), _relPc2(new SiconosVector(3)),
    _Nc(new SiconosVector(3))
  {
    /*_ds1=NULL;_ds2=NULL;*/
  }

  /** destructor
  */
  virtual ~NewtonEulerFrom1DLocalFrameR() {};

  virtual void computeJachq(double time, Interaction& inter, SP::BlockVector q0);

  /* Default implementation consists in multiplying jachq and T (see NewtonEulerR::computeJachqT)
   * but here we compute the operator from the the contact point locations
   * and the local frame at contact
   *  \param inter interaction that owns the relation
   *  \param q0  the block vector to the dynamical system position
   */
  virtual void computeJachqT(Interaction& inter, SP::BlockVector q0);

  /* Default implementation of computeh updates contact points and
   * distance for q if different than qold. */
  virtual void computeh(double time, BlockVector& q0, SiconosVector &y);

  inline SP::SiconosVector pc1() const
  {
    return _Pc1;
  }
  inline SP::SiconosVector pc2() const
  {
    return _Pc2;
  }
  inline SP::SiconosVector nc() const
  {
    return _Nc;
  }

  /** set the coordinates of first contact point
  * \param npc new coordinates
  */
  void setpc1(SP::SiconosVector npc)
  {
    _Pc1 = npc;
  };

  /** set the coordinates of second contact point
  * \param npc new coordinates
  */
  void setpc2(SP::SiconosVector npc)
  {
    _Pc2 = npc;
  };

  /** set the coordinates of inside normal vector at the contact point
  * \param nnc new coordinates
  */
  void setnc(SP::SiconosVector nnc)
  {
    _Nc = nnc;
  };

  // visitors hook
  ACCEPT_STD_VISITORS();

};
#endif // NEWTONEULERRIMPACT_H
