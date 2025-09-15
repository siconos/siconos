/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
/*! \file Lagrangian2d#DR.hpp

 */
#ifndef Lagrangian2d3DR_H
#define Lagrangian2d3DR_H

#include "LagrangianDS.hpp"
#include "LagrangianScleronomousR.hpp"

using namespace RELATION;
/** Lagrangian2d3DR
 *
 * This class is an interface for a relation with impact.  It
 * implements the computation of the jacoboian of h from the points of
 * contacts and the normal.  Use this class consists in overloading
 * the method computeh, by setting the member pc1, pc2, nc and y.  The
 * matrix jachq is used both for the building of the OSNSP (with T)
 * and for the predictor of activation of deactivation of the Interaction.
 *
 */

class Lagrangian2d3DR : public LagrangianScleronomousR {
protected:

  ACCEPT_SERIALIZATION(Lagrangian2d3DR);

  /* Current Contact Points, may be updated within Newton loop based
   * on _relPc1, _relPc2. */
  SP::SiconosVector _Pc1;
  SP::SiconosVector _Pc2;

  /* Inward Normal at the contact.
   * \todo The meaning of "Inward" has to be explained carefully.
   */
  SP::SiconosVector _Nc;

  /* _Nc must be calculated relative to q2 */
  SP::SiconosVector _relNc;

  /* Rotation matrix converting the absolute coordinate to the contact frame
   * coordinate. This matrix contains the unit vector(s)of the contact frame in
   * row.
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

  /** Set the coordinates of first contact point.  Must only be done
   *  in a computeh() override.
   *
   *  \param npc new coordinates
   */
  void setpc1(SP::SiconosVector npc) { _Pc1 = npc; };

  /** Set the coordinates of second contact point.  Must only be done
   *  in a computeh() override.
   *
   *  \param npc new coordinates
   */
  void setpc2(SP::SiconosVector npc) { _Pc2 = npc; };

  /** Set the coordinates of inside normal vector at the contact point.
   *  Must only be done in a computeh() override.
   *
   *  \param nnc new coordinates
   */
  void setnc(SP::SiconosVector nnc) { _Nc = nnc; };

public:
  /** V.A. boolean _isOnCOntact ?? Why is it public members ?
   *  seems parametrize the projection algorithm
   *  the projection is done on the surface \f$ y=0 \f$ or on \f$ y \geq 0 \f$
   */
  bool _isOnContact = false;

  /** constructor
   */
  Lagrangian2d3DR()
      : LagrangianScleronomousR(), _Pc1(new SiconosVector(2)),
        _Pc2(new SiconosVector(2)), _Nc(new SiconosVector(2))
  {
    /*_ds1=nullptr;_ds2=nullptr;*/
  }

  /** destructor
   */
  virtual ~Lagrangian2d3DR() noexcept {};

  void initialize(Interaction &inter) override;

  /**
     to compute the output y = h(q,z) of the Relation
     
     \param q coordinates of the dynamical systems involved in the relation
     \param z user defined parameters (optional)
     \param y the resulting vector
  */
  void computeh(const BlockVector &q, BlockVector &z,
                SiconosVector &y) override;

  /**
     to compute the jacobian of h(...). Set attribute _jachq (access: jacqhq())
     
     \param q coordinates of the dynamical systems involved in the relation
     \param z user defined parameters (optional)
  */
  void computeJachq(const BlockVector &q, BlockVector &z) override;

  /** Return the distance between pc1 and pc, with sign according to normal */
  double distance() const;

  inline SP::SiconosVector pc1() const { return _Pc1; }
  inline SP::SiconosVector pc2() const { return _Pc2; }
  inline SP::SiconosVector nc() const { return _Nc; }

  inline SP::SiconosVector relNc() const { return _relNc; }

  /** Set the coordinates of inside normal vector at the contact point in ds2
   *  frame. It will be used to compute _Nc during computeh().
   *
   *  \param nnc new coordinates
   */
  void setRelNc(SP::SiconosVector nnc) { _relNc = nnc; };
  void display() const override;

  ACCEPT_STD_VISITORS();
};
TYPEDEF_SPTR(Lagrangian2d3DR)
#endif // NEWTONEULERRIMPACT_H
