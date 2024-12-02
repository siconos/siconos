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
/*! \file NewtonEuler1DR.hpp

 */
#ifndef NEWTONEULERIMPACT_H
#define NEWTONEULERIMPACT_H

#include "NewtonEulerR.hpp"

namespace siconos::modeling {
/**
   This class is an interface for a relation with impact.  It
   implements the computation of the jacoboian of h from the points of
   contacts and the normal.  Use this class consists in overloading
   the method computeh, by setting the member pc1, pc2, nc and y.  The
   matrix jachq is used both for the building of the OSNSP (with T)
   and for the predictor of activation of deactivation of the Interaction.

*/
class NewtonEuler1DR : public NewtonEulerR {
 protected:
  ACCEPT_SERIALIZATION(NewtonEuler1DR);

  /** Current Contact Points, may be updated within Newton loop based
   * on _relPc1, _relPc2. */
  std::shared_ptr<siconos::algebra::SiconosVector> _Pc1{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _Pc2{nullptr};

  /** Contact Points in coordinates relative to attached DS->q.  Set
   * these if _Pc1/_Pc2 are not calculated within the Newton loop. */
  std::shared_ptr<siconos::algebra::SiconosVector> _relPc1{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _relPc2{nullptr};

  /** Inward Normal at the contact.
   *  \todo The meaning of "Inward" has to be explained carefully.
   */
  std::shared_ptr<siconos::algebra::SiconosVector> _Nc{nullptr};

  /** _Nc must be calculated relative to q2 */
  std::shared_ptr<siconos::algebra::SiconosVector> _relNc{nullptr};

  /** Rotation matrix converting the absolute coordinate to the contact frame
   *  coordinate. This matrix contains the unit vector(s)of the contact frame in
   *  row.
   */
  std::shared_ptr<siconos::algebra::SimpleMatrix> _rotationAbsoluteToContactFrame{nullptr};

  /** Matrix converting */
  std::shared_ptr<siconos::algebra::SimpleMatrix> _rotationBodyToAbsoluteFrame{nullptr};

  /** Cross product matrices that correspond the lever arm from
   *  contact point to center of mass*/
  std::shared_ptr<siconos::algebra::SimpleMatrix> _NPG1{nullptr};
  std::shared_ptr<siconos::algebra::SimpleMatrix> _NPG2{nullptr};

  /*buffer matrices*/
  std::shared_ptr<siconos::algebra::SimpleMatrix> _AUX1{nullptr};
  std::shared_ptr<siconos::algebra::SimpleMatrix> _AUX2{nullptr};

  /** Set the coordinates of first contact point.  Must only be done
   *  in a computeh() override.
   *
   *  \param npc new coordinates
   */
  void setpc1(std::shared_ptr<siconos::algebra::SiconosVector> npc) { _Pc1 = npc; };

  /** Set the coordinates of second contact point.  Must only be done
   *  in a computeh() override.
   *
   *  \param npc new coordinates
   */
  void setpc2(std::shared_ptr<siconos::algebra::SiconosVector> npc) { _Pc2 = npc; };

  /** Set the coordinates of inside normal vector at the contact point.
   *  Must only be done in a computeh() override.
   *
   *  \param nnc new coordinates
   */
  void setnc(std::shared_ptr<siconos::algebra::SiconosVector> nnc) { _Nc = nnc; };

 private:
  void NIcomputeJachqTFromContacts(std::shared_ptr<siconos::algebra::SiconosVector> q1);
  void NIcomputeJachqTFromContacts(std::shared_ptr<siconos::algebra::SiconosVector> q1,
                                   std::shared_ptr<siconos::algebra::SiconosVector> q2);

 public:
  /** V.A. boolean _isOnCOntact ?? Why is it public members ?
   *  seems parametrize the projection algorithm
   *  the projection is done on the surface  \f$ y=0 \f$  or on  \f$ y \geq 0 \f$
   */
  bool _isOnContact = false;

  /** constructor */
  NewtonEuler1DR();

  /** destructor
   */
  virtual ~NewtonEuler1DR() noexcept = default;

  void computeJachq(double time, Interaction &inter,
                    std::shared_ptr<siconos::algebra::BlockVector> q0) override;

  void initialize(Interaction &inter) override;

  /** Default implementation consists in multiplying jachq and T (see
   *  NewtonEulerR::computeJachqT) but here we compute the operator from the the
   *  contact point locations and the local frame at contact
   *
   *  \param inter interaction that owns the relation
   *  \param q0 the block vector to the dynamical system position
   */
  void computeJachqT(Interaction &inter,
                     std::shared_ptr<siconos::algebra::BlockVector> q0) override;

  /**
      to compute the output y = h(t,q,z) of the Relation
      with the relative contact points
      \param time current time value
      \param q coordinates of the dynamical systems involved in the relation
      \param y the resulting vector
   */
  void computehFromRelativeContactPoints(double time, const siconos::algebra::BlockVector &q0,
                                         siconos::algebra::SiconosVector &y);

  /** Return the distance between pc1 and pc, with sign according to normal */
  double distance() const;

  inline std::shared_ptr<siconos::algebra::SiconosVector> pc1() const { return _Pc1; }
  inline std::shared_ptr<siconos::algebra::SiconosVector> pc2() const { return _Pc2; }
  inline std::shared_ptr<siconos::algebra::SiconosVector> nc() const { return _Nc; }

  inline std::shared_ptr<siconos::algebra::SiconosVector> relPc1() const { return _relPc1; }
  inline std::shared_ptr<siconos::algebra::SiconosVector> relPc2() const { return _relPc2; }
  inline std::shared_ptr<siconos::algebra::SiconosVector> relNc() const { return _relNc; }

  /** Set the coordinates of first contact point in ds1 frame.
   *  It will be used to compute _Pc1 during computeh().
   *
   *  \param npc new coordinates
   */
  void setRelPc1(std::shared_ptr<siconos::algebra::SiconosVector> npc) { _relPc1 = npc; };

  /** Set the coordinates of second contact point in ds2 frame
   *  It will be used to compute _Pc2 during computeh().
   *
   *  \param npc new coordinates
   */
  void setRelPc2(std::shared_ptr<siconos::algebra::SiconosVector> npc) { _relPc2 = npc; };

  /** Set the coordinates of inside normal vector at the contact point in ds2
   *  frame. It will be used to compute _Nc during computeh().
   *
   *  \param nnc new coordinates
   */
  void setRelNc(std::shared_ptr<siconos::algebra::SiconosVector> nnc) { _relNc = nnc; };
  void display() const override {}
};
}  // namespace siconos::modeling
#endif
