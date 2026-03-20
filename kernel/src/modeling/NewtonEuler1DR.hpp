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
#include "SiconosVector.hpp"

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
   * on relPc1_, relPc2_. */
  siconos::algebra::SiconosVector3 contactPoint1_ = siconos::algebra::SiconosVector3::Zero();
  siconos::algebra::SiconosVector3 contactPoint2_ = siconos::algebra::SiconosVector3::Zero();

  /** Contact Points in coordinates relative to attached DS->q.  Set
   * these if contactPoint1_/contactPoint2_ are not calculated within the Newton loop. */
  siconos::algebra::SiconosVector3 relPc1_ = siconos::algebra::SiconosVector3::Zero();
  siconos::algebra::SiconosVector3 relPc2_ = siconos::algebra::SiconosVector3::Zero();

  /** Inward Normal at the contact.
   *  \todo The meaning of "Inward" has to be explained carefully.
   */
  siconos::algebra::SiconosVector3 nc_ = siconos::algebra::SiconosVector3::Zero();

  /** _Nc must be calculated relative to q2 */
  siconos::algebra::SiconosVector3 relNc_ = siconos::algebra::SiconosVector3::Zero();

  /** compute the jacobian of h w.r.t. q
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  q states vectors of the related the dynamical systems
   */
  virtual void computeH_NE_(double time, siconos::modeling::Interaction &inter,
                            const siconos::algebra::BlockVector &q0) override;

  /** Cross product matrices that correspond the lever arm from
   *  contact point to center of mass  - Internal buffer*/
  siconos::algebra::SiconosMatrix33 NPG_buffer_ = siconos::algebra::SiconosMatrix33::Zero();

  /** Matrix converting - Internal buffer*/
  siconos::algebra::SiconosMatrix33 rotationBodyToAbsoluteFrame_ =
      siconos::algebra::SiconosMatrix33::Zero();

 private:
  /** Rotation matrix converting the absolute coordinate to the contact frame^
   *  coordinate. This matrix contains the unit vector(s)of the contact frame in
   *  row. Internal bufer.
   */
  Eigen::RowVector3d rotationAbsoluteToContactFrame_ = Eigen::RowVector3d::Zero();

  void NIcomputeJachqTFromContacts(
      const Eigen::Ref<const siconos::algebra::SiconosVector7> &q1);
  void NIcomputeJachqTFromContacts(
      const Eigen::Ref<const siconos::algebra::SiconosVector7> &q1,
      const Eigen::Ref<const siconos::algebra::SiconosVector7> &q2);

 public:
  using NewtonEulerR::NewtonEulerR;

  /** V.A. boolean _isOnCOntact ?? Why is it public members ?
   *  seems parametrize the projection algorithm
   *  the projection is done on the surface  \f$ y=0 \f$  or on  \f$ y \geq 0 \f$
   */
  bool _isOnContact = false;

  /** destructor
   */
  virtual ~NewtonEuler1DR() noexcept = default;

  void initialize(Interaction &inter) override;

  /** Default implementation consists in multiplying jachq and T (see
   *  NewtonEulerR::computeJachqT) but here we compute the operator from the the
   *  contact point locations and the local frame at contact
   *
   *  \param inter interaction that owns the relation
   *  \param q0 the block vector to the dynamical system position
   */
  void computeH_NE_prod_T(const Interaction &inter,
                          const siconos::algebra::BlockVector &q0) override;

  /**
      to compute the output y = h(t,q,z) of the Relation
      with the relative contact points
      \param time current time value
      \param[in] q1 generalized coordinates vector of the fist dynamical system involved
      in the relation
      \param[in] q2 generalized coordinates vector of the second dynamical system
      involved in the relation
      \param[in,out] y the resulting vector
   */
  void computehFromRelativeContactPoints(
      double time, const Eigen::Ref<const siconos::algebra::SiconosVector7> &q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>> &q2,
      Eigen::Ref<siconos::algebra::SiconosVector> y);

  /** Return the distance between pc1 and pc, with sign according to normal */
  double distance() const;

  inline auto pc1() const {
    return siconos::algebra::ConstMapVector3Type(contactPoint1_.data(), 3);
  }

  inline const siconos::algebra::SiconosVector3 &pc2() const { return contactPoint2_; }

  inline const siconos::algebra::SiconosVector3 &nc() const { return nc_; }

  inline const siconos::algebra::SiconosVector3 &relPc1() const { return relPc1_; }
  inline const siconos::algebra::SiconosVector3 &relPc2() const { return relPc2_; }
  inline const siconos::algebra::SiconosVector3 &relNc() const { return relNc_; }

  /** Set the coordinates of first contact point in ds1 frame.
   *  It will be used to compute contactPoint1_ during computeh().
   *
   *  \param npc new coordinates
   */
  void setRelPc1(const siconos::algebra::SiconosVector3 &npc) { relPc1_ = npc; };

  /** Set the coordinates of second contact point in ds2 frame
   *  It will be used to compute contactPoint2_ during computeh().
   *
   *  \param npc new coordinates
   */
  void setRelPc2(const siconos::algebra::SiconosVector3 &npc) { relPc2_ = npc; };

  /** Set the coordinates of inside normal vector at the contact point in ds2
   *  frame. It will be used to compute _Nc during computeh().
   *
   *  \param nnc new coordinates
   */
  void setRelNc(const siconos::algebra::SiconosVector3 &nnc) { relNc_ = nnc; };
  void display() const override {}
  virtual void accept(relations::Visitor &tourist) const override { tourist.visit(*this); }
};
}  // namespace siconos::modeling
#endif
