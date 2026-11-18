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
 *  \brief Contact relation between Newton-Euler bodies: normal relative velocity.
 */
#ifndef NEWTONEULER1DR_HPP
#define NEWTONEULER1DR_HPP

#include "NewtonEulerR.hpp"
#include "SiconosVector.hpp"

namespace siconos::modeling {
/** Unilateral contact between body 1 and body 2 (or body 1 and a fixed obstacle):
  normal component of the relative velocity at the contact point.

  This class is also the base of the contact relations NewtonEuler3DR (with friction)
  and NewtonEulerVelocityAngularVelocityR (with friction and rolling resistance), which
  share the geometry and the notation below. See NewtonEulerR for the notation of the
  bodies (\f$ q_i, v_i, T \f$).

  <b>Contact geometry.</b> The contact points \f$ P_1 \f$ (on body 1), \f$ P_2 \f$ (on
  body 2) and the unit normal \f$ n \f$ are either set directly by a derived class
  (contactPoint1_, contactPoint2_, nc_), or computed by computehFromRelativeContactPoints
  from data attached to the bodies (setRelPc1, setRelPc2, setRelNc):
  \f[
    P_1 = x_{G_1} + R_1 P_1^{\mathrm{loc}}, \qquad
    P_2 = x_{G_2} + R_2 P_2^{\mathrm{loc}}, \qquad
    n = R_2\, n^{\mathrm{loc}},
  \f]
  where \f$ R_i = R(p_i) \f$ rotates body-frame vectors into the inertial frame. With a
  single body, \f$ P_2 \f$ and \f$ n \f$ are given in the inertial frame. The normal
  points from body 2 towards body 1, and the gap (distance()) is
  \f[
    g = \operatorname{sign}\bigl( n^\top (P_1 - P_2) \bigr)\, \| P_1 - P_2 \|,
  \f]
  positive when the bodies are separated.

  <b>Relative velocity.</b> With \f$ \omega_i = R_i \Omega_i \f$ the angular velocities in
  the inertial frame and \f$ P = P_1 \f$, the relative velocity of body 1 with respect to
  body 2 at \f$ P \f$ is
  \f[
    u = v_{G_1} + \omega_1 \times (P - x_{G_1})
      - v_{G_2} - \omega_2 \times (P - x_{G_2}).
  \f]
  Let \f$ [a]_\times \f$ be the skew-symmetric matrix such that
  \f$ [a]_\times b = a \times b \f$, so that \f$ \omega \times (P - x_G) = [x_G - P]_\times\,
  \omega \f$, and let \f$ R_c = [\, n \;\; t_1 \;\; t_2 \,]^\top \f$ be the rotation from the
  inertial frame to the contact frame (\f$ t_1, t_2 \f$ from orthoBaseFromVector). Then
  \f[
    R_c\, u = \underbrace{R_c \begin{bmatrix}
        I_3 & [x_{G_1} - P]_\times R_1 & -I_3 & -[x_{G_2} - P]_\times R_2
      \end{bmatrix}}_{H_c(q)} v
  \f]
  (only the first two block columns with a single body).

  <b>This relation</b> (\f$ m = 1 \f$):
  \f[
    y = g, \qquad \dot y = n^\top u = (H_c)_{1,:}\, v,
  \f]
  i.e. \f$ H \f$ (H_NE_prod_T) is the first row of \f$ H_c \f$, computed directly from the
  geometry (computeH_NE_prod_T) rather than as \f$ J_h T \f$.

  <b>Position Jacobian.</b> \f$ J_h \f$ (H_NE, used for the projection onto the constraints)
  is computed by computeH_NE_ as
  \f$ J_h = \nabla_q^\top \bigl[ n^\top (P_1(q_1) - P_2(q_2)) \bigr] \f$ with \f$ n \f$
  held fixed. This is the exact gradient of \f$ g \f$ when \f$ P_1 - P_2 \f$ is parallel to
  \f$ n \f$ (which is the case for points given by a contact detection), and then
  \f$ J_h T v = n^\top u \f$. For the classes derived from NewtonEuler1DR, only the first
  row of \f$ J_h \f$ is filled.

  <b>Output.</b> computehFromRelativeContactPoints only updates the contact points and the
  normal; the value of \f$ y \f$ itself (e.g. \f$ y_1 = g \f$) must be set by the derived
  class, typically in computeh (see e.g. ContactR in mechanics).

  \see DevNotes, chapter "NewtonEulerR: computation of \f$ \nabla_q H \f$".
*/
class NewtonEuler1DR : public NewtonEulerR {
 protected:
  ACCEPT_SERIALIZATION(NewtonEuler1DR);

  /** Contact point \f$ P_1 \f$ on body 1, in the inertial frame. Set by the derived class,
   *  or computed from relPc1_ by computehFromRelativeContactPoints. */
  siconos::algebra::SiconosVector3 contactPoint1_ = siconos::algebra::SiconosVector3::Zero();

  /** Contact point \f$ P_2 \f$ on body 2 (or on the obstacle), in the inertial frame. */
  siconos::algebra::SiconosVector3 contactPoint2_ = siconos::algebra::SiconosVector3::Zero();

  /** \f$ P_1^{\mathrm{loc}} \f$: contact point 1 in the frame of body 1. */
  siconos::algebra::SiconosVector3 relPc1_ = siconos::algebra::SiconosVector3::Zero();

  /** \f$ P_2^{\mathrm{loc}} \f$: contact point 2 in the frame of body 2 (in the inertial
   *  frame if there is only one body). */
  siconos::algebra::SiconosVector3 relPc2_ = siconos::algebra::SiconosVector3::Zero();

  /** Unit normal \f$ n \f$ at the contact, in the inertial frame, pointing from body 2
   *  (or the obstacle) towards body 1. */
  siconos::algebra::SiconosVector3 nc_ = siconos::algebra::SiconosVector3::Zero();

  /** \f$ n^{\mathrm{loc}} \f$: the normal in the frame of body 2 (in the inertial frame if
   *  there is only one body). */
  siconos::algebra::SiconosVector3 relNc_ = siconos::algebra::SiconosVector3::Zero();

  /** Contact frame \f$ R_c \f$: rows \f$ n^\top, t_1^\top, t_2^\top \f$ (only the first row
   *  is used by NewtonEuler1DR). Internal buffer. */
  siconos::algebra::SiconosMatrix33 contactFrame_ = siconos::algebra::SiconosMatrix33::Zero();

  /** compute \f$ J_h \f$ (first row): gradient of \f$ n^\top (P_1 - P_2) \f$ with respect
   *  to \f$ q \f$, \f$ n \f$ being held fixed.
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  coordinates \f$ q \f$ of the bodies of the interaction
   */
  virtual void computeH_NE_(double time, siconos::modeling::Interaction& inter,
                            const siconos::algebra::BlockVector& q0) override;

  /** compute the contact frame \f$ R_c \f$ (contactFrame_) from the normal \f$ n \f$.
   *  Throws if \f$ n \f$ is zero.
   *
   *  \param withTangents if false, only the first row \f$ n^\top \f$ is set
   */
  void computeContactFrame_(bool withTangents);

  /** fill the first nRows rows of \f$ H \f$ (H_NE_prod_T_) with the first nRows rows of
   *  \f$ H_c(q) = R_c [\, I_3 \;\; [x_{G_1} - P]_\times R_1 \;\; -I_3 \;\;
   *  -[x_{G_2} - P]_\times R_2 \,] \f$, with \f$ P = P_1 \f$ and the current contactFrame_.
   *
   *  \param q0 coordinates \f$ q \f$ of the bodies of the interaction
   *  \param nRows number of rows (1: normal, 3: normal and tangential)
   */
  void computeContactVelocityJacobian_(const siconos::algebra::BlockVector& q0,
                                       Eigen::Index nRows);

 public:
  using NewtonEulerR::NewtonEulerR;

  /** Not used by the kernel (it was meant to select the projection onto
   *  \f$ y = 0 \f$ or \f$ y \geq 0 \f$). Kept for compatibility. */
  bool _isOnContact = false;

  /** destructor
   */
  virtual ~NewtonEuler1DR() noexcept = default;

  /** allocate \f$ J_h \f$ (1 row) and the operators of NewtonEulerR.
   *
   *  \param inter the interaction using this relation
   */
  void initialize(Interaction& inter) override;

  /** compute \f$ H = n^\top [\, I_3 \;\; [x_{G_1} - P]_\times R_1 \;\; -I_3 \;\;
   *  -[x_{G_2} - P]_\times R_2 \,] \f$ directly from the contact point and the normal
   *  (instead of the product \f$ J_h T \f$).
   *
   *  \param inter interaction that owns the relation
   *  \param q0 coordinates \f$ q \f$ of the bodies of the interaction
   */
  void computeH_NE_prod_T(const Interaction& inter,
                          const siconos::algebra::BlockVector& q0) override;

  /** update the contact points and the normal from the data attached to the bodies,
   *  \f$ P_1 = x_{G_1} + R_1 P_1^{\mathrm{loc}} \f$, \f$ P_2 = x_{G_2} + R_2
   *  P_2^{\mathrm{loc}} \f$, \f$ n = R_2 n^{\mathrm{loc}} \f$ (with one body, \f$ P_2 =
   *  P_2^{\mathrm{loc}} \f$ and \f$ n = n^{\mathrm{loc}} \f$), then call
   *  NewtonEulerR::computeh.
   *
   *  \warning NewtonEulerR::computeh computes \f$ y = J_h q + e \f$ with the current
   *  \f$ J_h \f$; the caller must set the meaningful components of \f$ y \f$ (e.g.
   *  \f$ y_1 = g \f$ with distance()).
   *
   *  \param time current time
   *  \param[in] q1 coordinates of body 1
   *  \param[in] q2 coordinates of body 2 (std::nullopt if there is only one body)
   *  \param[in,out] y the output vector
   */
  void computehFromRelativeContactPoints(
      double time, const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
      Eigen::Ref<siconos::algebra::SiconosVector> y);

  /** \return the gap \f$ g = \operatorname{sign}(n^\top (P_1 - P_2))\, \|P_1 - P_2\| \f$,
   *  positive when the bodies are separated */
  double distance() const;

  /** \return the contact point \f$ P_1 \f$ on body 1 (inertial frame) */
  inline auto pc1() const {
    return siconos::algebra::ConstMapVector3Type(contactPoint1_.data(), 3);
  }

  /** \return the contact point \f$ P_2 \f$ on body 2 (inertial frame) */
  inline const siconos::algebra::SiconosVector3& pc2() const { return contactPoint2_; }

  /** \return the unit normal \f$ n \f$ (inertial frame), from body 2 towards body 1 */
  inline const siconos::algebra::SiconosVector3& nc() const { return nc_; }

  /** \return \f$ P_1^{\mathrm{loc}} \f$, contact point 1 in the frame of body 1 */
  inline const siconos::algebra::SiconosVector3& relPc1() const { return relPc1_; }
  /** \return \f$ P_2^{\mathrm{loc}} \f$, contact point 2 in the frame of body 2 */
  inline const siconos::algebra::SiconosVector3& relPc2() const { return relPc2_; }
  /** \return \f$ n^{\mathrm{loc}} \f$, the normal in the frame of body 2 */
  inline const siconos::algebra::SiconosVector3& relNc() const { return relNc_; }

  /** Set \f$ P_1^{\mathrm{loc}} \f$, the contact point 1 in the frame of body 1.
   *  Used by computehFromRelativeContactPoints.
   *
   *  \param npc new coordinates
   */
  void setRelPc1(const siconos::algebra::SiconosVector3& npc) { relPc1_ = npc; };

  /** Set \f$ P_2^{\mathrm{loc}} \f$, the contact point 2 in the frame of body 2 (in the
   *  inertial frame if there is only one body). Used by computehFromRelativeContactPoints.
   *
   *  \param npc new coordinates
   */
  void setRelPc2(const siconos::algebra::SiconosVector3& npc) { relPc2_ = npc; };

  /** Set \f$ n^{\mathrm{loc}} \f$, the unit normal (from body 2 towards body 1) in the
   *  frame of body 2 (in the inertial frame if there is only one body). Used by
   *  computehFromRelativeContactPoints.
   *
   *  \param nnc new coordinates
   */
  void setRelNc(const siconos::algebra::SiconosVector3& nnc) { relNc_ = nnc; };

  /** Compute current contact points from relative contact points stored in cohesive zone model
   *  This is used for cohesive zone models to track the evolution of contact points.
   *
   *  \param q0 block vector of dynamical system positions
   *  \param r_pc1 relative contact point 1 (input)
   *  \param r_pc2 relative contact point 2 (input)
   *  \param r_nc relative normal (input)
   *  \param r_t1 relative tangent 1 (input)
   *  \param r_t2 relative tangent 2 (input)
   *  \param pc1 computed absolute contact point 1 (output)
   *  \param pc2 computed absolute contact point 2 (output)
   *  \param nc computed absolute normal (output)
   *  \param t1 computed absolute tangent 1 (output)
   *  \param t2 computed absolute tangent 2 (output)
   */
  void computeContactPointsFromRelativeContactPoints(
      const siconos::algebra::BlockVector& q0, const siconos::algebra::SiconosVector3& r_pc1,
      const siconos::algebra::SiconosVector3& r_pc2,
      const siconos::algebra::SiconosVector3& r_nc,
      const siconos::algebra::SiconosVector3& r_t1,
      const siconos::algebra::SiconosVector3& r_t2, siconos::algebra::SiconosVector3& pc1,
      siconos::algebra::SiconosVector3& pc2, siconos::algebra::SiconosVector3& nc,
      siconos::algebra::SiconosVector3& t1, siconos::algebra::SiconosVector3& t2);

  void display() const override {}
  virtual void accept(relations::Visitor& tourist) const override { tourist.visit(*this); }
};
}  // namespace siconos::modeling
#endif  // NEWTONEULER1DR_HPP
