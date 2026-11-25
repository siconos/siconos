// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file Lagrangian2dNormalVelocityR.hpp

 */
#ifndef LAGRANGIAN2DNORMALVELOCITYR_HPP
#define LAGRANGIAN2DNORMALVELOCITYR_HPP

#include "LagrangianScleronomousR.hpp"

namespace siconos::modeling {
/** Unilateral contact between two 2D rigid bodies (or a body and a fixed obstacle), as
 *  Lagrangian systems: normal component of the relative velocity at the contact point.
 *
 *  This class is also the model of Lagrangian2dVelocityR (with friction) and
 *  Lagrangian2dVelocityAngularVelocityR (with friction and rolling resistance), which share
 *  the notation below.
 *
 *  Body \f$ i \f$ has coordinates \f$ q_i = (x_{G_i}, \theta_i) \in \mathbb{R}^3 \f$ and
 *  velocities \f$ v_i = \dot q_i = (v_{G_i}, \omega_i) \f$. The contact points
 *  \f$ P_1, P_2 \f$ and the unit normal \f$ n \f$ (from body 2 towards body 1) are set by the
 *  user or a derived class (updateContactPoints); the gap is
 *  \f$ g = \operatorname{sign}(n^\top (P_1 - P_2))\, \|P_1 - P_2\| \f$ (distance()). With
 *  \f$ a^\perp = (-a_y, a_x) \f$, the tangent \f$ t = n^\perp \f$ and \f$ P = P_1 \f$, the
 *  relative velocity of body 1 with respect to body 2 at \f$ P \f$ is
 *  \f[
 *    u = v_{G_1} + \omega_1 (P - x_{G_1})^\perp - v_{G_2} - \omega_2 (P - x_{G_2})^\perp .
 *  \f]
 *  The contact configuration matrix \f$ H \f$ (contactMatrix(), computed in
 * computeContactMatrix) gives
 *  \f$ \dot y = H v \f$:
 *  \f[
 *    H = \begin{bmatrix}
 *      n^\top & n^\top (P - x_{G_1})^\perp & -n^\top & -n^\top (P - x_{G_2})^\perp \\
 *      t^\top & t^\top (P - x_{G_1})^\perp & -t^\top & -t^\top (P - x_{G_2})^\perp \\
 *      0 & 1 & 0 & -1
 *    \end{bmatrix}
 *  \f]
 *  (only the first two block columns with a single body):
 *  - Lagrangian2dNormalVelocityR: first row, \f$ y = g \f$, \f$ \dot y = u_N = n^\top u \f$;
 *  - Lagrangian2dVelocityR: first two rows, \f$ \dot y = (u_N, u_T) \f$;
 *  - Lagrangian2dVelocityAngularVelocityR: the three rows, \f$ \dot y = (u_N, u_T,
 *    \omega_1 - \omega_2) \f$.
 *
 *  The first row is the gradient of the gap when \f$ P_1 - P_2 \f$ is parallel to \f$ n \f$;
 *  the other rows are not the Jacobian of a function.
 *
 *  Former name: Lagrangian2d1DR (deprecated alias).
 */
class Lagrangian2dNormalVelocityR : public LagrangianScleronomousR {
 protected:
  ACCEPT_SERIALIZATION(Lagrangian2dNormalVelocityR);

  /* Current Contact Points, may be updated within Newton loop based
   * on _relPc1, _relPc2. */
  siconos::algebra::SiconosVector2 contactPoint1_ = siconos::algebra::SiconosVector2::Zero();
  siconos::algebra::SiconosVector2 contactPoint2_ = siconos::algebra::SiconosVector2::Zero();

  /* Inward Normal at the contact.
   * \todo The meaning of "Inward" has to be explained carefully.
   */
  siconos::algebra::SiconosVector2 nc_ = siconos::algebra::SiconosVector2::Zero();

 public:
  /** Not used by the kernel (it was meant to select the projection onto
   *  \f$ y = 0 \f$ or \f$ y \geq 0 \f$). Kept for compatibility. */
  bool _isOnContact = false;

  using LagrangianScleronomousR::LagrangianScleronomousR;

  /** destructor */
  virtual ~Lagrangian2dNormalVelocityR() noexcept {};

  void initialize(Interaction& inter) override;

  /** to compute the output y = h(q,z) of the Relation
      \param q coordinates of the dynamical systems involved in the relation
      \param z user defined parameters (optional)
      \param y the resulting vector
  */
  void computeh(const siconos::algebra::BlockVector& q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  /** to compute the jacobian of h(...). Set attribute jacobianhOver_q_ (access: jacqhq())
      \param q coordinates of the dynamical systems involved in the relation
      \param z user defined parameters (optional)
  */
  void computeContactMatrix(const siconos::algebra::BlockVector& q) override;

  /** Return the distance between pc1 and pc, with sign according to normal */
  double distance() const;

  inline const siconos::algebra::SiconosVector2& pc1() const { return contactPoint1_; }

  inline const siconos::algebra::SiconosVector2& pc2() const { return contactPoint2_; }

  inline const siconos::algebra::SiconosVector2& nc() const { return nc_; }

  /** @brief Update contact point info.
   *
   *  @param pos1 Position on ds1 in ds1 frame.
   *  @param pos2 Position on ds2 in ds2 frame (or world frame if ds2=null).
   *  @param normal Normal in ds2 frame (or world frame if ds2=null).
   */
  inline void updateContactPoints(const siconos::algebra::SiconosVector2& pos1,
                                  const siconos::algebra::SiconosVector2& pos2,
                                  const siconos::algebra::SiconosVector2& normal) {
    assert(normal.norm() != 0);
    contactPoint1_ = pos1;
    contactPoint2_ = pos2;
    nc_ = normal;
  }

  void display() const override;
};
/** \deprecated former name of Lagrangian2dNormalVelocityR, kept for one release. */
using Lagrangian2d1DR [[deprecated("use Lagrangian2dNormalVelocityR")]] =
    Lagrangian2dNormalVelocityR;

}  // namespace siconos::modeling
#endif  // NEWTONEULERRIMPACT_H
