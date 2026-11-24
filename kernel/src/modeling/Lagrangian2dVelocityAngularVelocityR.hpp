/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
/*! \file Lagrangian2d#DR.hpp

 */
#ifndef LAGRANGIAN2DVELOCITYANGULARVELOCITYR_HPP
#define LAGRANGIAN2DVELOCITYANGULARVELOCITYR_HPP

#include "LagrangianScleronomousR.hpp"

namespace siconos::modeling {
/** Unilateral contact with friction and rolling resistance between two 2D rigid bodies (or a
 *  body and a fixed obstacle), as Lagrangian systems: normal and tangential components of the
 *  relative velocity at the contact point and relative angular velocity,
 *  \f$ \dot y = (u_N, u_T, \omega_1 - \omega_2) = H v \f$, for a
 *  NewtonImpactRollingFrictionNSL of size 3. The first component of \f$ y \f$ is the gap.
 *
 *  See Lagrangian2dNormalVelocityR for the notation and \f$ H \f$ (three rows).
 *
 *  Former name: Lagrangian2d3DR (deprecated alias).
 */
class Lagrangian2dVelocityAngularVelocityR : public LagrangianScleronomousR {
 protected:
  ACCEPT_SERIALIZATION(Lagrangian2dVelocityAngularVelocityR);

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
  virtual ~Lagrangian2dVelocityAngularVelocityR() noexcept {};

  void initialize(Interaction& inter) override;

  /**
     to compute the output y = h(q,z) of the Relation

     \param q coordinates of the dynamical systems involved in the relation
     \param z user defined parameters (optional)
     \param y the resulting vector
  */
  void computeh(const siconos::algebra::BlockVector& q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  /**
     to compute the jacobian of h(...). Set attribute jacobianhOver_q_ (access: jacqhq())

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
  virtual void accept(relations::Visitor& tourist) const override { tourist.visit(*this); }
};
/** \deprecated former name of Lagrangian2dVelocityAngularVelocityR, kept for one release. */
using Lagrangian2d3DR [[deprecated("use Lagrangian2dVelocityAngularVelocityR")]] =
    Lagrangian2dVelocityAngularVelocityR;

}  // namespace siconos::modeling
#endif  // NEWTONEULERRIMPACT_H
