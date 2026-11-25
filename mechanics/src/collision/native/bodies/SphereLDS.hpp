// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file SphereLDS.hpp
  \brief Definition of a 3D sphere as a LagrangianDS (with Euler
         Angles)
*/

#ifndef SphereLDS_h
#define SphereLDS_h

#include "LagrangianDS.hpp"

namespace siconos::collision::native::bodies {

class SphereLDS : public siconos::modeling::LagrangianDS,
                  public std::enable_shared_from_this<SphereLDS> {
 protected:
  ACCEPT_SERIALIZATION(SphereLDS);

  double radius_{0.};
  double massValue_{0.};
  double inertia_{0.};

 public:
  /** constructor from initial state and velocity
   *
   *  \param R radius of the sphere
   *  \param m mass of the sphere
   *  \param position initial coordinates
   *  \param velocity initial velocity
   */

  SphereLDS(double, double, const siconos::algebra::SiconosVector6& position,
            const siconos::algebra::SiconosVector6& velocity);

  ~SphereLDS() noexcept = default;

  /** to compute the mass matrix operator \f$ M(q) \f$
   *
   *  \param position q vector
   */
  void computeMass(const Eigen::Ref<const siconos::algebra::SiconosVector>& position) override;

  double getQ(unsigned int pos);

  double getVelocity(unsigned int pos);

  inline double getMassValue() const { return massValue_; };

  inline double getRadius() const { return radius_; };
  virtual void acceptSP(modeling::dynamical_systems::Visitor& tourist) override {
    tourist.visit(shared_from_this());
  }
  virtual void accept(modeling::dynamical_systems::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};

}  // namespace siconos::collision::native::bodies
#endif /* SphereLDS_h */
