// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file DiskDiskR.hpp
 */
#ifndef DiskDiskR_h
#define DiskDiskR_h

#include "CircularR.hpp"

namespace siconos::collision::native::bodies {

/** \class DiskDiskR
 *  \brief Two disks relation - Inherits from LagrangianScleronomousR
 */
class DiskDiskR : public CircularR, public std::enable_shared_from_this<DiskDiskR> {
 private:
  ACCEPT_SERIALIZATION(DiskDiskR);

  double r1pr2{0.};

 public:
  /** Constructor
   *
   *  \param disk1 radius
   *  \param disk2 radius
   */
  DiskDiskR(double disk1, double disk2);

  ~DiskDiskR() noexcept = default;

  double distance(double, double, double, double, double, double) override;

  /**
    to compute the output y = h(q) of the Relation

    \param q coordinates of the dynamical systems involved in the relation
    \param y the resulting vector
  */
  void computeh(const siconos::algebra::BlockVector& q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  /** Computes \f$ \nabla^\top_q h(q) \f$
   * \param q coordinates of the dynamical systems involved in the relation
   */
  void computeContactMatrix(const siconos::algebra::BlockVector& q) override;
  virtual void accept(modeling::relations::Visitor& tourist) const override {
    tourist.visit(*this);
  }
  void display() const override;
};
}  // namespace siconos::collision::native::bodies
#endif /* DiskDiskR_h */
