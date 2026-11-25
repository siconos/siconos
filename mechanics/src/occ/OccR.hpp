/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
/** \file OccR.hpp
    \brief A Siconos Newton Euler 3D friction relation between
    two BRep contact points.
 */

#ifndef OccR_hpp
#define OccR_hpp

#include "Geometer.hpp"              // For  DistanceCalculator
#include "NewtonEulerVelocityR.hpp"  // Base class

namespace siconos::mechanics::occ {

struct ContactPoint;

class OccR : public siconos::modeling::NewtonEulerVelocityR {
 protected:
  const ContactPoint& _contact1;
  const ContactPoint& _contact2;

  DistanceCalculator dt{CadmbtbDistanceType{}};

  //  std::shared_ptr<Geometer> _geometer{nullptr};

  double _offset1{0.};
  double _offset2{0.};

 public:
  /** Constructor from contact points.
   *
   *  \param contact1 : the first contact.
   *  \param contact2 : the second contact.
   */
  OccR(const ContactPoint& contact1, const ContactPoint& contact2,
       const DistanceCalculator& distance_calculator);

  // Rule of 5
  OccR() = delete;
  OccR(const OccR&) = delete;
  OccR(OccR&&) = delete;
  OccR& operator=(const OccR&) = delete;
  OccR& operator=(OccR&&) = delete;
  ~OccR() noexcept = default;

  /**
     to compute the output y = h(q) of the Relation

      \param[in] q1 generalized coordinates vector of the fist dynamical system involved
      in the relation
      \param[in] q2 generalized coordinates vector of the second dynamical system
      involved in the relation
      \param[in,out] y the resulting vector
  */
  virtual void computeh(
      const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
      Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  /** Set offset1, offset from first contact.
   *
   *  \param val : the new value.
   */
  void setOffset1(double val) { _offset1 = val; };

  /** Set offset2, offset from second contact.
   *
   *  \param val : the new value.
   */
  void setOffset2(double val) { _offset2 = val; };

  // /** Get geometer.
  //  *
  //  *  \return a std::shared_ptr<Geometer> object.
  //  */
  // std::shared_ptr<Geometer> geometer() { return _geometer; }

  // /** Set geometer.
  //  *
  //  *  \param geometer the new geometer
  //  */
  // void setGeometer(std::shared_ptr<Geometer> geometer) { _geometer = geometer; }
};
}  // namespace siconos::mechanics::occ
#endif
