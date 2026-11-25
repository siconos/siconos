/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
/*! \file NewtonImpactNSL.hpp

*/
#ifndef NEWTONIMPACTNSL_H
#define NEWTONIMPACTNSL_H

#include "NonSmoothLaw.hpp"

namespace siconos::modeling {
/**
   Newton impact Non Smooth Law

   This class formalizes the Newton Impact law together with a complementarity
   condition. i.e.

   \f[
   \left\{\begin{array}{l}
   y \geq 0, \lambda \geq 0, y^{T} \lambda=0\				\
   if y \leq 0 \quad \mbox{then} \quad \dot y(t^{+}) - e \dot y(t^{-}) \geq
   0, \quad  \lambda \geq 0, (\dot y(t^{+}) - e \dot y(t^{-}))^{T} \lambda=0
   \end{array}\right.
   \f]


   nsLawSize is equal to 1.

*/
class NewtonImpactNSL : public NonSmoothLaw {
 private:
  ACCEPT_SERIALIZATION(NewtonImpactNSL);

  /** The Newton normal coefficient of restitution  */
  double _e{0.};

  /** default constructor
   */
  NewtonImpactNSL() = delete;

 public:
  /** constructor with the value of the NewtonImpactNSL attributes
   *
   *  \param e the value of the coefficient of restitution
   */
  NewtonImpactNSL(double e) : NonSmoothLaw(1), _e(e) {};

  /** Apply multiple-axis impact */
  NewtonImpactNSL(siconos::algebra::Index size, double e) : NonSmoothLaw(size), _e(e) {};

  /** destructor
   */
  ~NewtonImpactNSL() noexcept = default;

  /** check the ns law to see if it is verified
   *
   *  \return a boolean value whioch determines if the NS Law is verified
   */
  bool isVerified() const override;

  /** \return the value of e
   */
  inline double e() const { return _e; };

  /** setter of e
   *
   *  \param newVal a double to set e
   */
  inline void setE(double newVal) { _e = newVal; };

  /** print the data to the screen
   */
  void display() const override;

  // visitors hook
  virtual void accept(nonsmooth_laws::Visitor& tourist) const override {
    tourist.visit(*this);
  }

  Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling

#endif  // NewtonImpactNSL_H
