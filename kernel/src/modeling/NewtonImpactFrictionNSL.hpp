// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file NewtonImpactFrictionNSL.hpp
  Newton-Impact Non-Smooth Law
*/

#ifndef NEWTONIMPACTFRICTIONNSLAW_H
#define NEWTONIMPACTFRICTIONNSLAW_H

#include "NonSmoothLaw.hpp"

namespace siconos::modeling {

/** Newton Impact-Friction Non Smooth Law
 *
 */
class NewtonImpactFrictionNSL : public NonSmoothLaw {
 private:
  ACCEPT_SERIALIZATION(NewtonImpactFrictionNSL);

  /** The Newton coefficient of restitution
   */
  double _en{0.};
  double _et{0.};
  /** friction coefficient */
  double _mu{0.};

  /** default constructor
   */
  NewtonImpactFrictionNSL() = default;

 public:
  /** basic constructor
   *
   *  \param size size of the ns law
   */
  NewtonImpactFrictionNSL(siconos::algebra::Index size) : NonSmoothLaw(size) {};

  /** constructor with the value of the NewtonImpactFrictionNSL attributes
   *
   *  \param en double : normal e coefficient
   *  \param et double : tangent e coefficient
   *  \param mu double : friction coefficient
   *  \param size size of the ns law
   */
  NewtonImpactFrictionNSL(double en, double et, double mu, siconos::algebra::Index size)
      : NonSmoothLaw(size), _en(en), _et(et), _mu(mu) {};

  /** Destructor */
  ~NewtonImpactFrictionNSL() noexcept = default;

  /** check the ns law to see if it is verified
   *
   *  \return a boolean value whioch determines if the NS Law is verified
   */
  bool isVerified(void) const override;

  // GETTERS/SETTERS

  /** \return the value of en
   */
  inline double en() const { return _en; };

  /** setter of en
   *
   *  \param newVal a double to set en
   */
  inline void setEn(double newVal) { _en = newVal; };

  /** \return the value of et
   */
  inline double et() const { return _et; };

  /** setter of et
   *
   *  \param newVal a double to set et
   */
  inline void setEt(double newVal) { _et = newVal; };

  /** \return the value of mu
   */
  inline double mu() const { return _mu; };

  /** setter of mu
   *
   *  \param newVal a double to set mu
   */
  inline void setMu(double newVal) { _mu = newVal; };

  // OTHER FUNCTIONS

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
#endif  // NewtonImpactFrictionNSL_H
