/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/*! \file MultipleImpactNSL.hpp
 */

#ifndef _MULTIPLEIMPACTNSL_
#define _MULTIPLEIMPACTNSL_

#include "NonSmoothLaw.hpp"

namespace siconos::modeling {

class MultipleImpactNSL : public NonSmoothLaw {
 private:
  ACCEPT_SERIALIZATION(MultipleImpactNSL);

  // Energytical restitution coefficient
  double _ResCof{0.};
  // Normal stiffness at contact
  double _Stiff{0.};
  // Elasticity coefficient
  double _ElasCof{0.};

 public:
  // Default Constructor
  MultipleImpactNSL() : NonSmoothLaw(1) {};
  // Constructor with parameters
  MultipleImpactNSL(double, double, double, siconos::algebra::Index _dim = 1);
  // Destructor
  ~MultipleImpactNSL() noexcept = default;
  // Get the value of the energytical restitution coefficientx
  inline double ResCof() const { return _ResCof; };
  // Get the value of the stiffness
  inline double Stiff() const { return _Stiff; };
  // Get the value of the elasticity coefficient
  inline double ElasCof() const { return _ElasCof; }
  // Set the value to the restitution coefficient
  void setResCof(double newResCof);
  // Set the value to the stiffness
  void setStiff(double newStiff);
  // Set the value to the elasticity cofficient
  void setElasCoeff(double _newElasCoef);
  //
  bool isVerified() const override;
  // Display the information about the multiple impact law
  void display() const override;
  // visitors hook
  virtual void accept(nonsmooth_laws::Visitor& tourist) const override {
    tourist.visit(*this);
  }

  Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling
#endif
