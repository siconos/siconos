// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file EqualityConditionNSL.hpp

*/
#ifndef EQUALITYCONDITIONNSLAW_H
#define EQUALITYCONDITIONNSLAW_H

#include "NonSmoothLaw.hpp"
namespace siconos::modeling {
/** Equality NonSmoothLaw
 *
 **/
class EqualityConditionNSL : public NonSmoothLaw {
 private:
  ACCEPT_SERIALIZATION(EqualityConditionNSL);

  /** default constructor
   */
  EqualityConditionNSL() = default;

 public:
  /** basic constructor
   *
   *  \param size of the non smooth law
   */
  EqualityConditionNSL(siconos::algebra::Index size) : NonSmoothLaw(size) {};

  /** Destructor */
  ~EqualityConditionNSL() noexcept = default;

  /** display the data of the NonSmoothLaw on the standard output
   *
   */
  void display() const override {};

  // visitors hook
  virtual void accept(nonsmooth_laws::Visitor& tourist) const override {
    tourist.visit(*this);
  }

  Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling

#endif  // EQUALITYCONDITIONNSLAW_H
