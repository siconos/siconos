// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file ComplementarityConditionNSL.hpp

*/
#ifndef COMPLEMENTARITYCONDITIONNSLAW_H
#define COMPLEMENTARITYCONDITIONNSLAW_H

#include "NonSmoothLaw.hpp"
namespace siconos::modeling {
/** Complementarity NonSmoothLaw
 *
 **/
class ComplementarityConditionNSL : public NonSmoothLaw {
 private:
  ACCEPT_SERIALIZATION(ComplementarityConditionNSL);

  /** default constructor
   */
  ComplementarityConditionNSL() {};

 public:
  /** basic constructor
   *
   *  \param size of the non smooth law
   */
  ComplementarityConditionNSL(siconos::algebra::Index size) : NonSmoothLaw(size) {};

  /** Destructor */
  ~ComplementarityConditionNSL() noexcept = default;

  /** print the data to the screen
   */
  inline void display() const override {};

  virtual void accept(nonsmooth_laws::Visitor& tourist) const override {
    tourist.visit(*this);
  }

  Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling

#endif  // COMPLEMENTARITYCONDITIONNSLAW_H
