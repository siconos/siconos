// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! @file FEM-GlobalFrictionContact.hpp
  specific implementation of Primal Fricton-Contact Non-Smooth Problem
  with SolidLinearTIDS
*/
#ifndef FEMGlobalFrictionContact_H
#define FEMGlobalFrictionContact_H

#include "GlobalFrictionContact.hpp"

namespace siconos::nonsmooth_formulations {
class GlobalFrictionContact;
}

namespace siconos::mechanics::fem::nonsmooth_formulations {

/**
   Formalization and Resolution of a Friction-Contact Problem
   when SolidLinearTIDS are used in the NSDS

   See details in GlobalFrictionContact base class
*/
class GlobalFrictionContact : public siconos::nonsmooth_formulations::GlobalFrictionContact {
 protected:
  ACCEPT_SERIALIZATION(GlobalFrictionContact);

  /** update the rhs of the ns problem with nonsmooth laws contributions
   *
   *  @param indexSet set of active interactions
   */
  virtual void compute_nslaw_contribution(
      siconos::graphs::InteractionsGraph& indexSet) override;

 public:
  using siconos::nonsmooth_formulations::GlobalFrictionContact::GlobalFrictionContact;

  /** destructor
   */
  virtual ~GlobalFrictionContact() noexcept = default;

  /** compute vector q */
  virtual void compute_q() override;

  /** update dynamical systems state with the current global velocity values */
  virtual void update_dynamicalsystems_state() override;
};
}  // namespace siconos::mechanics::fem::nonsmooth_formulations

#endif
