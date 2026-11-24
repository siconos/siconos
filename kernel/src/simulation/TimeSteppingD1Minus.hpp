// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*!\file
 * TimeSteppingD1Minus simulation
 */

#ifndef TIMESTEPPINGD1MINUS_H
#define TIMESTEPPINGD1MINUS_H

#include "Simulation.hpp"
#include "Topology.hpp"

namespace siconos::simulation {

/** TimeSteppingD1Minus Timestepping Strategy
 *
 *  see Schindler/Acary : Timestepping Schemes for Nonsmooth Dynamics Based
 *  on Discontinuous Galerkin Methods: Definition and Outlook
 */
class TimeSteppingD1Minus : public Simulation {
 private:
  ACCEPT_SERIALIZATION(TimeSteppingD1Minus);

  // /** default constructor */
  // TimeSteppingD1Minus() = default;

 protected:
  /** initialisation specific to TimeSteppingD1Minus for OneStepNSProblem */
  void initializeOneStepNSProblem() override;

 public:
  /** constructor with the time-discretisation
   * \param nsds the current nonsmooth dynamical system
   * \param td pointer to a TimeDiscretisation
   * \param nb number of non smooth problem
   */
  TimeSteppingD1Minus(std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds,
                      std::shared_ptr<TimeDiscretisation> td, int nb);

  /** destructor */
  ~TimeSteppingD1Minus() noexcept = default;

  /** updateIndexSet using current y and lambda values of interactions
   *  \param i the  number of the set to be updated
   *  0 : ALL interactions (NEVER)
   *  1 : ACTIVE interactions for IMPACTS
   *  2 : ACTIVE interactions for CONTACTS
   */
  void updateIndexSet(siconos::simulation::Topology::size_type i) override;

  /** run the simulation, from t0 to T */
  void run() override;

  /** step from current event to next event of EventsManager */
  void advanceToEvent() override;

  /** compute residu */
  void computeResidu();

  /** integrate DynamicalSystems taking not into account non-smooth part */
  void computeFreeState();
};
}  // namespace siconos::simulation

#endif  // TIMESTEPPINGD1MINUS_H
