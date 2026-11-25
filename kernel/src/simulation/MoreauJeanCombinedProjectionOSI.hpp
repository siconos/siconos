// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file
  MoreauJeanOSI Time-Integrator for Dynamical Systems for Combined Projection
  Algorithm
*/

#ifndef MOREAUCOMBINEDPROJECTIONOSI_H
#define MOREAUCOMBINEDPROJECTIONOSI_H

#include "MoreauJeanOSI.hpp"

namespace siconos::integrators {
/**
    One Step time Integrator for First Order Dynamical Systems  for
    mechanical Systems (LagrangianDS and NewtonEulerDS) with  Combined
    Projection Algorithm

    This class reimplement a special activation of constraints
    in the MoreauJeanOSI for the Combined Projection Algorithm
    (Algorithm 2, used with siconos::simulation::TimeSteppingCombinedProjection).

    Index sets (numberOfIndexSets() = 3):
    - indexSet 1, active set \f$ \mathcal{A} \f$ of the velocity problem: cleared at the
      beginning of each step; a constraint is added when it is violated at the end of an
      index set iteration, \f$ g_\alpha(q) < -\epsilon \f$, with
      \f$ \epsilon \f$ = TimeSteppingCombinedProjection::constraintTolUnilateral()
      (0 with another simulation). Closed contacts (\f$ g_\alpha = 0 \f$) are not added:
      they would be held closed by the equality projection. Constraints are never removed
      within a step (removeInteractionFromIndexSet() is not used).
    - indexSet 2, set of the position projection: all of indexSet 1.

    The input of level 0 (p(0), the projection impulse) is also allocated
    (_levelMinForInput = 0). See docs/devel/PROJECTION_SCHEMES.md.

    References :

    V. Acary. Projected event-capturing time-stepping schemes for nonsmooth
    mechanical systems with unilateral contact and coulomb’s friction. Computer
    Methods in Applied Mechanics and Engineering, 256:224 – 250, 2013. ISSN
    0045-7825. URL
    http://www.sciencedirect.com/science/article/pii/S0045782512003829.

*/
class MoreauJeanCombinedProjectionOSI : public MoreauJeanOSI {
 protected:
  ACCEPT_SERIALIZATION(MoreauJeanCombinedProjectionOSI);

 public:
  /** constructor from theta value only
   *
   *  \param theta value for all these DS.
   */
  explicit MoreauJeanCombinedProjectionOSI(double theta);

  /** destructor
   */
  virtual ~MoreauJeanCombinedProjectionOSI() noexcept = default;

  // --- OTHER FUNCTIONS ---

  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi
   *
   *  \param t time of initialization
   *  \param ds the dynamical system
   */
  void initializeWorkVectorsForDS(
      double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) override;

  /** initialization of the work vectors and matrices (properties) related to
   *  one interaction on the graph and needed by the osi
   *
   *  \param inter the interaction
   *  \param interProp the properties on the graph
   *  \param DSG the dynamical systems graph
   */
  void initializeWorkVectorsForInteraction(
      siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
      siconos::graphs::DynamicalSystemsGraph& DSG) override;

  /** get the number of index sets required for the simulation
   *
   *  \return unsigned int
   */
  unsigned int numberOfIndexSets() const override { return 3; };

  /** Apply the rule to one Interaction to known if is it should be included
   *  in the IndexSet of level i
   *
   *  - i = 1: \f$ g_\alpha(q) < -\epsilon \f$ (violated constraint, \f$ g = \f$ y(0)),
   *    \f$ \epsilon \f$ = TimeSteppingCombinedProjection::constraintTolUnilateral();
   *  - i = 2: always true (every active contact is projected).
   *
   *  \param inter concerned interaction
   *  \param i level
   *  \return bool
   */
  bool addInteractionInIndexSet(std::shared_ptr<siconos::modeling::Interaction> inter,
                                siconos::graphs::InteractionsGraph::size_type i) override;

  /** Apply the rule to one Interaction to known if is it should be removed
   *  in the IndexSet of level i
   *
   *  \param inter concerned interaction
   *  \param i level
   *  \return bool
   */
  bool removeInteractionFromIndexSet(std::shared_ptr<siconos::modeling::Interaction> inter,
                                     siconos::graphs::InteractionsGraph::size_type i) override;
};
}  // namespace siconos::integrators

#endif  // MOREAUCOMBINEDPROJECTIONOSI_H
