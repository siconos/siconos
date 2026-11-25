/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
/*! \file
  MoreauJeanOSI Time-Integrator for Dynamical Systems for Combined Projection
  Algorithm
*/

#ifndef MOREAUPROJECTONCONSTRAINTSOSI_H
#define MOREAUPROJECTONCONSTRAINTSOSI_H

#include "MoreauJeanOSI.hpp"

namespace siconos::integrators {

/**
    One Step time Integrator for First Order Dynamical Systems  for
    mechanical Systems (LagrangianDS and NewtonEulerDS) with  Direct
    Projection Algorithm

    This class reimplement a special activation of constraints
    in the MoreauJeanOSI for the Direct Projection Algorithm
    (Algorithm 1, used with siconos::simulation::TimeSteppingDirectProjection).

    Activation of indexSet 1 (contacts of the velocity problem), rule
    FIRSTWAY_ACTIVATION of MoreauJeanDirectProjectionOSI.cpp, with the forecast gap
    \f$ \tilde g_\alpha = g_\alpha + \gamma h \dot g_\alpha \f$ evaluated at the
    beginning of the step (outputs y(0), y(1)):
    - add if \f$ \tilde g_\alpha \leq \f$ activateYPosThreshold() (1e-7 by default);
    - remove if \f$ \tilde g_\alpha > \f$ deactivateYPosThreshold() (1e-7) and
      \f$ \dot g_\alpha \geq \f$ deactivateYVelThreshold() (0).

    \f$ \gamma = 1/2 \f$ unless set by MoreauJeanOSI::setGamma() or the constructor
    (gamma). The scheme of Acary (CMAME 2013) and the Scilab code Hots use
    \f$ \gamma = 1 \f$ and a threshold of 1e-10: see docs/devel/PROJECTION_SCHEMES.md.

    numberOfIndexSets() = 2: indexSet 0 (all interactions, used for the projection when
    TimeSteppingDirectProjection is built with level 0) and indexSet 1.

    References :

    V. Acary. Projected event-capturing time-stepping schemes for nonsmooth
    mechanical systems with unilateral contact and coulomb’s friction. Computer
    Methods in Applied Mechanics and Engineering, 256:224 – 250, 2013. ISSN
    0045-7825. URL
    http://www.sciencedirect.com/science/article/pii/S0045782512003829.

*/
class MoreauJeanDirectProjectionOSI : public MoreauJeanOSI {
 protected:
  ACCEPT_SERIALIZATION(MoreauJeanDirectProjectionOSI);

  static constexpr unsigned int MOREAUPROJECTONCONSTRAINTSOSISTEPSINMEMORY = 1;
  static constexpr double SICONOS_MPC_DEFAULT_ACTIVATION_POS_THRESHOLD = 1.e-7;
  static constexpr double SICONOS_MPC_DEFAULT_ACTIVATION_VEL_THRESHOLD = 0.0;
  static constexpr double SICONOS_MPC_DEFAULT_DEACTIVATION_POS_THRESHOLD = 1.e-7;
  static constexpr double SICONOS_MPC_DEFAULT_DEACTIVATION_VEL_THRESHOLD = 0.0;

  double _deactivateYPosThreshold{SICONOS_MPC_DEFAULT_DEACTIVATION_POS_THRESHOLD};
  double _deactivateYVelThreshold{SICONOS_MPC_DEFAULT_DEACTIVATION_VEL_THRESHOLD};
  double _activateYPosThreshold{SICONOS_MPC_DEFAULT_ACTIVATION_POS_THRESHOLD};
  double _activateYVelThreshold{SICONOS_MPC_DEFAULT_ACTIVATION_VEL_THRESHOLD};

 public:
  /** constructor from theta value only
   *
   *  \param theta value for all these DS.
   */
  explicit MoreauJeanDirectProjectionOSI(double theta);

  /** constructor from theta value only
   *
   *  \param theta value for all these DS.
   *  \param gamma value for all these DS.
   */
  explicit MoreauJeanDirectProjectionOSI(double theta, double gamma);

  /** destructor
   */
  virtual ~MoreauJeanDirectProjectionOSI() noexcept = default;

  inline double deactivateYPosThreshold() { return _deactivateYPosThreshold; };

  inline void setDeactivateYPosThreshold(double newValue) {
    _deactivateYPosThreshold = newValue;
  };

  inline double deactivateYVelThreshold() { return _deactivateYVelThreshold; };

  inline void setDeactivateYVelThreshold(double newValue) {
    _deactivateYVelThreshold = newValue;
  };

  inline double activateYPosThreshold() { return _activateYPosThreshold; };

  inline void setActivateYPosThreshold(double newValue) { _activateYPosThreshold = newValue; };

  inline double activateYVelThreshold() { return _activateYVelThreshold; };

  inline void setActivateYVelThreshold(double newValue) { _activateYVelThreshold = newValue; };

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
  unsigned int numberOfIndexSets() const override { return 2; };

  /** Apply the rule to one Interaction to known if is it should be included
   *  in the IndexSet of level i
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

  /** Perform the integration of the dynamical systems linked to this integrator
   *  without taking into account the nonsmooth input (_r or _p)
   */
  void computeFreeState() override;
};
}  // namespace siconos::integrators
#endif  // MOREAUPROJECTONCONSTRAINTSOSI_H
