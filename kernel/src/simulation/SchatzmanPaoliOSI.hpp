// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file
  SchatzmanPaoliOSI Time-Integrator for Dynamical Systems
*/

#ifndef SCHATZMANPAOLIOSI_H
#define SCHATZMANPAOLIOSI_H

#include "OneStepIntegrator.hpp"

namespace siconos::integrators {

/**
   Schatzman-Paoli time-integrator (two-step, position-based scheme for impacts).

   Only implemented for siconos::modeling::LagrangianLinearTIDS:
   \f$ M\ddot q + C\dot q + Kq = F_{ext}(t) + G^T\lambda \f$, with unilateral constraints
   \f$ g(q) = Hq + b \ge 0 \f$ and a Newton impact law of restitution coefficient e.

   With \f$ \Delta^2 q = q_{k+1} - 2q_k + q_{k-1} \f$, the step \f$ k \to k+1 \f$ solves

   \f[
   \begin{array}{l}
   M\Delta^2 q
   + h^2 C\left(\theta\,\frac{q_{k+1}-q_{k-1}}{2h} + (1-\theta)\,v_k\right)
   + h^2 K\left(\theta^2 q_{k+1} + (1-\theta^2)\,q_k\right)
   - h^2\left(\theta F_{ext}(t_{k+1}) + (1-\theta)F_{ext}(t_k)\right) = G^T r_{k+1},\\[2pt]
   0 \le g(q_{k+1}) + e\,g(q_{k-1}) \perp r_{k+1} \ge 0,\qquad
   v_{k+1} = \frac{q_{k+1}-q_{k-1}}{2h}.
   \end{array}
   \f]

   In the code:
   - \f$ W = M + \frac{h\theta}{2}C + h^2\theta^2 K \f$ (initializeIterationMatrix());
   - computeResidu() evaluates the left-hand side at \f$ q_{k+1} = q_k \f$ (RESIDU_FREE),
     computeFreeState() gives \f$ q_{free} = q_k - W^{-1}\,\mathrm{RESIDU\_FREE} \f$;
   - the LCP (index set 0: all the constraints) is
     \f$ w = H W^{-1} H^T r + g(q_{free}) + e\,g(q_{k-1}) \f$, where \f$ g(q_{k-1}) \f$ is
     read in the memory of the output, Interaction::yMemory(0) at index 1;
   - updateState(): \f$ q_{k+1} = q_{free} + W^{-1}H^T r \f$ (r is stored in p[0]) and the
     velocity is the centred difference quotient above.

   Two-step start: at \f$ t_0 \f$ the memories hold \f$ (q_0, q_{-1}) \f$ with
   \f$ q_{-1} = q_0 - h v_0 + \frac{h^2}{2}M^{-1}(F_{ext}(t_0) - Kq_0 - Cv_0) \f$
   (initializeWorkVectorsForDS()), so that \f$ q_1 \f$ is a second-order Taylor
   expansion. The memory of the output is filled with \f$ g(q_0) \f$ twice.

   The output memory must keep \f$ g(q_k) \f$ and \f$ g(q_{k-1}) \f$: at the beginning
   of each step, OneStepIntegrator::updateAndSwapAllOutput() only refreshes the latest
   entry (Interaction::updateLatestInMemory()) instead of pushing it again.

   See docs/devel/SCHATZMAN_PAOLI_SCHEME.md for the derivation, the history of the
   fixes and the validation with hots-python.

   Main functions:

   - computeFreeState(): computes xfree (or vfree), dynamical systems
   state without taking non-smooth part into account \n

   - updateState(): computes x (q,v), the complete dynamical systems
   states.

*/
class SchatzmanPaoliOSI : public OneStepIntegrator {
 public:
  enum SchatzmanPaoliOSI_ds_workVector_id { RESIDU_FREE, FREE, LOCAL_BUFFER, WORK_LENGTH };

  enum SchatzmanPaoliOSI_interaction_workVector_id { OSNSP_RHS, WORK_INTERACTION_LENGTH };

  enum SchatzmanPaoliOSI_workBlockVector_id { xfree, BLOCK_WORK_LENGTH };

  static constexpr unsigned int SCHATZMANPAOLISTEPSINMEMORY = 2;

 protected:
  ACCEPT_SERIALIZATION(SchatzmanPaoliOSI);

  /** Stl map that associates a theta parameter for the integration
   *  scheme to each DynamicalSystem of the OSI */
  double _theta{0.};

  /** A gamma parameter for the integration scheme to each DynamicalSystem of
   * the OSI This parameter is used to apply a theta-method to the input $r$
   */
  double _gamma{1.};

  /** a boolean to known is the parameter must be used or not
   */
  bool _useGamma{false};

  /** a boolean to known is the parameter must be used or not
   */
  bool _useGammaForRelation{false};

  /** nslaw effects
   */
  struct _NSLEffectOnFreeOutput;

 public:
  /** constructor from theta value only
   *
   *  \param theta value for all these DS.
   */
  SchatzmanPaoliOSI(double theta);

  /** constructor from theta value only
   *
   *  \param theta value for all these DS.
   *  \param gamma value for all these DS.
   */
  SchatzmanPaoliOSI(double theta, double gamma);

  /** destructor
   */
  virtual ~SchatzmanPaoliOSI() noexcept = default;

  /** get IterationMatrixBoundaryConditions corresponding to DynamicalSystem ds
   *
   *  \param ds a pointer to DynamicalSystem, optional, default =
   *  nullptr. get IterationMatrixBoundaryConditions[0] in that case
   * \return pointer to a SiconosMatrix
   */
  std::shared_ptr<siconos::algebra::SiconosMatrix> IterationMatrixBoundaryConditions(
      std::shared_ptr<siconos::modeling::DynamicalSystem> ds);

  // -- theta --

  /** get theta
   *
   *  \return a double
   */
  inline double theta() { return _theta; };

  /** set the value of theta
   *
   *  \param newTheta a double
   */
  inline void setTheta(double newTheta) { _theta = newTheta; };

  // -- gamma --

  /** get gamma
   *
   *  \return a double
   */
  inline double gamma() { return _gamma; };

  /** set the value of gamma
   *
   *  \param newGamma a double
   */
  inline void setGamma(double newGamma) {
    _gamma = newGamma;
    _useGamma = true;
  };

  // -- useGamma --

  /** get bool useGamma
   *
   *  \return a bool
   */
  inline bool useGamma() { return _useGamma; };

  /** set the boolean to indicate that we use gamma
   *
   *  \param newUseGamma a bool
   */
  inline void setUseGamma(bool newUseGamma) { _useGamma = newUseGamma; };

  /** get bool gammaForRelation for the relation
   *
   *  \return a
   */
  inline bool useGammaForRelation() { return _useGammaForRelation; };

  /** set the boolean to indicate that we use gamma for the relation
   *
   *  \param newUseGammaForRelation a bool
   */
  inline void setUseGammaForRelation(bool newUseGammaForRelation) {
    _useGammaForRelation = newUseGammaForRelation;
    if (_useGammaForRelation) _useGamma = false;
  };

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
  unsigned int numberOfIndexSets() const override { return 1; };

  /** initialize iteration matrix W SchatzmanPaoliOSI matrix at time t
   *
   *  \param time (double)
   *  \param ds a pointer to DynamicalSystem
   */
  void initializeIterationMatrix(double time,
                                 std::shared_ptr<siconos::modeling::DynamicalSystem> ds) const;

  /** return the maximum of all norms for the "SchatzmanPaoliOSI-discretized"
   *  residus of DS \return a double
   */
  double computeResidu() override;

  /** integrates the Dynamical System linked to this integrator
   *  without boring the constraints
   */
  void computeFreeState() override;

  /** integrates the Interaction linked to this integrator, without taking
   *  non-smooth effects into account \param vertex_inter of the interaction
   *  graph \param osnsp pointer to siconos::nonsmooth_formulations::OneStepNSProblem
   */
  void computeFreeOutput(siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
                         siconos::nonsmooth_formulations::OneStepNSProblem* osnsp) override;

  /** return the workVector corresponding to the right hand side of the OneStepNonsmooth
   * problem
   */
  siconos::algebra::SiconosVector& osnsp_rhs(
      siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
      siconos::graphs::InteractionsGraph& indexSet) override {
    return *(*indexSet.properties(vertex_inter).workVectors)[SchatzmanPaoliOSI::OSNSP_RHS];
  };

  void prepareNewtonIteration(double time) override;

  /** integrate the system, between tinit and tend (->iout=true), with possible
   *  stop at tout (->iout=false) \param tinit initial time \param tend end time
   *
   *  \param tout real end time
   *  \param idid useless flag (for SchatzmanPaoliOSI, used in LsodarOSI)
   */
  void integrate(double& tinit, double& tend, double& tout, int& idid) override;

  /** updates the state of the Dynamical Systems
   *
   *  \param level level of interest for the dynamics: not used at the time
   */
  void updateState(const unsigned int level) override {};

  /** compute the current iteration
   *
   */
  void computeIteration() override;

  /** Displays the data of the SchatzmanPaoliOSI's integrator
   */
  void display() const override;
};
}  // namespace siconos::integrators
#endif  // SCHATZMANPAOLIOSI_H
