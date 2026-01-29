/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/*! \file
  SchatzmanPaoliOSI Time-Integrator for Dynamical Systems
*/

#ifndef SCHATZMANPAOLIOSI_H
#define SCHATZMANPAOLIOSI_H

#include "OneStepIntegrator.hpp"

namespace siconos::integrators {

/**
   SchatzmanPaoliOSI Time-Integrator for Dynamical Systems

   SchatzmanPaoliOSI class is used to define some time-integrators methods for a
   list of dynamical systems.

   A SchatzmanPaoliOSI instance is defined by the value of theta and the list of
   concerned dynamical systems.  Each DynamicalSystem is associated to
   a SiconosMatrix, named "W"

   W matrices are initialized and computed in initializeIterationMatrix and
   computeIterationMatrix. Depending on the DS type, they may depend on time and DS
   state (x).

   For Lagrangian systems, the implementation uses p[0] for storing the
   discrete multiplier.

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
  void updateState(const unsigned int level) override;

  /** Displays the data of the SchatzmanPaoliOSI's integrator
   */
  void display() const override;
};
}  // namespace siconos::integrators
#endif  // SCHATZMANPAOLIOSI_H
