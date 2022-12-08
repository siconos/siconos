/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

/*! \file  MoreauJeanBilbaoOSI.hpp */

#ifndef MoreauJeanBilbaoOSI_H
#define MoreauJeanBilbaoOSI_H

#include "OneStepIntegrator.hpp"

#include <limits>

/**
   One-step integrator for event-capturing simulation combining Moreau-Jean and
   Bilbao numerical scheme.

   Numerical scheme which combines an exact method (Bilbao) for the linear
   (non-contacting) part of the equations of motion with a Moreau-Jean
   time-stepping for the nonsmooth part.

   Check details in users'guide.


*/
class MoreauJeanBilbaoOSI : public OneStepIntegrator {
protected:

  ACCEPT_SERIALIZATION(MoreauJeanBilbaoOSI);

  /** nslaw effects */
  struct _NSLEffectOnFreeOutput;

  friend struct _NSLEffectOnFreeOutput;

public:
  enum MoreauJeanBilbaoOSI_ds_workVector_id {
    TWO_DT_SIGMA_STAR,
    ONE_MINUS_THETA,
    VFREE,
    WORK_LENGTH
  };

  enum MoreauJeanBilbaoOSI_interaction_workVector_id {
    OSNSP_RHS,
    WORK_INTERACTION_LENGTH
  };

  enum MoreauJeanBilbaoOSI_interaction_workBlockVector_id {
    xfree,
    BLOCK_WORK_LENGTH
  };

  /** Constructor - No extra parameters: depends only on connected ds and
   *  simulation time step
   */
  MoreauJeanBilbaoOSI();

  /** destructor */
  virtual ~MoreauJeanBilbaoOSI(){};

  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi
   *
   *  \param t time of initialization
   *  \param ds the dynamical system
   */
  void initializeWorkVectorsForDS(double t, SP::DynamicalSystem ds) override;

  /** initialization of the work vectors and matrices (properties) related to
   *  one interaction on the graph and needed by the osi
   *
   *  \param inter the interaction
   *  \param interProp the properties on the graph
   *  \param DSG the dynamical systems graph
   */
  void initializeWorkVectorsForInteraction(Interaction &inter,
                                           InteractionProperties &interProp,
                                           DynamicalSystemsGraph &DSG) override;

  void _initialize_iteration_matrix(SP::DynamicalSystem ds);

  /** Initialization process of the nonsmooth problems
      linked to this OSI*/
  void initialize_nonsmooth_problems() override;

  /** get the number of index sets required for the simulation
   *
   *  \return unsigned int
   */
  unsigned int numberOfIndexSets() const override { return 2; };

  void compute_parameters(double time_step, double omega, double sigma,
                          double &theta, double &sigma_star);

  /** get iteration_matrix (pointer link) corresponding to DynamicalSystem ds
   *
   *  \param ds a pointer to DynamicalSystem
   *  \return pointer to a SiconosMatrix
   */
  inline SP::SimpleMatrix iteration_matrix(SP::DynamicalSystem ds)
  {
    return _dynamicalSystemsGraph
        ->properties(_dynamicalSystemsGraph->descriptor(ds))
        .W;
  }

  /** integrate the system, between tinit and tend (->iout=true), with possible
   *  stop at tout (->iout=false)
   *
   *  \param tinit initial time
   *  \param tend end time
   *  \param tout real end time
   *  \param idid flag used in EventDriven schemes
   */
  void integrate(double &tinit, double &tend, double &tout, int &idid) override;

  /** return the maximum of all norms for the "MoreauJeanOSI-discretized"
      residus of DS

      \return a double
   */
  double computeResidu() override;

  /** Perform the integration of the dynamical systems linked to this integrator
   *  without taking into account the nonsmooth input ( _p)
   */
  void computeFreeState() override;

  /** integrates the Interaction linked to this integrator, without taking
   *  non-smooth effects into account
   *
   *  \param vertex_inter vertex of the interaction graph
   *  \param osnsp pointer to OneStepNSProblem
   */
  void computeFreeOutput(InteractionsGraph::VDescriptor &vertex_inter,
                         OneStepNSProblem *osnsp) override;

  /** return the workVector corresponding to the right hand side of the OneStepNonsmooth problem
   */
  SiconosVector& osnsp_rhs(InteractionsGraph::VDescriptor& vertex_inter, InteractionsGraph& indexSet) override
  {
    return *(*indexSet.properties(vertex_inter).workVectors)[MoreauJeanBilbaoOSI::OSNSP_RHS];
  };

  /** update the state of the dynamical systems
      \param ds the dynamical to update
   */
  void updatePosition(DynamicalSystem &ds);

  /** update the state of the DynamicalSystem attached to this Integrator
   *
   *  \param level level of interest for the dynamics
   *  level is set to 0 by default since in all time-stepping schemes we update
   *  all the state whatever the value of level is.
   */
  void updateState(const unsigned int level) override;

  /** print the data to the screen
   */
  void display() override;

  void prepareNewtonIteration(double time) override;

  /** Apply the rule to one Interaction to know if it should be included in the
   *  IndexSet of level i
   *
   *  \param inter the Interaction to test
   *  \param i level of the IndexSet
   *  \return Boolean
   */
  bool addInteractionInIndexSet(SP::Interaction inter, unsigned int i) override;

  /** Apply the rule to one Interaction to know if it should be removed from the
   *  IndexSet of level i
   *
   *  \param inter the Interaction to test
   *  \param i level of the IndexSet
   *  \return Boolean
   */
  bool removeInteractionFromIndexSet(SP::Interaction inter,
                                     unsigned int i) override;

  ACCEPT_STD_VISITORS();
};

#endif // MoreauJeanBilbaoOSI_H
