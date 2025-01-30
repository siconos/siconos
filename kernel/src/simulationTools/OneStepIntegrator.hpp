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
/*! \file OneStepIntegrator.hpp

  Base class (i.e. common interface) for dynamical system integration over a time step.
*/

#ifndef ONESTEPINTEGRATOR_H
#define ONESTEPINTEGRATOR_H

#include "ExtraAdditionalTerms.hpp"
#include "OneStepIntegratorTypes.hpp"  // IntegratorType
#include "SimulationGraphs.hpp"

namespace siconos::simulation {
class Simulation;
}  // namespace siconos::simulation

namespace siconos::nonsmooth_formulations {
class OneStepNSProblem;
}

namespace siconos::integrators {

/**
   Generic class to manage DynamicalSystem(s) time-integration

   This is a virtual class, interface for some specific integrators.

   At the time, available integrators are:

   - EulerMoreauOSI
   - MoreauJeanOSI
   - MoreauJeanCombinedProjectionOSI
   - MoreauJeanDirectProjectionOSI
   - MoreauJeanBilbaoOSI
   - D1MinusLinearOSI
   - SchatzmanPaoliOSI
   - LsodarOSI
   - Hem5OSI
   - NewMarkAlphaOSI
   - ZeroOrderHoldOSI
*/
class OneStepIntegrator : public std::enable_shared_from_this<OneStepIntegrator> {
 protected:
  ACCEPT_SERIALIZATION(OneStepIntegrator);

  /** type/name of the Integrator */
  IntegratorType _integratorType;

  /** a graph of dynamical systems to integrate
   *  For the moment, we point to the graph of dynamical systems in
   *  the topology. We use the properties "osi" to check if the dynamical
   *  system is integrated by this osi. It has to be improved by using a subgraph
   *  to avoid the use of checkOSI
   */
  std::shared_ptr<siconos::graphs::DynamicalSystemsGraph> _dynamicalSystemsGraph{nullptr};

  /** size of the memory for the integrator */
  unsigned int _sizeMem{1};

  /** steps of the integrator */
  unsigned int _steps{0};

  /** _levelMinForOutput is the minimum level for the output
   *  needed by the OneStepIntegrator
   */
  unsigned int _levelMinForOutput{0};

  /** _levelMaxForOutput is the maximum level for the output
   *  needed by the OneStepIntegrator
   */
  unsigned int _levelMaxForOutput{0};

  /** _levelMinForInput is the minimum level for the input
   *  needed by the OneStepIntegrator
   */
  unsigned int _levelMinForInput{0};

  /** _levelMaxForInput is the maximum level for the input
   *  needed by the OneStepIntegrator
   */
  unsigned int _levelMaxForInput{0};

  bool _isInitialized{false};

  /** boolean variable to force an explicit evaluation of the Jacobians
   *  mapping of relations only at the beginning of the time--step and
   *  not in the Newton iteration
   */
  bool _explicitJacobiansOfRelation{false};

  /** A link to the simulation that owns this OSI */
  std::shared_ptr<siconos::simulation::Simulation> _simulation{nullptr};

  /** struct to add terms in the integration. Useful for Control */
  std::shared_ptr<ExtraAdditionalTerms> _extraAdditionalTerms{nullptr};

  /** basic constructor with OSI Id and default values for the other variables
   *
   *  \param type integrator type/name
   *  \param steps number of steps of the integrator
   *  \param lmin_output minimum level for the output
   *  \param lmax_output minimum level for the output
   *  \param lmin_input minimum level for the input
   *  \param lmax_input minimum level for the input
   */
  OneStepIntegrator(const IntegratorType& type, unsigned int steps, unsigned int lmin_output,
                    unsigned int lmax_output, unsigned int lmin_input, unsigned int lmax_input)
      : _integratorType(type),
        _steps{steps},
        _levelMinForOutput{lmin_output},
        _levelMaxForOutput{lmax_output},
        _levelMinForInput{lmin_input},
        _levelMaxForInput{lmax_input} {
            // Set levels. This may depend on the nonsmooth law and will be updated during
            // initializeWorkVectorsForInteraction(...) call.
        };

  /**
     Compare interaction and current OSI levels for input and output.
     Reset interaction if they are not compliant.

     \param inter a reference to an Interaction
  */
  void _check_and_update_interaction_levels(siconos::modeling::Interaction& inter);

  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi -- common code.
   *
   *  \param ds the dynamical system
   */
  std::shared_ptr<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>
  _initializeDSWorkVectors(std::shared_ptr<siconos::modeling::DynamicalSystem> ds);

 private:
  // Rule of five
  OneStepIntegrator() = delete;

  OneStepIntegrator(const OneStepIntegrator&) = delete;
  OneStepIntegrator(OneStepIntegrator&&) = delete;
  OneStepIntegrator& operator=(const OneStepIntegrator&) = delete;
  OneStepIntegrator& operator=(OneStepIntegrator&&) = delete;

 public:
  /** destructor
   */
  virtual ~OneStepIntegrator() noexcept = default;

  /** \return the id of the integrator (see list in IntegratorType enum)
   */
  inline IntegratorType getType() const { return _integratorType; }

  /**
     \return a pointer to the graph of dynamical systems associated with the Integrator
     warning: returns the whole ds graph, not only ds integrated by the present osi.
  */
  inline std::shared_ptr<siconos::graphs::DynamicalSystemsGraph> dynamicalSystemsGraph()
      const {
    return _dynamicalSystemsGraph;
  };

  /**
     set the graph of dynamical systems associated with the Integrator
  */
  inline void setDynamicalSystemsGraph(
      std::shared_ptr<siconos::graphs::DynamicalSystemsGraph> dsg) {
    _dynamicalSystemsGraph = dsg;
  };

  /** \return the number of internal memory vectors needed in dynamical systems integrated with
   * this osi.
   */
  inline unsigned int getSizeMem() const { return _sizeMem; };

  /** get the Simulation that owns the OneStepIntegrator (pointer link)
   *
   *  \return a pointer to Simulation
   */
  inline std::shared_ptr<siconos::simulation::Simulation> simulation() const {
    return _simulation;
  }

  /** set the Simulation of the OneStepIntegrator
   *
   *  \param newS a pointer to Simulation
   */
  inline void setSimulationPtr(std::shared_ptr<siconos::simulation::Simulation> newS) {
    _simulation = newS;
  }

  /**
     minimal level required for output var used with this integration scheme.
     var[level] is the derivative of order 'level' of var.
  */
  virtual unsigned int levelMinForOutput() { return _levelMinForOutput; }

  /**
      maximal level required for output var used with this integration scheme.
      var[level] is the derivative of order 'level' of var.
   */
  virtual unsigned int levelMaxForOutput() { return _levelMaxForOutput; }

  /**
      minimal level required for input var used with this integration scheme.
      var[level] is the derivative of order 'level' of var.
   */
  virtual unsigned int levelMinForInput() { return _levelMinForInput; }

  /**
      maximal level required for input var used with this integration scheme.
      var[level] is the derivative of order 'level' of var.
   */
  virtual unsigned int levelMaxForInput() { return _levelMaxForInput; }

  /** get the number of index sets required for the simulation
   *
   *  \return unsigned int
   */
  virtual unsigned int numberOfIndexSets() const = 0;

  inline bool isInitialized() { return _isInitialized; };

  inline void setIsInitialized(bool value) { _isInitialized = value; };

  bool explicitJacobiansOfRelation() { return _explicitJacobiansOfRelation; }

  void setExplicitJacobiansOfRelation(bool newval) { _explicitJacobiansOfRelation = newval; };

  /** initialise the integrator
   */
  virtual void initialize();

  /**
     Initialization process of the nonsmooth problems
     linked to this OSI*/
  virtual void initialize_nonsmooth_problems(){};

  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi
   *
   *  \param t time of initialization
   *  \param ds the dynamical system
   */
  virtual void initializeWorkVectorsForDS(
      double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) = 0;

  /** initialization of the work vectors and matrices (properties) related to
   *  one interaction on the graph and needed by the osi
   *
   *  \param inter the interaction
   *  \param interProp the properties on the graph
   *  \param DSG the dynamical systems graph
   */
  virtual void initializeWorkVectorsForInteraction(
      siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
      siconos::graphs::DynamicalSystemsGraph& DSG) = 0;

  /**
     compute interaction output (y) for all levels and swaps in memory
     \param time value for output computation
     \param interaction_properties properties of the interaction, in the Interaction Graph I0
  */
  void updateAndSwapAllOutput(double time);

  /**
     compute interaction output (y) for all levels and swaps in memory

     \param inter the interaction to update
     \param time value for output computation
  */
  void updateAndSwapAllOutput(siconos::modeling::Interaction& inter, double time);

  /** compute the initial state (for dynamical system variables) of the Newton loop. */
  virtual void computeInitialNewtonState() {
    // Default behavior :  does nothing and used the current state as starting state of the
    // Newton iteration
  }

  /** return the maximum of all norms for the discretized residus of DS
   *
   *  \return a double
   */
  virtual double computeResidu();

  /**
      integrates the Dynamical System linked to this integrator, without taking constraints
      into account.
  */
  virtual void computeFreeState();

  /** integrates the Interaction linked to this integrator, without taking non-smooth effects
   * into account
   *
   *  \param vertex_inter of the interaction graph
   *  \param osnsp pointer to siconos::nonsmooth_formulations::OneStepNSProblem
   */
  virtual void computeFreeOutput(siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
                                 siconos::nonsmooth_formulations::OneStepNSProblem* osnsp);

  /** compute the residu of the output of the relation (y)
   *  This computation depends on the type of OSI
   *
   *  \param time time of computation
   *  \param indexSet the index set of the interaction that are concerned
   */
  virtual double computeResiduOutput(
      double time, std::shared_ptr<siconos::graphs::InteractionsGraph> indexSet);
  /** compute the residu of the input of the relation (R or p)
   *  This computation depends on the type of OSI
   *
   *  \param time time of computation
   *  \param indexSet the index set of the interaction that are concerned
   */
  virtual double computeResiduInput(
      double time, std::shared_ptr<siconos::graphs::InteractionsGraph> indexSet);

  /** integrate the system, between tinit and tend, with possible stop at tout
   *
   *  \param tinit start time
   *  \param tend expected end time
   *  \param tout real end time
   *  \param idid extra flag, meaningful only for OSI used in EventDriven schemes
   */
  virtual void integrate(double& tinit, double& tend, double& tout, int& idid) = 0;

  /** set to zero all the r vectors of the DynamicalSystems integrated by this OSI
   */
  void resetAllNonSmoothParts();

  /** set to zero all the r vectors of the DynamicalSystems of the present OSI for a given
   *  level
   *
   *  \param level
   */
  void resetNonSmoothPart(unsigned int level);

  /** \return the iteration matrix corresponding to a given dynamical system
   *
   *  \param ds a dynamical system
   */
  virtual std::shared_ptr<siconos::algebra::SiconosMatrix> iterationMatrix(
      std::shared_ptr<siconos::modeling::DynamicalSystem> ds);

  /** \return the LU factorization of the iteration matrix corresponding to a given dynamical
   * system
   *
   *  \param ds a dynamical system
   */
  virtual std::shared_ptr<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>> LUiterationMatrix(
      std::shared_ptr<siconos::modeling::DynamicalSystem> ds);

  /** update the state of the DynamicalSystem attached to this Integrator
   *
   *  \param level level of interest for the dynamics
   *  level is set to 0 by default since in all time-stepping schemes we update all the state
   *  whatever the value of level
   */
  virtual void updateState(const unsigned int level) = 0;

  /** update the state of the DynamicalSystem attached to this Integrator
   *  level is set to 0 by default since in all time-stepping schemes we update all the state
   *  whatever the value of level
   */
  void updateState() { updateState(0); }

  /** update the output of the Interaction attached to this Integrator
   */
  virtual void updateOutput(double time);

  /** update the input of the Interaction attached to this Integrator
   */
  virtual void updateInput(double time);

  /** update the output of the Interaction attached to this Integrator
   *
   *  \param time current time
   *  \param level level of interest for the dynamics
   */
  virtual void updateOutput(double time, unsigned int level);

  /** update the input of the Interaction attached to this Integrator
   *
   *  \param time current time
   *  \param level level of interest for the dynamics
   *  \warning VA: 27/10/2022 Whatever the level, the updateInput method loops over indexSet0
   *  This is sometimes necessary for some OSI but for some others it may burden the
   *  computational time for nothing. For instance, in standard MoreauJEANOSI, p[1] is only
   *  defined on indexSet1. we should go towards virtual void updateInput(double time, unsigned
   *  int pLevel, unsigned int indexSetLevel );
   */
  virtual void updateInput(double time, unsigned int level);

  virtual void prepareNewtonIteration(double time) = 0;

  /** print the data to the screen */
  virtual void display() const = 0;

  /** Apply the rule to one Interaction to known if is it should be included
   *  in the IndexSet of level i
   *
   *  \param inter
   *  \param i
   *  \return bool
   */
  virtual bool addInteractionInIndexSet(std::shared_ptr<siconos::modeling::Interaction> inter,
                                        unsigned int i);

  /** Apply the rule to one Interaction to know if is it should be removed
   *  from the IndexSet of level i
   *
   *  \param inter
   *  \param i
   *  \return bool
   */
  virtual bool removeInteractionFromIndexSet(
      std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i);

  /** get the ExtraAdditionalTerms.
   *
   *  \return the ExtraAdditionalTerms
   */
  inline std::shared_ptr<ExtraAdditionalTerms> extraAdditionalTerms() {
    return _extraAdditionalTerms;
  }

  /** set the ExtraAdditionalTerms to add smooth terms for the integration process.
   *  Useful when a control loop is added to a DynamicalSystem.
   *
   *  \param eat the ExtraAdditionalTerms to use
   */
  inline void setExtraAdditionalTerms(std::shared_ptr<ExtraAdditionalTerms> eat) {
    _extraAdditionalTerms = eat;
  }

  /**
      True if the dynamical system (a vertex in the ds graph) is integrated by this osi.

      \param dsi the iterator on the node of the graph corresponding to the dynamical system
      of interest.
  */
  inline bool checkOSI(siconos::graphs::DynamicalSystemsGraph::VIterator dsi) const {
    return (_dynamicalSystemsGraph->properties(*dsi).osi.get()) == this;
  };

  /**
      True if the dynamical system (a vertex in the ds graph) is integrated by this osi.

      \param dsgv the descriptor of the node in the graph corresponding to the dynamical
      system of interest.
  */
  inline bool checkOSI(siconos::graphs::DynamicalSystemsGraph::VDescriptor dsgv) const {
    return (_dynamicalSystemsGraph->properties(dsgv).osi.get()) == this;
  };

  /**
     True if the dynamical system (a vertex in the ds graph) is integrated by this osi.

     \param dsi the iterator on the node of the graph corresponding to the dynamical system of
     interest.
   */
  inline bool checkInteractionOSI(siconos::graphs::InteractionsGraph& indexSet0,
                                  siconos::graphs::InteractionsGraph::VIterator ui) {
    return (indexSet0.properties(*ui).osi1.get()) == this;
  };

  virtual siconos::algebra::SiconosVector& osnsp_rhs(
      siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
      siconos::graphs::InteractionsGraph& indexSet) = 0;
};
}  // namespace siconos::integrators

#endif  // ONESTEPINTEGRATOR_H
