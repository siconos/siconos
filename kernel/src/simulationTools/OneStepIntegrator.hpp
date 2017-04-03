/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

  Basic class to handle with dynamical system integrators over a time step.
*/

#ifndef ONESTEPINTEGRATOR_H
#define ONESTEPINTEGRATOR_H

#include "SiconosConst.hpp"
#include "SimulationTypeDef.hpp"
#include "DynamicalSystemsSet.hpp"
#include "OneStepIntegratorTypes.hpp"
#include "SimulationGraphs.hpp"

// work around gccxml bug that does not accept an argument ot the deprecated attribute
#ifndef __GCCXML__

#ifdef __GNUC__
#define DEPRECATED_OSI_API(func) func __attribute__ ((deprecated ("This constructor of function is deprecrated and will be removed in the next major Siconos release! Use the insertDynamicalSystem method in NonSmoothDynamicalSystem !")))
#elif defined(_MSC_VER)
#define DEPRECATED_OSI_API(func) __declspec(deprecated("This constructor will be removed in the next major Siconos release and does not work with MSVC 2013 ! Use the insertDynamicalSystem method in NonSmoothDynamicalSystem !")) func
#else
#define DEPRECATED_OSI_API(func) func
#endif

#else
#define DEPRECATED_OSI_API(func) func
#endif

/**  Generic class to manage DynamicalSystem(s) time-integration
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 * !!! This is a virtual class, interface for some specific integrators !!!
 *
 * At the time, available integrators are:
 * <ul>
 * <li> EulerMoreauOSI </li>
 * <li> MoreauJeanOSI </li>
 * <li> MoreauJeanCombinedProjectionOSI </li>
 * <li> MoreauJeanDirectProjectionOSI </li>
 * <li> D1MinusLinearOSI </li>
 * <li> D1MinusLinearOSIHalfExplicitAccelerationLevelOSI </li>
 * <li> D1MinusLinearOSIHalfExplicitVelocityLevelOSI </li>
 * <li> SchatzmanPaoliOSI </li>
 * <li> LsodarOSI </li>
 * <li> Hem5OSI </li>
 * <li> NewMarkAlphaOSI </li>
 * <li> ZeroOrderHoldOSI </li>
 * </ul>
 *
 */
class OneStepIntegrator :public std11::enable_shared_from_this<OneStepIntegrator>
{

public :
  /** List of indices used to save tmp work matrices and vectors (last one is the size of the present list) */
  enum OSI_DSWorkVectorId {local_buffer, residu, residu_free, free, free_tdg, x_partial_ns, delta_x_for_relation,
                           qtmp, acce_memory, acce_like, work_vector_of_vector_size};

  enum OSI_DSWorkMatrixId {dense_output_coefficients, work_vector_of_matrix_size};

protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(OneStepIntegrator);


  /** type/name of the Integrator */
  OSI::TYPES _integratorType;

  /** a graph of dynamical systems to integrate
   * For the moment, we point to the graph of dynamical systems in
   * in the topology. We use the properties "osi" to check if the dynamical
   * system is integrated by this osi. It has to be improved by using a subgraph
   * to avoid the use of checkOSI
   */
  SP::DynamicalSystemsGraph _dynamicalSystemsGraph;

  /** size of the memory for the integrator */
  unsigned int _sizeMem;

  /** steps of the integrator
   */
  unsigned int _steps;
  
  /** _levelMinForOutput is the minimum level for the output
   * needed by the OneStepIntegrator
   */
  unsigned int _levelMinForOutput;


  /** _levelMaxForOutput is the maximum level for the output
   * needed by the OneStepIntegrator
   */
  unsigned int _levelMaxForOutput;


  /** _levelMinForInput is the minimum level for the input
   * needed by the OneStepIntegrator
   */
  unsigned int _levelMinForInput;


  /** _levelMaxForInput is the maximum level for the input
   * needed by the OneStepIntegrator
   */
  unsigned int _levelMaxForInput;

/** A link to the simulation that owns this OSI */
  SP::Simulation _simulation;

  /** basic constructor with Id
   *  \param type integrator type/name
   */
  OneStepIntegrator(const OSI::TYPES& type): _integratorType(type), _sizeMem(1) {};

  /** default constructor
   */
  OneStepIntegrator() {};

  /** struct to add terms in the integration. Useful for Control */
  SP::ExtraAdditionalTerms _extraAdditionalTerms;

  /** Compare interaction and current OSI levels for input and output.
      Reset interaction if they are not compliant.
      \param inter a reference to an Interaction
  */
  void _check_and_update_interaction_levels(Interaction& inter);

  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi -- common code.
   * \param ds the dynamical system
   */
  SP::VectorOfVectors initializeDynamicalSystemWorkVectors(SP::DynamicalSystem ds);

private:


  /** copy constructor, private, no copy nor pass-by value allowed */
  OneStepIntegrator(const OneStepIntegrator&);

  /** assignment (private => forbidden)
   * \param  OSI
   * \return OneStepIntegrator&
   */
  OneStepIntegrator& operator=(const OneStepIntegrator& OSI);


public:

  /** destructor
   */
  virtual ~OneStepIntegrator() {};

  // --- GETTERS/SETTERS ---

  /** get the type of the OneStepIntegrator
   *  \return std::string : the type of the OneStepIntegrator
   */
  inline OSI::TYPES getType() const
  {
    return _integratorType;
  }

  /** set the type of the OneStepIntegrator
   *  \param newType std::string : the type of the OneStepIntegrator
   */
  inline void setType(const OSI::TYPES& newType)
  {
    _integratorType = newType;
  };

  /** Check if the dynamical system bundle in the node of the
   * _dynamicalSystemGraph is integrated by this osi.
   * \param dsi the iterator on the node of the graph
   */
  inline bool checkOSI(DynamicalSystemsGraph::VIterator dsi)
  {
    return  (_dynamicalSystemsGraph->properties(*dsi).osi.get()) == this;
  };

  /** Check if the dynamical system bundle in the node of the
   * _dynamicalSystemGraph is integrated by this osi.
   * \param dsgv the descriptor on the node of the graph
   */
  inline bool checkOSI(DynamicalSystemsGraph::VDescriptor dsgv)
  {
    return  (_dynamicalSystemsGraph->properties(dsgv).osi.get()) == this;
  };

  /** get the set of DynamicalSystem associated with the Integrator
   *  \return a SP::DynamicalSystemsGraph
   */
  inline SP::DynamicalSystemsGraph dynamicalSystemsGraph() const
  {
    return _dynamicalSystemsGraph;
  };

  /** get number of internal memory vectors needed for the integration scheme
   *  \return an unsigned int
   */
  inline unsigned int getSizeMem() const
  {
    return _sizeMem;
  };

  /** set the number of internal memory vectors needed for the integration scheme
   *  \param newValue an unsigned int
   */
  inline void setSizeMem(unsigned int newValue)
  {
    _sizeMem = newValue;
  };

  /** get the Simulation that owns the OneStepIntegrator
   *  \return a pointer to Simulation
   */
  inline SP::Simulation simulation() const
  {
    return _simulation;
  }

  /** set the Simulation of the OneStepIntegrator
   *  \param newS a pointer to Simulation
   */
  inline void setSimulationPtr(SP::Simulation newS)
  {
    _simulation = newS;
  }

  // --- OTHERS ... ---

  /** initialise the integrator
   * \param m a Model
   */
  virtual void initialize(Model& m );

  /** Initialization process of the nonsmooth problems
   linked to this OSI*/
  virtual void initialize_nonsmooth_problems(){};

  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi
   * \param m the Model
   * \param t time of initialization
   * \param ds the dynamical system
   */
  virtual void initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds) = 0 ;

  /** initialization of the work vectors and matrices (properties) related to
   *  one interaction on the graph and needed by the osi
   * \param inter the interaction
   * \param interProp the properties on the graph
   * \param DSG the dynamical systems graph
   */
  virtual void fill_ds_links(Interaction &inter,
			     InteractionProperties& interProp,
			     DynamicalSystemsGraph & DSG) = 0 ;

  /** compute interaction output (y) for all levels and swaps in memory
      \param inter the interaction to update
      \param time value for output computation
      \param interaction_properties properties of the interaction, in the Interaction Graph I0
   */
  void update_interaction_output(Interaction& inter, double time, InteractionProperties& interaction_properties);
  
  /** get the number of index sets required for the simulation
   * \return unsigned int
   */
  virtual unsigned int numberOfIndexSets() const=0;

  virtual unsigned int levelMinForOutput()
  {
    return _levelMinForOutput;
  }
  virtual unsigned int levelMaxForOutput()
  {
    return _levelMaxForOutput;
  }
  virtual unsigned int levelMinForInput()
  {
    return _levelMinForInput;
  }
  virtual unsigned int levelMaxForInput()
  {
    return _levelMaxForInput;
  }

  /** compute the initial state of the Newton loop.
   */
  virtual void computeInitialNewtonState(){
    // Default behavior :  does nothing and used the current state as starting state of the Newton iteration
  }

  /** return the maximum of all norms for the discretized residus of DS
   *  \return a double
   */
  virtual double computeResidu(){
    // default : error
    RuntimeException::selfThrow("OneStepIntegrator::computeResidu not implemented for integrator of type " + _integratorType);
    return 0.0;
  }

  /** integrates the Dynamical System linked to this integrator, without taking constraints
   * into account.
   */
  virtual void computeFreeState(){
    // default : error
    RuntimeException::selfThrow("OneStepIntegrator::computeFreeState not implemented for integrator of type " + _integratorType);
  }


  /** integrates the Interaction linked to this integrator, without taking non-smooth effects into account
   * \param vertex_inter of the interaction graph
   * \param osnsp pointer to OneStepNSProblem
   */
  virtual void computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
  {
    // default : error
    RuntimeException::selfThrow("OneStepIntegrator::computeFreeOutput not implemented for integrator of type " + _integratorType);
  }

  /** integrate the system, between tinit and tend (->iout=true), with possible stop at tout (->iout=false)
   *  \param tinit initial time
   *  \param tend end time
   *  \param tout real end time
   *  \param idid flag used in EventDriven schemes
   */
  virtual void integrate(double& tinit, double& tend, double& tout, int& idid) = 0;

  /** set to zero all the r vectors of the DynamicalSystems of the present OSI
   */
  void resetNonSmoothPart();

  /** set to zero all the r vectors of the DynamicalSystems of the present OSI for a given level
   * \param level
   */
  void resetNonSmoothPart(unsigned int level);

  /** update the state of the DynamicalSystem attached to this Integrator
   *  \param level level of interest for the dynamics
   * level is set to 0 by default since in all time-stepping schemes we update all the state
   * whatever the value of level
   */
  virtual void updateState(const unsigned int level) = 0;

  /** update the state of the DynamicalSystem attached to this Integrator
   * level is set to 0 by default since in all time-stepping schemes we update all the state
   * whatever the value of level
   */
  void updateState() { updateState(0); }

  /** update the output of the Interaction attached to this Integrator
   */
  void updateOutput(double time);

  /** update the input of the Interaction attached to this Integrator
   */
  void updateInput(double time);

  /** update the output of the Interaction attached to this Integrator
   *  \param level level of interest for the dynamics
   */
  void updateOutput(double time, unsigned int level);

  /** update the input of the Interaction attached to this Integrator
   *  \param level level of interest for the dynamics
   */
  void updateInput(double time, unsigned int level);

  /** print the data to the screen
   */
  virtual void display() = 0;

  /**
   */
  virtual void prepareNewtonIteration(double time) = 0;

  /** Apply the rule to one Interaction to known if is it should be included
   * in the IndexSet of level i
   * \param inter
   * \param i
   * \return bool
   */
  virtual bool addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
  {
    RuntimeException::selfThrow("OneStepIntegrator::addInteractionInIndexSet - Should be called at this level");
    return 0;
  }
  ;

  /** Apply the rule to one Interaction to know if is it should be removed
   * from the IndexSet of level i
   * \param inter
   * \param i
   * \return bool
   */
  virtual bool removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)
  {
    RuntimeException::selfThrow("OneStepIntegrator::removeInteractionInIndexSet - Should not be called at this level");
    return 0;
  };

  /** get the ExtraAdditionalTerms.
   * \return the ExtraAdditionalTerms
   */
  inline SP::ExtraAdditionalTerms extraAdditionalTerms()
  {
    return _extraAdditionalTerms;
  }

  /** set the ExtraAdditionalTerms to add smooth terms for the integration process.
   * Useful when a control loop is added to a DynamicalSystem.
   * \param eat the ExtraAdditionalTerms to use
   */
  inline void setExtraAdditionalTerms(SP::ExtraAdditionalTerms eat)
  {
    _extraAdditionalTerms = eat;
  }

  /** visitors hook
   */
  VIRTUAL_ACCEPT_VISITORS(OneStepIntegrator);

};

#endif // ONESTEPINTEGRATOR_H
