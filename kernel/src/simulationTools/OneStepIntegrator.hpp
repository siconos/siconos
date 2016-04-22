/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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


/**  Generic object to manage DynamicalSystem(s) time-integration
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 * !!! This is a virtual class, interface for some specific integrators !!!
 *
 * At the time, available integrators are: MoreauJeanOSI, EulerMoreauOSI, LsodarOSI, Hem5OSI
 *
 */
class OneStepIntegrator :public std11::enable_shared_from_this<OneStepIntegrator>
{
protected:
/** serialization hooks
 */
  ACCEPT_SERIALIZATION(OneStepIntegrator);


/** type/name of the Integrator */
  OSI::TYPES _integratorType;

/** a set of DynamicalSystem to integrate */
  SP::DynamicalSystemsSet OSIDynamicalSystems;

/** a graph of dynamical to integrate */
  SP::DynamicalSystemsGraph _dynamicalSystemsGraph; 

/** size of the memory for the integrator */
  unsigned int _sizeMem;

/** A link to the simulation that owns this OSI */
  SP::Simulation _simulation;

/** basic constructor with Id
 *  \param type integrator type/name
 */
  OneStepIntegrator(const OSI::TYPES& type);

/** default constructor
 */
  OneStepIntegrator() {};

  /** struct to add terms in the integration. Useful for Control */
  SP::ExtraAdditionalTerms _extraAdditionalTerms;

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

  inline bool checkOSI(DynamicalSystemsGraph::VIterator dsi)
  {
    return  (_dynamicalSystemsGraph->osi[*dsi].get()) == this;
  };
  inline bool checkOSI(DynamicalSystemsGraph::VDescriptor dsgv)
  {
    return  (_dynamicalSystemsGraph->osi[dsgv].get()) == this;
  };

/** get the set of DynamicalSystem associated with the Integrator
 *  \return a SP::DynamicalSystemsGraph
 */
  inline SP::DynamicalSystemsGraph dynamicalSystemsGraph() const
  {
    return _dynamicalSystemsGraph;
  };

/** get the set of DynamicalSystem associated with the Integrator
 *  \return a DynamicalSystemsSet
 */
  inline SP::DynamicalSystemsSet dynamicalSystems() const
  {
    std::cout << "call to dynamicalSystems()" << std::endl;
    return OSIDynamicalSystems;
  };

/** gets an iterator to the first element of the OSIDynamicalSystems set.
 *  \return a DSIterator.
 */
  inline DSIterator dynamicalSystemsBegin()
  {
    std::cout << "call to dynamicalSystemsBegin()" << std::endl;
    return OSIDynamicalSystems->begin();
  };

/** gets an iterator equal to OSIDynamicalSystems->end().
 *  \return a DSIterator.
 */
  inline DSIterator dynamicalSystemsEnd()
  {
    std::cout << "call to dynamicalSystemsEnd()" << std::endl;
    return OSIDynamicalSystems->end();
  };

/** gets a const iterator to the first element of the OSIDynamicalSystems set.
 *  \return a ConstDSIterator.
 */
  inline ConstDSIterator dynamicalSystemsBegin() const
  {
    std::cout << "call to dynamicalSystemsBegin()" << std::endl;
    return OSIDynamicalSystems->begin();
  };

/** gets a const iterator equal to OSIDynamicalSystems->end().
 *  \return a ConstDSIterator.
 */
  inline ConstDSIterator dynamicalSystemsEnd() const
  {
    std::cout << "call to dynamicalSystemsEnd()" << std::endl;
    return OSIDynamicalSystems->end();
  };

/** get _sizeMem value
 *  \return an unsigned int
 */
  inline unsigned int getSizeMem() const
  {
    return _sizeMem;
  };

/** set _sizeMem
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
   */
  virtual void initialize() = 0;

  /** Save Dynamical Systems data into memory.
   */
  void saveInMemory();

  /** compute the initial state of the Newton loop.
   */
  virtual void computeInitialNewtonState();

  /** return the maximum of all norms for the discretized residus of DS
      \return a double
  */
  virtual double computeResidu();

  /** integrates the Dynamical System linked to this integrator, without taking constraints
   * into account.
   */
  virtual void computeFreeState();

  /** integrates the Interaction linked to this integrator, without taking non-smooth effects into account
   * \param vertex_inter of the interaction graph
   * \param osnsp pointer to OneStepNSProblem
   */
  virtual void computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp);

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
   */
  virtual void updateState(const unsigned int level) = 0;

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
