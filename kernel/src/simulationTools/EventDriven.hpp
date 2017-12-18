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
  Event Driven Simulation
  */
#ifndef EventDriven_H
#define EventDriven_H

#include "Simulation.hpp"
#include "SiconosFwd.hpp"               // for OneStepIntegrator, etc

/** Simulation based on event driven method, ie events detection (see theoretical manual for more details).
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 * WARNING: at the time only written for Lagrangian systems !!!
 *
 */
class EventDriven : public Simulation
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(EventDriven);


  /** flag used in DLSODAR -
   *  As input: 1 if first call, else 2
   *  As output: 2 if no root found, else 3
   */
  int _istate;

  void initOSNS();

  /** Initialize OneStepIntergrators */
  void initOSIs();

  /** Initialize the Rhs of the OSI */
  void initOSIRhs();

protected:
  /** an epsilon to define the contraint g for Interaction in IndexSet[1]
   */
  static double _TOL_ED;

  /** boolean variable to known whether Newton iterations converges or not */
  bool _isNewtonConverge;

  /** Default maximum number of Newton iteration*/
  unsigned int _newtonMaxIteration;

  /** Number of steps performed is the Newton Loop */
  unsigned int _newtonNbIterations;

  /** Maximum Residual for the Dynamical system */
  double _newtonResiduDSMax;

  /** Maximum Residual for the output of the relation */
  double _newtonResiduYMax;

  /** Tolerance to check Newton iteration convergence */
  double _newtonTolerance;

  /** Maximum number of iterations to localize events */
  unsigned int _localizeEventMaxIter;

  /**  number of OneStepNSProblems considered in the simulation */
  unsigned int _numberOfOneStepNSproblems;

  /** store the indexSet0 for performance reason (Lsodar)*/
  SP::InteractionsGraph _indexSet0;

  /** store the graph of DynamicalSystems for performance reason (Lsodar)*/
  SP::DynamicalSystemsGraph _DSG0;

public:
  /** defaut constructor
   *  \param td time discretisation
   */
  EventDriven(SP::TimeDiscretisation td);

  /** constructor with data
   *  \param td time discretisation
   *  \param nb number of NSProblem
   */
  EventDriven(SP::TimeDiscretisation td, int nb);

  /** defaut constructor (needed for serialization)
   */
  EventDriven() : _isNewtonConverge(false) {};

  /** destructor
  */
  ~EventDriven() {};

  /** Overload Simulation::initialize
   * \param m the Model
   * \param withOSI
   */
  void initialize(SP::Model m, bool withOSI = true);

  /* Getters and setters */

  
  /** Set value to _istate 
   * \param newValue
   */
  inline void setIstate(int newValue)
  {
    _istate = newValue;
  }
  /** Get value of _istate
   * \return _istate a double
   */
  inline double istate()
  {
    return _istate;
  }
  
  /** Set value to _TOL_ED
   * \param var the new tolerance
   */
  inline void setToleranceED(double var)
  {
    _TOL_ED = var;
  }
  /** Get value of _epsilon
   * \return double
   */
  inline double toleranceED()
  {
    return _TOL_ED;
  }

  /** To know if Newton Iteratio is convergent
   * \return _isNewtonConverge
   */
  bool isNewtonConverge()
  {
    return _isNewtonConverge;
  };

  /** To known the number of steps performed by the Newton algorithm 
   * \return _newtonNbIterations
   */
  unsigned int getNewtonNbIterations()
  {
    return _newtonNbIterations;
  }

  /** Set value to the maximum number of iterations
   * \param maxStep maximum number of step
   */
  void setNewtonMaxIteration(unsigned int maxStep)
  {
    _newtonMaxIteration = maxStep;
  };

  /** To set maximum number of iterations to localize events
   *\param maxIter maximum number of iterations
   */
  void setLocalizeEventsMaxIteration(unsigned int maxIter)
  {
    _localizeEventMaxIter = maxIter;
  }

  /** get the maximum number of Newton iteration
   *  \return unsigned int
   */
  double newtonMaxIteration()
  {
    return _newtonMaxIteration;
  };

  /** get the maximum number of iterations to localize events
   * \return unsigned int: maximum number of iterations
   */
  unsigned int LocalizeEventsMaxIteration()
  {
    return _localizeEventMaxIter;
  }

  /** accessor to _newtonResiduDSMax
   * \return double _newtonResiduDSMax
   */
  double newtonResiduDSMax()
  {
    return _newtonResiduDSMax;
  };

  /** accessor to _newtonResiduYMax 
   * \return double _newtonResiduYMax
   */
  double newtonResiduYMax()
  {
    return _newtonResiduYMax;
  };

  /** set the Default Newton tolerance
   *  \param tol new tolerance
   */
  void setNewtonTolerance(double tol)
  {
    _newtonTolerance = tol;
  };

  /** get the Newton tolerance
   *  \return tolerance
   */
  double newtonTolerance()
  {
    return _newtonTolerance;
  };

  /** Redefine method insertIntegrator of the class Simulation */
  void insertIntegrator(SP::OneStepIntegrator);

  /** update indexSets[i] of the topology, using current y and lambda values of Interactions.
   *  \param i the number of the set to be updated
   */

  void updateIndexSet(unsigned int i);

  /** update indexSets[1] and [2] (using current y and lambda values of Interactions) with conditions on y[2] AND lambda[2].
  */
  void updateIndexSetsWithDoubleCondition();

  /** compute right-hand side of xdot = f(x,t), for the integrator osi.
   *  \param osi the integrator (Lsodar)
   *  \param sizeOfX size of vector x
   *  \param time current time given by the integrator
   *  \param x state vector
   *  \param xdot derivative of x
   */
  void computef(OneStepIntegrator& osi, integer* sizeOfX, doublereal* time, doublereal* x, doublereal* xdot);

  /** compute jacobian of the right-hand side
   *  \param osi the integrator (Lsodar)
   *  \param sizeOfX size of vector x
   *  \param time current time given by the integrator
   *  \param x state vector
   *  \param jacob jacobian of f according to x
   */
  void computeJacobianfx(OneStepIntegrator& osi, integer* sizeOfX, doublereal* time, doublereal* x,  doublereal* jacob);

  /** compute the size of constraint function g(x,t,...) for osi 
   * \return unsigned int 
   */
  virtual unsigned int computeSizeOfg();

  /** compute constraint function g(x,t,...) for osi.
   *  \param osi pointer to OneStepIntegrator.
   *  \param sizeX integer*, size of vector x
   *  \param time doublereal*, time
   *  \param x doublereal*, x:array of double
   *  \param sizeG integer*, size of vector g (ie number of constraints)
   *  \param g doublereal*, g (in-out parameter)
   */
  virtual void computeg(SP::OneStepIntegrator osi, integer* sizeX,
                        doublereal* time, doublereal* x,
                        integer* sizeG, doublereal* g);

  /** update input for impact case (ie compute p[1])
  */
  void updateImpactState();

  /** update state for smooth dynamic case (i.e. compute p[2] and update acceleration) */

  void updateSmoothState();

  /** update input
   *  \param level of lambda  used to compute input
   */
  void updateInput(unsigned int level) {}

  /** update state.
   *  \param level of lambda  used to compute input
   */
  void updateState(unsigned int level);

  /** update output and indexSets.
   *  \param level of lambda  used to compute input
   */
  void updateOutput(unsigned int level);

  /** Initialize EventDriven **/

  /** run simulation from one Event to the next, according to events manager settings.
   *
   */
  void advanceToEvent();

  /** visitors hook
  */

  /** Methods for NewMarkAlphaOSI scheme */


  /** compute maximum residu over all gap functions of Index2 contacts
   * \return double: maximum residu for all gap functions
   */
  double computeResiduConstraints();

  /** prepare for Newton iterations for all OSIs
   *\return maximum residu over all DSs
  */
  void prepareNewtonIteration();

  /** Check convergence of Newton iteration
   *\return bool: true if convergent, false otherwise
   */
  bool newtonCheckConvergence(double);

  /** Predict the state of all Dynamical Systems before Newton iterations */
  void predictionNewtonIteration();

  /** Correct the state of all Dynamical Systems during Newton iterations */
  void correctionNewtonIteration();

  /** Newton iteration to get the state of all Dynamical Systems at the end of step
   *\param criterion tolerance to check convergence
   *\param maxStep maximum number of steps
   */
  void newtonSolve(double criterion, unsigned int maxStep);

  /** Detect whether or not events occur during each integration step
   *\param updateIstate true if we need to update the flag _istate, false otherwise
   *\return double, maximum of absolute values of constraint fonctions over all activated ot deactivated contacts
   */
  double detectEvents(bool updateIstate = true);

  /** Localize time of the first event */
  void LocalizeFirstEvent();


  ACCEPT_STD_VISITORS();

};
#endif // EventDriven_H
