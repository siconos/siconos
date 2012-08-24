/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
  Event Driven Simulation
  */
#ifndef EVENTDRIVEN_H
#define EVENTDRIVEN_H

#include "Simulation.hpp"
#include "SiconosNumerics.h"
const double DEFAULT_TOL_ED  = 1000 * DEFAULT_TOLERANCE;

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

  /** an epsilon to define the contraint g for Urs in IndexSet[1]
  */
  static double TOL_ED;

  /** initialisation specific to EventDriven for OneStepNSProblem.
  */
  void initOSNS();

  /** initialisation specific to EventDriven for Interactions.
   */
  void initializeInteraction(SP::Interaction inter);


  // /** compute LevelMin */
  // void initLevelMin();

  // /** compute LevelMax */
  // void initLevelMax();


public:



  /** defaut constructor
   *  \param a SP::timeDiscretisation
   */
  EventDriven(SP::TimeDiscretisation td);

  /** constructor with data
   *  \param a SP::timeDiscretisation
   *  \param an int (number of NSProblem)
   */
  EventDriven(SP::TimeDiscretisation, int);

  /** constructor with XML object of the EventDriven
    \param SimulationXML* : the XML object corresponding
    \param initial time
    \param final time
    \param the set of all DS in the NSDS
    \param the set of all interactions in the NSDS
    */
  EventDriven(SP::SimulationXML, double, double, SP::DynamicalSystemsSet , SP::InteractionsSet);

  /** defaut constructor (needed for serialization)
  */
  EventDriven() {};

  /** destructor
  */
  ~EventDriven() {};

  /* type name because parent class needs it */
  inline std::string typeName()
  {
    return Type::name(*this);
  };

  /* Getters and setters */
  /** Set value to _istate */
  inline void setIstate(int newValue)
  {
    _istate = newValue;
  }
  /** Get value of _istate */
  inline double istate()
  {
    return _istate;
  }
  /** Set value to _epsilon */
  inline void setToleranceED(double var)
  {
    TOL_ED = var;
  }
  /** Get value of _epsilon */
  inline double toleranceED()
  {
    return TOL_ED;
  }
  /** update indexSets[i] of the topology, using current y and lambda values of Interactions.
   *  \param unsigned int: the number of the set to be updated
   */
  void updateIndexSet(unsigned int);

  /** update indexSets[1] and [2] (using current y and lambda values of Interactions) with conditions on y[2] AND lambda[2].
  */
  void updateIndexSetsWithDoubleCondition();

  /** compute right-hand side of xdot = f(x,t), for the integrator osi.
   *  \param pointer to OneStepIntegrator.
   *  \param integer*, size of vector x
   *  \param doublereal*, time
   *  \param doublereal*, x:array of double
   *  \param doublereal*, derivative of x (in-out parameter)
   */
  void computef(SP::OneStepIntegrator, integer*, doublereal*, doublereal*, doublereal*);

  /** compute jacobian of the right-hand side
   *  \param pointer to OneStepIntegrator.
   *  \param integer*, size of vector x
   *  \param doublereal*, time
   *  \param doublereal*, x:array of double
   *  \param doublereal*, jacobian of f according to x (in-out parameter)
   */
  void computeJacobianfx(SP::OneStepIntegrator, integer*, doublereal*, doublereal*,  doublereal*);

  /** compute constraint function g(x,t,...) for osi.
   *  \param pointer to OneStepIntegrator.
   *  \param integer*, size of vector x
   *  \param doublereal*, time
   *  \param doublereal*, x:array of double
   *  \param integer*, size of vector g (ie number of constraints)
   *  \param doublereal*, g (in-out parameter)
   */
  void computeg(SP::OneStepIntegrator, integer*, doublereal*, doublereal*, integer*, doublereal*);

  /** update input for impact case (ie compute p[1])
  */
  void updateImpactState();

  using Simulation::update;

  /** update input, output and indexSets.
   *  \param lambda order used to compute input
   */
  void update(unsigned int);

  /** Initialize EventDriven **/

  /** run simulation from one Event to the next, according to events manager settings.
   *
   */
  void advanceToEvent();

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param SP::Simulation : the Simulation which must be converted
   * \return a pointer on the Simulation if it is of the right type, NULL otherwise
   */
  static EventDriven* convert(Simulation* str);

  /** visitors hook
  */
  ACCEPT_STD_VISITORS();

};
DEFINE_SPTR(EventDriven);
#endif // EVENTDRIVEN_H
