/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
  Time-Stepping simulation
*/
#ifndef TIMESTEPPING_H
#define TIMESTEPPING_H

#include "Simulation.hpp"

#include "SiconosPointers.hpp"

/** type of function used to post-treat output info from solver. */
typedef void (*CheckSolverFPtr)(int, Simulation*);

/** Time-Stepping scheme
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 */
class TimeStepping : public Simulation
{
private:

  /** Default Constructor
   */
  TimeStepping() {};

  /** initialisation specific to TimeStepping for OneStepNSProblem.
   */
  void initOSNS();

  /** compute LevelMax */
  void initLevelMax();

  bool mComputeResiduY;

  unsigned int mNbNewtonSteps;


public:

  /** Constructor with the time-discretisation.
  *  \param a pointer to a timeDiscretisation (linked to the model
  *  that owns this simulation)
     \param a one step integrator (default none)
     \param a one step non smooth problem (default none)
  */
  TimeStepping(SP::TimeDiscretisation,
               SP::OneStepIntegrator = SP::OneStepIntegrator(),
               SP::OneStepNSProblem = SP::OneStepNSProblem());

  /** constructor with XML object for TimeStepping
      \param SimulationXML* : the XML object corresponding
      \param initial time
      \param final time
      \param the set of all DS in the NSDS
      \param the set of all interactions in the NSDS
  */
  TimeStepping(SP::SimulationXML, double, double, SP::DynamicalSystemsSet , SP::InteractionsSet);

  /** Destructor.
   */
  ~TimeStepping();

  /* type name because parent class needs it */
  inline std::string typeName()
  {
    return Type::name(*this);
  };

  /** add a OneStepNSProblem of the Simulation (if its not the first, it needs to have an id clearly defined)
   *  \param a pointer to OneStepNSProblem

  void insertNonSmoothProblem(SP::OneStepNSProblem);
  */
  /** update indexSets[i] of the topology, using current y and lambda values of Interactions.
  *  \param unsigned int: the number of the set to be updated
  */
  void updateIndexSet(unsigned int);

  /** increment model current time according to User TimeDiscretisation and call SaveInMemory.
   */
  void nextStep();

  /** update input, state of each dynamical system and output
   *  \param lambda order used to compute input
   */
  void update(unsigned int);

  /** integrates all the DynamicalSystems taking not into account nslaw, reactions (ie non-smooth part) ...
   */
  void computeFreeState();

  /** step from current event to next event of EventsManager
   */
  void advanceToEvent();

  /** run one step of the simulation
   */
  void computeOneStep();

  /** newton algorithm
   * \param double, convergence criterion
   * \param unsigned int: maximum number of Newton steps
   */
  void newtonSolve(double, unsigned int);

  /*
   * To known the number of steps performed by the Newton algorithm.
   *
   */
  unsigned int getNewtonNbSteps()
  {
    return mNbNewtonSteps;
  }

  /** compute initial residu
   * It computes the initial residu to start the newton algorithm.
   *
   */
  void computeInitialResidu();


  void   prepareNewtonIteration();

  /** check the convergence of Newton algorithm according to criterion
   * \param double, convergence criterion
   */
  bool newtonCheckConvergence(double);
  void saveYandLambdaInMemory();
  /** run the simulation, from t0 to T
   * \param: simulation option. Default = "linear", else "Newton", used a Newton algorithm.
   * \param: convergence criterion for Newton. Default = 0.005.
   * \param: maximum iteration number for Newton. Default = 500.
   */
  void run(const std::string& = "linear", double = 0.005, unsigned int = 500);

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param SP::Simulation : the Simulation which must be converted
   * \return a pointer on the Simulation if it is of the right type, NULL otherwise
   */
  static TimeStepping* convert(Simulation* str);

  /** check returning value from computeOneStepNSProblem and process
   *  \param: an int
   */
  void DefaultCheckSolverOutput(int);

  /** Set CheckSolverOutput function */
  void setCheckSolverFunction(CheckSolverFPtr);

  /**
   *To specify if the interaction residu must be computed.
   * \param: a bool
   */
  void setComputeResiduY(bool v)
  {
    mComputeResiduY = v;
  };
  /**
   *To known if the interaction residu must be computed.
   */
  bool getComputeResiduY(bool v)
  {
    return mComputeResiduY;
  };

  /** visitors hook
   */
  ACCEPT_STD_VISITORS();

};

DEFINE_SPTR(TimeStepping);

#endif // TIMESTEPPING_H





