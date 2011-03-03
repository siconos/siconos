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
#define SICONOS_TS_LINEAR 1
#define SICONOS_TS_LINEAR_IMPLICIT 2
#define SICONOS_TS_NONLINEAR 3
class TimeStepping : public Simulation
{
private:

  /** Default Constructor
   */
  TimeStepping() {};

  /** compute LevelMax */
  void initLevelMax();

  /** boolean variable to known whether the ResiduY has to be computed or not
   *  if true, the ResiduY is computed and the convergence is checked
   */
  bool _computeResiduY;

  /** boolean variable to known whether the ResiduR has to be computed or not
   *  if true, the ResiduR is computed and the convergence is checked
   */
  bool _computeResiduR;

  /** Default Newton tolerance used in call of run() of ComputeOneStep() */
  double _newtonTolerance;

  /** Default maximum number of Newton iteration*/
  unsigned int _newtonMaxIteration;

  /** Number of steps perfomed is the Newton Loop */
  unsigned int _newtonNbSteps;

  /** Maximum Residual for the Dynamical system */
  double _newtonResiduDSMax;

  /** Maximum Residual for the output of the relation */
  double _newtonResiduYMax;

  /** Maximum Residual for the input of the relation */
  double _newtonResiduRMax;


  /** unsigned int  _newtonOptions
   *  option in the Newon iteration
   *  SICONOS_TS_LINEAR or SICONOS_TS_LINEAR_IMPLICIT SICONOS_TS_NONLINEAR will force a single iteration of the Newton Solver
   * SICONOS_TS_NONLINEAR (default) will perform the newton iteration up to convergence
   */
  unsigned int _newtonOptions;
protected:
  /** initialisation specific to TimeStepping for OneStepNSProblem.
   */
  virtual void initOSNS();
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
  /** Used by the updateIndexSet function in order to deactivate an UR.
      \param SP::UnitaryRelation ur: an UR activated.
      \param unsigned int: the number of the set to be updated
      Return true iff the ur must be deactivate.
   */
  virtual bool predictorDeactivate(SP::UnitaryRelation ur, unsigned int i);
  /** Used by the updateIndexSet function in order to activate an UR.
      \param SP::UnitaryRelation ur: an UR deactivated.
      \param unsigned int: the number of the set to be updated
      Return true iff the ur must be activate.
   */
  virtual bool predictorActivate(SP::UnitaryRelation ur, unsigned int i);
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

  /** run one time--step of the simulation
   */
  void computeOneStep();

  /** newton algorithm
   * \param double, convergence criterion
   * \param unsigned int: maximum number of Newton steps
   */
  virtual void newtonSolve(double, unsigned int);

  /** To known the number of steps performed by the Newton algorithm.
   *
   */
  unsigned int getNewtonNbSteps()
  {
    return _newtonNbSteps;
  }

  /** compute initial residu
   * It computes the initial residu to start the newton algorithm.
   *
   */
  void computeInitialResidu();


  void prepareNewtonIteration();

  /** check the convergence of Newton algorithm according to criterion
   * \param double, convergence criterion
   */
  bool newtonCheckConvergence(double);

  /*save y_k^p, the current Newton iteration*/
  void saveYandLambdaInMemory();

  /** run the simulation, from t0 to T
   * with default parameters if any setting has been done
   */
  void run();

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

  /** To specify if the output interaction residu must be computed.
   *  \param: a bool
   */
  void setComputeResiduY(bool v)
  {
    _computeResiduY = v;
  };

  /** To known if the output interaction residu must be computed. */
  bool computeResiduY(bool v)
  {
    return _computeResiduY;
  };


  /** To specify if the input interaction residu must be computed.
   *  \param: a bool
   */
  void setComputeResiduR(bool v)
  {
    _computeResiduR = v;
  };

  /** To known if the input interaction residu must be computed. */
  bool computeResiduR(bool v)
  {
    return _computeResiduR;
  };


  /** set the Default Newton tolerance
   *  \param: double
   */
  void setNewtonTolerance(double tol)
  {
    _newtonTolerance = tol;
  };

  /** get the Newton tolerance
   *  \return double
   */
  double newtonTolerance()
  {
    return   _newtonTolerance;
  };

  /** set the maximum number of Newton iteration
   *  \param: unsigned int
   */
  void setNewtonMaxIteration(unsigned int maxStep)
  {
    _newtonMaxIteration = maxStep;
  };

  /** get the maximum number of Newton iteration
   *  \return unsigned int
   */
  double newtonMaxIteration()
  {
    return _newtonMaxIteration;
  };

  /** set the NewtonOptions
   *  \param: std::string
   */
  void setNewtonOptions(unsigned int v)
  {
    _newtonOptions = v;
  };

  /** get the NewtonOptions
   *  \return unsigned int
   */
  unsigned int newtonOptions()
  {
    return _newtonOptions;
  };


  /** accessor to _newtonResiduDSMax
   */
  double newtonResiduDSMax()
  {
    return _newtonResiduDSMax;
  };

  /** accessor to _newtonResiduYMax
   */
  double newtonResiduYMax()
  {
    return _newtonResiduYMax;
  };

  /** accessor to _newtonResiduRMax
   */
  double newtonResiduRMax()
  {
    return _newtonResiduRMax;
  };


  /*TS set the ds->q memory, the world (CAD model for example) must be updated.
    Overload this method to update user model.*/
  virtual void updateWorldFromDS()
  {
    ;
  };



  /** visitors hook
   */
  ACCEPT_STD_VISITORS();

};

DEFINE_SPTR(TimeStepping);

#endif // TIMESTEPPING_H





