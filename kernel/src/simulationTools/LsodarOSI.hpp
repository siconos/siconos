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
  LsodarOSI solver (from odepack)
*/
#ifndef LsodarOSI_H
#define LsodarOSI_H

#include "OneStepIntegrator.hpp"

#include <vector>

#define ATOL_DEFAULT 100 * MACHINE_PREC;
#define RTOL_DEFAULT 10 * MACHINE_PREC;

/** LsodarOSI solver (odepack)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 * Many parameters are required as input/output for LSODAR. See the documentation of this function
 * in externals/odepack/opkdmain.f to have a full description of these parameters.  \n
 * Most of them are read-only parameters (ie can not be set by user). \n
 *  Except: \n
 *  - jt: Jacobian type indicator (1 means a user-supplied full Jacobian, 2 means an internally generated full Jacobian). \n
 *    Default = 2.
 *  - itol, rtol and atol \n
 *    ITOL   = an indicator for the type of error control. \n
 *    RTOL   = a relative error tolerance parameter, either a scalar or array of length NEQ. \n
 *    ATOL   = an absolute error tolerance parameter, either a scalar or an array of length NEQ.  Input only.
 */
class LsodarOSI : public OneStepIntegrator
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LsodarOSI);


  /** neq, ng, itol, itask, istate, iopt, lrw, liw, jt
   * See opkdmain.f and lsodar routine for details on those variables.
   */
  std::vector<integer> _intData;
  /** relative tolerance */
  SA::doublereal rtol;
  /** absolute tolerance */
  SA::doublereal atol;
  /** real work array */
  SA::doublereal rwork;
  /** integer work array */
  SA::integer iwork;
  /** integer array used for output of root information */
  SA::integer jroot;
  /** temporary vector to save x values */
  SP::BlockVector _xWork;

  SP::SiconosVector _xtmp;
  /** nslaw effects
   */
  struct _NSLEffectOnFreeOutput;
  friend struct _NSLEffectOnFreeOutput;

public:
  
  enum {OSNSP_RHS,WORK_INTERACTION_LENGTH};


  /** Lsodar counter : Number of steps taken for the problem so far. */
  static int count_NST;
  /** Number of RHS evaluations for the problem so far. */
  static int count_NFE;

  /** Default constructor */
  LsodarOSI();

  /** destructor
   */
  ~LsodarOSI() {};

  /** get vector of integer parameters for lsodar
   *  \return a vector<integer>
   */
  inline const std::vector<integer> intData() const
  {
    return _intData;
  }

  /** get _intData[i]
   * \param i index number (starting from 0)
   *  \return an integer
   */
  inline integer intData(unsigned int i) const
  {
    return _intData[i];
  }

  /** set _intData[i]
   * \param i index number (starting from 0)
   * \param newValue the new value
   */
  inline void setIntData(unsigned int i, int newValue)
  {
    _intData[i] = newValue;
  }

  /** get relative tolerance parameter for lsodar
   *  \return a doublereal*
   */
  inline const SA::doublereal getRtol() const
  {
    return rtol;
  }

  /** get absolute tolerance parameter for lsodar
   *  \return a doublereal*
   */
  inline const SA::doublereal getAtol() const
  {
    return atol;
  }

  /** get the maximum number of steps for one call
  *\return an interger
  */
  inline int getMaxNstep() const
  {
    return iwork[5];
  }

  /** get real work vector parameter for lsodar
   *  \return a doublereal*
   */
  inline const SA::doublereal getRwork() const
  {
    return rwork;
  }

  /** get iwork
   *  \return a pointer to integer
   */
  inline SA::integer getIwork() const
  {
    return iwork;
  }

  /** get output of root information
   *  \return a pointer to integer
   */
  inline SA::integer getJroot() const
  {
    return jroot;
  }

  /** set Jt value, Jacobian type indicator. Excerpts from the lsodar documentation.
   *    1 means a user-supplied full (neq by neq) jacobian.
   *    2 means an internally generated (difference quotient) full
   *      jacobian (using neq extra calls to f per df/dy value).
   *    4 means a user-supplied banded jacobian.
   *    5 means an internally generated banded jacobian (using
   *      ml+mu+1 extra calls to f per df/dy evaluation).
   *  if jt = 1 or 4, the user must supply a subroutine jac
   *  (the name is arbitrary) as described above under jac.
   *  if jt = 2 or 5, a dummy argument can be used.
   *  \param newJT new value for the jt parameter.
   */
  inline void setJT(integer newJT)
  {
    _intData[8] = newJT;
  };

  /** set itol, rtol and atol (tolerance parameters for lsodar)
   *  \param newItol itol value
   *  \param newRtol rtol value
   *  \param newAtol atol value
   */
  void setTol(integer newItol, SA::doublereal newRtol, SA::doublereal newAtol);

  /** set itol, rtol and atol (scalar tolerance parameters for lsodar)
   *  \param newItol itol value
   *  \param newRtol rtol value
   *  \param newAtol atol value
   */
  void setTol(integer newItol, doublereal newRtol, doublereal newAtol);

  /** set the maximum number of steps for one call of Lsodar
   * \param maxNumberSteps the maximum number of steps
   */
  void setMaxNstep(integer maxNumberSteps);

  /** set the minimum and maximum step sizes
   *\param minStep minimum step size
   *\param maxStep maximum step size
   */
  void setMinMaxStepSizes(doublereal minStep, doublereal maxStep);

  /** set maximum method order
   *\param maxorderNonStiff maximum order for nonstiff methods
   *\param maxorderStiff maximum order for stiff methods
   */
  void setMaxOrder(integer maxorderNonStiff, integer maxorderStiff);

  /** update doubleData and iwork memory size, when changes occur in _intData.
   */
  void updateData();

  /** fill xWork with a doublereal
   *  \param size size of x array
   *  \param array x array of double
   */
  void fillXWork(integer* size, doublereal* array);

  /** compute rhs(t) for all dynamical systems in the set
   * \param t current time of simulation
   * \param DSG0 the graph of DynamicalSystem
   */
  void computeRhs(double t, DynamicalSystemsGraph& DSG0);

  /** compute jacobian of the rhs at time t for all dynamical systems in the set
   * \param t current time of simulation
   * \param DSG0 the graph of DynamicalSystem
   */
  void computeJacobianRhs(double t, DynamicalSystemsGraph& DSG0);

  void f(integer* sizeOfX, doublereal* time, doublereal* x, doublereal* xdot);

  void g(integer* nEq, doublereal* time, doublereal* x, integer* ng, doublereal* gOut);

  void jacobianfx(integer*, doublereal*, doublereal*, integer*, integer*,  doublereal*, integer*);

  /** initialization of the integrator
   */
  void initialize(Model& m);

  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi
   * \param m the Model
   * \param t time of initialization
   * \param ds the dynamical system
   */
  void initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds);

  /** initialization of the work vectors and matrices (properties) related to
   *  one interaction on the graph and needed by the osi
   * \param inter the interaction
   * \param interProp the properties on the graph
   * \param DSG the dynamical systems graph
   */
  void fillDSLinks(Interaction &inter,
		     InteractionProperties& interProp,
		     DynamicalSystemsGraph & DSG);

  /** get the number of index sets required for the simulation
   * \return unsigned int
   */
  unsigned int numberOfIndexSets() const {return 3;};
  
  /** integrate the system, between tinit and tend (->iout=true), with possible stop at tout (->iout=false)
   *  \param tinit initial time
   *  \param tend end time
   *  \param tout real end time
   *  \param ioparam in-out parameter, input: 1 for first call, else 2. Output: 2 if no root was found, else 3.
   */
  void integrate(double& tinit, double& tend, double& tout, int& ioparam);

  /** update the state of the DynamicalSystems attached to this Integrator
   *  \param level level of interest for the dynamics
   */
  void updateState(const unsigned int level);


  void prepareNewtonIteration(double time)
  {
    assert(0);
  };

  /** integrates the Interaction linked to this integrator, without taking non-smooth effects into account
   * \param vertex_descr descriptor vertex of the interaction graph
   * \param osnsp pointer to OneStepNSProblem
   */
  virtual void computeFreeOutput(InteractionsGraph::VDescriptor& vertex_descr, OneStepNSProblem* osnsp);

  /** print the data to the screen
   */
  void display();

  /** Return current number of rhs call (for all lsodar-like OSIs!)
   * \return int
   */
  static int count_rhs_call()
  {
    return count_NFE;
  }

  /** Return the number of lsodar steps already done (for all lsodar-like OSIs!)
   * \return int
   */
  static int count_steps()
  {
    return count_NST;
  }

  /** visitors hook
  */
  ACCEPT_STD_VISITORS();
};

#endif // LsodarOSI_H
