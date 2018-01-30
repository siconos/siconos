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
/** \file Hem5OSI.hpp
 * Hem5OSI solver (from E. Hairer software lists)
 */
#ifndef Hem5OSI_H
#define Hem5OSI_H

#include "OneStepIntegrator.hpp"
#include <vector>

#define HEM5_ATOL_DEFAULT 100 * MACHINE_PREC;
#define HEM5_RTOL_DEFAULT 10 * MACHINE_PREC;

class Hem5OSI_impl;
TYPEDEF_SPTR(Hem5OSI_impl);

/** Hem5OSI solver (odepack)
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
class Hem5OSI : public OneStepIntegrator
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Hem5OSI);

  /** vector of integer data for the integrator
   * _intData[0] NQ size of the position vector q
   * _intData[1] NV size of the velocity  vector v NQ >= NQ
   * _intData[2] NU size of the external dynamic vector u
   * _intData[3] NL size of the Lagrange multiplier vector lambda
   * _intData[4] ITOL indicates whether RTOL and ATOL are scalar (ITOL=0), or array of
   *             dimension NQ + NV + NU (ITOL=1)
   * _intData[5] IOUT selects the dense output formula
   * _intData[6] LWK length of real array rwork
   * _intData[7] LIWK length of integer array iwork
   * See hem5.f
   */
  std::vector<integer> _intData;


  integer _idid;

  /** relative tolerance */
  SA::doublereal rtol;
  /** absolute tolerance */
  SA::doublereal atol;
  /** real work array */
  SA::doublereal rwork;
  /** integer work array */
  SA::integer iwork;


  doublereal _timeStep; // initial step size guess

  /** temporary vector to save q values */
  SP::BlockVector _qWork;
  /** temporary vector to save v values */
  SP::BlockVector _vWork;
  /** temporary vector to save v values */
  SP::BlockVector _uWork;
  /** temporary vector to save a values */
  SP::BlockVector _aWork;
  /** temporary vector to save lambda values */
  SP::BlockVector _lambdaWork;
  /** temporary vector to save forces values */
  SP::BlockVector _forcesWork;

  SP::SiconosVector _qtmp;
  SP::SiconosVector _vtmp;
  SP::SiconosVector _utmp;
  SP::SiconosVector _atmp;
  SP::SiconosVector _lambdatmp;
  SP::SiconosVector _forcestmp;

  /** nslaw effects
   */
  struct _NSLEffectOnFreeOutput;
  friend struct _NSLEffectOnFreeOutput;


public:
  SP::Hem5OSI_impl _impl;
  friend class Hem5OSI_impl;

  enum {OSNSP_RHS,WORK_INTERACTION_LENGTH};

  /** constructor from a minimum set of data
   */
  Hem5OSI();

  /** destructor
   */
  ~Hem5OSI() {};

  /** get vector of integer parameters for lsodar
   *  \return a vector<integer>
   */
  inline const std::vector<integer> intData() const
  {
    return _intData;
  }

  /** get _intData[i]
   * \param i index
   * \return an integer
   */
  inline integer intData(unsigned int i) const
  {
    return _intData[i];
  }
  /** set _intData[i]
   * \param i index
   * \param newValue
   */
  inline void setIntData(unsigned int i, int newValue)
  {
    _intData[i] = newValue;
  }

  /** get relative tolerance parameter for Hem5
   *  \return a doublereal*
   */
  inline const SA::doublereal getRtol() const
  {
    return rtol;
  }

  /** get absolute tolerance parameter for Hem5
   *  \return a doublereal*
   */
  inline const SA::doublereal getAtol() const
  {
    return atol;
  }

  /** get the maximum number of steps for one call
  *\return an interger
  */
  inline  int getMaxNstep()const
  {
    return iwork[11];
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

  /** set itol, rtol and atol (tolerance parameters for Hem5)
   *  \param itol integer (itol value)
   *  \param rtol doublereal * (rtol)
   *  \param atol doublereal * (atol)
   */
  void setTol(integer itol, SA::doublereal rtol, SA::doublereal atol);

  /** set itol, rtol and atol (scalar tolerance parameters for Hem5)
   *  \param itol integer (itol value)
   *  \param rtol double (rtol)
   *  \param atol double (atol)
   */
  void setTol(integer itol, doublereal rtol, doublereal atol);

  /** set the maximul number of steps for one call of Hem5OSI
   *\param nstepmax an integer
   */
  void setMaxNstep(integer nstepmax);

  /** set the minimum and maximum step sizes
   * \param maxstepsize double (maximul step size)
   */
  void setMaxStepSize(doublereal maxstepsize);

  /** update _intData
   */
  void updateIntData();

  /** update doubleData and iwork memory size, when changes occur in _intData.
   */
  void updateData();

  /** fill qWork with a doublereal
   *  \param sizex integer*, size of x array
   *  \param x doublereal* x:array of double
   */
  void fillqWork(integer* sizex, doublereal* x) ;

  /** fill vWork with a doublereal
   *  \param sizex integer*, size of x array
   *  \param x doublereal* x:array of double
   */
  void fillvWork(integer* sizex, doublereal* x) ;

  /** compute rhs(t) for all dynamical systems in the set
   */
  void computeRhs(double) ;

  /** compute jacobian of the rhs at time t for all dynamical systems in the set
   */
  void computeJacobianRhs(double) ;

  unsigned int numberOfConstraints();

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
   *  \param idid in-out parameter, input: 1 for first call, else 2. Output: 2 if no root was found, else 3.
   */
  void integrate(double& tinit, double& tend, double& tout, int& idid);

  /** update the state of the DynamicalSystems attached to this Integrator
   *  \param level level of interest for the dynamics
   */
  void updateState(const unsigned int level);

  void prepareNewtonIteration(double time)
  {
    assert(0);
  };

  /** integrates the Interaction linked to this integrator, without taking non-smooth effects into account
   * \param vertex_inter of the interaction graph
   * \param osnsp pointer to OneStepNSProblem
   */
  virtual void computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp);

  /** print the data to the screen
   */
  void display();

  /** visitors hook
  */
  ACCEPT_STD_VISITORS();
};

#endif // Hem5OSI_H
