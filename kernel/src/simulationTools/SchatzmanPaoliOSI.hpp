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
  SchatzmanPaoliOSI Time-Integrator for Dynamical Systems
*/

#ifndef SCHATZMANPAOLIOSI_H
#define SCHATZMANPAOLIOSI_H

#include "OneStepIntegrator.hpp"
#include "SimpleMatrix.hpp"


const unsigned int SCHATZMANPAOLISTEPSINMEMORY = 2;

/**  SchatzmanPaoliOSI Time-Integrator for Dynamical Systems
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 *
 * SchatzmanPaoliOSI class is used to define some time-integrators methods for a
 * list of dynamical systems.

 * A SchatzmanPaoliOSI instance is defined by the value of theta and the list of
 * concerned dynamical systems.  Each DynamicalSystem is associated to
 * a SiconosMatrix, named "W"
 *
 * W matrices are initialized and computed in initializeIterationMatrixW and
 * computeW. Depending on the DS type, they may depend on time and DS
 * state (x).
 *
 *
 * For Lagrangian systems, the implementation uses _p[0] for storing the
 * discrete multiplier.
 *
 * Main functions:
 *
 * - computeFreeState(): computes xfree (or vfree), dynamical systems
 *   state without taking non-smooth part into account \n
 *
 * - updateState(): computes x (q,v), the complete dynamical systems
 *    states.
 *
 */
class SchatzmanPaoliOSI : public OneStepIntegrator
{
public:
   enum {OSNSP_RHS,WORK_INTERACTION_LENGTH};

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SchatzmanPaoliOSI);

 

  /** Stl map that associates a theta parameter for the integration
  *  scheme to each DynamicalSystem of the OSI */
  double _theta;

  /** A gamma parameter for the integration scheme to each DynamicalSystem of the OSI
   * This parameter is used to apply a theta-method to the input $r$
   */
  double _gamma;

  /** a boolean to known is the parameter must be used or not
   */
  bool _useGamma;

  /** a boolean to known is the parameter must be used or not
   */
  bool _useGammaForRelation;

  /** nslaw effects
   */
  struct _NSLEffectOnFreeOutput;
  friend struct _NSLEffectOnFreeOutput;



  /** Default constructor
   */
  SchatzmanPaoliOSI() {};

public:

  /** constructor from theta value only
   *  \param theta value for all these DS.
   */
  SchatzmanPaoliOSI(double theta);

  /** constructor from theta value only
   *  \param theta value for all these DS.
   *  \param gamma value for all these DS.
   */
  SchatzmanPaoliOSI(double theta, double gamma);

  /** destructor
   */
  virtual ~SchatzmanPaoliOSI() {};

  // --- GETTERS/SETTERS ---

  /** get the value of W corresponding to DynamicalSystem ds
   * \param ds a pointer to DynamicalSystem, optional, default =
   * NULL. get W[0] in that case
   *  \return SimpleMatrix
   */
  const SimpleMatrix getW(SP::DynamicalSystem ds = SP::DynamicalSystem());

  /** get W corresponding to DynamicalSystem ds
   * \param ds a pointer to DynamicalSystem, optional, default =
   * NULL. get W[0] in that case
   * \return pointer to a SiconosMatrix
   */
  SP::SimpleMatrix W(SP::DynamicalSystem ds);

  /** set the value of W[ds] to newValue
   * \param newValue SiconosMatrix
   * \param ds  a pointer to DynamicalSystem,
   */
  void setW(const SiconosMatrix& newValue, SP::DynamicalSystem ds);

  /** set W[ds] to pointer newPtr
   * \param newPtr SP::SiconosMatrix
   * \param ds a pointer to DynamicalSystem
   */
  void setWPtr(SP::SimpleMatrix newPtr, SP::DynamicalSystem ds);

  // -- WBoundaryConditions --

  /** get the value of WBoundaryConditions corresponding to DynamicalSystem ds
   * \param ds a pointer to DynamicalSystem, optional, default =
   * NULL. get WBoundaryConditions[0] in that case
   *  \return SimpleMatrix
   */
  const SimpleMatrix getWBoundaryConditions(SP::DynamicalSystem ds = SP::DynamicalSystem());

  /** get WBoundaryConditions corresponding to DynamicalSystem ds
   * \param ds a pointer to DynamicalSystem, optional, default =
   * NULL. get WBoundaryConditions[0] in that case
   * \return pointer to a SiconosMatrix
   */
  SP::SiconosMatrix WBoundaryConditions(SP::DynamicalSystem ds);

  // -- theta --

  /** get theta
   *  \return a double
   */
  inline double theta()
  {
    return _theta;
  };

  /** set the value of theta
   *  \param newTheta a double
   */
  inline void setTheta(double newTheta)
  {
    _theta = newTheta;
  };

  // -- gamma --

  /** get gamma
   *  \return a double
   */
  inline double gamma()
  {
    return _gamma;
  };

  /** set the value of gamma
   *  \param newGamma a double
   */
  inline void setGamma(double newGamma)
  {
    _gamma = newGamma;
    _useGamma = true;
  };

  // -- useGamma --

  /** get bool useGamma
   *  \return a bool
   */
  inline bool useGamma()
  {
    return _useGamma;
  };

  /** set the boolean to indicate that we use gamma
   *  \param newUseGamma a bool
   */
  inline void setUseGamma(bool newUseGamma)
  {
    _useGamma = newUseGamma;
  };

  /** get bool gammaForRelation for the relation
   *  \return a
   */
  inline bool useGammaForRelation()
  {
    return _useGammaForRelation;
  };


  /** set the boolean to indicate that we use gamma for the relation
   *  \param newUseGammaForRelation a bool
   */
  inline void setUseGammaForRelation(bool newUseGammaForRelation)
  {
    _useGammaForRelation = newUseGammaForRelation;
    if(_useGammaForRelation) _useGamma = false;
  };


  // --- OTHER FUNCTIONS ---

  /** initialization of the SchatzmanPaoliOSI integrator; for linear time
      invariant systems, we compute time invariant operator (example :
      W)
   */
  //void initialize(Model& m);

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
  unsigned int numberOfIndexSets() const {return 1;};
  /** initialize iteration matrix W SchatzmanPaoliOSI matrix at time t
   *  \param time (double)
   *  \param ds a pointer to DynamicalSystem
   */
  void initializeIterationMatrixW(double time, SP::DynamicalSystem ds);

  /** compute W SchatzmanPaoliOSI matrix at time t
   *  \param time the time (double)
   *  \param ds a pointer to DynamicalSystem
   *  \param W the matrix to compute
   */
  void computeW(double time, SP::DynamicalSystem ds, SiconosMatrix& W);

  /** compute WBoundaryConditionsMap[ds] SchatzmanPaoliOSI matrix at time t
   *  \param ds a pointer to DynamicalSystem
   *  \param WBoundaryConditions write the result in WBoundaryConditions
   */
  void computeWBoundaryConditions(SP::DynamicalSystem ds, SiconosMatrix& WBoundaryConditions);

  /** initialize iteration matrix WBoundaryConditionsMap[ds] SchatzmanPaoliOSI
   *  \param ds a pointer to DynamicalSystem
   *  \param dsv a descriptor of the ds on the graph (redundant to avoid invocation)
   */
  void initializeIterationMatrixWBoundaryConditions(SP::DynamicalSystem ds, const DynamicalSystemsGraph::VDescriptor& dsv);

  /** return the maximum of all norms for the "SchatzmanPaoliOSI-discretized" residus of DS
   *  \return a double
   */
  double computeResidu();

  /** integrates the Dynamical System linked to this integrator
   *  without boring the constraints
   */
  virtual void computeFreeState();

  /** integrates the Interaction linked to this integrator, without taking non-smooth effects into account
   * \param vertex_inter of the interaction graph
   * \param osnsp pointer to OneStepNSProblem
   */
  virtual void computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp);

  void prepareNewtonIteration(double time);

  /** integrate the system, between tinit and tend (->iout=true), with possible stop at tout (->iout=false)
   *  \param tinit initial time
   *  \param tend end time
   *  \param tout real end time
   *  \param idid useless flag (for SchatzmanPaoliOSI, used in LsodarOSI)
   */
  void integrate(double& tinit, double& tend, double& tout, int& idid);

  /** updates the state of the Dynamical Systems
   *  \param level level of interest for the dynamics: not used at the time
   */
  virtual void updateState(const unsigned int level);

  /** Displays the data of the SchatzmanPaoliOSI's integrator
   */
  void display();

  /** visitors hook
  */
  ACCEPT_STD_VISITORS();

};

#endif // SCHATZMANPAOLIOSI_H
