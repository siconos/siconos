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

/** \file EulerMoreauOSI.hpp */

#ifndef EulerMoreauOSI_H
#define EulerMoreauOSI_H

#include "OneStepIntegrator.hpp"


const unsigned int EULERMOREAUSTEPSINMEMORY = 1;

/** \class EulerMoreauOSI
 * Time-Integrator for Dynamical Systems
 *  \brief One Step time Integrator for First Order Dynamical Systems.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.7.0.
 *  \date (Creation) Dec 26, 2013
 *
 * This integrator is the work horse of the event--capturing time stepping schemes
 * for first order systems.
 * It is mainly based on some extensions of the Backward Euler and \f$\theta-\gamma\f$
 * schemes proposed in the pionnering work of J.J. Moreau for the sweeping process
 *
 * J.J. Moreau. Evolution problem associated with a moving convex set in a Hilbert space.
 * Journal of Differential Equations, 26, pp 347--374, 1977.
 *
 * Variants are now used to integrate LCS, Relay systems, Higher order sweeping process see
 * for instance
 *
 *
 * Consistency of a time-stepping method for a class of piecewise linear networks
 *
 * M.K. Camlibel, W.P.M.H. Heemels, and J.M. Schumacher
 * IEEE Transactions on Circuits and Systems I, 2002, 49(3):349--357
 *
 * Numerical methods for nonsmooth dynamical systems: applications in mechanics and electronics
 *
 * V Acary, B Brogliato Springer Verlag 2008
 *
 * Convergence of time-stepping schemes for passive and extended linear complementarity systems
 * L. Han, A. Tiwari, M.K. Camlibel, and J.-S. Pang SIAM Journal on Numerical Analysis 2009, 47(5):3768-3796
 *
 * On preserving dissipativity properties of linear complementarity dynamical systems with the &theta-method
 *
 * Greenhalgh Scott, Acary Vincent, Brogliato Bernard Numer. Math., , 2013.
 *
 * Main time--integration schemes are based on the following \f$\theta-\gamma\f$ scheme
 *
 * \f{equation}{
 *  \begin{cases}
 *   \label{eq:toto1}
 *     M x_{k+1} = M x_{k} +h\theta f(x_{k+1},t_{k+1})+h(1-\theta) f(x_k,t_k) + h \gamma r(t_{k+1})
 *   + h(1-\gamma)r(t_k)  \\[2mm]
 *   y_{k+1} =  h(t_{k+1},x_{k+1},\lambda _{k+1}) \\[2mm]
 *   r_{k+1} = g(x_{k+1},\lambda_{k+1},t_{k+1})\\[2mm]
 *    \mbox{nslaw} ( y_{k+1} , \lambda_{k+1})
 * \end{cases}
 * \f}
 * where \f$\theta = [0,1]\f$ and \f$\gamma \in [0,1]\f$.
 * As in Acary & Brogliato 2008, we call the previous problem  the ``one--step nonsmooth problem''.
 *
 * Another variant can also be used (FullThetaGamma scheme)
 *  \f{equation}{
 *   \begin{cases}
 *     M x_{k+1} = M x_{k} +h f(x_{k+\theta},t_{k+1}) + h r(t_{k+\gamma}) \\[2mm]
 *     y_{k+\gamma} =  h(t_{k+\gamma},x_{k+\gamma},\lambda _{k+\gamma}) \\[2mm]
 *     r_{k+\gamma} = g(x_{k+\gamma},\lambda_{k+\gamma},t_{k+\gamma})\\[2mm]
 *     \mbox{nslaw} ( y_{k+\gamma} , \lambda_{k+\gamma})
 *   \end{cases}
 * \f}
 *
 *
 * EulerMoreauOSI class is used to define some time-integrators methods for a
 * list of first order dynamical systems. A EulerMoreauOSI instance is defined by
 * the value of theta  and possibly gamma and the list of
 * concerned dynamical systems.
 *
 * Each DynamicalSystem is associated to a SiconosMatrix, named "W", which is the "iteration" matrix.
 * W matrices are initialized and computed in initializeIterationMatrixW and computeW. Depending on the DS type, they may
 * depend on time t and DS state x.
 *
 * For first order systems, the implementation uses _r for storing the
 * the input due to the nonsmooth law. This EulerMoreauOSI scheme assumes that the
 * relative degree is zero or one and one level for _r is sufficient
 *
 * Main functions:
 *
 * - computeFreeState(): computes xfree (or vfree), dynamical systems
 *   state without taking non-smooth part into account \n
 *
 * - updateState(): computes x (q,v), the complete dynamical systems
 *    states.
 *
 * See User's guide, for details.
 *
 */

class EulerMoreauOSI : public OneStepIntegrator
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(EulerMoreauOSI);


  /** Stl map that associates a theta parameter for the integration
  *  scheme to each DynamicalSystem of the OSI */
  double _theta;

  /** A gamma parameter for the integration scheme to each DynamicalSystem of the OSI
   * This parameter is used to apply a theta-method to the input $r$
   */
  double _gamma;

  /** a boolean to know if the parameter must be used or not
   */
  bool _useGamma;

  /** a boolean to know if the parameter must be used or not
   */
  bool _useGammaForRelation;

  /** nslaw effects
   */
  struct _NSLEffectOnFreeOutput;
  friend struct _NSLEffectOnFreeOutput;


  /** Default constructor
   */
  EulerMoreauOSI() {};

public:
  /** constructor from theta value only
   *  \param theta value for all DS.
   */
  EulerMoreauOSI(double theta);

  /** constructor from theta value only
   *  \param theta value for all linked DS.
   *  \param gamma value for all linked DS.
   */
  EulerMoreauOSI(double theta, double gamma);

  /** destructor
   */
  virtual ~EulerMoreauOSI() {};

  // --- GETTERS/SETTERS ---

  /** get the value of W corresponding to DynamicalSystem ds
   * \param ds a pointer to DynamicalSystem, optional, default =
   * NULL. get W[0] in that case
   *  \return SimpleMatrix
   */
  const SimpleMatrix getW(SP::DynamicalSystem ds = SP::DynamicalSystem());

  /** get W corresponding to DynamicalSystem ds
   * \param ds a pointer to DynamicalSystem
   * \return pointer to a SiconosMatrix
   */
  SP::SimpleMatrix W(SP::DynamicalSystem ds);

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
   *  \param b true if gamma has to be used, false otherwise
   */
  inline void setUseGamma(bool b)
  {
    _useGamma = b;
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

  /** initialization of the EulerMoreauOSI integrator; for linear time
      invariant systems, we compute time invariant operator (example :
      W)
   */
  //virtual void initialize(Model& m);
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

  /** initialize iteration matrix W EulerMoreauOSI matrix at time t
   *  \param time the time (double)
   *  \param ds a pointer to DynamicalSystem
   */
  void initializeIterationMatrixW(double time, SP::DynamicalSystem ds);

  /** compute W EulerMoreauOSI matrix at time t
   *  \param time the current time
   *  \param ds the DynamicalSystem
   *  \param dsv a descriptor of the ds on the graph (redundant to avoid invocation)
   *  \param W the matrix to compute
   */
  void computeW(double time, DynamicalSystem& ds, DynamicalSystemsGraph::VDescriptor& dsv, SiconosMatrix& W);

  /** compute WBoundaryConditionsMap[ds] EulerMoreauOSI matrix at time t
   *  \param ds a pointer to DynamicalSystem
   */
  void computeWBoundaryConditions(SP::DynamicalSystem ds);

  /** initialize iteration matrix WBoundaryConditionsMap[ds] EulerMoreauOSI
   *  \param ds a pointer to DynamicalSystem
   */
  void initializeIterationMatrixWBoundaryConditions(SP::DynamicalSystem ds);

  /** Computes the residuFree and residu of all the DynamicalSystems
   *  \return the maximum of the 2-norm over all the residu
   */
  double computeResidu();

  /** Perform the integration of the dynamical systems linked to this integrator
   *  without taking into account the nonsmooth input r
   */
  virtual void computeFreeState();

  
  double computeResiduOutput(double time, SP::InteractionsGraph indexSet);
  
  double computeResiduInput(double time, SP::InteractionsGraph indexSet);

  /** integrates the Interaction linked to this integrator, without taking non-smooth effects into account
   * \param vertex_inter of the interaction graph
   * \param osnsp pointer to OneStepNSProblem
   */
  virtual void computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp);

  /** computes all the W matrices
   * \param time current time
   */
  void prepareNewtonIteration(double time);


  /** integrate the system, between tinit and tend (->iout=true), with possible stop at tout (->iout=false)
   *  \param tinit initial time
   *  \param tend end time
   *  \param tout real end time
   *  \param useless flag (for EulerMoreauOSI, used in LsodarOSI)
   */
  void integrate(double& tinit, double& tend, double& tout, int& useless);

  /** updates the state of the Dynamical Systems
   *  \param level the level of interest for the dynamics: not used at the time
   */
  virtual void updateState(const unsigned int level);

  /** Displays the data of the EulerMoreauOSI's integrator
   */
  void display();

  /** visitors hook
  */
  ACCEPT_STD_VISITORS();

};

#endif // EulerMoreauOSI_H
