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

/*! \file  MoreauJeanOSI.hpp */

#ifndef MoreauJeanOSI_H
#define MoreauJeanOSI_H

#include "OneStepIntegrator.hpp"

#include <limits>

const unsigned int MOREAUSTEPSINMEMORY = 1;

/**  \class MoreauJeanOSI 
 *   \brief One Step time Integrator for First Order Dynamical Systems  for
 *    mechanical Systems (LagrangianDS and NewtonEulerDS)
 *   \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date (Creation) Apr 26, 2004
 *
 * This integrator is the work horse of the event--capturing time stepping schemes
 * for mechanical systems.  It is mainly based on the pioneering works of M. Jean and
 * J.J. Moreau for the time integration of mechanical systems
 * with unilateral contact, impact and Coulomb's friction with \f$\theta\f$ scheme
 *
 * For the linear Lagrangina system, the scheme reads as
 * \f{equation}{
 * \begin{cases}
 *  M (v_{k+1}-v_k)
 *  + h K q_{k+\theta} + h C v_{k+\theta}     -   h F_{k+\theta} = p_{k+1} = G P_{k+1},\label{eq:MoreauTS-motion}\\[1mm] 
 *  q_{k+1} = q_{k} + h v_{k+\theta}, \quad \\[1mm]
 *  U_{k+1} = G^\top\, v_{k+1}, \\[1mm]
 *  \begin{array}{lcl}
 *    0 \leq U^\alpha_{k+1} + e  U^\alpha_{k} \perp P^\alpha_{k+1}  \geq 0,& \quad&\alpha \in \mathcal I_1, \\[1mm]
 *    P^\alpha_{k+1}  =0,&\quad& \alpha \in \mathcal I \setminus \mathcal I_1,
 * \end{array}
 * \end{cases} \f}
 * with  \f$\theta \in [0,1]\f$. The index set \f$\mathcal I_1\f$ is the discrete equivalent
 * to the rule that allows us to apply the Signorini  condition at the velocity level. 
 * In the numerical practice, we choose to define this set by
 * \f{equation}{
 *   \label{eq:index-set1}
 *  \mathcal I_1 = \{\alpha \in \mathcal I \mid G^\top (q_{k} + h v_{k}) + w \leq 0\text{ and } U_k \leq 0 \}.
 * \f}.
 *
 * For more details, we refer to
 *
 * M. Jean and J.J. Moreau. Dynamics in the presence of unilateral contacts and dry friction: a numerical approach.
 * In G. Del Pietro and F. Maceri, editors, Unilateral problems in structural analysis. 
 * II, pages 151–196. CISM 304, Spinger Verlag, 1987.
 *
 * J.J. Moreau. Unilateral contact and dry friction in finite freedom dynamics. 
 * In J.J. Moreau and Panagiotopoulos P.D., editors, Nonsmooth Mechanics and Applications, 
 * number 302 in CISM, Courses and lectures, pages 1–82. CISM 302, Spinger Verlag, Wien- New York, 1988a.
 *
 * J.J. Moreau. Numerical aspects of the sweeping process. 
 * Computer Methods in Applied Mechanics and Engineering, 177:329–349, 1999. 
 *
 * M. Jean. The non smooth contact dynamics method. 
 * Computer Methods in Applied Mechanics and Engineering, 177:235–257, 1999.
 *
 * and for a review :
 *
 * V. Acary and B. Brogliato. Numerical Methods for Nonsmooth Dynamical Systems: 
 * Applications in Mechanics and Electronics, volume 35 of Lecture Notes in
 * Applied and Computational Mechanics. Springer Verlag, 2008.
 *
 * 
 * MoreauJeanOSI class is used to define some time-integrators methods for a
 * list of dynamical systems. A MoreauJeanOSI instance is defined by the value
 * of theta and the list of concerned dynamical systems.
 *
 * Each DynamicalSystem is associated to a SiconosMatrix, named "W", the "iteration" matrix"
 * W matrices are initialized and computed in initW and computeW. Depending on the DS type, 
 * they may depend on time t and DS state x.
 *
 * For mechanical systems, the implementation uses _p for storing the
 * the input due to the nonsmooth law. This MoreauJeanOSI scheme assumes that the
 * relative degree is two. 
 *
 * For Lagrangian systems, the implementation uses _p[1] for storing the
 * discrete impulse.
 *
 * Main functions:
 *
 * - computeFreeState(): computes xfree (or vfree), dynamical systems
 *   state without taking non-smooth part into account \n
 *
 * - updateState(): computes x (q,v), the complete dynamical systems
 *    states.
 * See User's guide for details.
 *
 */
class MoreauJeanOSI : public OneStepIntegrator
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(MoreauJeanOSI);


  /** Stl map that associates a W MoreauJeanOSI matrix to each DynamicalSystem of the OSI */
  MapOfDSMatrices WMap;

  /** Stl map that associates the columns of  W MoreauJeanOSI matrix to each DynamicalSystem of the OSI if it has some boundary conditions */
  MapOfDSMatrices _WBoundaryConditionsMap;

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

  /** a boolean to force the evaluation of T in an explicit way
   */
  bool _explicitNewtonEulerDSOperators;

  /** nslaw effects
   */
  struct _NSLEffectOnFreeOutput;
  friend struct _NSLEffectOnFreeOutput;

public:

  /** constructor from theta value only
   *  \param theta value for all linked DS (default = 0.5).
   *  \param gamma value for all linked DS (default = NaN and gamma is not used).
   */
  MoreauJeanOSI(double theta = 0.5, double gamma = std::numeric_limits<double>::quiet_NaN());

  /** destructor
   */
  virtual ~MoreauJeanOSI() {};

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

  /** set the value of W[ds] to newValue
   * \param newValue SiconosMatrix
   * \param ds a pointer to DynamicalSystem,
   */
  void setW(const SiconosMatrix& newValue, SP::DynamicalSystem ds);

  /** set W[ds] to pointer newPtr
   * \param newPtr
   * \param ds a pointer to DynamicalSystem
   */
  void setWPtr(SP::SimpleMatrix newPtr, SP::DynamicalSystem ds);

  // -- WBoundaryConditions --

  /** get WBoundaryConditions map
   *  \return a MapOfDSMatrices
   */
  inline MapOfDSMatrices getWBoundaryConditionsMap() const
  {
    return _WBoundaryConditionsMap;
  };

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

  /** set the Boolean to indicate that we use gamma
   *  \param newUseGamma  a  Boolean variable
   */
  inline void setUseGamma(bool newUseGamma)
  {
    _useGamma = newUseGamma;
  };

  /** get bool gammaForRelation for the relation
   *  \return a Boolean
   */
  inline bool useGammaForRelation()
  {
    return _useGammaForRelation;
  };

  /** set the boolean to indicate that we use gamma for the relation
   *  \param newUseGammaForRelation a Boolean
   */
  inline void setUseGammaForRelation(bool newUseGammaForRelation)
  {
    _useGammaForRelation = newUseGammaForRelation;
    if (_useGammaForRelation) _useGamma = false;
  };

  /** get boolean _explicitNewtonEulerDSOperators for the relation
   *  \return a Boolean
   */
  inline bool explicitNewtonEulerDSOperators()
  {
    return _explicitNewtonEulerDSOperators;
  };

  /** set the boolean to indicate that we use gamma for the relation
   *  \param newUseGammaForRelation a Boolean
   */
  inline void setExplicitNewtonEulerDSOperators(bool newExplicitNewtonEulerDSOperators)
  {
    _explicitNewtonEulerDSOperators = newExplicitNewtonEulerDSOperators;
  };

  // --- OTHER FUNCTIONS ---

  /** initialization of the MoreauJeanOSI integrator; for linear time
      invariant systems, we compute time invariant operator (example :
      W)
   */
  virtual void initialize(Model& m);

  /** init WMap[ds] MoreauJeanOSI matrix at time t
   *  \param time
   *  \param ds a pointer to DynamicalSystem
   */
  void initW(double time, SP::DynamicalSystem ds );

  /** compute WMap[ds] MoreauJeanOSI matrix at time t
   *  \param time (double)
   *  \param ds a pointer to DynamicalSystem
   */
  void computeW(double time , SP::DynamicalSystem ds);

  /** compute WBoundaryConditionsMap[ds] MoreauJeanOSI matrix at time t
   *  \param ds a pointer to DynamicalSystem
   */
  void computeWBoundaryConditions(SP::DynamicalSystem ds);

  /** init WBoundaryConditionsMap[ds] MoreauJeanOSI
   *  \param ds a pointer to DynamicalSystem
   */
  void initWBoundaryConditions(SP::DynamicalSystem ds);


  /** compute the initial state of the Newton loop.
   */
  void computeInitialNewtonState();


  /** return the maximum of all norms for the "MoreauJeanOSI-discretized" residus of DS
      \return a double
   */
  double computeResidu();

  /** Perform the integration of the dynamical systems linked to this integrator
   *  without taking into account the nonsmooth input (_r or _p)
   */
  virtual void computeFreeState();

  /** integrates the Interaction linked to this integrator, without taking non-smooth effects into account
   * \param vertex_inter vertex of the interaction graph
   * \param osnsp pointer to OneStepNSProblem
   */
  virtual void computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp);

  /** Apply the rule to one Interaction to known if is it should be included
   * in the IndexSet of level i
   * \param inter the Interaction to test
   * \param i level of the IndexSet
   * \return Boolean
   */
  virtual bool addInteractionInIndexSet(SP::Interaction inter, unsigned int i);

  /** Apply the rule to one Interaction to known if is it should be removed
   * in the IndexSet of level i
   * \param inter the Interaction to test
   * \param i level of the IndexSet
   * \return Boolean
   */
  virtual bool removeInteractionInIndexSet(SP::Interaction inter, unsigned int i);


  /** method to prepare the fist Newton iteration
   *   \param time
   */
  void prepareNewtonIteration(double time);


  /** integrate the system, between tinit and tend (->iout=true), with possible stop at tout (->iout=false)
   *  \param tinit the initial time
   *  \param tend the end time
   *  \param tout the real end time
   *  \param notUsed useless flag (for MoreauJeanOSI, used in LsodarOSI)
   */
  void integrate(double& tinit, double& tend, double& tout, int& notUsed);

  /** update the state of the dynamical systems
      \param ds the dynamical to update
   */
  virtual void updatePosition(SP::DynamicalSystem ds);

  /** update the state of the dynamical systems
   *  \param level the level of interest for the dynamics: not used at the time
   */
  virtual void updateState(const unsigned int level);

  /** Displays the data of the MoreauJeanOSI's integrator
   */
  void display();

  /** visitors hook
  */
  ACCEPT_STD_VISITORS();

};

#endif // MoreauJeanOSI_H
