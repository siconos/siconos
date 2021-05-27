/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

/*! \file  MoreauJeanGOSI.hpp */

#ifndef MoreauJeanGOSI_H
#define MoreauJeanGOSI_H

#include "MoreauJeanOSI.hpp"
#include "OneStepNSProblem.hpp"

#include <limits>

/**  \class MoreauJeanGOSI
 *   \brief One Step time Integrator for First Order Dynamical Systems  for
 *    mechanical Systems (LagrangianDS and NewtonEulerDS)
 *
 * This integrator is the work horse of the event--capturing time stepping schemes
 * for mechanical systems.  It is mainly based on the pioneering works of M. Jean and
 * J.J. Moreau for the time integration of mechanical systems
 * with unilateral contact, impact and Coulomb's friction with \f$\theta\f$ scheme
 *
 * For the linear Lagrangina system, the scheme reads as
 * \rststar
 *
 * .. math::
 *   :nowrap:
 * 
 *   \begin{cases}
 *    M (v_{k+1}-v_k)
 *     + h K q_{k+\theta} + h C v_{k+\theta}     -   h F_{k+\theta} = p_{k+1} = G P_{k+1},\label{eq:MoreauTS-motion}\\[1mm]
 *     q_{k+1} = q_{k} + h v_{k+\theta}, \quad \\[1mm]
 *     U_{k+1} = G^\top\, v_{k+1}, \\[1mm]
 *     \begin{array}{lcl}
 *      0 \leq U^\alpha_{k+1} + e  U^\alpha_{k} \perp P^\alpha_{k+1}  \geq 0,& \quad&\alpha \in \mathcal I_1, \\[1mm]
 *     P^\alpha_{k+1}  =0,&\quad& \alpha \in \mathcal I \setminus \mathcal I_1,
 *    \end{array}
 *     \end{cases}
 *
 * \endrststar
 *
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
 * MoreauJeanGOSI class is used to define some time-integrators methods for a
 * list of dynamical systems. A MoreauJeanGOSI instance is defined by the value
 * of theta and the list of concerned dynamical systems.
 *
 * Each DynamicalSystem is associated to a SiconosMatrix, named "W", the "iteration" matrix"
 * W matrices are initialized and computed in initializeIterationMatrixW and computeW. Depending on the DS type,
 * they may depend on time t and DS state x.
 *
 * For mechanical systems, the implementation uses _p for storing the
 * the input due to the nonsmooth law. This MoreauJeanGOSI scheme assumes that the
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
class MoreauJeanGOSI : public MoreauJeanOSI
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(MoreauJeanGOSI);

  /** Stl map that associates the columns of  W MoreauJeanGOSI matrix to each DynamicalSystem of the OSI if it has some boundary conditions */
  std::map<unsigned int, SP::SimpleMatrix> _WBoundaryConditionsMap;

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


public:
  
  enum MoreauJeanGOSI_ds_workVector_id {RESIDU_FREE, FREE, LOCAL_BUFFER, WORK_LENGTH};

  enum MoreauJeanGOSI_interaction_workVector_id{OSNSP_RHS, WORK_INTERACTION_LENGTH};

  enum MoreauJeanGOSI_workBlockVector_id{xfree, BLOCK_WORK_LENGTH};

  /** constructor from theta value only
   *  \param theta value for all linked DS (default = 0.5).
   *  \param gamma value for all linked DS (default = NaN and gamma is not used).
   */
  MoreauJeanGOSI(double theta = 0.5, double gamma = std::numeric_limits<double>::quiet_NaN()): MoreauJeanOSI(theta,gamma){};

  /** destructor
   */
  virtual ~MoreauJeanGOSI() {};

  // --- GETTERS/SETTERS ---
  
  // --- OTHER FUNCTIONS ---
  
  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi
   * \param t time of initialization
   * \param ds the dynamical system
   */
  void initializeWorkVectorsForDS( double t, SP::DynamicalSystem ds);

  /** initialization of the work vectors and matrices (properties) related to
   *  one interaction on the graph and needed by the osi
   * \param inter the interaction
   * \param interProp the properties on the graph
   * \param DSG the dynamical systems graph
   */
  void initializeWorkVectorsForInteraction(Interaction &inter,
		     InteractionProperties& interProp,
		     DynamicalSystemsGraph & DSG);

  /** return the maximum of all norms for the "MoreauJeanGOSI-discretized" residus of DS
      \return a double
   */
  double computeResidu();

  /** Perform the integration of the dynamical systems linked to this integrator
   *  without taking into account the nonsmooth input (_r or _p)
   */
  virtual void computeFreeState();

  /** integrate the system, between tinit and tend (->iout=true), with possible stop at tout (->iout=false)
   *  \param tinit the initial time
   *  \param tend the end time
   *  \param tout the real end time
   *  \param notUsed useless flag (for MoreauJeanGOSI, used in LsodarOSI)
   */
  void integrate(double& tinit, double& tend, double& tout, int& notUsed);

  /** update the state of the dynamical systems
      \param ds the dynamical to update
   */
  virtual void updatePosition(DynamicalSystem &ds);

  /** update the state of the dynamical systems
   *  \param level the level of interest for the dynamics: not used at the time
   */
  virtual void updateState(const unsigned int level);

  /** Compute the nonsmooth law contribution to the output
   * \param inter the interaction (for y_k)
   * \param osnsp the non-smooth integrator
   */
  void NonSmoothLawContributionToOutput(SP::Interaction inter, OneStepNSProblem& osnsp);

  /** Displays the data of the MoreauJeanGOSI's integrator
   */
  void display();

  /** visitors hook
  */
  ACCEPT_STD_VISITORS();

};

#endif // MoreauJeanGOSI_H
