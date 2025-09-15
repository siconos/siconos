/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*! \file  MoreauJeanOSI.hpp */

#ifndef MoreauJeanOSI_H
#define MoreauJeanOSI_H

#include "OneStepIntegrator.hpp"

#include <limits>

const unsigned int MOREAUSTEPSINMEMORY = 1;

/**
   One Step time Integrator, Moreau-Jean algorithm.
   This integrator is the work horse of the event--capturing time stepping
   schemes
   for mechanical systems.  It is mainly based on the pioneering works of M.
   Jean and
   J.J. Moreau for the time integration of mechanical systems
   with unilateral contact, impact and Coulomb's friction with \f$ \theta \f$
   scheme

   For the linear Lagrangian system, the scheme reads as

   \f[

   \begin{cases}
   M (v_{k+1}-v_k)
   + h K q_{k+\theta} + h C v_{k+\theta}     -   h F_{k+\theta} = p_{k+1} = G
   P_{k+1},\label{eq:MoreauTS-motion}\\[1mm] q_{k+1} = q_{k} + h v_{k+\theta},
   \quad \\[1mm] U_{k+1} = G^\top\, v_{k+1}, \\[1mm] \begin{array}{lcl} 0 \leq
   U^\alpha_{k+1} + e  U^\alpha_{k} \perp P^\alpha_{k+1}  \geq 0,& \quad&\alpha
   \in \mathcal I_1, \\[1mm] P^\alpha_{k+1}  =0,&\quad& \alpha \in \mathcal I
   \setminus \mathcal I_1, \end{array} \end{cases}

   \f]

   with  \f$ \theta \in [0,1] \f$. The index set \f$ \mathcal I_1 \f$ is the
   discrete equivalent
   to the rule that allows us to apply the Signorini  condition at the velocity
   level.
   In the numerical practice, we choose to define this set by

   \f[
   \mathcal I_1 = \{\alpha \in \mathcal I \mid G^\top (q_{k} + h v_{k}) + w \leq 0\text{ and } U_k \leq 0 \}.
   \f]

   For more details, we refer to

   M. Jean and J.J. Moreau. Dynamics in the presence of unilateral contacts and
   dry friction: a numerical approach.
   In G. Del Pietro and F. Maceri, editors, Unilateral problems in structural
   analysis.
   II, pages 151–196. CISM 304, Spinger Verlag, 1987.

   J.J. Moreau. Unilateral contact and dry friction in finite freedom dynamics.
   In J.J. Moreau and Panagiotopoulos P.D., editors, Nonsmooth Mechanics and
   Applications,
   number 302 in CISM, Courses and lectures, pages 1–82. CISM 302, Spinger
   Verlag, Wien- New York, 1988a.

   J.J. Moreau. Numerical aspects of the sweeping process.
   Computer Methods in Applied Mechanics and Engineering, 177:329–349, 1999.

   M. Jean. The non smooth contact dynamics method.
   Computer Methods in Applied Mechanics and Engineering, 177:235–257, 1999.

   and for a review :

   V. Acary and B. Brogliato. Numerical Methods for Nonsmooth Dynamical Systems:
   Applications in Mechanics and Electronics, volume 35 of Lecture Notes in
   Applied and Computational Mechanics. Springer Verlag, 2008.

   MoreauJeanOSI class is used to define some time-integrators methods for a
   list of dynamical systems. A MoreauJeanOSI instance is defined by the value
   of theta and the list of concerned dynamical systems.

   Each DynamicalSystem is associated to a SiconosMatrix, named "W", the
   "teration" matrix"
   W matrices are initialized and computed in initializeIterationMatrixW and
   computeW. Depending on the DS type, they may depend on time t and DS state x.

   For mechanical systems, the implementation uses _p for storing the
   the input due to the nonsmooth law. This MoreauJeanOSI scheme assumes that
   the relative degree is two.

   For Lagrangian systems, the implementation uses _p[1] for storing the discrete impulse.

   Main functions:

   - computeFreeState(): computes xfree (or vfree), dynamical systems
   state without taking non-smooth part into account \n

   - updateState(): computes x (q,v), the complete dynamical systems
   states.
   See User's guide for details.

*/
class MoreauJeanOSI : public OneStepIntegrator {
protected:
  ACCEPT_SERIALIZATION(MoreauJeanOSI);

  /** theta-scheme parameter */
  double _theta;

  /** A gamma parameter for the forecast of activation of constraints
   *  leap-frog estimation of the constraints
   *  \f$ \tilde y_k =  y_k + \gamma * h * ydot \f$
   */
  double _gamma;

  /** a boolean to know if the gamma-parameter must be used or not
   */
  bool _useGamma;

  /** Constraint activation threshold
   *
   */
  double _constraintActivationThreshold;

  /** a boolean to know if the parameter must be used or not
   */
  bool _useGammaForRelation;

  /** a boolean to force the evaluation of T in an explicit way
   */
  bool _explicitNewtonEulerDSOperators;

  /** a boolean to know if the matrix W is symmetric definite positive
   */
  bool _isWSymmetricDefinitePositive;

  /** a boolean to perform activation with negative relative velocity
   */
  bool _activateWithNegativeRelativeVelocity;

  /** Constraint activation threshold
   *
   */
  double _constraintActivationThresholdVelocity;

  /**
      A set of work indices for the selected coordinates when
      we subprod in computeFreeOuput
  */
  std::vector<std::size_t> _selected_coordinates;

  /** nslaw effects
   */
  // struct _NSLEffectOnFreeOutput;
  struct _NSLEffectOnFreeOutput : public SiconosVisitor {
    using SiconosVisitor::visit;

    OneStepNSProblem &_osnsp;
    Interaction &_inter;
    InteractionProperties &_interProp;
    double _theta;

    _NSLEffectOnFreeOutput(OneStepNSProblem &p, Interaction &inter,
                           InteractionProperties &interProp,
			   double theta)
      : _osnsp(p), _inter(inter), _interProp(interProp), _theta(theta) {};

    void visit(const NewtonImpactNSL &nslaw);
    void visit(const RelayNSL &nslaw);
    void visit(const NewtonImpactFrictionNSL &nslaw);
    void visit(const FremondImpactFrictionNSL &nslaw);
    void visit(const NewtonImpactRollingFrictionNSL &nslaw);
    void visit(const EqualityConditionNSL &nslaw);
    void visit(const MixedComplementarityConditionNSL &nslaw);
    void visit(const ComplementarityConditionNSL &nslaw);
  };

  friend struct _NSLEffectOnFreeOutput;

public:
  enum MoreauJeanOSI_ds_workVector_id {
    RESIDU_FREE,
    VFREE,
    BUFFER,
    QTMP,
    WORK_LENGTH
  };

  enum MoreauJeanOSI_interaction_workVector_id {
    OSNSP_RHS,
    WORK_INTERACTION_LENGTH
  };

  enum MoreauJeanOSI_interaction_workBlockVector_id {
    xfree,
    BLOCK_WORK_LENGTH
  };

  /** constructor from theta value only
   *
   *  \param theta value for all linked DS (default = 0.5).
   *  \param gamma value for all linked DS (default = NaN and gamma is not
   *  used).
   */
  MoreauJeanOSI(double theta = 0.5,
                double gamma = std::numeric_limits<double>::quiet_NaN());

  /** destructor
   */
  virtual ~MoreauJeanOSI(){};

  // --- GETTERS/SETTERS ---

  /** get the value of W corresponding to DynamicalSystem ds
   *
   *  \param ds a pointer to DynamicalSystem, optional, default =
   *  nullptr. get W[0] in that case
   *  \return SimpleMatrix
   */
  const SimpleMatrix getW(SP::DynamicalSystem ds = SP::DynamicalSystem());

  /** get W corresponding to DynamicalSystem ds
   *
   *  \param ds a pointer to DynamicalSystem
   *  \return pointer to a SiconosMatrix
   */
  SP::SimpleMatrix W(SP::DynamicalSystem ds);

  inline bool isWSymmetricDefinitePositive() const
  {
    return _isWSymmetricDefinitePositive;
  };

  inline void setIsWSymmetricDefinitePositive(bool b)
  {
    _isWSymmetricDefinitePositive = b;
  };

  // -- WBoundaryConditions --

  /** Get the value of WBoundaryConditions corresponding to DynamicalSystem ds
   *
   *  \param ds a pointer to DynamicalSystem, optional, default =
   *  nullptr. get WBoundaryConditions[0] in that case
   *  \return SimpleMatrix
   */
  const SimpleMatrix
  getWBoundaryConditions(SP::DynamicalSystem ds = SP::DynamicalSystem());

  /** get WBoundaryConditions corresponding to DynamicalSystem ds
   *
   *  \param ds a pointer to DynamicalSystem, optional, default =
   *  nullptr. get WBoundaryConditions[0] in that case
   *  \return pointer to a SiconosMatrix
   */
  SP::SiconosMatrix WBoundaryConditions(SP::DynamicalSystem ds);

  // -- theta --

  /** get theta
   *
   *  \return a double
   */
  inline double theta() { return _theta; };

  /** set the value of theta
   *
   *  \param newTheta a double
   */
  inline void setTheta(double newTheta) { _theta = newTheta; };

  // -- gamma --

  /** get gamma
   *
   *  \return a double
   */
  inline double gamma() { return _gamma; };

  /** set the value of gamma
   *
   *  \param newGamma a double
   */
  inline void setGamma(double newGamma)
  {
    _gamma = newGamma;
    _useGamma = true;
  };

  // -- useGamma --

  /** get bool useGamma
   *
   *  \return a bool
   */
  inline bool useGamma() { return _useGamma; };

  /** set the Boolean to indicate that we use gamma
   *
   *  \param newUseGamma  a  Boolean variable
   */
  inline void setUseGamma(bool newUseGamma) { _useGamma = newUseGamma; };

  /** get bool gammaForRelation for the relation
   *
   *  \return a Boolean
   */
  inline bool useGammaForRelation() { return _useGammaForRelation; };

  /** set the boolean to indicate that we use gamma for the relation
   *
   *  \param newUseGammaForRelation a Boolean
   */
  inline void setUseGammaForRelation(bool newUseGammaForRelation)
  {
    _useGammaForRelation = newUseGammaForRelation;
    if (_useGammaForRelation)
      _useGamma = false;
  };
  /** set the constraint activation threshold */
  inline void setConstraintActivationThreshold(double v)
  {
    _constraintActivationThreshold = v;
  }

  /** get the constraint activation threshold */
  inline double constraintActivationThreshold()
  {
    return _constraintActivationThreshold;
  }
  /** set the constraint activation threshold */
  inline void setConstraintActivationThresholdVelocity(double v)
  {
    _constraintActivationThresholdVelocity = v;
  }

  /** get the constraint activation threshold */
  inline double constraintActivationThresholdVelocity()
  {
    return _constraintActivationThresholdVelocity;
  }

  /** get boolean _explicitNewtonEulerDSOperators for the relation
   *
   *  \return a Boolean
   */
  inline bool explicitNewtonEulerDSOperators()
  {
    return _explicitNewtonEulerDSOperators;
  };

  /** set the boolean to indicate that we use gamma for the relation
   *
   *  \param newExplicitNewtonEulerDSOperators a Boolean
   */
  inline void
  setExplicitNewtonEulerDSOperators(bool newExplicitNewtonEulerDSOperators)
  {
    _explicitNewtonEulerDSOperators = newExplicitNewtonEulerDSOperators;
  };

  /** get boolean _activateWithNegativeRelativeVelocity
   *
   *  \return a Boolean
   */
  inline bool activateWithNegativeRelativeVelocity()
  {
    return _activateWithNegativeRelativeVelocity;
  };

  /** set the boolean to perform activation with negative relative velocity
   *
   *  \param newActivateWithNegativeRealtiveVelocity a Boolean
   */
  inline void
  setActivateWithNegativeRelativeVelocity(bool newActivateWithNegativeRelativeVelocity)
  {
    _activateWithNegativeRelativeVelocity = newActivateWithNegativeRelativeVelocity;
  };


  // --- OTHER FUNCTIONS ---

  /**
     initialization of the MoreauJeanOSI integrator; for linear time
     invariant systems, we compute time invariant operator (example :
     W)
   */
  // virtual void initialize(Model& m);

  /**
     Initialization process of the nonsmooth problems
      linked to this OSI*/
  void initialize_nonsmooth_problems() override;

  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi
   *
   *  \param t time of initialization
   *  \param ds the dynamical system
   */
  void initializeWorkVectorsForDS(double t, SP::DynamicalSystem ds) override;

  /** initialization of the work vectors and matrices (properties) related to
   *  one interaction on the graph and needed by the osi
   *
   *  \param inter the interaction
   *  \param interProp the properties on the graph
   *  \param DSG the dynamical systems graph
   */
  void initializeWorkVectorsForInteraction(Interaction &inter,
                                           InteractionProperties &interProp,
                                           DynamicalSystemsGraph &DSG) override;

  /** get the number of index sets required for the simulation
   *
   *  \return unsigned int
   */
  unsigned int numberOfIndexSets() const override { return 2; };

  /** initialize iteration matrix W MoreauJeanOSI matrix at time t
   *
   *  \param time
   *  \param ds a pointer to DynamicalSystem
   */
  void initializeIterationMatrixW(double time, SP::SecondOrderDS ds);

  /** compute W MoreauJeanOSI matrix at time t
   *
   *  \param time (double)
   *  \param ds a  DynamicalSystem
   *  \param W the result in W
   */
  void computeW(double time, SecondOrderDS &ds, SiconosMatrix &W);

  /** get and compute if needed W MoreauJeanOSI matrix at time t
   *
   *  \param time (double)
   *  \param ds a  DynamicalSystem
   *  \param Winverse the result in Winverse
   *  \param keepW
   */
  SP::SimpleMatrix Winverse(SP::SecondOrderDS ds, bool keepW = false);

  /** compute WBoundaryConditionsMap[ds] MoreauJeanOSI matrix at time t
   *
   *  \param ds a pointer to DynamicalSystem
   *  \param WBoundaryConditions write the result in WBoundaryConditions
   *  \param iteration_matrix the OSI iteration matrix (W)
   */
  void _computeWBoundaryConditions(SecondOrderDS &ds,
                                   SiconosMatrix &WBoundaryConditions,
                                   SiconosMatrix &iteration_matrix);

  /** initialize iteration matrix WBoundaryConditionsMap[ds] MoreauJeanOSI
   *
   *  \param ds a pointer to DynamicalSystem
   *  \param dsv a descriptor of the ds on the graph (redundant)
   */
  void _initializeIterationMatrixWBoundaryConditions(
      SecondOrderDS &ds, const DynamicalSystemsGraph::VDescriptor &dsv);

  void applyBoundaryConditions(SecondOrderDS &d, SiconosVector &residu,
                               DynamicalSystemsGraph::VIterator dsi, double t,
                               const SiconosVector &v);

  /** compute the initial state of the Newton loop.
   */
  void computeInitialNewtonState() override;

  /**
     return the maximum of all norms for the "MoreauJeanOSI-discretized"
     residus of DS

     \return a double
   */
  double computeResidu() override;

  /** Perform the integration of the dynamical systems linked to this integrator
   *  without taking into account the nonsmooth input (_r or _p)
   */
  void computeFreeState() override;

  /** integrates the Interaction linked to this integrator, without taking
   *  non-smooth effects into account
   *
   *  \param vertex_inter vertex of the interaction graph
   *  \param osnsp pointer to OneStepNSProblem
   */
  void computeFreeOutput(InteractionsGraph::VDescriptor &vertex_inter,
                         OneStepNSProblem *osnsp) override;

  /** return the workVector corresponding to the right hand side of the OneStepNonsmooth problem
   */
  SiconosVector& osnsp_rhs(InteractionsGraph::VDescriptor& vertex_inter,   InteractionsGraph& indexSet) override
  {
    return *(*indexSet.properties(vertex_inter).workVectors)[MoreauJeanOSI::OSNSP_RHS];
  };

  /** Apply the rule to one Interaction to know if it should be included in the
   *  IndexSet of level i
   *
   *  \param inter the Interaction to test
   *  \param i level of the IndexSet
   *  \return Boolean
   */
  bool addInteractionInIndexSet(SP::Interaction inter, unsigned int i) override;

  /** Apply the rule to one Interaction to know if it should be removed from the
   *  IndexSet of level i
   *
   *  \param inter the Interaction to test
   *  \param i level of the IndexSet
   *  \return Boolean
   */
  bool removeInteractionFromIndexSet(SP::Interaction inter,
                                     unsigned int i) override;

  /** method to prepare the fist Newton iteration
   *   \param time
   */
  void prepareNewtonIteration(double time) override;

  /** integrate the system, between tinit and tend (->iout=true), with possible
   *  stop at tout (->iout=false)
   *
   *  \param tinit the initial time
   *  \param tend the end time
   *  \param tout the real end time
   *  \param notUsed useless flag (for MoreauJeanOSI, used in LsodarOSI)
   */
  void integrate(double &tinit, double &tend, double &tout,
                 int &notUsed) override;

  /**
      update the state of the dynamical systems

      \param ds the dynamical to update
   */
  virtual void updatePosition(DynamicalSystem &ds);

  /** update the state of the dynamical systems
   *
   *  \param level the level of interest for the dynamics: not used at the time
   */
  void updateState(const unsigned int level) override;

  /** Compute the matrix of work of forces by ds
     \return SP::Siconosmatrix
   */
  SP::SimpleMatrix computeWorkForces();

  /** Displays the data of the MoreauJeanOSI's integrator
   */
  void display() override;

  ACCEPT_STD_VISITORS();
};

#endif // MoreauJeanOSI_H
