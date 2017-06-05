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

/*!\file
 * D1MinusLinearOSI Time-Integrator for Dynamical Systems
 * T. Schindler, V. Acary
 */

#ifndef D1MINUSLINEAR_H
#define D1MINUSLINEAR_H

#ifdef DEBUG_D1MINUSLINEAR
//#define DEBUG_MESSAGES
#endif

#include "OneStepIntegrator.hpp"
#include "SimpleMatrix.hpp"


/** \class D1MinusLinearOSI D1MinusLinearOSI.hpp Time-Integrator for Dynamical Systems
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.6.0 -- 3.7.x
 *  \date (Creation) September 01, 2011
 *
 * Reference:
 *
 * Timestepping schemes for nonsmooth dynamics based on discontinuous Galerkin methods:
 * Definition and outlook
 * Thorsten Schindler; Vincent Acary
 * Mathematics and Computers in Simulation,
 * Elsevier, 2014, 95, pp. 180-199
 *
 *  A D1MinusLinearOSI instance is defined by the list of concerned dynamical systems.
 *
 *  Main functions:
 *
 *  - computeFreeState(): computes xfree (or vfree), dynamical systems
 *    state without taking non-smooth part into account \n
 *
 *  - updateState(): computes x(q,v), the complete dynamical systems
 *    states.
 *

 *
 * \section sec_D1MinusLinear_implementation Some details on the implementation
 *
 * \subsection sec_D1MinusLinear_implementation_various Various implementations at various levels.
 *
 * To be completed
 *
 * The (lower-case) lambda vector that corresponds to the (non-impulsive) contact forces are computed
 * for the closed contact with vanishing relative velocity (the number of the index set depends on the type of D1minulinear)
 * The result, stored in lambda(2) and p(2) is computed by solving
 *       (*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]
 * at different specific time
 *
 * The impact equation are solved at the velocity level as in MoreauJeanOSI using
 * lambda(1) and p(1). The number of the index set depends on the type of D1minulinear
 *
 *
 * \subsection sec_D1MinusLinear_implementation_halfAcc  HalfExplicitAccelerationLevel
 *
 * The (lower-case) lambda \f$\lambda\f$ that corresponds to the (non-impulsive) contact
 * forces are computed at the acceleration level for the closed contact with vanishing
 * relative velocity (indexSet2)
 *
 * The problem, solved for  lambda(2) and p(2)   (*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]
 * is solved at time told for \f$\lambda^+_{k}\f$ and time t for \f$\lambda^-_{k+1}\f$
 *
 * The problem, solved for  lambda(1) and p(1)   (*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]
 * is solved at time t for \f$Lambda_{k+1}\f$.
 *
 * We use four index sets to compute the multipliers:
 * <ul>
 * <li> indexSet0 is the set of all constraints </li>
 * <li> indexSet1 is the set of closed constraints (computed with \f$ q_{k,1} \f$) </li>
 * <li> indexSet2 is the set of closed constraints and vanishing velocities (computed with \f$q_{k,1}\f$ and \f$v_{k,1}\f$  )  </li>
 * <li> indexSet3 is the set of the constraints that have been closed within the time--step (computed with \f$q_{k,0}\f$ and \f$q_{k,1}\f$  )  </li>
 * </ul>
 *
 * \f[
 * \begin{cases}
 * v_{k,0} = \mbox{vold} \\
 * q_{k,0} = \mbox{qold} \\
 * F^+_{k} = \mbox{F(told,qold,vold)} \\
 * v_{k,1} = v_{k,0} + h M^{-1}_k (P^+_{2,k}+F^+_{k}) \\[2mm]
 * q_{k,1} = q_{k,0} + \frac{h}{2} (v_{k,0} + v_{k,1})  \\[2mm]
 * F^-_{k+1} = F(t^{-}_{k+1},q_{k,1},v_{k,1}) \\[2mm]
 * R_{free} = - \frac{h}{2}  M^{-1}_k (P^+_{2,k}+F^+_{k}) -  \frac{h}{2}  M^{-1}_{k+1} (P^-_{2,k+1}+F^-_{k+1}) \\[2mm]
 * \end{cases}
 * \f]
 *
 * \subsection sec_D1MinusLinear_implementation_halfVel  HalfExplicitVelocityLevel
 *
 *
 */


const double DEFAULT_TOL_D1MINUS  = 1e-8;



class D1MinusLinearOSI : public OneStepIntegrator
{
public :

protected:

  /** nslaw effects */
  struct _NSLEffectOnFreeOutput : public SiconosVisitor
  {

    using SiconosVisitor::visit;

    OneStepNSProblem* _osnsp;
    SP::Interaction _inter;
    InteractionProperties& _interProp;
    _NSLEffectOnFreeOutput(OneStepNSProblem *p, SP::Interaction inter, InteractionProperties& interProp) :
      _osnsp(p), _inter(inter), _interProp(interProp) {};

    void visit(const NewtonImpactNSL& nslaw);
    void visit(const EqualityConditionNSL& nslaw)
    {
      ;
    }
  };

  //friend struct _NSLEffectOnFreeOutput;
  bool _isThereImpactInTheTimeStep ;
  unsigned int _typeOfD1MinusLinearOSI;


public:

  /** Switching variable for various versions of D1MinusLinear
   */
  enum ListOfTypeOfD1MinusLinearOSI {halfexplicit_acceleration_level,
                                     halfexplicit_acceleration_level_full,
                                     explicit_velocity_level,
                                     halfexplicit_velocity_level,
                                     numberOfTypeOfD1MinusLinearOSI
                                    };

  enum {OSNSP_RHS,WORK_INTERACTION_LENGTH};

  /** basic constructor
   */
  D1MinusLinearOSI();

  /** Constructor with type of D1MinusLinear
   * \param type unsigned int that specifies the type of D1MinusLinear
   * D1MinusLinearOSI::halfexplicit_acceleration_level,
   * D1MinusLinearOSI::halfexplicit_acceleration_level_full,
   * D1MinusLinearOSI::explicit_velocity_level,
   * D1MinusLinearOSI::halfexplicit_velocity_level,
   * D1MinusLinearOSI::numberOfTypeOfD1MinusLinearOSI
   */
  D1MinusLinearOSI(unsigned int type);

  /** destructor */
  virtual ~D1MinusLinearOSI() {};

  /** Set the type of the D1MinusLinear
   * \param type  the type to set
   * D1MinusLinearOSI::halfexplicit_acceleration_level,
   * D1MinusLinearOSI::halfexplicit_acceleration_level_full,
   * D1MinusLinearOSI::explicit_velocity_level,
   * D1MinusLinearOSI::halfexplicit_velocity_level,
   * D1MinusLinearOSI::numberOfTypeOfD1MinusLinearOSI
   */
  void setTypeOfD1MinusLinearOSI(unsigned int type);

  /** get the type of the D1MinusLinear
   * \return unsigned int type  the type to set
   * D1MinusLinearOSI::halfexplicit_acceleration_level,
   * D1MinusLinearOSI::halfexplicit_acceleration_level_full,
   * D1MinusLinearOSI::explicit_velocity_level,
   * D1MinusLinearOSI::halfexplicit_velocity_level,
   * D1MinusLinearOSI::numberOfTypeOfD1MinusLinearOSI
   */
  unsigned int typeOfD1MinusLinearOSI()
  {
    return _typeOfD1MinusLinearOSI;
  };

  /** get the number of index sets required for the simulation
   * \return unsigned int
   */
  unsigned int numberOfIndexSets() const;

  /** initialization of the D1MinusLinearOSI integrator; for linear time
   *  invariant systems, we compute time invariant operator
   */
  // virtual void initialize(Model& m);
  virtual void initialize_nonsmooth_problems();

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
  void fillDSLinks(Interaction &inter, InteractionProperties& interProp,
			     DynamicalSystemsGraph & DSG);
  
  /** return the maximum of all norms for the residus of DS
   *  \post{ds->residuFree will be calculated, ds->q() contains new position, ds->velocity contains predicted velocity}
   *  \return double
   */
  virtual double computeResidu();

  /** return the maximum of all norms for the residus of DS for the type explicit_acceleration_level
   *  \post{ds->residuFree will be calculated, ds->q() contains new position, ds->velocity contains predicted velocity}
   *  \return double
   */
  virtual double computeResiduHalfExplicitAccelerationLevel();

  /** return the maximum of all norms for the residus of DS for the type explicit_acceleration_level_full
   *  \post{ds->residuFree will be calculated, ds->q() contains new position, ds->velocity contains predicted velocity}
   *  \return double
   */
  virtual double computeResiduHalfExplicitAccelerationLevelFull();

  /** return the maximum of all norms for the residus of DS for the type explicit_acceleration_level_full
   *  \post{ds->residuFree will be calculated, ds->q() contains new position, ds->velocity contains predicted velocity}
   *  \return double
   */
  virtual double computeResiduHalfExplicitVelocityLevel();


  /** integrates the Dynamical System linked to this integrator without taking non-smooth effects into account
   *  \post{ds->velocity contains predicted velocity}
   */
  virtual void computeFreeState();

  /** integrates the Interaction linked to this integrator, without taking non-smooth effects into account
   * \param vertex_inter of the interaction graph
   * \param osnsp pointer to OneStepNSProblem
   */
  virtual void computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp);

  /** integrates the Interaction linked to this integrator, without taking non-smooth effects into account
   * \param vertex_inter of the interaction graph
   * \param osnsp pointer to OneStepNSProblem
   */
  virtual void computeFreeOutputHalfExplicitAccelerationLevel(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp);

  /** integrates the Interaction linked to this integrator, without taking non-smooth effects into account
   * \param vertex_inter of the interaction graph
   * \param osnsp pointer to OneStepNSProblem
   */
  virtual void computeFreeOutputHalfExplicitVelocityLevel(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp);



  /** integrate the system, between tinit and tend (->iout=true), with possible stop at tout (->iout=false)
   *  \param ti initial time
   *  \param tf end time
   *  \param t real end time
   *  \param flag useless flag (for D1MinusLinearOSI, used in LsodarOSI)
   */
  virtual void integrate(double& ti, double& tf, double& t , int& flag)
  {
    RuntimeException::selfThrow("D1MinusLinearOSI::integrate - not implemented!");
  }

  /** updates the state of the Dynamical Systems
   *  \param level level of interest for the dynamics: not used at the time
   *  \post{ds->velocity contains new velocity}
   */
  virtual void updateState(const unsigned int level);


  /** Apply the rule to one Interaction to known if is it should be included
   * in the IndexSet of level i
   * \param inter the involved interaction
   * \param i the index set level
   * \return a boolean if it needs to be added or not
   */
  virtual bool addInteractionInIndexSet(SP::Interaction inter, unsigned int i);

  /** Apply the rule to one Interaction to known if is it should be removed
   * in the IndexSet of level i
   * \param inter the involved interaction
   * \param i the index set level
   * \return a boolean if it needs to be removed or not
  */
  virtual bool removeInteractionInIndexSet(SP::Interaction inter, unsigned int i);

  /** Apply the rule to one Interaction to known if is it should be included
   * in the IndexSet of level i
   * \param inter the involved interaction
   * \param i the index set level
   * \return a boolean if it needs to be added or not
   */
  virtual bool addInteractionInIndexSetHalfExplicitAccelerationLevel(SP::Interaction inter, unsigned int i);

  /** Apply the rule to one Interaction to known if is it should be removed
   * in the IndexSet of level i
   * \param inter the involved interaction
   * \param i the index set level
   * \return a boolean if it needs to be removed or not
  */
  virtual bool removeInteractionInIndexSetHalfExplicitAccelerationLevel(SP::Interaction inter, unsigned int i);

  /** Apply the rule to one Interaction to known if is it should be included
    * in the IndexSet of level i
    * \param inter the involved interaction
    * \param i the index set level
    * \return a boolean if it needs to be added or not
    */
  virtual bool addInteractionInIndexSetHalfExplicitVelocityLevel(SP::Interaction inter, unsigned int i);

  /** Apply the rule to one Interaction to known if is it should be removed
   * in the IndexSet of level i
   * \param inter the involved interaction
   * \param i the index set level
   * \return a boolean if it needs to be removed or not
  */
  virtual bool removeInteractionInIndexSetHalfExplicitVelocityLevel(SP::Interaction inter, unsigned int i);

  /** displays the data of the D1MinusLinearOSI's integrator */
  virtual void display()
  {
    RuntimeException::selfThrow("D1MinusLinearOSI::display - not implemented!");
  }

  /** preparations for Newton iteration
   *  \param time time
   */
  virtual void prepareNewtonIteration(double time)
  {
    RuntimeException::selfThrow("D1MinusLinearOSI::prepareNewtonIteration - not implemented!");
  }

  /** visitors hook */
  ACCEPT_STD_VISITORS();
};

#endif // D1MINUSLINEAR_H
