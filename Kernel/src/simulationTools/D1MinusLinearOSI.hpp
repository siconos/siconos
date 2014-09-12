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

/** D1MinusLinearOSI Time-Integrator for Dynamical Systems
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.6.0.
 *  \date (Creation) September 01, 2011
 *
 *  see Schindler/Acary : Timestepping Schemes for Nonsmooth Dynamics Based
 *  on Discontinuous Galerkin Methods: Definition and Outlook
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
 */
class D1MinusLinearOSI : public OneStepIntegrator
{
protected:

  /** nslaw effects */
  struct _NSLEffectOnFreeOutput;
  friend struct _NSLEffectOnFreeOutput;
  bool _isThereImpactInTheTimeStep ;
  unsigned int _typeOfD1MinusLinearOSI;

  /** Switching variable for various versions of D1MinusLinear
   */
  enum ListOfTypeOfD1MinusLinearOSI {explicit_acceleration_level,
                                     explicit_acceleration_level_full,
                                     explicit_velocity_level,
                                     halfexplicit_velocity_level,
                                     numberOfTypeOfD1MinusLinearOSI };


public:

  /** basic constructor
   */
  D1MinusLinearOSI();

  /** constructor from one dynamical system
   *  \param newDS DynamicalSystem to be integrated
   */
  DEPRECATED_OSI_API(D1MinusLinearOSI(SP::DynamicalSystem newDS));

  /** destructor */
  virtual ~D1MinusLinearOSI() {};

  /** initialization of the D1MinusLinearOSI integrator; for linear time
  c   *  invariant systems, we compute time invariant operator
   */
  virtual void initialize();

  /** return the maximum of all norms for the residus of DS
   *  \post{ds->residuFree will be calculated, ds->q() contains new position, ds->velocity contains predicted velocity}
   *  \return double
   */
  virtual double computeResidu();

  /** return the maximum of all norms for the residus of DS for the type explicit_acceleration_level
   *  \post{ds->residuFree will be calculated, ds->q() contains new position, ds->velocity contains predicted velocity}
   *  \return double
   */
  virtual double computeResidu_explicit_acceleration_level();

  /** return the maximum of all norms for the residus of DS for the type explicit_acceleration_level_full
   *  \post{ds->residuFree will be calculated, ds->q() contains new position, ds->velocity contains predicted velocity}
   *  \return double
   */
  virtual double computeResidu_explicit_acceleration_level_full();


  /** integrates the Dynamical System linked to this integrator without taking non-smooth effects into account
   *  \post{ds->velocity contains predicted velocity}
   */
  virtual void computeFreeState();

  /** integrates the Interaction linked to this integrator, without taking non-smooth effects into account
   * \param vertex_inter of the interaction graph
   * \param osnsp pointer to OneStepNSProblem
   */
  virtual void computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp);

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
   */
  virtual bool addInteractionInIndexSet(SP::Interaction inter, unsigned int i);

  /** Apply the rule to one Interaction to known if is it should be removed
   * in the IndexSet of level i
   */
  virtual bool removeInteractionInIndexSet(SP::Interaction inter, unsigned int i);




  /** displays the data of the D1MinusLinearOSI's integrator */
  virtual void display()
  {
    RuntimeException::selfThrow("D1MinusLinearOSI::display - not implemented!");
  }

  /** preparations for Newton iteration */
  virtual void prepareNewtonIteration(double time)
  {
    RuntimeException::selfThrow("D1MinusLinearOSI::prepareNewtonIteration - not implemented!");
  }

  /** visitors hook */
  ACCEPT_STD_VISITORS();
};

#endif // D1MINUSLINEAR_H
