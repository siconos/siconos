/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
 * D1MinusLinear Time-Integrator for Dynamical Systems
 * T. Schindler, V. Acary
 */

#ifndef D1MINUSLINEAR_H
#define D1MINUSLINEAR_H

#ifdef DEBUG_D1MINUSLINEAR
#define DEBUG_MESSAGES
#endif

#include "OneStepIntegrator.hpp"
#include "SimpleMatrix.hpp"

class SiconosMatrix;

/** D1MinusLinear Time-Integrator for Dynamical Systems
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.3.0.
 *  \date (Creation) September 01, 2011
 *
 *  see Schindler/Acary : Timestepping Schemes for Nonsmooth Dynamics Based
 *  on Discontinuous Galerkin Methods: Definition and Outlook
 *
 *  A D1MinusLinear instance is defined by the list of concerned dynamical systems.
 *
 *  Main functions:
 *
 *  - computeFreeState(): computes xfree (or vfree), dynamical systems
 *    state without taking non-smooth part into account \n
 *
 *  - updateState(): computes x(q,v), the complete dynamical systems
 *    states.
 */
class D1MinusLinear : public OneStepIntegrator
{
protected:

  /** nslaw effects */
  struct _NSLEffectOnFreeOutput;
  friend class _NSLEffectOnFreeOutput;

  /** default constructor */
  D1MinusLinear() {};

public:

  /** constructor from one dynamical system
   *  \param DynamicalSystem to be integrated
   */
  D1MinusLinear(SP::DynamicalSystem);

  /** constructor from a list of dynamical systems
   *  \param list of DynamicalSystems to be integrated
   */
  D1MinusLinear(DynamicalSystemsSet&);

  /** destructor */
  virtual ~D1MinusLinear() {};

  /** initialization of the D1MinusLinear integrator; for linear time
   *  invariant systems, we compute time invariant operator
   */
  virtual void initialize();

  /** return the maximum of all norms for the residus of DS
   *  \post{ds->residuFree will be calculated, ds->p() contains new position, ds->velocity contains predicted velocity}
   *  \return double
   */
  virtual double computeResidu();

  /** integrates the Dynamical System linked to this integrator without taking non-smooth effects into account
   *  \post{ds->velocity contains free velocity}
   */
  virtual void computeFreeState();

  /** integrates the UnitaryRelation linked to this integrator, without taking non-smooth effects into account
   * \param pointer to UnitaryRelation
   * \param pointer to OneStepNSProblem
   */
  virtual void computeFreeOutput(SP::UnitaryRelation UR, OneStepNSProblem* osnsp);

  /** integrate the system, between tinit and tend (->iout=true), with possible stop at tout (->iout=false)
   *  \param initial time
   *  \param end time
   *  \param real end time
   *  \param useless flag (for D1MinusLinear, used in Lsodar)
   */
  virtual void integrate(double&, double&, double&, int&)
  {
    RuntimeException::selfThrow("D1MinusLinear::integrate - not implemented!");
  }

  /** updates the state of the Dynamical Systems
   *  \param level of interest for the dynamics: not used at the time
   *  \post{ds->velocity contains new velocity}
   */
  virtual void updateState(unsigned int);

  /** displays the data of the D1MinusLinear's integrator */
  virtual void display()
  {
    RuntimeException::selfThrow("D1MinusLinear::display - not implemented!");
  }

  /** preparations for Newton iteration */
  virtual void prepareNewtonIteration(double time)
  {
    RuntimeException::selfThrow("D1MinusLinear::prepareNewtonIteration - not implemented!");
  }

  /** insert a dynamical system in this Integrator
   *  \param pointer to DynamicalSystem
   */
  virtual void insertDynamicalSystem(SP::DynamicalSystem ds);

  /** encapsulates an operation of dynamic casting
   *  needed by Python interface
   *  \param integrator which must be converted
   *  \return pointer to the integrator if it is of the right type, NULL otherwise
   */
  static D1MinusLinear* convert(OneStepIntegrator* osi);

  /** visitors hook */
  ACCEPT_STD_VISITORS();
};

#endif // D1MINUSLINEAR_H
