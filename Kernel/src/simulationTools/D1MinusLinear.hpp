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
 *  Each DynamicalSystem is associated to a SiconosMatrix, named "W"
 *
 *  W matrices are initialized and computed in initW and
 *  computeW. Depending on the DS type, they may depend on time and DS
 *  state (x).
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

  /** Stl map that associates a W D1MinusLinear matrix to each DynamicalSystem of the OSI */
  MapOfDSMatrices WMap;

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

  // --- GETTERS/SETTERS ---

  /** get value of W corresponding to DynamicalSystem ds
   *  \param pointer to DynamicalSystem, optional, default = NULL. get W[0] in that case
   *  \return SimpleMatrix
   */
  const SimpleMatrix getW(SP::DynamicalSystem = SP::DynamicalSystem());

  /** get W corresponding to DynamicalSystem ds
   * \param pointer to DynamicalSystem, optional, default = NULL. get W[0] in that case
   * \return pointer to SiconosMatrix
   * \todo identifier W is in general not a good choice, e.g. it is the Delassus operator in Acary2008 or global-local velocity projection matrix in Pfeiffer2008; I suggest e.g. IterMat (Thorsten Schindler, 01.09.2011)
   */
  SP::SimpleMatrix W(SP::DynamicalSystem ds);

  /** set the value of W[ds] to newValue
   * \param SiconosMatrix newValue
   * \param pointer to DynamicalSystem,
   */
  void setW(const SiconosMatrix&, SP::DynamicalSystem);

  /** set W[ds] to pointer newPtr
   * \param newPtr
   * \param pointer to DynamicalSystem
   */
  void setWPtr(SP::SimpleMatrix newPtr, SP::DynamicalSystem);

  // --- OTHER FUNCTIONS ---

  /** initialization of the D1MinusLinear integrator; for linear time
   *  invariant systems, we compute time invariant operator (example : W)
   */
  virtual void initialize();

  /** init WMap[ds] D1MinusLinear matrix at time t
   *  allocate memory for W and insert into WMap with ds as key
   *  \param time
   *  \param pointer to DynamicalSystem
   */
  void initW(double, SP::DynamicalSystem);

  /** compute WMap[ds] D1MinusLinear matrix at time t
   *  \param time
   *  \param pointer to DynamicalSystem
   */
  void computeW(double, SP::DynamicalSystem);

  /** return the maximum of all norms for the residus of DS
   *  \post{ds->residuFree will be calculated, ds->workFree contains ds->residuFree-p, ds->p() contains new position}
   *  \todo T is used just explicitely in position update of Newton Euler; it has to be included locally in fixed point iteration considering the norm of the position residuum
   *  \return double
   */
  virtual double computeResidu();

  /** integrates the Dynamical System linked to this integrator without boring the constraints
   *  \post{ds->workFree contains free velocity}
   */
  virtual void computeFreeState();

  /** integrates the UnitaryRelation linked to this integrator, without taking constraints into account
   * \param pointer to UnitaryRelation
   * \param pointer to OneStepNSProblem
   */
  virtual void computeFreeOutput(SP::UnitaryRelation UR, OneStepNSProblem* osnsp);

  void prepareNewtonIteration(double time);


  /** integrate the system, between tinit and tend (->iout=true), with possible stop at tout (->iout=false)
   *  \param initial time
   *  \param end time
   *  \param real end time
   *  \param useless flag (for D1MinusLinear, used in Lsodar)
   */
  virtual void integrate(double&, double&, double&, int&)
  {
    RuntimeException::selfThrow("D1MinusLinear::integrate - not yet implemented!");
  }

  /** updates the state of the Dynamical Systems
   *  \param level of interest for the dynamics: not used at the time
   *  \post{ds->velocity contains new velocity}
   */
  virtual void updateState(unsigned int);

  /** displays the data of the D1MinusLinear's integrator */
  virtual void display();

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

  /** visitors hook
   */
  ACCEPT_STD_VISITORS();


};

#endif // D1MINUSLINEAR_H
