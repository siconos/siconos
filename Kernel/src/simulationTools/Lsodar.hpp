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
 *
 * WARNING: at the time mainly written for Lagrangian systems !!!
 */
/*! \file
  Lsodar solver (from odepack)
*/
#ifndef Lsodar_H
#define Lsodar_H

#include"OneStepIntegrator.hpp"
#include"SiconosNumerics.h"
#include<vector>
const doublereal ATOL_DEFAUTL = 100 * MACHINE_PREC;
const doublereal RTOL_DEFAULT = 10 * MACHINE_PREC;
class BlockVector;

/** Lsodar solver (odepack)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 * Many parameters are required as input/output for LSODAR. See the documentation of this function
 * in Numerics/src/odepack/opkdmain.f to have a full description of these parameters.  \n
 * Most of them are read-only parameters (ie can not be set by user). \n
 *  Except: \n
 *  - jt: Jacobian type indicator (1 means a user-supplied full Jacobian, 2 means an internally generated full Jacobian). \n
 *    Default = 2.
 *  - itol, rtol and atol \n
 *    ITOL   = an indicator for the type of error control. \n
 *    RTOL   = a relative error tolerance parameter, either a scalar or array of length NEQ. \n
 *    ATOL   = an absolute error tolerance parameter, either a scalar or an array of length NEQ.  Input only.
 */
class Lsodar : public OneStepIntegrator, public boost::enable_shared_from_this<Lsodar>
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Lsodar);


  /** neq, ng, itol, itask, istate, iopt, lrw, liw, jt
   * See opkdmain.f and lsodar routine for details on those variables.
   */
  std::vector<integer> intData;
  /** relative tolerance */
  SA::doublereal rtol;
  /** absolute tolerance */
  SA::doublereal atol;
  /** real work array */
  SA::doublereal rwork;
  /** integer work array */
  SA::integer iwork;
  /** integer array used for output of root information */
  SA::integer jroot;
  /** temporary vector to save x values */
  SP::BlockVector xWork;
  Lsodar() {};

  /** nslaw effects
   */
  struct _NSLEffectOnFreeOutput;
  friend class _NSLEffectOnFreeOutput;


public:

  /** constructor from xml file
      \param OneStepIntegratorXML* : the XML object
      \param the set of all DS in the NSDS
      \param the set of all interactions in the NSDS
  */
  Lsodar(SP::OneStepIntegratorXML, SP::DynamicalSystemsSet , SP::InteractionsSet);

  /** constructor from a minimum set of data
   *  \param SP::DynamicalSystem : the DynamicalSystem linked to the OneStepIntegrator
   */
  Lsodar(SP::DynamicalSystem);

  /** constructor from a list of Dynamical Systems
   *  \param DynamicalSystemsSet : the list of DynamicalSystems to be integrated
   */
  Lsodar(DynamicalSystemsSet&);

  /** destructor
   */
  ~Lsodar() {};

  /** get vector of integer parameters for lsodar
   *  \return a vector<integer>
   */
  inline const std::vector<integer> getIntData() const
  {
    return intData;
  }

  /** get intData[i]
   *  \return an integer
   */
  inline integer getIntData(unsigned int i) const
  {
    return intData[i];
  }

  /** get relative tolerance parameter for lsodar
   *  \return a doublereal*
   */
  inline const SA::doublereal getRtol() const
  {
    return rtol;
  }

  /** get absolute tolerance parameter for lsodar
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
    return iwork[5];
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

  /** get output of root information
   *  \return a pointer to integer
   */
  inline SA::integer getJroot() const
  {
    return jroot;
  }

  /** set Jt value, Jacobian type indicator
   *  \param pointer to integer
   */
  inline void setJT(integer newValue)
  {
    intData[8] = newValue;
  };

  /** set itol, rtol and atol (tolerance parameters for lsodar)
   *  \param integer (itol value)
   *  \param doublereal * (rtol)
   *  \param doublereal * (atol)
   */
  void setTol(integer, SA::doublereal, SA::doublereal);

  /** set itol, rtol and atol (scalar tolerance parameters for lsodar)
   *  \param integer (itol value)
   *  \param double (rtol)
   *  \param double (atol)
   */
  void setTol(integer, doublereal, doublereal);

  /** set the maximul number of steps for one call of Lsodar
   *\param an integer
   */
  void setMaxNstep(int);

  /** update doubleData and iwork memory size, when changes occur in intData.
   */
  void updateData();

  /** fill xWork with a doublereal
   *  \param integer*, size of x array
   *  \param doublereal* x:array of double
   */
  void fillXWork(integer*, doublereal*) ;

  /** compute rhs(t) for all dynamical systems in the set
   */
  void computeRhs(double) ;

  /** compute jacobian of the rhs at time t for all dynamical systems in the set
   */
  void computeJacobianRhs(double) ;

  void f(integer* sizeOfX, doublereal* time, doublereal* x, doublereal* xdot);

  void g(integer* nEq, doublereal* time, doublereal* x, integer* ng, doublereal* gOut);

  void jacobianfx(integer*, doublereal*, doublereal*, integer*, integer*,  doublereal*, integer*);

  /** initialization of the integrator
   */
  void initialize();

  /** integrate the system, between tinit and tend (->iout=true), with possible stop at tout (->iout=false)
   *  \param double: tinit, initial time
   *  \param double: tend, end time
   *  \param double: tout, real end time
   *  \param int&: in-out parameter, input: 1 for first call, else 2. Output: 2 if no root was found, else 3.
   */
  void integrate(double&, double&, double&, int&);

  /** update the state of the DynamicalSystems attached to this Integrator
   *  \param unsigned int: level of interest for the dynamics
   */
  void updateState(unsigned int);

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, NULL otherwise
   */
  //static Lsodar* convert (OneStepIntegrator* osi);

  /** integrates the UnitaryRelation linked to this integrator, without taking constraints
     * into account.
     */
  virtual void computeFreeOutput(SP::UnitaryRelation UR, OneStepNSProblem * osnsp);

  /** print the data to the screen
   */
  void display();
};

#endif // Lsodar_H
