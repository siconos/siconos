/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 *
 * WARNING: at the time mainly written for Lagrangian systems !!!
 */
/*! \file
 Lsodar solver (from odepack)
*/
#ifndef Lsodar_H
#define Lsodar_H

#include"OneStepIntegrator.h"
#include"SiconosNumerics.h"
#include<vector>

class BlockVector;

/** Lsodar solver (odepack)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.0.
 *  \date (Creation) Apr 26, 2004
 */
class Lsodar : public OneStepIntegrator
{
private:

  /** The time discretisation, specific to Lsodar
   *  This discretisation must include each time step of the simulation time
   *  discretisation, plus all event detection time steps.
   */
  TimeDiscretisation *localTimeDiscretisation;

  /** a bool to check whether localTimeDiscretisation
   * has been allocated inside the current class or not
   */
  bool isLocalTimeDiscretisationAllocatedIn;

  /** neq, ng, itol, itask, istate, iopt, lrw, liw, jt
   * See opkdmain.f and lsodar routine for details on those variables.
   */
  std::vector<integer> intData;
  /** rtol, atol, rwork, jroot */
  std::vector<doublereal*> doubleData;
  /** iwork */
  integer * iwork;

  /** temporary vector to save x values */
  BlockVector* xWork;

  /** default constructor
  */
  Lsodar();

public:

  /** constructor from xml file
  *  \param OneStepIntegratorXML* : the XML object
  *  \param Simulation * : the simulation that owns the osi
  */
  Lsodar(OneStepIntegratorXML*, Simulation*);

  /** constructor from a minimum set of data
  *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
  *  \param Simulation * : the simulation that owns the osi
  */
  Lsodar(DynamicalSystem* , Simulation*);

  /** constructor from a list of Dynamical Systems
  *  \param DynamicalSystemsSet : the list of DynamicalSystems to be integrated
  *  \param Simulation * : the simulation that owns the osi
  */
  Lsodar(DynamicalSystemsSet&, Simulation*);

  /** destructor
  */
  ~Lsodar();

  /** get the TimeDiscretisation of lsodar
  *  \return the TimeDiscretisation
  */
  inline TimeDiscretisation* getTimeDiscretisationPtr() const
  {
    return localTimeDiscretisation;
  };

  /** set timeDiscretisation of lsodar
  *  \param the TimeDiscretisation to set
  */
  void setTimeDiscretisationPtr(TimeDiscretisation*);

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
  inline const integer getIntData(const unsigned int i) const
  {
    return intData[i];
  }

  /** set vector intData to newVector with a copy.
  *  \param std::vector<integer>
  */
  void setIntData(const std::vector<integer>&);

  /** set intData[i] to newValue
  *  \param a unsigned int (index) and an integer (value)
  */
  inline void setIntData(const unsigned int  i, const integer  newValue)
  {
    intData[i] = newValue;
  }

  /** get vector of doublereal* parameters for lsodar
  *  \return a vector<doublereal*>
  */
  inline const std::vector<doublereal*> getDoubleData() const
  {
    return doubleData;
  }

  /** get doubleData[i]
  *  \return a pointer on doublereal.
  */
  inline doublereal* getDoubleData(const unsigned int i) const
  {
    return doubleData[i];
  }

  /** set vector doubleData to newVector with a copy -> memory allocation
  *  \param std::vector<doublereal*>
  */
  void setDoubleData(const std::vector<doublereal*>&);

  /** set doubleData[i] to newPtr
  *  \param a unsigned int (index) and a pointer to doublereal
  */
  // void setDoubleData(const unsigned int , doublereal*);

  /** get iwork
  *  \return a pointer to integer
  */
  inline integer* getIwork() const
  {
    return iwork;
  }

  /** set iwork to newValue with a copy.
  *  \param pointer to integer
  */
  void setIwork(integer*);

  /** update doubleData and iwork memory size, when changes occur in intData.
  */
  void updateData();

  /** fill xWork with a doublereal
  *  \param integer*, size of x array
  *  \param doublereal* x:array of double
  */
  void fillXWork(integer*, doublereal *) ;

  /** compute rhs(t) for all dynamical systems in the set
  */
  void computeRhs(const double) ;

  /** compute jacobian of the rhs at time t for all dynamical systems in the set
  */
  void computeJacobianRhs(const double) ;

  void f(integer * sizeOfX, doublereal * time, doublereal * x, doublereal * xdot);

  void g(integer * nEq, doublereal * time, doublereal* x, integer * ng, doublereal * gOut);

  void jacobianF(integer *, doublereal *, doublereal *, integer *, integer *,  doublereal *, integer *);

  /** initialise the integrator
  */
  void initialize();

  /** compute the free state of the dynamical system
  */
  void computeFreeState();

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
  void updateState(const unsigned int);

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
  *  \param OneStepIntegrator* : the integrator which must be converted
  * \return a pointer on the integrator if it is of the right type, NULL otherwise
  */
  //static Lsodar* convert (OneStepIntegrator* osi);

  /** print the data to the screen
  */
  void display();


};

#endif // Lsodar_H
