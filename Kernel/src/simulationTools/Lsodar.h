/* Siconos version 1.0, Copyright INRIA 2005.
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
 */
#ifndef Lsodar_H
#define Lsodar_H

#include "OneStepIntegrator.h"
#include "LsodarXML.h"

/** \class Lsodar
 *  \brief It's a kind of single-step Integrator
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 */

// ===== Lsodar class =====
class Lsodar : public OneStepIntegrator
{
private:

  /** The time discretisation, specific to Lsodar
   *  This discretisation must include each time step of the strategy time
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

  /** \fn Lsodar()
   *  \brief default constructor
   */
  Lsodar();

public:

  /** \fn Lsodar(OneStepIntegratorXML*)
   *  \brief constructor from xml file
   *  \param OneStepIntegratorXML* : the XML object
   */
  Lsodar(OneStepIntegratorXML*);

  /** \fn Lsodar(TimeDiscretisation*, DynamicalSystem* )
   *  \brief constructor from a minimum set of data
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   */
  Lsodar(TimeDiscretisation*, DynamicalSystem*);

  ~Lsodar();

  /** \fn TimeDiscretisation* getTimeDiscretisationPtr()
   *  \brief get the TimeDiscretisation of lsodar
   *  \return the TimeDiscretisation
   */
  inline TimeDiscretisation* getTimeDiscretisationPtr() const
  {
    return localTimeDiscretisation;
  };

  /** \fn void setTimeDiscretisationPtr(TimeDiscretisation*)
   *  \brief set timeDiscretisation of lsodar
   *  \param the TimeDiscretisation to set
   */
  void setTimeDiscretisationPtr(TimeDiscretisation*);

  /** \fn  const vector<integer> getIntData() const
   *  \brief get vector of integer parameters for lsodar
   *  \return a vector<integer>
   */
  inline const std::vector<integer> getIntData() const
  {
    return intData;
  }

  /** \fn  const integer getIntData(const unsigned int & i) const
   *  \brief get intData[i]
   *  \return an integer
   */
  inline const integer getIntData(const unsigned int& i) const
  {
    return intData[i];
  }

  /** \fn void setIntData (const vector<integer>& newVector)
   *  \brief set vector intData to newVector with a copy.
   *  \param std::vector<integer>
   */
  void setIntData(const std::vector<integer>&);

  /** \fn void setIntData (const unsigned int & i, const integer& newValue);
   *  \brief set intData[i] to newValue
   *  \param a unsigned int (index) and an integer (value)
   */
  inline void setIntData(const unsigned int & i, const integer & newValue)
  {
    intData[i] = newValue;
  }

  /** \fn  const vector<doublereal*> getDoubleData() const
   *  \brief get vector of doublereal* parameters for lsodar
   *  \return a vector<doublereal*>
   */
  inline const std::vector<doublereal*> getDoubleData() const
  {
    return doubleData;
  }

  /** \fn  const doublereal getDoubleData(const unsigned int & i) const
   *  \brief get doubleData[i]
   *  \return a pointer on doublereal.
   */
  inline doublereal* getDoubleData(const unsigned int& i) const
  {
    return doubleData[i];
  }

  /** \fn void setDoubleData (const vector<doublereal*>& newVector)
   *  \brief set vector doubleData to newVector with a copy -> memory allocation
   *  \param std::vector<doublereal*>
   */
  void setDoubleData(const std::vector<doublereal*>&);

  /** \fn void setDoubleData (const unsigned int & i, doublereal* newPtr);
   *  \brief set doubleData[i] to newPtr
   *  \param a unsigned int (index) and a pointer to doublereal
   */
  // void setDoubleData(const unsigned int &, doublereal*);

  /** \fn integer* getIwork() const
   *  \brief get iwork
   *  \return a pointer to integer
   */
  inline integer* getIwork() const
  {
    return iwork;
  }

  /** \fn void setIwork (integer*)
   *  \brief set iwork to newValue with a copy.
   *  \param pointer to integer
   */
  void setIwork(integer*);

  /** \fn void updateData()
   *  \brief update doubleData and iwork memory size, when changes occur in intData.
   */
  void updateData();

  void f(integer * sizeOfX, doublereal * time, doublereal * x, doublereal * xdot);

  void g(integer * nEq, doublereal * time, doublereal* x, integer * ng, doublereal * gOut);

  void jacobianF(integer *, doublereal *, doublereal *, integer *, integer *,  doublereal *, integer *);

  /** \fn void initialize()
   *  \brief initialise the integrator
   */
  void initialize();

  /** \fn void computeFreeState()
   *   \brief compute the free state of the dynamical system
   */
  void computeFreeState();

  /** \fn void integrate()
   *   \brief integrates the dynamical system
   */
  void integrate();

  /** \fn void updateState()
   *  \brief update the state of the DynamicalSystem attached to this Integrator
   */
  void updateState();

  /** \fn Lsodar* convert (OneStepIntegrator* osi)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, NULL otherwise
   */
  //static Lsodar* convert (OneStepIntegrator* osi);

};

#endif // Lsodar_H
