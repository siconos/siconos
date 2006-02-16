/* Siconos-Kernel version 1.1.1, Copyright INRIA 2005-2006.
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
#ifndef MOREAU_H
#define MOREAU_H

#include "OneStepIntegrator.h"
#include "MoreauXML.h"

const unsigned int MOREAUSTEPSINMEMORY = 1;

/** \class Moreau
 *  \brief It's a kind of single-step Integrator
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.1.
 *  \date (Creation) Apr 26, 2004
 *
 *
 * \todo Add the LU Factorization of W in the initialization,  and adapt the resolution in the iteration
 */
class Moreau : public OneStepIntegrator
{
private:

  /** \fn Moreau()
   *  \brief Default constructor
   */
  Moreau();

  /** a specific matrix of the Moreau Integrator */
  SiconosMatrix *W;

  /** boolean indicator, to check whether W has been allocated in the present class or not */
  bool isWAllocatedIn;

  /** parameter of the theta-method */
  double theta;

public:

  /** \fn Moreau(OneStepIntegratorXML*,TimeDiscretisation*, DynamicalSystem* )
   *  \brief constructor from xml file
   *  \param OneStepIntegratorXML* : the XML object corresponding
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   */
  Moreau(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem*);

  /** \fn Moreau(TimeDiscretisation*, DynamicalSystem* , const double& theta)
   *  \brief constructor from a minimum set of data
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   *  \param Theta value
   */
  Moreau(TimeDiscretisation*, DynamicalSystem*, const double& theta);

  /** \fn ~Moreau()
   *  \brief destructor
   */
  ~Moreau();

  // --- GETTERS/SETTERS ---
  // -- W --

  /** \fn  const SiconosMatrix getW(void) const
   *  \brief get the value of W
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getW() const
  {
    return *W;
  }

  /** \fn SiconosMatrix* getWPtr(void) const
   *  \brief get W
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getWPtr() const
  {
    return W;
  }

  /** \fn void setW (const SiconosMatrix& newValue)
   *  \brief set the value of W to newValue
   *  \param SiconosMatrix newValue
   */
  void setW(const SiconosMatrix& newValue);

  /** \fn void setWPtr(SiconosMatrix* newPtr)
   *  \brief set W to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setWPtr(SiconosMatrix *newPtr);

  // -- theta --

  /** \fn const double getTheta() const
   *  \brief allows to get the double value theta of the Moreau's Integrator
   *  \return double value theta
   */
  inline const double getTheta() const
  {
    return theta;
  }

  /** \fn double setTheta(const double&)
   *  \brief set the value of theta
   *  \param ref on a double
   */
  inline void setTheta(const double& newTheta)
  {
    theta = newTheta;
  }

  // --- OTHER FUNCTIONS ---

  /** \fn void initialize()
   *  \brief initialization of the Moreau integrator; for linear time invariant systems, we compute time invariant operator (example : W)
   *  \todo LU factorization of time invariant operator (example : W)
   */
  void initialize();

  /** \fn void computeW(const double& t)
   *  \brief compute W Moreau matrix at time t
   *  \param the time (double)
   */
  void computeW(const double&);

  /** \fn void computeFreeState()
   *  \brief integrates the Dynamical System linked to this integrator without boring the constraints
   */
  void computeFreeState();

  /** \fn void integrate()
   *  \brief makes computations to integrate the data of a Dynamical System with the Moreau Integrator
   */
  void integrate();

  /** \fn void updateState()
   *  \brief updates the state of the Dynamical System
   */
  void updateState();

  /** \fn void saveIntegratorToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveIntegratorToXML();

  /** \fn void saveWToXML()
   *  \brief copy the matrix W of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveWToXML();

  /** \fn display()
   *  \brief Displays the data of the Moreau's integrator
   */
  void display() const;

  /** \fn Moreau* convert (OneStepIntegrator* osi)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, 0 otherwise
   */
  static Moreau* convert(OneStepIntegrator* osi);

};

#endif // MOREAU_H
