/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
 *  \version 1.2.0.
 *  \date (Creation) Apr 26, 2004
 *
 *
 * \todo Add the LU Factorization of W in the initialization,  and adapt the resolution in the iteration
 */

/** map of SiconosMatrix; key = the related DS*/
typedef std::map<DynamicalSystem*, SiconosMatrix*> mapOfMatrices;

/** map of SiconosMatrix; key = the related DS*/
typedef std::map<DynamicalSystem*, bool> mapOfBool;

/** map of double; key = the related DS */
typedef std::map<DynamicalSystem*, double> mapOfDouble;

/** iterator through a map of matrices */
typedef mapOfMatrices::iterator matIterator;

/** iterator through a map of double */
typedef mapOfDouble::iterator doubleIterator;


class Moreau : public OneStepIntegrator
{
private:

  /** Stl map that associates a W Moreau matrix to each DynamicalSystem of the OSI */
  mapOfMatrices WMap;

  /** Stl map that associates a theta parameter for the integration scheme to each DynamicalSystem of the OSI */
  mapOfDouble thetaMap;

  /** Stl map that associates a bool to each DynamicalSystem of the OSI - If true, then W has been allocated inside the class.*/
  mapOfBool isWAllocatedInMap;

  /** \fn Moreau(Strategy* = NULL)
   *  \brief Default constructor
   *  \param Strategy * : the strategy that owns the osi, default = NULL
   */
  Moreau(Strategy*);

public:

  /** \fn Moreau(OneStepIntegratorXML*,Strategy* = NULL);
   *  \brief constructor from xml file
   *  \param OneStepIntegratorXML* : the XML object corresponding
   *  \param Strategy * : the strategy that owns the osi, default = NULL
   */
  Moreau(OneStepIntegratorXML*, Strategy* = NULL);

  /** \fn Moreau(DynamicalSystem* , const double&, Strategy* = NULL)
   *  \brief constructor from a minimum set of data: one DS and its theta
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   *  \param Theta value
   *  \param Strategy * : the strategy that owns the osi, default = NULL
   */
  Moreau(DynamicalSystem*, const double&, Strategy* = NULL);

  /** \fn ~Moreau()
   *  \brief destructor
   */
  ~Moreau();

  // --- GETTERS/SETTERS ---
  // -- W --

  /** \fn mapOfMatrices getWMap()
    *  \brief get W map
    *  \return a mapOfMatrices
    */
  inline mapOfMatrices getWMap() const
  {
    return WMap;
  };

  /** \fn mapOfMatrices getIsWAllocatedInMap()
    *  \brief get isWAllocatedIn map
    *  \return a mapOfBool
    */
  inline mapOfBool getIsWAllocatedInMap() const
  {
    return isWAllocatedInMap;
  };

  /** \fn void setWMap(const mapOfMatrices& newMap)
   *  \brief set W map to newMap
   *  \param a mapOfMatrices
   */
  void setWMap(const mapOfMatrices&);

  /** \fn const SimpleMatrix getW(DynamicalSystem* ds = NULL) const
   *  \brief get the value of W corresponding to DynamicalSystem ds
   * \param a pointer to DynamicalSystem, optional, default = NULL. get W[0] in that case
   *  \return SimpleMatrix
   */
  const SimpleMatrix getW(DynamicalSystem* = NULL);

  /**\fn SiconosMatrix* getWPtr(DynamicalSystem * ds) const
   * \brief get W corresponding to DynamicalSystem ds
   * \param a pointer to DynamicalSystem, optional, default = NULL. get W[0] in that case
   * \return pointer to a SiconosMatrix
   */
  SiconosMatrix* getWPtr(DynamicalSystem* ds);

  /**\fn void setW (const SiconosMatrix& newValue, DynamicalSystem* ds)
   * \brief set the value of W[ds] to newValue
   * \param SiconosMatrix newValue
   * \param a pointer to DynamicalSystem,
   */
  void setW(const SiconosMatrix&, DynamicalSystem*);

  /**\fn void setWPtr(SiconosMatrix* newPtr, DynamicalSystem* ds)
   * \brief set W[ds] to pointer newPtr
   * \param SiconosMatrix * newPtr
   * \param a pointer to DynamicalSystem
   */
  void setWPtr(SiconosMatrix *newPtr, DynamicalSystem*);

  // -- theta --

  /** \fn mapOfDouble getThetaMap()
    *  \brief get theta map
    *  \return a mapOfDouble
    */
  inline mapOfDouble getThetaMap() const
  {
    return thetaMap;
  };

  /** \fn void setThetaMap(const mapOfDouble&)
   *  \brief set theta map
   *  \param a mapOfDouble
   */
  void setThetaMap(const mapOfDouble&);

  /** \fn const double getTheta(DynamicalSystem* ds)
   *  \brief get thetaMap[ds]
   *  \param a DynamicalSystem
   *  \return a double
   */
  const double getTheta(DynamicalSystem*);

  /** \fn double setTheta(const double&, DynamicalSystem* ds)
   *  \brief set the value of thetaMap[ds]
   *  \param a double
   *  \param a DynamicalSystem
   */
  void setTheta(const double&, DynamicalSystem*);

  // --- OTHER FUNCTIONS ---

  /** \fn void initialize()
   *  \brief initialization of the Moreau integrator; for linear time invariant systems, we compute time invariant operator (example : W)
   *  \todo LU factorization of time invariant operator (example : W)
   */
  void initialize();

  /** \fn void computeW(const double& t, DynamicalSystem* ds)
   *  \brief compute WMap[ds] Moreau matrix at time t
   *  \param the time (double)
   *  \param a pointer to DynamicalSystem
   */
  void computeW(const double&, DynamicalSystem*);

  /** \fn void computeFreeState()
   *  \brief integrates the Dynamical System linked to this integrator without boring the constraints
   */
  void computeFreeState();

  /** \fn void integrate(const double&, const double&, double&, bool&)
   *  \brief integrate the system, between tinit and tend (->iout=true), with possible stop at tout (->iout=false)
   *  \param double: tinit, initial time
   *  \param double: tend, end time
   *  \param double: tout, real end time
   *  \param bool: true if tend is reached, else false.
   */
  void integrate(const double&, const double&, double&, bool&);

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
  void display();

  /** \fn Moreau* convert (OneStepIntegrator* osi)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, 0 otherwise
   */
  static Moreau* convert(OneStepIntegrator* osi);

};

#endif // MOREAU_H
