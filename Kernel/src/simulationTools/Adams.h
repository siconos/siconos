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
#ifndef ADAMS_H
#define ADAMS_H

#include "OneStepIntegrator.h"
#include "AdamsXML.h"

/** \class Adams
 *  \brief Adams is a kind of multi-step integrator.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 *
 *
 */
class Adams : public OneStepIntegrator
{
private:
  /** \fn Adams()
   *  \brief default constructor
   */
  Adams();

  /**  */
  int r;

public:

  /** \fn Adams(OneStepIntegratorXML*)
   *  \brief constructor from xml file
   *  \param OneStepIntegratorXML* : the XML object corresponding
   */
  Adams(OneStepIntegratorXML*);
  /** \fn Adams(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem* )
   *  \brief constructor from a minimum set of data
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   */
  Adams(TimeDiscretisation*, DynamicalSystem*);

  ~Adams();

  /** \fn double const getR() const
   *   \brief Return the r of the OneStepIntegrator
   *   \return int : the value of r
   */
  inline const int getR() const
  {
    return r;
  }

  /** \fn void setR(const int& r)
   *   \brief Return the r of OneStepIntegrator
   *   \param double : the value to set r
   */
  inline void setR(const int& newR)
  {
    r = newR;
  }

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

  /** \fn Adams* convert (OneStepIntegrator* osi)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, NULL otherwise
   */
  static Adams* convert(OneStepIntegrator* osi);

};

#endif // ADAMS_H
