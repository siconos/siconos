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
class Lsodar : public OneStepIntegrator
{
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
   Integrator
   */
  Lsodar(TimeDiscretisation*, DynamicalSystem*);

  ~Lsodar();

  /** \fn void computeFreeState()
  *   \brief compute the free state of the dynamical system
  */
  void computeFreeState();

  /** \fn void integrate()
  *   \brief integrates the dynamical system
  */
  void integrate();

  /** \fn Lsodar* convert (OneStepIntegrator* osi)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, NULL otherwise
   */
  static Lsodar* convert(OneStepIntegrator* osi);

private:
  /** \fn Lsodar()
   *  \brief default constructor
   */
  Lsodar();
};


typedef void (* fctPtr)(int *sizeOfX, double *time, double *x, double *xdot);
extern "C" void tryfunction(fctPtr);

#endif // Lsodar_H
