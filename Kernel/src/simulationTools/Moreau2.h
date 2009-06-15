/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
/*! \file
  Moreau Time-Integrator for Dynamical Systems
*/

#ifndef MOREAU2_H
#define MOREAU2_H

#include "Moreau.h"
#include "SimpleMatrix.h"
#include "SiconosVector.h"
#include "FirstOrderLinearDS.h"


/**  Moreau Time-Integrator for Dynamical Systems
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 * See User's guide, \ref docSimuMoreauTS for details.
 *
 * Moreau class is used to define some time-integrators methods for a list of dynamical systems.
 * Each DynamicalSystem is associated to a SiconosMatrix, named "W", and a double, "theta", through two
 * STL maps:
 * - WMap, with WMap[ds] = a pointer to a SiconosMatrix
 * - thetaMap, thetaMap[ds] = a double
 * ds being a SP::DynamicalSystem
 *
 * W matrices are initialized and computed in initW and computeW. Depending on the DS type, they
 * may depend on time and DS state (x).
 *
 * Main functions:
 *
 * - computeFreeState(): computes Ffree (or vfree), dynamical systems state without taking non-smooth part into account \n
 * - updateState():
 *
 */
class Moreau2 : public Moreau
{
private:

  /** Default constructor
   */
  Moreau2() {};

public:

  /** constructor from a minimum set of data: one DS and its theta
   *  \param SP::DynamicalSystem : the DynamicalSystem linked to the OneStepIntegrator
   *  \param Theta value
   */
  Moreau2(SP::DynamicalSystem, double);

  /** constructor from a minimum set of data: one DS and its theta
   *  \param DynamicalSystemSet* : the DynamicalSystem set linked to the OneStepIntegrator
   *  \param Theta value
   */
  Moreau2(DynamicalSystemsSet&, double);

  ~Moreau2();

  //  SP::SiconosVector  getFfree(FirstOrderLinearDS *d);

  /** integrates the Dynamical System linked to this integrator without boring the constraints
   */
  void computeFreeState();


  /** updates the state of the Dynamical Systems
   *  \param unsigned int: level of interest for the dynamics: not used at the time
   */
  void updateState(unsigned int);
  /*  virtual SP::SiconosVector  getWorkX(SP::DynamicalSystem d);*/

  static Moreau2* convert(OneStepIntegrator* osi);


};

#endif // MOREAU2_H
