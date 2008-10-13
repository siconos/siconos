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

#ifndef MOREAU_H
#define MOREAU_H

#include "OneStepIntegrator.h"
#include "SimpleMatrix.h"

class Simulation;
class SiconosMatrix;

const unsigned int MOREAUSTEPSINMEMORY = 1;

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
 * - computeFreeState(): computes xfree (or vfree), dynamical systems state without taking non-smooth part into account \n
 * - updateState(): computes x (q,v), the complete dynamical systems states.
 *
 */
class Moreau : public OneStepIntegrator
{
protected:

  /** Stl map that associates a W Moreau matrix to each DynamicalSystem of the OSI */
  MapOfDSMatrices WMap;

  /** Stl map that associates a theta parameter for the integration scheme to each DynamicalSystem of the OSI */
  MapOfDouble thetaMap;

private:
  /** Default constructor
   */
  Moreau();

public:

  /** constructor from xml file
   *  \param OneStepIntegratorXML* : the XML object corresponding
   *  \param DynamicalSystemsSet: set of all DS in the NSDS
   */
  Moreau(SP::OneStepIntegratorXML, SP::DynamicalSystemsSet);

  /** constructor from a minimum set of data: one DS and its theta
   *  \param SP::DynamicalSystem : the DynamicalSystem linked to the OneStepIntegrator
   *  \param Theta value
   */
  Moreau(SP::DynamicalSystem, double);

  /** constructor from a minimum set of data
   *  \param DynamicalSystemsSet : the list of DynamicalSystems to be integrated
   *  \param theta value for all these DS.
   */
  Moreau(DynamicalSystemsSet&, double);

  /** constructor from a minimum set of data
   *  \param DynamicalSystemsSet : the list of DynamicalSystems to be integrated
   *  \param Map of theta values for the DS.
   */
  Moreau(DynamicalSystemsSet&, const MapOfDouble&);

  /** destructor
   */
  virtual ~Moreau() {};

  // --- GETTERS/SETTERS ---
  // -- W --

  /** get W map
   *  \return a MapOfDSMatrices
   */
  inline MapOfDSMatrices getWMap() const
  {
    return WMap;
  };

  /** set W map to newMap
   *  \param a MapOfDSMatrices
   */
  void setWMap(const MapOfDSMatrices&);

  /** get the value of W corresponding to DynamicalSystem ds
   * \param a pointer to DynamicalSystem, optional, default = NULL. get W[0] in that case
   *  \return SimpleMatrix
   */
  const SimpleMatrix getW(SP::DynamicalSystem = SP::DynamicalSystem());

  /** get W corresponding to DynamicalSystem ds
   * \param a pointer to DynamicalSystem, optional, default = NULL. get W[0] in that case
   * \return pointer to a SiconosMatrix
   */
  SP::SiconosMatrix getWPtr(SP::DynamicalSystem ds);

  /** set the value of W[ds] to newValue
   * \param SiconosMatrix newValue
   * \param a pointer to DynamicalSystem,
   */
  void setW(const SiconosMatrix&, SP::DynamicalSystem);

  /** set W[ds] to pointer newPtr
   * \param SP::SiconosMatrix  newPtr
   * \param a pointer to DynamicalSystem
   */
  void setWPtr(SP::SiconosMatrix newPtr, SP::DynamicalSystem);

  // -- theta --

  /** get theta map
   *  \return a MapOfDouble
   */
  inline MapOfDouble getThetaMap() const
  {
    return thetaMap;
  };

  /** set theta map
   *  \param a MapOfDouble
   */
  void setThetaMap(const MapOfDouble&);

  /** get thetaMap[ds]
   *  \param a DynamicalSystem
   *  \return a double
   */
  const double getTheta(SP::DynamicalSystem);

  /** set the value of thetaMap[ds]
   *  \param a double
   *  \param a DynamicalSystem
   */
  void setTheta(double, SP::DynamicalSystem);

  // --- OTHER FUNCTIONS ---

  /** initialization of the Moreau integrator; for linear time invariant systems, we compute time invariant operator (example : W)
      \param the simulation, owner of this OSI
   */
  void initialize(SP::Simulation);

  /** init WMap[ds] Moreau matrix at time t
   *  \param the time (double)
   *  \param a pointer to DynamicalSystem
   */
  void initW(double, SP::DynamicalSystem);

  /** compute WMap[ds] Moreau matrix at time t
   *  \param the time (double)
   *  \param a pointer to DynamicalSystem
   */
  void computeW(double, SP::DynamicalSystem);

  /** return the maximum of all norms for the "Moreau-discretized" residus of DS
      \return a double
   */
  double computeResidu();

  /** integrates the Dynamical System linked to this integrator without boring the constraints
   */
  virtual void computeFreeState();

  /** integrate the system, between tinit and tend (->iout=true), with possible stop at tout (->iout=false)
   *  \param double: tinit, initial time
   *  \param double: tend, end time
   *  \param double: tout, real end time
   *  \param int: useless flag (for Moreau, used in Lsodar)
   */
  void integrate(double&, double&, double&, int&);

  /** updates the state of the Dynamical Systems
   *  \param unsigned int: level of interest for the dynamics: not used at the time
   */
  virtual void updateState(unsigned int);

  /** copy the matrix W of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveWToXML();

  /** Displays the data of the Moreau's integrator
   */
  void display();

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, 0 otherwise
   */
  static Moreau* convert(OneStepIntegrator* osi);

};

#endif // MOREAU_H
