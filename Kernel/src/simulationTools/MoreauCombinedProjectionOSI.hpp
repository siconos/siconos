/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
/*! \file
  Moreau Time-Integrator for Dynamical Systems for Combined Projection Algorithm
*/

#ifndef MOREAUCOMBINEDPROJECTIONOSI_H
#define MOREAUCOMBINEDPROJECTIONOSI_H

#include "OneStepIntegrator.hpp"
#include "Moreau.hpp"
#include "SimpleMatrix.hpp"

class Simulation;
class SiconosMatrix;

const unsigned int MOREAUCOMBINEDPROJECTIONOSISTEPSINMEMORY = 1;

/**  Moreau Time-Integrator for Dynamical Systems for Combined Projection Algorithm
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.4.0.
 *  \date (Creation) May 02, 2012
 *
 * See User's guide, \ref docSimuMoreauTS for details.
 *
 * Moreau class is used to define some time-integrators methods for a
 * list of dynamical systems.
 *
 * This class reimplement a special activation of constraints
 * for the Combined Projection Algorithm
 *
 *
 */
class MoreauCombinedProjectionOSI : public Moreau
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(MoreauCombinedProjectionOSI);

  /** Default constructor
   */
  MoreauCombinedProjectionOSI() {};




public:

  /** constructor from theta value only
   *  \param theta value for all these DS.
   */
  explicit MoreauCombinedProjectionOSI(double theta) : Moreau(theta) {}  ;

  /** constructor from a minimum set of data: one DS and its theta
   *  \param ds SP::DynamicalSystem : the DynamicalSystem linked to the OneStepIntegrator
   *  \param theta value of the parameter
   */
  MoreauCombinedProjectionOSI(SP::DynamicalSystem ds, double theta) : Moreau(ds, theta) {}  ;

  /** constructor from a minimum set of data
   *  \param DynamicalSystemsSet : the list of DynamicalSystems to be integrated
   *  \param theta value for all these DS.
   */
  MoreauCombinedProjectionOSI(DynamicalSystemsSet& dsSet, double theta) :     Moreau(dsSet, theta) {}   ;

  /** destructor
   */
  virtual ~MoreauCombinedProjectionOSI() {};

  // --- OTHER FUNCTIONS ---

  /** initialization of the integrator; for linear time
      invariant systems, we compute time invariant operator (example :
      W)
  */
  void initialize();





  /** Apply the rule to one Interaction to known if is it should be included
   * in the IndexSet of level i
   */
  virtual bool addInteractionInIndexSet(SP::Interaction inter, unsigned int i);

  /** Apply the rule to one Interaction to known if is it should be removed
   * in the IndexSet of level i
   */
  virtual bool removeInteractionInIndexSet(SP::Interaction inter, unsigned int i);

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, 0 otherwise
   */
  static MoreauCombinedProjectionOSI* convert(OneStepIntegrator* osi);

  /** visitors hook
  */
  ACCEPT_STD_VISITORS();

};

#endif // MOREAUCOMBINEDPROJECTIONOSI_H
