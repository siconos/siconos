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
  MoreauJeanOSI Time-Integrator for Dynamical Systems for Combined Projection Algorithm
*/

#ifndef MOREAUCOMBINEDPROJECTIONOSI_H
#define MOREAUCOMBINEDPROJECTIONOSI_H

#include "OneStepIntegrator.hpp"
#include "MoreauJeanOSI.hpp"
#include "SimpleMatrix.hpp"


const unsigned int MOREAUCOMBINEDPROJECTIONOSISTEPSINMEMORY = 1;

/**  \class MoreauJeanCombinedProjectionOSI 
 *   \brief One Step time Integrator for First Order Dynamical Systems  for
 *    mechanical Systems (LagrangianDS and NewtonEulerDS) with  Combined Projection Algorithm
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.4.0.
 *  \date (Creation) May 02, 2012
 *
 * This class reimplement a special activation of constraints
 * in the MoreauJeanOSI for the Combined Projection Algorithm
 *
 * References :
 *
 * V. Acary. Projected event-capturing time-stepping schemes for nonsmooth mechanical systems with unilateral contact 
 * and coulomb’s friction. 
 * Computer Methods in Applied Mechanics and Engineering, 256:224 – 250, 2013. ISSN 0045-7825. 
 * URL http://www.sciencedirect.com/science/article/pii/S0045782512003829.
 *
 */
class MoreauJeanCombinedProjectionOSI : public MoreauJeanOSI
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(MoreauJeanCombinedProjectionOSI);

  /** Default constructor
   */
  MoreauJeanCombinedProjectionOSI() {};




public:

  /** constructor from theta value only
   *  \param theta value for all these DS.
   */
  explicit MoreauJeanCombinedProjectionOSI(double theta) : MoreauJeanOSI(theta) {}  ;

  /** destructor
   */
  virtual ~MoreauJeanCombinedProjectionOSI() {};

  // --- OTHER FUNCTIONS ---

  /** initialization of the integrator; for linear time
      invariant systems, we compute time invariant operator (example :
      W)
  */
  void initialize(Model& m);





  /** Apply the rule to one Interaction to known if is it should be included
   * in the IndexSet of level i
   * \param inter concerned interaction 
   * \param i level
   * \return bool
   */
  virtual bool addInteractionInIndexSet(SP::Interaction inter, unsigned int i);

  /** Apply the rule to one Interaction to known if is it should be removed
   * in the IndexSet of level i
   * \param inter concerned interaction 
   * \param i level
   * \return bool
   */
  virtual bool removeInteractionInIndexSet(SP::Interaction inter, unsigned int i);

  /** visitors hook
  */
  ACCEPT_STD_VISITORS();

};

#endif // MOREAUCOMBINEDPROJECTIONOSI_H
