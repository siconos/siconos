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
#ifndef UNITARYRELATION_H
#define UNITARYRELATION_H

#include "Interaction.h"

// tools
#include "SimpleVector.h"

// const
#include "SiconosConst.h"

/** \class UnitaryRelation
 *  \brief this class is an interface to single relations from Interactions
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.2.0.
 *  \date (Creation) June 06, 2006
 *
 * Remind that each Interaction is composed with one (pointer to) relation and a non smooth law (size nsLawSize). The number of relations
 * is sizeOfInteraction/nsLawSize. A UnitaryRelation is one of the relations of the Interaction. Actually, this class only provides an interface
 * to handle single relations, this for IndexSets used in Topology.
 * Each UnitaryRelation has a pointer to its "mother" Interaction  and methods to compute y, lambda and so on.
 *
 */

class UnitaryRelation
{

  // === PRIVATE MEMBERS ===

private:

  /** link to Interaction that owns this relation **/
  Interaction * mainInteraction;

  /** relative position of the present relation in the Interaction - For example if the present relation takes place from index 2 to 4 in y vector
   of mainInteraction, thus the relative position is equal to 2. */
  unsigned int relativePosition;

public:

  /** \fn UnitaryRelation()
   *  \brief default constructor
   */
  UnitaryRelation();

  // === PUBLIC FUNCTIONS ===

  // === CONSTRUCTORS/DESTRUCTOR ===

  /** \fn UnitaryRelation(const UnitaryRelation&)
   *  \brief copy constructor
   *  \param UnitaryRelation* : the object to copy
   */
  UnitaryRelation(const UnitaryRelation& inter);

  /** \fn ~UnitaryRelation()
   * \brief destructor
   */
  ~UnitaryRelation();

  /** \fn SimpleVector* getYPtr(const unsigned int& i) const
   *  \brief get y[i], derivative number i of output for the present relation
   *  \return pointer on a SimpleVector
   */
  SimpleVector* getYPtr(const unsigned int& i) const ;

  /** \fn SimpleVector* getLambdaPtr(const unsigned int& i) const
   *  \brief get lambda[i], derivative number i of input for the present relation
   *  \return pointer on a SimpleVector
   */
  SimpleVector* getLambdaPtr(const unsigned int& i) const ;


};

#endif // UNITARYRELATION_H
