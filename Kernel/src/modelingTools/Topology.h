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
#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "NonSmoothDynamicalSystem.h"
#include "Interaction.h"

// const
#include "SiconosConst.h"

#include "InteractionsSet.h"
#include "UnitaryRelationsSet.h"
#include <vector>
#include <string>
#include <map>
#include <set>

class NonSmoothDynamicalSystem;
class InteractionsSet;
class Interaction;
class SiconosMatrix;
class UnitaryRelationsSet;
class UnitaryRelation;

/** \class Topology
 *  \brief this class provides maps to describe the topology of interactions of a NonSmoothDynamicalSystem
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.2.0.
 *  \date (Creation) July 20, 2005
 *
 *
 * Rule: no Unitary Relation can be created outside of the present class (ie the only calls to new UnitaryRelation(...) are part of Topology methods)
 * Thus, no flags for inside-class allocation - All UnitaryRelation pointers are cleared during Topology destructor call.
 *
 */

/** map that links each unitary relation with its relative degree */
typedef std::map< UnitaryRelation*, unsigned int > UnitaryRelationsIntMap;

/** and the corresponding iterators */
typedef UnitaryRelationsIntMap::iterator IteratorForRelativeDegrees;
typedef UnitaryRelationsIntMap::const_iterator ConstIteratorForRelativeDegrees;
class Topology
{

private:

  // --- MEMBERS ---

  /** the set of all the interactions of the system */
  InteractionsSet allInteractions;

  /** index set I0, ie a set of all the Unitary Relations - This corresponds to indexSets[0] of the Simulation */
  UnitaryRelationsSet indexSet0;

  /** map that links UnitaryRelations with their relative degrees */
  UnitaryRelationsIntMap relativeDegrees;

  /** check if topology has been updated since nsds modifications occur */
  bool isTopologyUpToDate;

  /** check if topology is static (all relative degrees = 0 or 1) or not */
  bool isTopologyTimeInvariant;

  /** the NonSmoothDynamicalSystem that owns this topology */
  NonSmoothDynamicalSystem * nsds;

  // === PRIVATE FUNCTIONS ===

  /** \fn const bool addInteractionInIndexSet(Interaction* inter)
   *  \brief schedules the relations of Interaction inter in IndexSet0 (ie creates UnitaryRelations)
   * \param: a pointer to Interaction
   */
  const bool addInteractionInIndexSet(Interaction*);

  /** \fn void computeRelativeDegrees()
   *   \brief compute the  RelativeDegrees Map
   */
  void computeRelativeDegrees();

  /** \fn Topology()
   *  \brief default constructor
   */
  Topology();

public:

  // --- CONSTRUCTORS/DESTRUCTOR ---

  /** \fn Topology(NonSmoothDynamicalSystem*)
   *  \brief constructor from nsds that owns that topology
   * \param: a NonSmoothDynamicalSystem*
   */
  Topology(NonSmoothDynamicalSystem*) ;

  /** \fn ~Topology()
   *  \brief destructor */
  ~Topology();

  // === GETTERS/SETTERS ===

  /** \fn const InteractionsSet getInteractions()
   *  \brief get all the Interactions of the Topology problem (saved in a set)
   *  \return an InteractionsSet
   */
  inline const InteractionsSet getInteractions() const
  {
    return allInteractions;
  }

  /** \fn void setInteractions(const InteractionsSet&)
   *  \brief to set allInteractions
   *  \param an InteractionsSet
   */
  void setInteractions(const InteractionsSet&) ;

  /** \fn const bool hasInteraction(Interaction * inter)
   *  \brief check if Interaction inter is in the set
   *  \param a pointer to Interaction
   *  \return a bool
   */
  const bool hasInteraction(Interaction*) const;

  /** \fn  const UnitaryRelationsSet getIndexSet0()
   *  \brief get the index set of all Unitary Relations
   *  \return a UnitaryRelationsSet
   */
  inline const UnitaryRelationsSet getIndexSet0() const
  {
    return indexSet0;
  }

  // --- relativeDegreesMap ---

  /** \fn  UnitaryRelationsIntMap getRelativeDegrees()
   *  \brief get the relativeDegrees Map of this topology
   *  \return a UnitaryRelationsIntMap
   */
  inline const UnitaryRelationsIntMap getRelativeDegrees() const
  {
    return relativeDegrees;
  }

  /** \fn const unsigned int getRelativeDegree(UnitaryRelation*) const
   *  \brief get the relativeDegree vector of a specific UnitaryRelation
   *  \param a pointer to UnitaryRelation
   *  \return an unsigned int
   */
  inline const unsigned int getRelativeDegree(UnitaryRelation* UR)
  {
    return relativeDegrees[UR];
  }

  /** \fn const unsigned int getMaxRelativeDegree() const
   *  \brief for all relative degrees (one per Unitary Relation), find the maximum value.
   *  \return an unsigned int
   */
  const unsigned int getMaxRelativeDegree();

  /** \fn const unsigned int getMinRelativeDegree() const
   *  \brief for all relative degrees (one per Unitary Relation), find the minimum value.
   *  \return an unsigned int
   */
  const unsigned int getMinRelativeDegree();

  // --- isTopologyUpToDate ---

  /** \fn void  setUpToDate(const bool val)
   *  \brief set isTopologyUpToDate to val
   *  \param a bool
   */
  inline void setUpToDate(const bool val)
  {
    isTopologyUpToDate = val;
  }

  /** \fn bool isUpToDate()
   *  \brief check if topology has been updated since modifications occurs on nsds
   *  \return a bool
   */
  inline bool isUpToDate()
  {
    return isTopologyUpToDate;
  }

  // --- isTopologyTimeInvariant ---

  /** \fn void  setTimeInvariant(const bool val)
   *  \brief set isTopologyTimeInvariant to val
   *  \param a bool
   */
  inline void setTimeInvariant(const bool val)
  {
    isTopologyTimeInvariant = val;
  }

  /** \fn bool isTimeInvariant()
   *  \brief check if all relative degrees are equal to 0 or 1
   *  \return a bool
   */
  inline bool isTimeInvariant()
  {
    return isTopologyTimeInvariant;
  }

  /** \fn void initialize();
   *   \brief initializes the topology
   */
  void initialize();
};

#endif // TOPOLOGY_H
