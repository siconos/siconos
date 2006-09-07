/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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
#ifndef UnitaryRelationsSET_H
#define UnitaryRelationsSET_H

#include "RuntimeCmp.h"
#include "check.h"
#include "UnitaryRelation.h"
#include <iostream>
#include <set>
#include <map>

class UnitaryRelation;

/** \class UnitaryRelationsSet
 *  \brief set (stl) of pointers to Unitary Relation - Used in Index sets in Topology.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
 *  \date (Creation) June 8, 2006
 *
 * Unitary Relations pointers are sorted in according to their address
 * Only one occurence of a Unitary Relation can be present in the set.
 * Possible operations are insert, erase, get or check presence of an UR.
 *
 * Warning: any call to operator = or copy constructor results in a false copy: (if set1=set2, pointers to ds of set1 are equal to those of set2)
 *
 */

// Structure used for sorting in Unitary Relation set. The address is used to compare two UR.
struct compareUR
{
  bool operator()(const UnitaryRelation* r1, const UnitaryRelation* r2) const
  {
    return (r1 < r2);
  }
};

/** set of Unitary Relations */
typedef std::set<UnitaryRelation*, compareUR > URSet;

/** iterator through a set of Unitary Relations */
typedef URSet::iterator UnitaryRelationIterator;

/** const iterator through a set of Unitary Relations */
typedef URSet::const_iterator ConstUnitaryRelationIterator;

/** return type value for insert function - bool = false if insertion failed. */
typedef std::pair<URSet::iterator, bool> CheckInsertUnitaryRelation;

class UnitaryRelationsSet
{
protected:

  /** a set of UnitaryRelation, sorted thanks to their id number */
  URSet setOfUnitaryRelations;

  /** a map of bool to check inside-class allocation.
   *  isURAllocatedIn[ur] = true if ur has been allocated in a method of the present class.
   */
  std::map<UnitaryRelation*, bool > isURAllocatedIn;

public:

  /** \fn UnitaryRelationsSet()
   *  \brief default constructor
   */
  UnitaryRelationsSet();

  /** \fn UnitaryRelationsSet(const UnitaryRelationsSet&)
   *  \brief copy constructor
   *  \param a UnitaryRelationsSet to be copied
   */
  UnitaryRelationsSet(const UnitaryRelationsSet&);

  /** \fn ~UnitaryRelationsSet()
   *  \brief destructor
   */
  ~UnitaryRelationsSet();

  /** \fn   UnitaryRelationsSet& operator=( const UnitaryRelationsSet& );
   *  \brief assignment
   */
  UnitaryRelationsSet& operator=(const UnitaryRelationsSet&);

  /** \fn const unsigned int size() const
   *  \brief return the number of UnitaryRelations in the set
   *  \return an unsigned int
   */
  inline const unsigned int size() const
  {
    return setOfUnitaryRelations.size();
  };

  /** \fn const bool isEmpty() const
   *  \brief return true if the set is empty, else false
   *  \return a bool
   */
  inline const bool isEmpty() const
  {
    return setOfUnitaryRelations.empty();
  };

  /** \fn UnitaryRelationIterator begin() const
   *  \brief return iterator on the first element of setOfUnitaryRelations
   *  \return a UnitaryRelationIterator
   */
  inline UnitaryRelationIterator begin() const
  {
    return setOfUnitaryRelations.begin();
  };

  /** \fn UnitaryRelationIterator end() const
   *  \brief return iterator on setOfUnitaryRelations.end()
   *  \return a UnitaryRelationIterator
   */
  inline UnitaryRelationIterator end() const
  {
    return setOfUnitaryRelations.end();
  }

  /** \fn URSet getSetOfUnitaryRelations()
   *  \brief return setOfUnitaryRelations
   *  \return a UnitaryRelationSet
   */
  inline const URSet getSetOfUnitaryRelations() const
  {
    return setOfUnitaryRelations;
  }

  /** \fn bool isUnitaryRelationIn(UnitaryRelation* ur)
   *  \brief return true if ur is in the set
   *  \param a pointer to UnitaryRelation
   *  \return a bool
   */
  const bool isUnitaryRelationIn(UnitaryRelation*) const;

  /** \fn UnitaryRelationIterator find(UnitaryRelation* ur)
   *  \brief same as find function of stl set
   *  \param a pointer to UnitaryRelation
   *  \param a UnitaryRelationIterator
   */
  UnitaryRelationIterator find(UnitaryRelation*);

  /** \fn CheckInsertUnitaryRelation insert(UnitaryRelation* ur)
   *  \brief insert Unitary Relation ur into the set
   *  \param a pointer to UnitaryRelation
   *  \return a CheckInsertUnitaryRelation ( boolean type information)
   */
  CheckInsertUnitaryRelation insert(UnitaryRelation*);

  /** \fn void erase(UnitaryRelation* ur)
   *  \brief remove Unitary Relation ur from the set
   *  \param a pointer to UnitaryRelation
   */
  void erase(UnitaryRelation*);

  /** \fn clear() {setOfUnitaryRelations.clear();};
   *  \brief remove all Unitary Relations from the set
   */
  void clear();

  /** \fn void display() const
   *  \brief screen-display of the numbers of the Unitary Relations present in the set.
   */
  void display() const;

  /** \fn const UnitaryRelationsSet intersection(const UnitaryRelationsSet& s1, const UnitaryRelationsSet& s2) const
   *  \brief return the intersection of s1 and s2 (-> set_intersection stl function)
   */
  friend const UnitaryRelationsSet intersection(const UnitaryRelationsSet& s1, const UnitaryRelationsSet& s2);

  /** \fn const UnitaryRelationsSet operator-(const UnitaryRelationsSet& s1, const UnitaryRelationsSet& s2) const
   *  \brief return the difference betwee s1 and s2 (-> set_difference stl function)
   */
  friend const UnitaryRelationsSet operator-(const UnitaryRelationsSet& s1, const UnitaryRelationsSet& s2);
};

#endif
