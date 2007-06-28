/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
/*! \file DynamicalSystemsSet.h

*/
#ifndef DSSET_H
#define DSSET_H

#include "RuntimeCmp.h"
#include <set>
#include <map>

class DynamicalSystem;

/** set of Dynamical Systems */
typedef std::set<DynamicalSystem*, RuntimeCmp<DynamicalSystem> > DSSet;

/** iterator through a set of Dynamical Systems */
typedef DSSet::iterator DSIterator;

/** const iterator through a set of Dynamical Systems */
typedef DSSet::const_iterator ConstDSIterator;

/** return type value for insert function - bool = false if insertion failed. */
typedef std::pair<DSSet::iterator, bool> CheckInsertDS;

/** Set (STL) of pointers to Dynamical Systems - Useful for NSDS, Interaction, OSI ...
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.0.
 *  \date (Creation) May 17, 2006
 *
 * DS are sorted in according to their id-number in a growing order.
 * Only one occurence of a DS can be present in the set.
 * Possible operations are insert, erase, get or check presence of Dynamical Systems.
 *
 * Warning: any call to operator = or copy constructor results in a false copy: (if set1=set2, pointers to ds of set1 are equal to those of set2)
 *
 */
class DynamicalSystemsSet
{
protected:

  /** a set of DynamicalSystem, sorted thanks to their id number */
  DSSet setOfDS;

  /** a map of bool to check inside-class allocation.
   *  isDSAllocatedIn[ds] = true if ds has been allocated in a method of the present class.
   */
  std::map<DynamicalSystem*, bool > isDSAllocatedIn;

public:

  /** default constructor
   */
  DynamicalSystemsSet();

  /** copy constructor
   *  \param a DynamicalSystemsSet to be copied
   */
  DynamicalSystemsSet(const DynamicalSystemsSet&);

  /** destructor
   */
  ~DynamicalSystemsSet();

  /** assignment
   */
  DynamicalSystemsSet& operator=(const DynamicalSystemsSet&);

  /** return the number of DS in the set
   *  \return an unsigned int
   */
  inline const unsigned int size() const
  {
    return setOfDS.size();
  };

  /** iterator equal to the first element of setOfDS
   *  \return a DSIterator
   */
  inline DSIterator begin()
  {
    return setOfDS.begin();
  };

  /** iterator equal to setOfDS.end()
   *  \return a DSIterator
   */
  inline DSIterator end()
  {
    return setOfDS.end();
  }

  /** const iterator equal to the first element of setOfDS
   *  \return a ConstDSIterator
   */
  inline ConstDSIterator begin() const
  {
    return setOfDS.begin();
  };

  /** const iterator equal to setOfDS.end()
   *  \return a ConstDSIterator
   */
  inline ConstDSIterator end() const
  {
    return setOfDS.end();
  }

  /** return setOfDS
   *  \return a DSSet
   */
  inline const DSSet getSetOfDS() const
  {
    return setOfDS;
  }

  /** get Dynamical System number num, if it is present in the set (else, exception)
   *  \return a pointer to DynamicalSystem
   */
  DynamicalSystem* getDynamicalSystemPtr(const int) const;

  /** return true if ds is in the set
   *  \param a pointer to DynamicalSystem
   *  \return a bool
   */
  const bool isDynamicalSystemIn(DynamicalSystem*) const;

  /** return true if DynamicalSystem number num is in the set
   *  \param an int
   *  \return a bool
   */
  const bool isDynamicalSystemIn(const int) const;

  /** same as find function of stl set
   *  \param a pointer to DynamicalSystem
   *  \param a DSIterator
   */
  DSIterator find(DynamicalSystem*);

  /** same as find function of stl set
   *  \param an int
   *  \return a DSIterator
   */
  DSIterator find(const int);

  /** insert Dynamical System ds into the set
   *  \param a pointer to DynamicalSystem
   *  \return a CheckInsertDS ( boolean type information)
   */
  CheckInsertDS insert(DynamicalSystem*);

  /** true if the DS set is empty
   *  \return a bool
   */
  inline const bool isEmpty() const
  {
    return setOfDS.empty();
  };

  /** remove Dynamical System ds from the set
   *  \param a pointer to DynamicalSystem
   */
  void erase(DynamicalSystem*);

  /** remove all Dynamical Systems from the set
   */
  void clear();

  /** screen-display of the numbers of the Dynamical Systems present in the set.
   */
  void display() const;

  /** return the intersection of s1 and s2 (-> set_intersection stl function)
   */
  friend const DynamicalSystemsSet intersection(const DynamicalSystemsSet& s1, const DynamicalSystemsSet& s2);

  /** return the difference betwee s1 and s2 (-> set_difference stl function)
   */
  friend const DynamicalSystemsSet operator-(const DynamicalSystemsSet& s1, const DynamicalSystemsSet& s2);
};

#endif
