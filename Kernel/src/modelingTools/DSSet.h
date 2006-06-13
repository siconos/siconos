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
#ifndef DSSET_H
#define DSSET_H

#include "RuntimeCmp.h"
#include "check.h"
#include <iostream>
#include <set>
#include <map>

class DynamicalSystem;

/** \class DSSet
 *  \brief set (stl) of pointers to Dynamical Systems - Useful for NSDS, Interaction, OSI ...
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.2.0.
 *  \date (Creation) May 17, 2006
 *
 * DS are sorted in according to their id-number in a growing order.
 * Only one occurence of a DS can be present in the set.
 * Possible operations are insert, erase, get or check presence of Dynamical Systems.
 *
 * Warning: any call to operator = or copy constructor results in a false copy: (if set1=set2, pointers to ds of set1 are equal to those of set2)
 *
 */

/** set of Dynamical Systems */
typedef std::set<DynamicalSystem*, RuntimeCmp<DynamicalSystem> > DynamicalSystemSet;

/** iterator through a set of Dynamical Systems */
typedef DynamicalSystemSet::iterator DSIterator;

/** const iterator through a set of Dynamical Systems */
typedef DynamicalSystemSet::const_iterator ConstDSIterator;

/** return type value for insert function - bool = false if insertion failed. */
typedef std::pair<DynamicalSystemSet::iterator, bool> CheckInsertDS;

class DSSet
{
protected:

  /** a set of DynamicalSystem, sorted thanks to their id number */
  DynamicalSystemSet setOfDS;

  /** a map of bool to check inside-class allocation.
   *  isDSAllocatedIn[ds] = true if ds has been allocated in a method of the present class.
   */
  std::map<DynamicalSystem*, bool > isDSAllocatedIn;

public:

  /** \fn DSSet()
   *  \brief default constructor
   */
  DSSet();

  /** \fn DSSet(const DSSet&)
   *  \brief copy constructor
   *  \param a DSSet to be copied
   */
  DSSet(const DSSet&);

  /** \fn ~DSSet()
   *  \brief destructor
   */
  ~DSSet();

  /** \fn   DSSet& operator=( const DSSet& );
   *  \brief assignment
   */
  DSSet& operator=(const DSSet&);

  /** \fn const unsigned int size() const
   *  \brief return the number of DS in the set
   *  \return an unsigned int
   */
  inline const unsigned int size() const
  {
    return setOfDS.size();
  };

  /** \fn DSIterator begin() const
   *  \brief return iterator on the first element of setOfDS
   *  \return a DSIterator
   */
  inline DSIterator begin() const
  {
    return setOfDS.begin();
  };

  /** \fn DSIterator end() const
   *  \brief return iterator on setOfDS.end()
   *  \return a DSIterator
   */
  inline DSIterator end() const
  {
    return setOfDS.end();
  }

  /** \fn DynamicalSystemSet getSetOfDS()
   *  \brief return setOfDS
   *  \return a DynamicalSystemSet
   */
  inline const DynamicalSystemSet getSetOfDS() const
  {
    return setOfDS;
  }

  /** \fn DynamicalSystem* getDynamicalSystem(const int& num)
   *  \brief get Dynamical System number num, if it is present in the set (else, exception)
   *  \return a pointer to DynamicalSystem
   */
  DynamicalSystem* getDynamicalSystem(const int&) const;

  /** \fn bool isDSIn(DynamicalSystem* ds)
   *  \brief return true if ds is in the set
   *  \param a pointer to DynamicalSystem
   *  \return a bool
   */
  const bool isDSIn(DynamicalSystem*) const;

  /** \fn bool isDSIn(const int& num)
   *  \brief return true if DynamicalSystem number num is in the set
   *  \param an int
   *  \return a bool
   */
  const bool isDSIn(const int&) const;

  /** \fn DSIterator find(DynamicalSystem* ds)
   *  \brief same as find function of stl set
   *  \param a pointer to DynamicalSystem
   *  \param a DSIterator
   */
  DSIterator find(DynamicalSystem*);

  /** \fn DSIterator find(const int& num)
   *  \brief same as find function of stl set
   *  \param an int
   *  \return a DSIterator
   */
  DSIterator find(const int&);

  /** \fn CheckInsertDS insert(DynamicalSystem* ds)
   *  \brief insert Dynamical System ds into the set
   *  \param a pointer to DynamicalSystem
   *  \return a CheckInsertDS ( boolean type information)
   */
  CheckInsertDS insert(DynamicalSystem*);

  /** \fn inline const bool isEmpty() const
   *  \brief true if the DS set is empty
   *  \return a bool
   */
  inline const bool isEmpty() const
  {
    return setOfDS.empty();
  };

  /** \fn void erase(DynamicalSystem* ds)
   *  \brief remove Dynamical System ds from the set
   *  \param a pointer to DynamicalSystem
   */
  void erase(DynamicalSystem*);

  /** \fn clear() {setOfDS.clear();};
   *  \brief remove all Dynamical Systems from the set
   */
  void clear();

  /** \fn void display() const
   *  \brief screen-display of the numbers of the Dynamical Systems present in the set.
   */
  void display() const;

  /** \fn const DSSet intersection(const DSSet& s1, const DSSet& s2) const
   *  \brief return the intersection of s1 and s2 (-> set_intersection stl function)
   */
  friend const DSSet intersection(const DSSet& s1, const DSSet& s2);

  /** \fn const DSSet operator-(const DSSet& s1, const DSSet& s2) const
   *  \brief return the difference betwee s1 and s2 (-> set_difference stl function)
   */
  friend const DSSet operator-(const DSSet& s1, const DSSet& s2);
};

#endif
