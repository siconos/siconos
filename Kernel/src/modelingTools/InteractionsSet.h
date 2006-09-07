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
#ifndef InteractionsSET_H
#define InteractionsSET_H

#include "RuntimeException.h"
#include "RuntimeCmp.h"
#include "check.h"
#include <iostream>
#include <set>
#include <map>

class Interaction;

/** \class InteractionsSet
 *  \brief set (stl) of pointers to Interactions - Useful for NSDS, Interaction, OSI ...
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
 *  \date (Creation) May 17, 2006
 *
 * Interactions are sorted according to their id-number in a growing order.
 * Only one occurence of an Interaction can be present in the set.
 * Possible operations are insert, erase, get or check presence of Interactions.
 *
 * Warning: any call to operator = or copy constructor results in a false copy: (if set1=set2, pointers to interaction of set1 are equal to those of set2)
 *
 */

// Structure used for Interactions sorting. The address is used to compare two Interactions.
struct compareInter
{
  bool operator()(const Interaction* I1, const Interaction* I2) const
  {
    return (I1 < I2);
  }
};

/** set of Interactions */
typedef std::set<Interaction*, compareInter > InterSet;

/** iterator through a set of Interactions */
typedef InterSet::iterator InteractionsIterator;

/** const iterator through a set of Interactions */
typedef InterSet::const_iterator ConstInteractionsIterator;

/** return type value for insert function - bool = false if insertion failed. */
typedef std::pair<InterSet::iterator, bool> CheckInsertInteraction;

class InteractionsSet
{
protected:

  /** a set of Interaction, sorted thanks to their id number */
  InterSet setOfInteractions;

  /** a map of bool to check inside-class allocation.
   *  isInteractionAllocatedIn[ds] = true if ds has been allocated in a method of the present class.
   */
  std::map<Interaction*, bool > isInteractionAllocatedIn;

public:

  /** \fn InteractionsSet()
   *  \brief default constructor
   */
  InteractionsSet();

  /** \fn InteractionsSet(const InteractionsSet&)
   *  \brief copy constructor
   *  \param a InteractionsSet to be copied
   */
  InteractionsSet(const InteractionsSet&);

  /** \fn ~InteractionsSet()
   *  \brief destructor
   */
  ~InteractionsSet();

  /** \fn   InteractionsSet& operator=( const InteractionsSet& );
   *  \brief assignment
   */
  InteractionsSet& operator=(const InteractionsSet&);

  /** \fn const unsigned int size() const
   *  \brief return the number of Interactions in the set
   *  \return an unsigned int
   */
  inline const unsigned int size() const
  {
    return setOfInteractions.size();
  };

  /** \fn const bool isEmpty() const
   *  \brief return true if the set is empty, else false
   *  \return a bool
   */
  inline const bool isEmpty() const
  {
    return setOfInteractions.empty();
  };

  /** \fn InteractionsIterator begin() const
   *  \brief return iterator on the first element of setOfInteractions
   *  \return a InteractionsIterator
   */
  inline InteractionsIterator begin() const
  {
    return setOfInteractions.begin();
  };

  /** \fn InteractionsIterator end() const
   *  \brief return iterator on setOfInteractions.end()
   *  \return a InteractionsIterator
   */
  inline InteractionsIterator end() const
  {
    return setOfInteractions.end();
  }

  /** \fn InterSet getSetOfInteractions()
   *  \brief return setOfInteractions
   *  \return an InterSet
   */
  inline const InterSet getSetOfInteractions() const
  {
    return setOfInteractions;
  }

  /** \fn Interaction* getInteraction(const int num)
   *  \brief get Dynamical System number num, if it is present in the set (else, exception)
   *  \return a pointer to Interaction
   */
  Interaction* getInteraction(const int) const;

  /** \fn bool isInteractionIn(Interaction* ds)
   *  \brief return true if ds is in the set
   *  \param a pointer to Interaction
   *  \return a bool
   */
  const bool isInteractionIn(Interaction*) const;

  /** \fn bool isInteractionIn(const int num)
   *  \brief return true if Interaction number num is in the set
   *  \param an int
   *  \return a bool
   */
  const bool isInteractionIn(const int) const;

  /** \fn InteractionsIterator find(Interaction* ds)
   *  \brief same as find function of stl set
   *  \param a pointer to Interaction
   *  \param a InteractionsIterator
   */
  InteractionsIterator find(Interaction*);

  /** \fn InteractionsIterator find(const int num)
   *  \brief same as find function of stl set
   *  \param an int
   *  \return a InteractionsIterator
   */
  InteractionsIterator find(const int);

  /** \fn CheckInsertInteraction insert(Interaction* ds)
   *  \brief insert Dynamical System ds into the set
   *  \param a pointer to Interaction
   *  \return a CheckInsertInteraction (boolean type information)
   */
  CheckInsertInteraction insert(Interaction*);

  /** \fn void erase(Interaction* ds)
   *  \brief remove Dynamical System ds from the set
   *  \param a pointer to Interaction
   */
  void erase(Interaction*);

  /** \fn clear() {setOfInteractions.clear();};
   *  \brief remove all Interactions from the set
   */
  void clear();

  /** \fn void display() const
   *  \brief screen-display of the numbers of the Interactions present in the set.
   */
  void display() const;

  /** \fn const InteractionsSet intersection(const InteractionsSet& s1, const InteractionsSet& s2) const
   *  \brief return the intersection of s1 and s2 (-> set_intersection stl function)
   */
  friend const InteractionsSet intersection(const InteractionsSet& s1, const InteractionsSet& s2);

  /** \fn const InteractionsSet operator-(const InteractionsSet& s1, const InteractionsSet& s2) const
   *  \brief return the difference betwee s1 and s2 (-> set_difference stl function)
   */
  friend const InteractionsSet operator-(const InteractionsSet& s1, const InteractionsSet& s2);
};

#endif
