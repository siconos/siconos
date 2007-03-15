/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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

/*! \file InteractionsSet.h

*/
#ifndef InteractionsSET_H
#define InteractionsSET_H

#include <set>
#include <map>

class Interaction;

/** Structure used for Interactions sorting. The address is used to compare two Interactions. */
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

/** Set (STL) of pointers to Interactions - Useful for NSDS, Interaction, OSI ...
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.0.1.
 *  \date (Creation) May 17, 2006
 *
 * Interactions are sorted according to their id-number in a growing order.
 * Only one occurence of an Interaction can be present in the set.
 * Possible operations are insert, erase, get or check presence of Interactions.
 *
 * Warning: any call to operator = or copy constructor results in a false copy: (if set1=set2, pointers to interaction of set1 are equal to those of set2)
 *
 */
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

  /** default constructor
  */
  InteractionsSet();

  /** copy constructor
  *  \param a InteractionsSet to be copied
  */
  InteractionsSet(const InteractionsSet&);

  /** destructor
  */
  ~InteractionsSet();

  /** assignment
  */
  InteractionsSet& operator=(const InteractionsSet&);

  /** return the number of Interactions in the set
  *  \return an unsigned int
  */
  inline const unsigned int size() const
  {
    return setOfInteractions.size();
  };

  /** return true if the set is empty, else false
  *  \return a bool
  */
  inline const bool isEmpty() const
  {
    return setOfInteractions.empty();
  };

  /** iterator equal to the first element of setOfInteractions
  *  \return a InteractionsIterator
  */
  inline InteractionsIterator begin()
  {
    return setOfInteractions.begin();
  };

  /** iterator equal to setOfInteractions.end()
  *  \return a InteractionsIterator
  */
  inline InteractionsIterator end()
  {
    return setOfInteractions.end();
  }

  /** const iterator equal to the first element of setOfInteractions
  *  \return a InteractionsIterator
  */
  inline ConstInteractionsIterator begin() const
  {
    return setOfInteractions.begin();
  };

  /** const iterator equal to setOfInteractions.end()
  *  \return a InteractionsIterator
  */
  inline ConstInteractionsIterator end() const
  {
    return setOfInteractions.end();
  }

  /** return setOfInteractions
  *  \return an InterSet
  */
  inline const InterSet getSetOfInteractions() const
  {
    return setOfInteractions;
  }

  /** get Dynamical System number num, if it is present in the set (else, exception)
  *  \return a pointer to Interaction
  */
  Interaction* getInteraction(const int) const;

  /** return true if ds is in the set
  *  \param a pointer to Interaction
  *  \return a bool
  */
  const bool isInteractionIn(Interaction*) const;

  /** return true if Interaction number num is in the set
  *  \param an int
  *  \return a bool
  */
  const bool isInteractionIn(const int) const;

  /** same as find function of stl set
  *  \param a pointer to Interaction
  *  \param a InteractionsIterator
  */
  InteractionsIterator find(Interaction*);

  /** same as find function of stl set
  *  \param an int
  *  \return a InteractionsIterator
  */
  InteractionsIterator find(const int);

  /** insert Dynamical System ds into the set
  *  \param a pointer to Interaction
  *  \return a CheckInsertInteraction (boolean type information)
  */
  CheckInsertInteraction insert(Interaction*);

  /** remove Dynamical System ds from the set
  *  \param a pointer to Interaction
  */
  void erase(Interaction*);

  /** remove all Interactions from the set
  */
  void clear();

  /** screen-display of the numbers of the Interactions present in the set.
  */
  void display() const;

  /** return the intersection of s1 and s2 (-> set_intersection stl function)
  */
  friend const InteractionsSet intersection(const InteractionsSet& s1, const InteractionsSet& s2);

  /** return the difference betwee s1 and s2 (-> set_difference stl function)
  */
  friend const InteractionsSet operator-(const InteractionsSet& s1, const InteractionsSet& s2);
};

#endif
