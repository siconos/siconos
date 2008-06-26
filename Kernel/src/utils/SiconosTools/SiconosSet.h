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

/*! \file SiconosSet.h
  Template class to define a set of Siconos object.

*/
#ifndef SiconosSet_H
#define SiconosSet_H

#include <set>
#include <map>
#include <iostream>
#include "Cmp.h"
#include "RuntimeException.h"

/** Set (STL) of pointers to T -
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) May 17, 2006
 *
 * A SiconosSet<T,U> handles a set of pointers to T, sorted in growing \n
 * order, according to the value (of type U) returned by a member \n
 * function "getSort" of class T.\n
 * Thus: T object must have a function named getSort:
 * \code
 * const U getSort() const
 * \endcode
 * "<" being defined for type-U objects. \n
 * See for example Interaction.h or DynamicalSystem.h
 *
 * Possible operations are insert, erase, get or check presence of T, \n
 * intersection or difference of sets (in a mathematical sense).
 *
 */
template <class T, class U> class SiconosSet
{
public:

  /** set of T */
  typedef std::set<T*, Cmp<T, U> > TSet;

  /** iterator through a set of T */
  typedef typename TSet::iterator TIterator;

  /** const iterator through a set of T */
  typedef typename TSet::const_iterator ConstTIterator;

  /** return type value for insert function - bool = false if insertion failed. */
  typedef typename std::pair<TIterator, bool> CheckInsertT;

protected:

  /** Pointer to function used in ordering relation */
  const U(T::*fpt)() const;

  /** a set of T, sorted thanks to their address */
  TSet * setOfT;

  /** a map of bool to check inside-class allocation.
   *  isTAllocatedIn[ds] = true if ds has been allocated in a method of the present class.
   */
  std::map<T*, bool > isTAllocatedIn;

private:

  /** copy constructor, private => forbidden
   *  \param a SiconosSet to be copied
   */
  SiconosSet(const SiconosSet&);

public:

  /** default constructor
   */
  SiconosSet(): fpt(NULL), setOfT(NULL)
  {
    fpt = &T::getSort;
    setOfT = new TSet(fpt);
  };

  /** destructor */
  ~SiconosSet()
  {
    TIterator it;
    for (it = setOfT->begin(); it != setOfT->end(); ++it)
    {
      if (isTAllocatedIn[*it]) delete(*it);
    }
    setOfT->clear();
    isTAllocatedIn.clear();
    delete setOfT;
  };

  /** return the number of Ts in the set
   *  \return an unsigned int
   */
  const unsigned int size() const
  {
    return setOfT->size();
  };

  /** return true if the set is empty, else false
   *  \return a bool
   */
  const bool isEmpty() const
  {
    return setOfT->empty();
  };

  /** iterator equal to the first element of setOfT
   *  \return a TIterator
   */
  TIterator begin()
  {
    return setOfT->begin();
  };

  /** iterator equal to setOfT->end()
   *  \return a TIterator
   */
  TIterator end()
  {
    return setOfT->end();
  }

  /** const iterator equal to the first element of setOfT
   *  \return a TIterator
   */
  ConstTIterator begin() const
  {
    return setOfT->begin();
  };

  /** const iterator equal to setOfT->end()
   *  \return a TIterator
   */
  ConstTIterator end() const
  {
    return setOfT->end();
  }

  /** return setOfT
   *  \return an InterSet
   */
  const TSet * getSet() const
  {
    return setOfT;
  }

  /** get T number num, if it is present in the set (else, exception)
   *  \return a pointer to T
   */
  T* getPtr(int num) const
  {
    ConstTIterator it;
    for (it = setOfT->begin(); it != setOfT->end(); ++it)
    {
      if (((*it)->getNumber()) == num)
        break;
    }
    if (it == setOfT->end())
      RuntimeException::selfThrow("SiconosSet::get(num): can not find an object number ""num"" in the set.");
    return *it;
  };

  /** return true if ds is in the set
   *  \param a pointer to T
   *  \return a bool
   */
  const bool isIn(T* t) const
  {
    TIterator it = setOfT->find(t);
    bool out = false;
    if (it != setOfT->end()) out = true;
    return out;
  };

  /** return true if T number num is in the set
   *  \param an int
   *  \return a bool
   */
  const bool isIn(int num) const
  {
    bool out = false;
    TIterator it;
    for (it = setOfT->begin(); it != setOfT->end(); ++it)
    {
      if (((*it)->getNumber()) == num)
      {
        out = true;
        break;
      }
    }
    return out;
  }


  /** same as find function of stl set
   *  \param a pointer to T
   *  \param a TIterator
   */
  TIterator find(T* t)
  {
    return setOfT->find(t);
  };

  /** same as find function of stl set
   *  \param an int
   *  \return a TIterator
   */
  TIterator find(int num)
  {
    TIterator it;
    for (it = setOfT->begin(); it != setOfT->end(); ++it)
    {
      if (((*it)->getNumber()) == num)
        break;
    }
    return it; // == this.end() if not found.
  };

  /** insert a T* into the set
   *  \param a pointer to T
   *  \return a CheckInsertT (boolean type information)
   */
  CheckInsertT insert(T* t)
  {
    return setOfT->insert(t);
    isTAllocatedIn[t] = false;
  };

  /** remove a T* from the set
   *  \param a pointer to T
   */
  void erase(T* t)
  {
    TIterator it = setOfT->find(t);
    if (it == setOfT->end())
      RuntimeException::selfThrow("SiconosSet::erase(t): t is not in the set!");

    // If ds has been allocated inside the class, it is first necessary to release memory.
    if (isTAllocatedIn[t])
      delete *it;
    isTAllocatedIn.erase(t);
    setOfT->erase(*it);
  };

  /** remove all Ts from the set
   */
  void clear()
  {
    TIterator it;
    for (it = setOfT->begin(); it != setOfT->end(); ++it)
    {
      if (isTAllocatedIn[*it]) delete *it;
    }
    setOfT->clear();
    isTAllocatedIn.clear();
  };

  /** screen-display of the numbers of the Ts present in the set.
   */
  void display() const
  {
    std::cout << "====> Set display - The following objects are present in the set (number):" << std::endl;
    TIterator it;
    for (it = setOfT->begin(); it != setOfT->end(); ++it)
      std::cout << "(" << (*it)->getNumber() << "), ";
    std::cout << std::endl;
    std::cout << "=============================================================================================" << std::endl;
  };

  /**   computes in s3 the intersection of s1 and s2 (-> set_intersection stl function) */
  /*       \param SiconosSet, s1 */
  /*       \param SiconosSet, s2 */
  /*       \param SiconosSet, s3 */
  friend void intersection(const SiconosSet& s1, const SiconosSet& s2, SiconosSet& commonT)
  {
    set_intersection(s1.setOfT->begin(), s1.setOfT->end(), s2.setOfT->begin(), s2.setOfT->end(),
                     inserter(*commonT.setOfT, commonT.setOfT->begin()), (commonT.setOfT)->value_comp());
  };

  /** computes in s3 the difference betwee s1 and s2 (-> set_difference stl function)
  \param SiconosSet, s1
  \param SiconosSet, s2
  \param SiconosSet, s3
  */
  friend void difference(const SiconosSet& s1, const SiconosSet& s2, SiconosSet& commonT)
  {
    set_difference(s1.setOfT->begin(), s1.setOfT->end(), s2.setOfT->begin(), s2.setOfT->end(),
                   inserter(*commonT.setOfT, commonT.setOfT->begin()), (commonT.setOfT)->value_comp());
  };
};

#endif


