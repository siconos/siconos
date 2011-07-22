/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
/*! \file DynamicalSystemsSet.hpp
Set of SP::DynamicalSystem
*/
#ifndef DSSET_H
#define DSSET_H

#include "DynamicalSystem.hpp"
#include "SiconosSet.hpp"

#include <boost/shared_ptr.hpp>

/** DEPRECATED : A set of pointers to dynamical systems, sorted in a growing order
    according to their numbers */
class DynamicalSystemsSet : public std::vector<SP::DynamicalSystem>
{
protected:
  ACCEPT_SERIALIZATION(DynamicalSystemsSet);

public:
  typedef std::vector<SP::DynamicalSystem>::iterator iterator;
  typedef std::vector<SP::DynamicalSystem>::const_iterator const_iterator;
  DynamicalSystemsSet() {};

  bool isEmpty() const
  {
    return size() == 0. ;
  };

  SP::DynamicalSystem getPtr(int num)
  {
    SP::DynamicalSystem ds;
    for (iterator it = begin(); it != end(); ++it)
    {
      if ((*it)->number() == num)
      {
        ds = *it;
        break;
      }
    }
    if (!ds) RuntimeException::selfThrow("getPtr: dynamicalSystem no found");
    return ds;
  };

  bool isIn(SP::DynamicalSystem& ds)
  {
    bool find = false;
    for (iterator it = begin(); it != end(); ++it)
    {
      if ((*it)->number() == ds->number())
      {
        find = true;
        break;
      }
    }
    return find;
  };

  bool isIn(int num)
  {
    bool find = false;
    for (iterator it = begin(); it != end(); ++it)
    {
      if ((*it)->number() == num)
      {
        find = true;
        break;
      }
    }
    return find;
  };

  iterator find(int num)
  {
    iterator it;
    for (it = begin(); it != end(); ++it)
    {
      if ((*it)->number() == num)
      {
        break;
      }
    }
    return it;
  };

  void insert(SP::DynamicalSystem s)
  {
    if (! isIn(s))
      this->push_back(s);
  };

  void insert(iterator begin,
              iterator end)
  {
    for (iterator it = begin; it != end; ++it)
    {
      insert(*it);
    }
  }

  void insert(const_iterator begin,
              const_iterator end)
  {
    for (const_iterator it = begin; it != end; ++it)
    {
      insert(*it);
    }
  }

};

/** Iterator through a set of DS */
typedef DynamicalSystemsSet::iterator DSIterator;

/** const Iterator through a set of DS */
typedef DynamicalSystemsSet::const_iterator ConstDSIterator;

/** return type value for insert function - bool = false if insertion
    failed. */
typedef std::pair<DSIterator, bool> CheckInsertDS;

TYPEDEF_SPTR(DynamicalSystemsSet);

#endif
