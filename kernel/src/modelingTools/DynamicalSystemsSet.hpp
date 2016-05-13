/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
/*! \file DynamicalSystemsSet.hpp
Set of SP::DynamicalSystem
*/
#ifndef DynamicalSystemsSet_H
#define DynamicalSystemsSet_H

#include "DynamicalSystem.hpp"

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

#endif
