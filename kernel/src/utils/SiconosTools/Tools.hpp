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

/*! \file Tools.hpp
  Various useful functions and typedef.
*/

#ifndef TOOLS_H
#define TOOLS_H

#include<string>
#include <sstream>
#include <vector>
#include <deque>

#include "SiconosPointers.hpp"

/** A vector of unsigned integers used in some cases in kernel */
TYPEDEF_TPL1_SPTR(UnsignedIntVector, std::vector, unsigned int)

/** Type used for inside-class allocation checking */
typedef std::deque<bool> AllocationFlags;


/** A function to convert any type to std::string*/
template <class T> std::string toString(const T& obj)
{
  static std::ostringstream o;
  o.str("");
  o << obj ;
  return o.str();
}

// Next two functions: from Thinking in C++ vol2, p 536.
/** to purge a STL container of pointers, assuming it owns all its pointers.
    \param a STL sequence container
*/
template<class Seq> void purge(Seq& c)
{
  typename Seq::iterator i;
  for (i = c.begin(); i != c.end(); ++i)
  {
    delete *i;
    *i = NULL;
  }
}

/** to purge a STL container of pointers, assuming it owns all its pointers.
    \param iterator equal to a container.begin()
    \param iterator equal to a container.end()
*/
template<class InpIt> void purge(InpIt begin, InpIt end)
{
  while (begin != end)
  {
    delete *begin;
    *begin = NULL;
    ++begin;
  }
}

/** to purge a STL container of pointers; only pointers owned by the container, ie for which second arg
    corresponding value is true, are deleted.
    \param a STL sequence container
    \param a std::vector<bool>
*/
template<class Seq> void purge(Seq& c, const std::vector<bool>& isAllocatedIn)
{
  typename Seq::iterator i;
  std::vector<bool>::const_iterator it = isAllocatedIn.begin();
  for (i = c.begin(); i != c.end(); ++i)
  {
    if (*it ++) delete *i;
    *i = NULL;
  }
}

/** to purge a STL container of pointers, assuming it owns all its pointers.
    \param iterator equal to a container.begin()
    \param iterator equal to a container.end()
*/
template<class InpIt> void purge(InpIt begin, InpIt end, const std::vector<bool>& isAllocatedIn)
{
  std::vector<bool>::const_iterator it = isAllocatedIn.begin();
  while (begin != end)
  {
    if (*it ++) delete *begin;
    *begin = NULL;
    ++begin;
  }
}

/** To copy a value into an object ( created if required)
    \param a smart pointer to T (SPT): the object to be filled - Must have op= and copy constructor from U.
    \param type U param, the value to be assigned.
 */
template <class T, class SPT, class U> void setObject(SPT& obj, const U& val)
{
  if (!obj)
    obj.reset(new T(val));
  else
    *obj = val;
}

#include "SiconosPointers.hpp"
/** Graph -> Set conversion */
template <class S, class G>
std11::shared_ptr<S> setOfGraph(std11::shared_ptr<G> g)
{
  std11::shared_ptr<S> r;
  r.reset(new S());
  for (typename G::VIterator vi = g->begin(); vi != g->end(); ++vi)
  {
    r->insert(g->bundle(*vi));
  }
  return(r);
}

// --- CONSTRUCTORS ---

#endif /* TOOLS_H */

#ifndef PRINTSEQUENCE_H
#define PRINTSEQUENCE_H
#include<algorithm>
#include<iostream>
#include<iterator>

/** Print the contents of any sequence - From Thinking in C++, vol 2 p365.
\param first, any iterator, beginning of the sequence
\param last, any iterator, end of the sequence
\param char*, optional message on top of output, default ""
\param char*, separator between sequence elements, default new line
\param ostream, output destination, default std::cout
*/
template<typename Iter>
void print(Iter first, Iter last, const char* nm = "", const char * sep = "\n", std::ostream& os =  std::cout)
{
  if (nm != NULL && *nm != '\0')
    os << nm << ": " << sep;
  typedef typename std::iterator_traits<Iter>::value_type T;
  std::copy(first, last, std::ostream_iterator<T>(os, sep));
  os << std::endl;
}
#endif // PRINTSEQUENCE_H


