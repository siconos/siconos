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

/*! \file Tools.hpp
  Various useful functions and typedef.
*/

#ifndef TOOLS
#define TOOLS

#include<string>
#include <sstream>
#include <vector>
#include <map>
#include <deque>

/** Type used for inside-class allocation checking */
typedef std::deque<bool> AllocationFlags;


/** A function to convert any type to string*/
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

#ifndef __APPLE__
#include <malloc.h>
static struct mallinfo malloc_info1, malloc_info2;
static int malldelta;
/** Tool to compute the increase of memory alloc */
static int TRM()
{
  malloc_info1 = malloc_info2;
  malloc_info2 = mallinfo();
  malldelta = malloc_info2.uordblks - malloc_info1.uordblks;
  return(malldelta);
}
#endif

#endif

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
\param ostream, output destination, default cout
*/
template<typename Iter>
void print(Iter first, Iter last, const char* nm = "", const char * sep = "\n", std::ostream& os = std::cout)
{
  if (nm != 0 && *nm != '\0')
    os << nm << ": " << sep;
  typedef typename std::iterator_traits<Iter>::value_type T;
  std::copy(first, last, std::ostream_iterator<T>(os, sep));
  os << std::endl;
}
#endif // PRINTSEQUENCE_H


