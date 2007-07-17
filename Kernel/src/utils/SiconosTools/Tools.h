/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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

/*! \file Tools.h
  Various useful functions and typedef.
*/

#ifndef TOOLS
#define TOOLS

#include<string>
#include <sstream>
#include <vector>
#include <map>
#include <deque>

/** A map to link string to bool (for plug-in flags)  */
typedef std::map<std::string, bool> BoolMap;

/** Type used for inside-class allocation checking */
typedef std::deque<bool> AllocationFlags;

/** Map used to save the list of plug-in names. */
typedef std::map<std::string, std::string> NamesList;

/** Iterator through a list of names. */
typedef NamesList::iterator NamesIterator;

/** const Iterator through a list of names. */
typedef NamesList::const_iterator NamesConstIterator;

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


#endif
