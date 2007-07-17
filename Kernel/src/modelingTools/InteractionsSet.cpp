/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
#include "InteractionsSet.h"
#include "RuntimeException.h"
#include "Interaction.h"

using namespace std;

//-- Default constructor --
InteractionsSet::InteractionsSet()
{}

//-- Copy constructor --
InteractionsSet::InteractionsSet(const InteractionsSet& newSet)
{

  // Warning: "false" copy since pointers links remain between Interactions of each set
  clear();
  ConstInteractionsIterator it;
  for (it = newSet.begin(); it != newSet.end(); ++it)
  {
    setOfInteractions.insert(*it);
    isInteractionAllocatedIn[*it] = false ;
  }
  // Warning: Copy of Interactions leads to Interactions with non id and -2 as a number.
  // Thus it is better to avoid InteractionsSet copy.
  //  RuntimeException::selfThrow("InteractionsSet::copy constructor, not implemented. ");
}

// --- Destructor ---
InteractionsSet::~InteractionsSet()
{
  InteractionsIterator it;
  for (it = setOfInteractions.begin(); it != setOfInteractions.end(); ++it)
  {
    if (isInteractionAllocatedIn[*it]) delete(*it);
  }
  setOfInteractions.clear();
  isInteractionAllocatedIn.clear();
}

InteractionsSet& InteractionsSet::operator=(const InteractionsSet& newSet)
{
  // Warning: "false" copy since pointers links remain between Interactions of each set
  clear();
  ConstInteractionsIterator it;
  for (it = newSet.begin(); it != newSet.end(); ++it)
  {
    setOfInteractions.insert(*it);
    isInteractionAllocatedIn[*it] = false ;
  }
  return *this;
}

Interaction* InteractionsSet::getInteraction(const int num) const
{
  ConstInteractionsIterator it;
  for (it = setOfInteractions.begin(); it != setOfInteractions.end(); ++it)
  {
    if (((*it)->getNumber()) == num)
      break;
  }
  if (it == setOfInteractions.end())
    RuntimeException::selfThrow("InteractionsSet::getInteraction(num): can not find this Dynamical System in the set.");

  return *it;
}

const bool InteractionsSet::isInteractionIn(Interaction* ds) const
{
  InteractionsIterator it = setOfInteractions.find(ds);
  bool out = false;
  if (it != setOfInteractions.end()) out = true;
  return out;
}

const bool InteractionsSet::isInteractionIn(const int num) const
{
  bool out = false;
  InteractionsIterator it;
  for (it = setOfInteractions.begin(); it != setOfInteractions.end(); ++it)
  {
    if (((*it)->getNumber()) == num)
    {
      out = true;
      break;
    }
  }

  return out;
}

InteractionsIterator InteractionsSet::find(Interaction* ds)
{
  return setOfInteractions.find(ds);
}

InteractionsIterator InteractionsSet::find(const int num)
{
  InteractionsIterator it;
  for (it = setOfInteractions.begin(); it != setOfInteractions.end(); ++it)
  {
    if (((*it)->getNumber()) == num)
      break;
  }
  return it; // == this.end() if not found.
}

CheckInsertInteraction InteractionsSet::insert(Interaction* ds)
{
  return setOfInteractions.insert(ds);
  isInteractionAllocatedIn[ds] = false;
}

void InteractionsSet::erase(Interaction* ds)
{
  InteractionsIterator it = setOfInteractions.find(ds);
  if (it == setOfInteractions.end())
    RuntimeException::selfThrow("InteractionsSet::erase(ds): ds is not in the set!");

  // If ds has been allocated inside the class, it is first necessary to release memory.
  if (isInteractionAllocatedIn[ds])
    delete *it;
  isInteractionAllocatedIn.erase(ds);
  setOfInteractions.erase(*it);
}

void InteractionsSet::clear()
{
  InteractionsIterator it;
  for (it = setOfInteractions.begin(); it != setOfInteractions.end(); ++it)
  {
    if (isInteractionAllocatedIn[*it]) delete *it;
  }
  setOfInteractions.clear();
  isInteractionAllocatedIn.clear();
}

void InteractionsSet::display() const
{
  cout << "====> InteractionsSet display - The following Interactions are present in the set ( id - number):" << endl;
  InteractionsIterator it;
  for (it = setOfInteractions.begin(); it != setOfInteractions.end(); ++it)
    cout << "(" << (*it)->getId() << "," << (*it)->getNumber() << "), ";
  cout << endl;
  cout << "=============================================================================================" << endl;
}

const InteractionsSet intersection(const InteractionsSet& s1, const InteractionsSet& s2)
{
  // output
  InteractionsSet commonInteractions;

  set_intersection(s1.setOfInteractions.begin(), s1.setOfInteractions.end(), s2.setOfInteractions.begin(), s2.setOfInteractions.end(),
                   inserter(commonInteractions.setOfInteractions, commonInteractions.setOfInteractions.begin()), compareInter());

  return commonInteractions;
}

const InteractionsSet operator - (const InteractionsSet& s1, const InteractionsSet& s2)
{
  // output
  InteractionsSet commonInteractions;

  set_difference(s1.setOfInteractions.begin(), s1.setOfInteractions.end(), s2.setOfInteractions.begin(), s2.setOfInteractions.end(),
                 inserter(commonInteractions.setOfInteractions, commonInteractions.setOfInteractions.begin()), compareInter());

  return commonInteractions;
}
