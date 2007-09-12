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
#include "DynamicalSystemsSet.h"
#include "DynamicalSystem.h"

using namespace std;

// RuntimeCmp object used to sort DynamicalSystems in the ds set
// !!! \todo Find a way to avoid global variable ... !!!
RuntimeCmp<DynamicalSystem> compareDS(&DynamicalSystem::getNumberForSorting);

//-- Default constructor --
DynamicalSystemsSet::DynamicalSystemsSet(): setOfDS(NULL)
{
  setOfDS = new DSSet(compareDS);
  //setOfDS = new DSSet();
}

// --- Destructor ---
DynamicalSystemsSet::~DynamicalSystemsSet()
{
  DSIterator it;
  for (it = setOfDS->begin(); it != setOfDS->end(); ++it)
  {
    if (isDSAllocatedIn[*it]) delete *it;
  }
  setOfDS->clear();
  isDSAllocatedIn.clear();
  delete setOfDS;
}

DynamicalSystem* DynamicalSystemsSet::getDynamicalSystemPtr(int num) const
{
  ConstDSIterator it;
  for (it = setOfDS->begin(); it != setOfDS->end(); ++it)
  {
    if (((*it)->getNumber()) == num)
      break;
  }
  if (it == setOfDS->end())
    RuntimeException::selfThrow("DynamicalSystemsSet::getDynamicalSystem(num): can not find this Dynamical System in the set.");

  return *it;
}

const bool DynamicalSystemsSet::isDynamicalSystemIn(DynamicalSystem* ds) const
{
  DSIterator it = setOfDS->find(ds);
  bool out = false;
  if (it != setOfDS->end()) out = true;
  return out;
}

const bool DynamicalSystemsSet::isDynamicalSystemIn(const int num) const
{
  bool out = false;
  DSIterator it;
  for (it = setOfDS->begin(); it != setOfDS->end(); ++it)
  {
    if (((*it)->getNumber()) == num)
    {
      out = true;
      break;
    }
  }

  return out;
}

DSIterator DynamicalSystemsSet::find(DynamicalSystem* ds)
{
  return setOfDS->find(ds);
}

DSIterator DynamicalSystemsSet::find(const int num)
{
  DSIterator it;
  for (it = setOfDS->begin(); it != setOfDS->end(); ++it)
  {
    if (((*it)->getNumber()) == num)
      break;
  }
  return it; // == this.end() if not found.
}

CheckInsertDS DynamicalSystemsSet::insert(DynamicalSystem* ds)
{
  return setOfDS->insert(ds);
  isDSAllocatedIn[ds] = false;
}

void DynamicalSystemsSet::erase(DynamicalSystem* ds)
{
  DSIterator it = setOfDS->find(ds);
  if (it == setOfDS->end())
    RuntimeException::selfThrow("DynamicalSystemsSet::erase(ds): ds is not in the set!");

  // If ds has been allocated inside the class, it is first necessary to release memory.
  if (isDSAllocatedIn[ds])
    delete *it;
  isDSAllocatedIn.erase(ds);
  setOfDS->erase(*it);
}

void DynamicalSystemsSet::clear()
{
  DSIterator it;
  for (it = setOfDS->begin(); it != setOfDS->end(); ++it)
  {
    if (isDSAllocatedIn[*it]) delete *it;
  }
  setOfDS->clear();
  isDSAllocatedIn.clear();
}

void DynamicalSystemsSet::display() const
{
  cout << "====> DynamicalSystemsSet display - The following Dynamical Systems are present in the set ( id - number):" << endl;
  DSIterator it;
  for (it = setOfDS->begin(); it != setOfDS->end(); ++it)
    cout << "(" << (*it)->getId() << "," << (*it)->getNumber() << "), ";
  cout << endl;
  cout << "=============================================================================================" << endl;
}

void intersection(const DynamicalSystemsSet& s1, const DynamicalSystemsSet& s2, DynamicalSystemsSet& commonDS)
{
  set_intersection(s1.setOfDS->begin(), s1.setOfDS->end(), s2.setOfDS->begin(), s2.setOfDS->end(),
                   inserter(*commonDS.setOfDS, commonDS.setOfDS->begin()), compareDS);
}

void difference(const DynamicalSystemsSet& s1, const DynamicalSystemsSet& s2, DynamicalSystemsSet& commonDS)
{
  set_difference(s1.setOfDS->begin(), s1.setOfDS->end(), s2.setOfDS->begin(), s2.setOfDS->end(),
                 inserter(*commonDS.setOfDS, commonDS.setOfDS->begin()), compareDS);
}
