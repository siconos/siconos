/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
#include "DSSet.h"
#include "LagrangianDS.h"
#include "LagrangianLinearTIDS.h"
#include "LinearDS.h"
#include "LinearTIDS.h"

using namespace std;

// RuntimeCmp object used to sort DynamicalSystems in the ds set
// !!! \todo Find a way to avoid global variable ... !!!
RuntimeCmp<DynamicalSystem> compareDS(&DynamicalSystem::getNumberForSorting);

//-- Default constructor --
DSSet::DSSet(): setOfDS(compareDS)
{}

//-- Copy constructor --
DSSet::DSSet(const DSSet& newSet): setOfDS(compareDS)
{

  // Warning: "false" copy since pointers links remain between Dynamical Systems of each set
  clear();
  ConstDSIterator it;
  for (it = newSet.begin(); it != newSet.end(); ++it)
  {
    setOfDS.insert(*it);
    isDSAllocatedIn[*it] = false ;
  }

  // Warning: Copy of DS leads to DS with non id and -2 as a number.
  // Thus it is better to avoid DSSet copy.
  //  RuntimeException::selfThrow("DSSet::copy constructor, not implemented. ");
  //   constDSIterator it;
  //   string type;
  //   // We sweep the DS set of newSet and, depending on the DS type, insert a new one in setOfDS, by using the DynamicalSystem copy constructor.
  //   for(it=newSet.begin();it!=newSet.end();++it)
  //     {
  //       // \todo use factories to improve this copy.
  //       type = (*it) ->getType();
  //       if( type ==  NLDS)
  //  setOfDS.insert( new DynamicalSystem( **it ));
  //       else if ( type ==  LDS)
  //  setOfDS.insert( new LinearDS( **it ));
  //       else if ( type ==  LITIDS)
  //  setOfDS.insert( new LinearTIDS( **it ));
  //       else if ( type ==  LNLDS)
  //  setOfDS.insert( new LagrangianDS( **it ));
  //       else if ( type ==  LTIDS)
  //  setOfDS.insert( new LagrangianLinearTIDS( **it ));
  //       else
  //  RuntimeException::selfThrow("DSSet::copy constructor, unknown Dynamical system type:"+type);
  //       isDSAllocatedIn[*it] = true ;
  //     }
}

// --- Destructor ---
DSSet::~DSSet()
{
  DSIterator it;
  for (it = setOfDS.begin(); it != setOfDS.end(); ++it)
  {
    if (isDSAllocatedIn[*it]) delete *it;
  }
  setOfDS.clear();
  isDSAllocatedIn.clear();
}

DSSet& DSSet::operator=(const DSSet& newSet)
{
  // Warning: "false" copy since pointers links remain between Dynamical Systems of each set
  clear();
  ConstDSIterator it;
  for (it = newSet.begin(); it != newSet.end(); ++it)
  {
    setOfDS.insert(*it);
    isDSAllocatedIn[*it] = false ;
  }
  return *this;
}

DynamicalSystem* DSSet::getDynamicalSystem(const int& num) const
{
  ConstDSIterator it;
  for (it = setOfDS.begin(); it != setOfDS.end(); ++it)
  {
    if (((*it)->getNumber()) == num)
      break;
  }
  if (it == setOfDS.end())
    RuntimeException::selfThrow("DSSet::getDynamicalSystem(num): can not find this Dynamical System in the set.");

  return *it;
}

const bool DSSet::isDSIn(DynamicalSystem* ds) const
{
  DSIterator it = setOfDS.find(ds);
  bool out = false;
  if (it != setOfDS.end()) out = true;
  return out;
}

const bool DSSet::isDSIn(const int& num) const
{
  bool out = false;
  DSIterator it;
  for (it = setOfDS.begin(); it != setOfDS.end(); ++it)
  {
    if (((*it)->getNumber()) == num)
    {
      out = true;
      break;
    }
  }

  return out;
}

DSIterator DSSet::find(DynamicalSystem* ds)
{
  return setOfDS.find(ds);
}

DSIterator DSSet::find(const int& num)
{
  DSIterator it;
  for (it = setOfDS.begin(); it != setOfDS.end(); ++it)
  {
    if (((*it)->getNumber()) == num)
      break;
  }
  return it; // == this.end() if not found.
}

CheckInsertDS DSSet::insert(DynamicalSystem* ds)
{
  return setOfDS.insert(ds);
  isDSAllocatedIn[ds] = false;
}

void DSSet::erase(DynamicalSystem* ds)
{
  DSIterator it = setOfDS.find(ds);
  if (it == setOfDS.end())
    RuntimeException::selfThrow("DSSet::erase(ds): ds is not in the set!");

  // If ds has been allocated inside the class, it is first necessary to release memory.
  if (isDSAllocatedIn[ds])
    delete *it;
  isDSAllocatedIn.erase(ds);
  setOfDS.erase(*it);
}

void DSSet::clear()
{
  DSIterator it;
  for (it = setOfDS.begin(); it != setOfDS.end(); ++it)
  {
    if (isDSAllocatedIn[*it]) delete *it;
  }
  setOfDS.clear();
  isDSAllocatedIn.clear();
}

void DSSet::display() const
{
  cout << "====> DSSet display - The following Dynamical Systems are present in the set ( id - number):" << endl;
  DSIterator it;
  for (it = setOfDS.begin(); it != setOfDS.end(); ++it)
    cout << "(" << (*it)->getId() << "," << (*it)->getNumber() << "), ";
  cout << endl;
  cout << "=============================================================================================" << endl;
}

const DSSet intersection(const DSSet& s1, const DSSet& s2)
{
  DSSet commonDS;
  DynamicalSystemSet tmp = commonDS.getSetOfDS();
  //  insert_iterator<DynamicalSystemSet> res_ins(commonDS.getSetOfDS(), commonDS.getSetOfDS().begin());
  set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(tmp, tmp.begin()));

  return commonDS;
}

// const DSSet difference(const DSSet& s1, const DSSet& s2)
// {
//   DSSet commonDS;

//   insert_iterator<DynamicalSystemSet> res_ins(commonDS.getSetOfDS(), commonDS.begin());
//   set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), res_ins,&DynamicalSystem::getNumberForSorting compareDS);

//   return commonDS;
// }
