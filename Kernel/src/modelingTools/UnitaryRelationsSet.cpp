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
#include "UnitaryRelationsSet.h"
#include "UnitaryRelation.h"
#include "Interaction.h"
#include "RuntimeException.h"
using namespace std;

//-- Default constructor --
UnitaryRelationsSet::UnitaryRelationsSet()
{}

//-- Copy constructor --
UnitaryRelationsSet::UnitaryRelationsSet(const UnitaryRelationsSet& newSet)
{

  // Warning: "false" copy since pointers links remain between Dynamical Systems of each set
  clear();
  ConstUnitaryRelationIterator it;
  for (it = newSet.begin(); it != newSet.end(); ++it)
  {
    setOfUnitaryRelations.insert(*it);
    isURAllocatedIn[*it] = false ;
  }
}

// --- Destructor ---
UnitaryRelationsSet::~UnitaryRelationsSet()
{
  UnitaryRelationIterator it;
  for (it = setOfUnitaryRelations.begin(); it != setOfUnitaryRelations.end(); ++it)
  {
    if (isURAllocatedIn[*it]) delete *it;
  }
  setOfUnitaryRelations.clear();
  isURAllocatedIn.clear();
}

UnitaryRelationsSet& UnitaryRelationsSet::operator=(const UnitaryRelationsSet& newSet)
{
  // Warning: "false" copy since pointers links remain between Dynamical Systems of each set
  clear();
  ConstUnitaryRelationIterator it;
  for (it = newSet.begin(); it != newSet.end(); ++it)
  {
    setOfUnitaryRelations.insert(*it);
    isURAllocatedIn[*it] = false ;
  }
  return *this;
}

const bool UnitaryRelationsSet::isUnitaryRelationIn(UnitaryRelation* ur) const
{
  UnitaryRelationIterator it = setOfUnitaryRelations.find(ur);
  bool out = false;
  if (it != setOfUnitaryRelations.end()) out = true;
  return out;
}

UnitaryRelationIterator UnitaryRelationsSet::find(UnitaryRelation* ur)
{
  return setOfUnitaryRelations.find(ur);
}

CheckInsertUnitaryRelation UnitaryRelationsSet::insert(UnitaryRelation* ur)
{
  return setOfUnitaryRelations.insert(ur);
  isURAllocatedIn[ur] = false;
}

void UnitaryRelationsSet::erase(UnitaryRelation* ur)
{
  UnitaryRelationIterator it = setOfUnitaryRelations.find(ur);
  if (it == setOfUnitaryRelations.end())
    RuntimeException::selfThrow("UnitaryRelationsSet::erase(ur): ur is not in the set!");

  // If ur has been allocated inside the class, it is first necessary to release memory.
  if (isURAllocatedIn[ur])
    delete *it;
  isURAllocatedIn.erase(ur);
  setOfUnitaryRelations.erase(*it);
}

void UnitaryRelationsSet::clear()
{
  UnitaryRelationIterator it;
  for (it = setOfUnitaryRelations.begin(); it != setOfUnitaryRelations.end(); ++it)
  {
    if (isURAllocatedIn[*it]) delete *it;
  }
  setOfUnitaryRelations.clear();
  isURAllocatedIn.clear();
}

void UnitaryRelationsSet::display() const
{
  cout << "====> UnitaryRelationsSet display: " << endl;
  cout << "There is(are) " << setOfUnitaryRelations.size() << " Unitary Relation(s) in the set (see the list below)." << endl;
  UnitaryRelationIterator it;
  for (it = setOfUnitaryRelations.begin(); it != setOfUnitaryRelations.end(); ++it)
    cout << "- UR belongs to Interaction named " << (*it)->getInteractionPtr()->getId() << ", for the relation number: " << (*it)->getRelativePosition() << endl;
  cout << endl;
  cout << "=============================================================================================" << endl;
}

const UnitaryRelationsSet intersection(const UnitaryRelationsSet& s1, const UnitaryRelationsSet& s2)
{
  // output
  UnitaryRelationsSet commonUnitaryRelations;

  //  insert_iterator<UnitaryRelationSet> res_ins(commonUnitaryRelation.getSetOfUnitaryRelations(), commonUnitaryRelation.getSetOfUnitaryRelations().begin());
  set_intersection(s1.setOfUnitaryRelations.begin(), s1.setOfUnitaryRelations.end(), s2.setOfUnitaryRelations.begin(), s2.setOfUnitaryRelations.end(),
                   inserter(commonUnitaryRelations.setOfUnitaryRelations, commonUnitaryRelations.setOfUnitaryRelations.begin()), compareUR());

  return commonUnitaryRelations;
}

const UnitaryRelationsSet operator - (const UnitaryRelationsSet& s1, const UnitaryRelationsSet& s2)
{
  // output
  UnitaryRelationsSet commonUnitaryRelations;

  set_difference(s1.setOfUnitaryRelations.begin(), s1.setOfUnitaryRelations.end(), s2.setOfUnitaryRelations.begin(), s2.setOfUnitaryRelations.end(),
                 inserter(commonUnitaryRelations.setOfUnitaryRelations, commonUnitaryRelations.setOfUnitaryRelations.begin()), compareUR());

  return commonUnitaryRelations;
}
