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
#include "Relation.h"
#include "RelationXML.h"
#include "Interaction.h"
#include "FirstOrderNonLinearDS.h"

using namespace std;

// Default constructor
Relation::Relation(RELATIONTYPES newType, RELATIONSUBTYPES newSub):
  relationType(newType), subType(newSub), interaction(NULL), relationxml(NULL), workX(NULL), workZ(NULL), workY(NULL), workL(NULL)
{}

// xml constructor
Relation::Relation(RelationXML* relxml, RELATIONTYPES newType, RELATIONSUBTYPES newSub):
  relationType(newType), subType(newSub), interaction(NULL), relationxml(relxml), workX(NULL), workZ(NULL), workY(NULL), workL(NULL)
{
  if (relationxml == NULL)
    RuntimeException::selfThrow("Relation::fillRelationWithRelationXML - object RelationXML does not exist");
}

Relation::~Relation()
{
  if (workX != NULL) delete workX;
  workX = NULL;
  if (workZ != NULL) delete workZ;
  workZ = NULL;
  if (workY != NULL) delete workY;
  workY = NULL;
  if (workL != NULL) delete workL;
  workL = NULL;
}

void Relation::display() const
{
  cout << "=====> Relation of type " << relationType << " and subtype " << subType << endl;
  NamesConstIterator it;
  cout << "The following operators are linked to plug-in: " << endl;
  for (it = pluginNames.begin(); it != pluginNames.end(); ++it)
    cout << (*it).first << " plugged to:" << (*it).second << endl;
  cout << "===================================== " << endl;
}

void Relation::saveRelationToXML() const
{
  RuntimeException::selfThrow("Relation - saveRelationToXML: not yet implemented for relation of type " + getType());
}
