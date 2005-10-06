/* Siconos version 1.0, Copyright INRIA 2005.
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

#include "LagrangianR.h"
using namespace std;



// Default constructor with optional interaction parameter
LagrangianR::LagrangianR(Interaction* inter): Relation(inter)
{
  relationType = LAGRANGIANRELATION;
}

// xml constructor
LagrangianR::LagrangianR(RelationXML* relxml, Interaction* inter):
  Relation(relxml, inter), computeJacobianPtr(NULL), computeHPtr(NULL)
{
  relationType = LAGRANGIANRELATION;
}

// constructor from a set of data
LagrangianR::LagrangianR(const string& computeInput, const string& computeOutput, Interaction* inter):
  Relation(inter), computeJacobianPtr(NULL), computeHPtr(NULL)
{
  relationType = LAGRANGIANRELATION;
  // computeInput
  setComputeInputFunction(cShared.getPluginName(computeInput), cShared.getPluginFunctionName(computeInput));
  // computeOutput
  setComputeOutputFunction(cShared.getPluginName(computeOutput), cShared.getPluginFunctionName(computeOutput));
}

// copy constructor (inter is optional)
LagrangianR::LagrangianR(const Relation & newLNLR, Interaction* inter):
  Relation(newLNLR, inter)
{
  if (relationType !=  LAGRANGIANRELATION || relationType !=  LAGRANGIANLINEARRELATION)
    RuntimeException::selfThrow("LagrangianR:: copy constructor, inconsistent relation types for copy");

  //const LagrangianR * lnlr = static_cast<const LagrangianR*>(&newLNLR);
  // \todo

}


LagrangianR::~LagrangianR()
{}

void LagrangianR::computeJacobian()
{
  if (computeJacobianPtr == NULL)
    RuntimeException::selfThrow("computeJacobian() is not linked to a plugin function");
  computeJacobianPtr(NULL, NULL, NULL, NULL);
}

void LagrangianR::saveRelationToXML()
{
  IN("LagrangianR::saveRelationToXML\n");
  if (relationxml != NULL)
  {
    relationxml->setComputeInputPlugin(computeInputName);
    relationxml->setComputeOutputPlugin(computeOutputName);
  }
  else RuntimeException::selfThrow("LagrangianR::saveRelationToXML - object RelationXML does not exist");
  OUT("LagrangianR::saveRelationToXML\n");
}

LagrangianR* LagrangianR::convert(Relation *r)
{
  cout << "LagrangianR::convert (Relation *r)" << endl;
  LagrangianR* lnlr = dynamic_cast<LagrangianR*>(r);
  return lnlr;
}

