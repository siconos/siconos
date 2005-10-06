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
#include "Relation.h"
using namespace std;


// Default constructor with optional interaction parameter
Relation::Relation(Interaction* inter): relationType("none"), interaction(inter), relationxml(NULL),
  computeInputName("none"), computeOutputName("none"),
  computeOutputPtr(NULL), computeInputPtr(NULL)
{}

// xml constructor
Relation::Relation(RelationXML* relxml, Interaction* inter):
  relationType("none"), interaction(inter), relationxml(relxml),
  computeInputName("none"), computeOutputName("none")
{
  if (relationxml != NULL)
  {
    string plugin;

    // computeInput
    if (relationxml->hasComputeInput())
    {
      plugin = (relationxml)->getComputeInputPlugin();
      setComputeInputFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    // computeOutput
    if (relationxml->hasComputeOutput())
    {
      plugin = (relationxml)->getComputeOutputPlugin();
      setComputeOutputFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
  }
  else RuntimeException::selfThrow("Relation::fillRelationWithRelationXML - object RelationXML does not exist");
}

// copy constructor (inter is optional)
Relation::Relation(const Relation& newRel, Interaction* inter):
  relationType(newRel.getType()), interaction(inter), relationxml(NULL),
  computeInputName(newRel.getComputeInputName()), computeOutputName(newRel.getComputeOutputName()),
  computeOutputPtr(NULL), computeInputPtr(NULL)
{
  // \warning:  interaction, relationxml and dsioVector are not copied !
  // Interaction can be set with optional parameter inter (default=NULL)
  // \todo: manage dsio copy when this class will be well implemented
}



Relation::~Relation()
{}

vector<DSInputOutput*> Relation::getDSInputOutputs(void)
{
  return dsioVector;
}

DSInputOutput* Relation::getDSInputOutput(const unsigned int& i)
{
  if (i >= dsioVector.size())
    RuntimeException::selfThrow("Relation - getDSInputOutput : \'i\' is out of range");
  return dsioVector[i];
}

void Relation::setDSInputOutputs(vector<DSInputOutput*> dsioVect)
{
  dsioVector = dsioVect;
}

void Relation::addDSInputOutput(DSInputOutput* dsio)
{
  /*
   *  in EqualityConstraint class, we don't create new objects in the DSInputOutput vector
   *    => we only save a link (pointer) on the DSInputOutputs of the DynamicalSystems !!
   */
  dsioVector.push_back(dsio);
}



void Relation::computeOutput(const double& time)
{
  //to do
  RuntimeException::selfThrow("Relation - computeOutput: not yet implemented for relation of type" + getType());
}

void Relation::computeFreeOutput(const double& time)
{
  RuntimeException::selfThrow("Relation - computeFreeOutput: not yet implemented for relation of type" + getType());
  //to do
}

void Relation::computeInput(const double& time)
{
  //to do
  RuntimeException::selfThrow("Relation - computeIntput: not yet implemented for relation of type" + getType());

}

void Relation::setComputeOutputFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&computeOutputPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  computeOutputName = plugin + ":" + functionName;
}

void Relation::setComputeInputFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&computeInputPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  computeInputName = plugin + ":" + functionName;
}

