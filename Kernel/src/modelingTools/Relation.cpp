/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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

void Relation::initParameter(const string id)
{
  if (parametersList[id] == NULL)
  {
    parametersList[id] = new SimpleVector(1);
    string alloc = "parameter_for_" + id;
    isAllocatedIn[alloc] = true;
  }
}

// Default constructor
Relation::Relation(const string newType):
  relationType(newType), interaction(NULL), relationxml(NULL), computeInputName("none"),
  computeOutputName("none"), isOutputPlugged(false), isInputPlugged(false),
  computeOutputPtr(NULL), computeInputPtr(NULL)
{

  // Set plug-in default functions for h and g.
  setComputeOutputFunction("DefaultPlugin.so", "computeOutput");
  setComputeInputFunction("DefaultPlugin.so", "computeInput");
  isOutputPlugged = false;
  isInputPlugged = false;
}

// xml constructor
Relation::Relation(RelationXML* relxml, const string newType):
  relationType(newType), interaction(NULL), relationxml(relxml),
  computeInputName("none"), computeOutputName("none"),
  isOutputPlugged(true), isInputPlugged(true),
  computeOutputPtr(NULL), computeInputPtr(NULL)
{
  if (relationxml == NULL)
    RuntimeException::selfThrow("Relation::fillRelationWithRelationXML - object RelationXML does not exist");

  string plugin;
  // computeInput
  if (relationxml->hasComputeInput())
  {
    plugin = (relationxml)->getComputeInputPlugin();
    setComputeInputFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
  }
  else
  {
    setComputeInputFunction("DefaultPlugin.so", "computeInput");
    isInputPlugged = false; //
  }

  // computeOutput
  if (relationxml->hasComputeOutput())
  {
    plugin = (relationxml)->getComputeOutputPlugin();
    setComputeOutputFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
  }
  else
  {
    setComputeOutputFunction("DefaultPlugin.so", "computeOutput");
    isOutputPlugged = false; //
  }
}

// copy constructor
Relation::Relation(const Relation& newRel):
  relationType(newRel.getType()), interaction(NULL), relationxml(NULL),
  computeInputName(newRel.getComputeInputName()), computeOutputName(newRel.getComputeOutputName()),
  isOutputPlugged(true), isInputPlugged(true),
  computeOutputPtr(NULL), computeInputPtr(NULL)
{
  // \warning:  interaction, relationxml and dsioVector are not copied !
  // \todo: manage dsio copy when this class will be well implemented

  setParameters(newRel.getParameters());   // Copy !!
  string plugin = computeInputName;
  setComputeInputFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
  if (computeInputName == "DefaultPlugin:computeInput")
    isInputPlugged = false;
  plugin = computeOutputName;
  setComputeOutputFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
  if (computeOutputName == "DefaultPlugin:computeOutput")
    isOutputPlugged = false;

  // Remark : set is...Plugged to false is useful for derived classes .
}

Relation::~Relation()
{
  std::map<string, SimpleVector*>::iterator it;
  for (it = parametersList.begin(); it != parametersList.end(); ++it)
  {
    string alloc = "parameter_for_" + it->first;
    if (isAllocatedIn[alloc]) delete it->second;
  }
  parametersList.clear();
  interaction = NULL;
}

void Relation::initialize()
{
  // Check connection with an Interaction
  if (interaction == NULL)
    RuntimeException::selfThrow("Relation - initialize: no interaction linked with the present relation.");

  if (interaction->getRelationPtr() != this)
    RuntimeException::selfThrow("Relation - initialize: inconsistent pointers links between the present relation and its interaction.");
}

void Relation::setParameters(const std::map<string, SimpleVector*>& newMap)
{
  // copy!!

  map<string, SimpleVector*>::const_iterator it;
  for (it = newMap.begin(); it != newMap.end(); ++it)
  {
    parametersList[it->first] = new SimpleVector(*(it->second));
    string alloc = "parameter_for_" + it->first;
    isAllocatedIn[alloc] = true;
  }
}

void Relation::setParameter(const SimpleVector& newValue, const string id)
{
  parametersList[id] = new SimpleVector(newValue);
  string alloc = "parameter_for_" + id;
  isAllocatedIn[alloc] = true;
}

void Relation::setParameterPtr(SimpleVector *newPtr, const string id)
{
  parametersList[id] = newPtr;
  string alloc = "parameter_for_" + id;
  isAllocatedIn[alloc] = false;
}

void Relation::computeOutput(const double time, const unsigned int)
{
  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0] (at the time).

  if (interaction == NULL)
    RuntimeException::selfThrow("Relation - computeOutput: relation is not connected to any interaction");

  if (computeOutputPtr == NULL)
    RuntimeException::selfThrow("computeOutput() is not linked to a plugin function");

  DynamicalSystemsSet vDS = interaction->getDynamicalSystems();
  BlockVector *xTmp = new BlockVector();
  BlockVector *uTmp = new BlockVector();
  DSIterator it;
  for (it = vDS.begin(); it != vDS.end(); it++)
  {
    // Put x and u of each DS into a block
    // Warning: use copy constructors, no link between pointers
    if (((*it)->getType() != LDS) && ((*it)->getType() != LITIDS))
      RuntimeException::selfThrow("LinearTIR - computeOutput: not yet implemented for DS type " + (*it)->getType());

    xTmp->add((*it)->getX());
    if ((*it)->getUPtr() != NULL)
      uTmp->add(*((*it)->getUPtr())) ;
  }
  unsigned int sizeU = uTmp->size();
  unsigned int sizeX = xTmp->size();

  SiconosVector *y = interaction->getYPtr(0);
  SiconosVector *lambda = interaction->getLambdaPtr(0);
  unsigned int sizeY = y->size();
  SiconosVector* param = parametersList["output"];


  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  SimpleVector * xTmp2 = new SimpleVector(*xTmp);
  SimpleVector * uTmp2 = new SimpleVector(*uTmp);

  SimpleVector * yTmp = new SimpleVector(*y);
  SimpleVector * lambdaTmp = new SimpleVector(*lambda);

  computeOutputPtr(sizeX, &(*xTmp2)(0), time, sizeY, &(*lambdaTmp)(0), sizeU,  &(*uTmp2)(0), &(*yTmp)(0), &(*param)(0));

  // Rebuilt lambda/y from Tmp
  *lambda = *lambdaTmp;
  *y = *yTmp;

  delete xTmp2;
  delete uTmp2;
  delete lambdaTmp;
  delete yTmp;
  delete xTmp;
  delete uTmp;
}

void Relation::computeFreeOutput(const double time, const unsigned int)
{
  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0] (at the time).

  if (interaction == NULL)
    RuntimeException::selfThrow("Relation - computeFreeOutput: relation is not connected to any interaction");

  if (computeOutputPtr == NULL)
    RuntimeException::selfThrow("computeOutput() is not linked to a plugin function");

  DynamicalSystemsSet vDS = interaction->getDynamicalSystems();
  BlockVector *xTmp = new BlockVector();
  BlockVector *uTmp = new BlockVector();
  DSIterator it;

  for (it = vDS.begin(); it != vDS.end(); it++)
  {
    // Put xFree and u of each DS into a block
    // Warning: use copy constructors, no link between pointers
    if (((*it)->getType() != LDS) && ((*it)->getType() != LITIDS))
      RuntimeException::selfThrow("LinearTIR - computeFreeOutput: not yet implemented for DS type " + (*it)->getType());

    xTmp->add((*it)->getXFree());
    if ((*it)->getUPtr() != NULL)
      uTmp->add(*((*it)->getUPtr())) ;
  }

  SiconosVector *yFree = interaction->getYPtr(0);
  // warning : yFree is saved in y !!
  unsigned int sizeU = uTmp->size();
  unsigned int sizeX = xTmp->size();

  SiconosVector *lambda = interaction->getLambdaPtr(0);
  unsigned int sizeY = yFree->size();
  SiconosVector* param = parametersList["output"];
  computeOutputPtr(sizeX, &(*xTmp)(0), time, sizeY, &(*lambda)(0), sizeU,  &(*uTmp)(0), &(*yFree)(0), &(*param)(0));
  delete xTmp;
  delete uTmp;

  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  SimpleVector * xTmp2 = new SimpleVector(*xTmp);
  SimpleVector * uTmp2 = new SimpleVector(*uTmp);

  SimpleVector * yFreeTmp = new SimpleVector(*yFree);
  SimpleVector * lambdaTmp = new SimpleVector(*lambda);

  computeOutputPtr(sizeX, &(*xTmp2)(0), time, sizeY, &(*lambdaTmp)(0), sizeU,  &(*uTmp2)(0), &(*yFreeTmp)(0), &(*param)(0));

  // Rebuilt lambda/y from Tmp
  *lambda = *lambdaTmp;
  *yFree = *yFreeTmp;

  delete xTmp2;
  delete uTmp2;
  delete lambdaTmp;
  delete yFreeTmp;
  delete xTmp;
  delete uTmp;

  // \todo update y, yDot ... depending on the relative degree.
}

void Relation::computeInput(const double time, const unsigned int level)
{
  if (interaction == NULL)
    RuntimeException::selfThrow("Relation - computeInput: relation is not connected to any interaction");

  if (computeInputPtr == NULL)
    RuntimeException::selfThrow("computeInput() is not linked to a plugin function");

  DynamicalSystemsSet vDS = interaction->getDynamicalSystems();
  DSIterator it;
  BlockVector *r = new BlockVector();
  for (it = vDS.begin(); it != vDS.end(); it++)
  {
    // Put r of each DS into a block
    // Warning: use addPtr -> link between pointers
    bool isComp = (*it)->getRPtr()->isBlock();
    if (isComp)
    {
      BlockVector * tmp = static_cast<BlockVector*>((*it)->getRPtr());
      r->addPtr(tmp->getVectorPtr(0));
      r->addPtr(tmp->getVectorPtr(1));
    }
    else
      r->addPtr(static_cast<SimpleVector*>((*it)->getRPtr()));
  }


  SiconosVector *lambda = interaction->getLambdaPtr(level);
  unsigned int sizeY = lambda->size();

  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  SimpleVector * rTmp = new SimpleVector(*r);
  SimpleVector * lambdaTmp = new SimpleVector(*lambda);

  SiconosVector* param = parametersList["input"];
  computeInputPtr(sizeY, &(*lambdaTmp)(0), time, &(*rTmp)(0), &(*param)(0));

  *r = *rTmp;

  delete rTmp;
  delete lambdaTmp;
  delete r;
}

void Relation::setComputeOutputFunction(const string pluginPath, const string functionName)
{
  cShared.setFunction(&computeOutputPtr, pluginPath, functionName);

  string plugin;
  initParameter("output");
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  computeOutputName = plugin + ":" + functionName;
  isOutputPlugged = true;
}

void Relation::setComputeInputFunction(const string pluginPath, const string functionName)
{
  cShared.setFunction(&computeInputPtr, pluginPath, functionName);

  string plugin;
  initParameter("input");
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  computeInputName = plugin + ":" + functionName;
  isInputPlugged = true;
}

void Relation::display() const
{
  cout << "===== Relation display ===== " << endl;
  cout << "- Relation type: " << relationType << endl;
  if (interaction != NULL) cout << "- Interaction id" << interaction->getId() << endl;
  else cout << "- Linked interaction -> NULL" << endl;
  cout << "- Input plug-in name: " << computeInputName << endl;
  cout << "- Output plug-in name: " << computeOutputName << endl;
  cout << "===================================== " << endl;
}

void Relation::saveRelationToXML() const
{
  //to do
  RuntimeException::selfThrow("Relation - saveRelationToXML: not yet implemented for relation of type" + getType());
}
