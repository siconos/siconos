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
LagrangianR::LagrangianR(Interaction* inter):
  Relation(inter), isHPlugin(true), hFunctionName("none"),
  computeHPtr(NULL), computeJacobianQHPtr(NULL)
{
  relationType = LAGRANGIANRELATION;
  jacobianH.resize(1, NULL);
  isJacobianHAllocatedIn.resize(1, false);
  isJacobianHPlugin.resize(1, false);
  jacobianHFunctionName.resize(1, "none");
  setComputeHFunction("DefaultPlugin.so", "h");
}


// xml constructor
LagrangianR::LagrangianR(RelationXML* relxml, Interaction* inter):
  Relation(relxml, inter), isHPlugin(true), hFunctionName("none"),
  computeHPtr(NULL), computeJacobianQHPtr(NULL)
{
  relationType = LAGRANGIANRELATION;
  jacobianH.resize(1, NULL);
  isJacobianHAllocatedIn.resize(1, false);
  isJacobianHPlugin.resize(1, false);
  jacobianHFunctionName.resize(1, "none");
  setComputeHFunction("DefaultPlugin.so", "h");
  // \todo implement xml loading for LagrangianR
}

// constructor from a set of data
LagrangianR::LagrangianR(const string& computeH, const string& computeJacobianQH, Interaction* inter):
  Relation(inter), isHPlugin(true), hFunctionName("none"), computeHPtr(NULL), computeJacobianQHPtr(NULL)
{
  relationType = LAGRANGIANRELATION;
  jacobianH.resize(1, NULL);
  isJacobianHAllocatedIn.resize(1, false);
  isJacobianHPlugin.resize(1, false);
  jacobianHFunctionName.resize(1, "none");

  // === Set plug-in for h and jacobianH functions
  // h
  setComputeHFunction(cShared.getPluginName(computeH), cShared.getPluginFunctionName(computeH));
  // jacobianQH
  setComputeJacobianQHFunction(cShared.getPluginName(computeJacobianQH), cShared.getPluginFunctionName(computeJacobianQH));
}

// copy constructor (inter is optional)
LagrangianR::LagrangianR(const Relation & newLNLR, Interaction* inter):
  Relation(newLNLR, inter), isHPlugin(true), hFunctionName("none"), computeHPtr(NULL), computeJacobianQHPtr(NULL)
{
  if (relationType !=  LAGRANGIANRELATION && relationType !=  LAGRANGIANLINEARRELATION)
    RuntimeException::selfThrow("LagrangianR:: copy constructor, inconsistent relation types for copy");

  jacobianH.resize(1, NULL);
  isJacobianHAllocatedIn.resize(1, false);
  isJacobianHPlugin.resize(1, false);
  jacobianHFunctionName.resize(1, "none");

  const LagrangianR * lnlr = static_cast<const LagrangianR*>(&newLNLR);

  // h plug-in
  if (lnlr->isHPlugged())
  {
    string hPluginName = lnlr->getHFunctionName();
    setComputeHFunction(cShared.getPluginName(hPluginName), cShared.getPluginFunctionName(hPluginName));
  }
  else
    setComputeHFunction("DefaultPlugin.so", "h");

  // jacobian H
  if (lnlr->getJacobianQHPtr() != NULL)
  {
    jacobianH[0] = new SiconosMatrix(lnlr->getJacobianQH()) ;
    isJacobianHAllocatedIn[0] = true;
  }
  if ((lnlr->isJacobianHPlugged())[0])
  {
    isJacobianHPlugin[0] = true;
    string jacoHPluginName = (lnlr->getJacobianHFunctionName())[0];
    setComputeJacobianQHFunction(cShared.getPluginName(jacoHPluginName), cShared.getPluginFunctionName(jacoHPluginName));
  }
}

LagrangianR::~LagrangianR()
{
  for (unsigned int i = 0; i < jacobianH.size(); i++)
  {
    if (isJacobianHAllocatedIn[i]) delete jacobianH[i];
    jacobianH[i] = NULL;
  }
  computeHPtr = NULL;
  computeJacobianQHPtr = NULL;
}

void LagrangianR::setJacobianQH(const SiconosMatrix& newValue)
{
  unsigned int sizeY = interaction->getNInteraction();
  unsigned int sizeQ = interaction->getSizeOfDS();
  if (newValue.size(0) != sizeQ || newValue.size(1) != sizeY)
    RuntimeException::selfThrow("LagrangianR - setJacobianQH: inconsistent input matrix size ");

  if (jacobianH[0] == NULL)
  {
    jacobianH[0] = new SiconosMatrix(newValue);
    isJacobianHAllocatedIn[0] = true;
  }
  else
    *(jacobianH[0]) = newValue;
  isJacobianHPlugin[0] = false;
}

void LagrangianR::setJacobianQHPtr(SiconosMatrix *newPtr)
{
  unsigned int sizeY = interaction->getNInteraction();
  unsigned int sizeQ = interaction->getSizeOfDS();
  if (newPtr->size(0) != sizeQ || newPtr->size(1) != sizeY)
    RuntimeException::selfThrow("LagrangianR - setJacobianQHPtr: inconsistent input matrix size ");

  if (isJacobianHAllocatedIn[0]) delete jacobianH[0];
  jacobianH[0] = newPtr;
  isJacobianHAllocatedIn[0] = false;
  isJacobianHPlugin[0] = false;
}

void LagrangianR::setComputeHFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianR::setComputeHFunction\n");

  cShared.setFunction(&computeHPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  hFunctionName = plugin + ":" + functionName;
  isHPlugin = true;
  OUT("LagrangianR::setComputeHFunction\n");
}

void LagrangianR::setComputeJacobianQHFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianR::setComputeJacobianQHFunction\n");

  if (interaction != NULL)
  {
    unsigned int sizeY = interaction->getNInteraction();
    unsigned int sizeQ = interaction->getSizeOfDS();
    if (jacobianH[0] == NULL)
    {
      jacobianH[0] = new SiconosMatrix(sizeQ, sizeY);
      isJacobianHAllocatedIn[0] = true;
    }
  }
  else
    cout << "LagrangianR::setComputeJacobianQHFunction: warning, no interaction linked with relation, can not set size of H matrix" << endl;

  cShared.setFunction(&computeJacobianQHPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  jacobianHFunctionName[0] = plugin + ":" + functionName;
  isJacobianHPlugin[0] = true;
  OUT("LagrangianR::setComputeJacobianQHFunction\n");
}

void LagrangianR::computeH(const double& time)
{
  if (computeHPtr == NULL)
    RuntimeException::selfThrow("LagrangianR:computeH() is not linked to a plugin function");

  if (interaction == NULL)
    RuntimeException::selfThrow("LagrangianR:computeH(): no interaction linked with relation");

  unsigned int sizeY = interaction->getNInteraction();
  unsigned int sizeQ = interaction->getSizeOfDS();

  // Get the DS concerned by the interaction of this relation
  vector<DynamicalSystem*> vDS = interaction->getDynamicalSystems();
  vector<LagrangianDS*> vLDS;

  unsigned int size = vDS.size(), i;
  CompositeVector *qTmp = new CompositeVector();
  for (i = 0; i < size; i++)
  {
    // check dynamical system type
    if (vDS[i]->getType() != LTIDS && vDS[i]->getType() != LNLDS)
      RuntimeException::selfThrow("LagrangianLinearR::computeOutput not yet implemented for dynamical system of type: " + vDS[i]->getType());

    // convert vDS systems into LagrangianDS and put them in vLDS
    vLDS.push_back(static_cast<LagrangianDS*>(vDS[i]));

    // Put q of each DS into a composite
    // Warning: use copy constructors (add function), no link between pointers
    qTmp->add(vLDS[i]->getQ());
  }

  // get vector y of the current interaction
  SimpleVector *y = interaction->getYPtr(0);

  computeHPtr(&sizeQ, &time, &sizeY, &(*qTmp)(0), &(*y)(0));
}

void LagrangianR::computeJacobianQH(const double & time)
{
  if (computeJacobianQHPtr == NULL)
    RuntimeException::selfThrow("computeJacobianQH() is not linked to a plugin function");

  if (interaction == NULL)
    RuntimeException::selfThrow("LagrangianR:computeJacobianH(): no interaction linked with relation");

  unsigned int sizeY = interaction->getNInteraction();
  unsigned int sizeQ = interaction->getSizeOfDS();

  // Get the DS concerned by the interaction of this relation
  vector<DynamicalSystem*> vDS = interaction->getDynamicalSystems();
  vector<LagrangianDS*> vLDS;

  unsigned int size = vDS.size(), i;
  CompositeVector *qTmp = new CompositeVector();
  for (i = 0; i < size; i++)
  {
    // check dynamical system type
    if (vDS[i]->getType() != LTIDS && vDS[i]->getType() != LNLDS)
      RuntimeException::selfThrow("LagrangianLinearR::computeOutput not yet implemented for dynamical system of type: " + vDS[i]->getType());

    // convert vDS systems into LagrangianDS and put them in vLDS
    vLDS.push_back(static_cast<LagrangianDS*>(vDS[i]));

    // Put q of each DS into a composite
    // Warning: use copy constructors (add function), no link between pointers
    qTmp->add(vLDS[i]->getQ());
  }

  SiconosMatrix * jacob = jacobianH[0];

  computeJacobianQHPtr(&sizeQ, &time, &sizeY, &(*qTmp)(0), &(*jacob)(0, 0));

}

void LagrangianR::computeOutput(const double&)
{
  cout << "Not yet implemented" << endl;
}
void LagrangianR::computeFreeOutput(const double&)
{
  cout << "Not yet implemented" << endl;
}

void LagrangianR::computeInput(const double&)
{
  cout << "Not yet implemented" << endl;
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

void LagrangianR::display() const
{
  IN("LagrangianR::display\n");
  cout << "===== Lagrangian Relation display ===== " << endl;
  cout << "- Relation type: " << relationType << endl;
  if (interaction != NULL) cout << "- Interaction id" << interaction->getId() << endl;
  else cout << "- Linked interaction -> NULL" << endl;
  if (isHPlugin) cout << "h is plugged on function: " << hFunctionName << endl;
  if (isJacobianHPlugin[0]) cout << "jacobianQH is plugged on function: " << jacobianHFunctionName[0] << endl;
  cout << " jacobianQH: " << endl;
  if (jacobianH[0] != NULL)
    jacobianH[0] -> display();
  else
    cout << " -> NULL " << endl;
  cout << "===================================== " << endl;
  OUT("LagrangianR::display\n");
}
