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
  Relation(inter), LagrangianRelationType("holonom"), isHPlugged(false), hFunctionName("none"),
  parametersList(NULL), isParametersListAllocatedIn(false),
  h0Ptr(NULL), G0Ptr(NULL), h1Ptr(NULL), G10Ptr(NULL), G11Ptr(NULL), h2Ptr(NULL), G20Ptr(NULL), G21Ptr(NULL)
{
  isOutputPlugged = false;
  isInputPlugged  = false;
  relationType = LAGRANGIANRELATION;
  isGPlugged.resize(1, false);
  G.resize(1, NULL);
  isGAllocatedIn.resize(1, false);
  GFunctionName.resize(1, "none");
  setComputeGFunction("DefaultPlugin.so", "G0");
  isGPlugged[0] = false;
  setComputeHFunction("DefaultPlugin.so", "h0");
  isHPlugged = false;
  parametersList = new SimpleVector(1);
  isParametersListAllocatedIn = true;
  (*parametersList)(0) = 1.0;
}

// xml constructor
LagrangianR::LagrangianR(RelationXML* relxml, Interaction* inter):
  Relation(relxml, inter), LagrangianRelationType("holonom"), isHPlugged(true), hFunctionName("none"),
  parametersList(NULL), isParametersListAllocatedIn(false),
  h0Ptr(NULL), G0Ptr(NULL), h1Ptr(NULL), G10Ptr(NULL), G11Ptr(NULL), h2Ptr(NULL), G20Ptr(NULL), G21Ptr(NULL)
{

  if (isOutputPlugged || isInputPlugged)
    RuntimeException::selfThrow("LagrangianR:: xml constructor, input and output plug-in must not be set for this kind of relation");

  relationType = LAGRANGIANRELATION;
  LagrangianRXML * LRxml = static_cast<LagrangianRXML *>(relationxml);
  string plugin;
  // problem type
  LagrangianRelationType = LRxml->getLagrangianType();

  // h plug-in
  if (LRxml->hasH())
  {
    plugin = LRxml->getHPlugin();
    setComputeHFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
  }
  else
  {
    // This case corresponds to LagrangianLinearR, else h is required !!
    setComputeHFunction("DefaultPlugin.so", "h0");
    isHPlugged = false ;
  }


  if (LRxml->hasG())
  {
    if (LagrangianRelationType == "holonom")
    {
      isGPlugged.resize(1, true);
      G.resize(1, NULL);
      isGAllocatedIn.resize(1, false);
      GFunctionName.resize(1, "none");
      if (LRxml->isGPlugin())
      {
        plugin = LRxml->getGPlugin();
        setComputeGFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      }
      else
      {
        G[0] = new SiconosMatrix(LRxml->getGMatrix());
        isGAllocatedIn[0] = true;
        isGPlugged[0] = false;
      }
    }
    else if (LagrangianRelationType == "non-holonom" || LagrangianRelationType == "holonom+lambda")
    {
      isGPlugged.resize(2, true);
      G.resize(2, NULL);
      isGAllocatedIn.resize(2, false);
      GFunctionName.resize(2, "none");
      if (LRxml->isGPlugin())
      {
        plugin = LRxml->getGPlugin(0);
        setComputeGFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin), 0);
        plugin = LRxml->getGPlugin(1);
        setComputeGFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin), 1);
      }
      else
      {
        G[0] = new SiconosMatrix(LRxml->getGMatrix(0));
        G[1] = new SiconosMatrix(LRxml->getGMatrix(1));
        isGAllocatedIn[0] = true;
        isGAllocatedIn[1] = true;
        isGPlugged[0] = false;
        isGPlugged[1] = false;
      }
    }
    else
      cout << "LagrangianR xml constructor: warning, xml loading not yet implemented for G for this kind of problem.\n You should fix it by yourself using set functions" << endl;
  }
  else
  {
    isGPlugged.resize(1, false);
    G.resize(1, NULL);
    isGAllocatedIn.resize(1, false);
    GFunctionName.resize(1, "none");
    setComputeGFunction("DefaultPlugin.so", "G0");
    isGPlugged[0] = false;
  }
  parametersList = new SimpleVector(1);
  isParametersListAllocatedIn = true;
  (*parametersList)(0) = 1.0;
}

// constructor from a set of data, holonom case
LagrangianR::LagrangianR(const string& computeH, const string& computeG, Interaction* inter):
  Relation(inter), LagrangianRelationType("holonom"), isHPlugged(true), hFunctionName("none"),
  parametersList(NULL), isParametersListAllocatedIn(false),
  h0Ptr(NULL), G0Ptr(NULL), h1Ptr(NULL), G10Ptr(NULL), G11Ptr(NULL), h2Ptr(NULL), G20Ptr(NULL), G21Ptr(NULL)
{
  relationType = LAGRANGIANRELATION;
  isOutputPlugged = false;
  isInputPlugged  = false;
  isGPlugged.resize(1, true);
  G.resize(1, NULL);
  isGAllocatedIn.resize(1, false);
  GFunctionName.resize(1, "none");
  parametersList = new SimpleVector(1);
  isParametersListAllocatedIn = true;
  (*parametersList)(0) = 1.0;

  // === Set plug-in for h and G functions
  // h
  setComputeHFunction(cShared.getPluginName(computeH), cShared.getPluginFunctionName(computeH));
  // jacobianQH
  setComputeGFunction(cShared.getPluginName(computeG), cShared.getPluginFunctionName(computeG));
}

// copy constructor (inter is optional)
LagrangianR::LagrangianR(const Relation & newLNLR, Interaction* inter):
  Relation(newLNLR, inter), LagrangianRelationType("none"), isHPlugged(true), hFunctionName("none"),
  parametersList(NULL), isParametersListAllocatedIn(false),
  h0Ptr(NULL), G0Ptr(NULL), h1Ptr(NULL), G10Ptr(NULL), G11Ptr(NULL)
{
  if (relationType !=  LAGRANGIANRELATION && relationType !=  LAGRANGIANLINEARRELATION)
    RuntimeException::selfThrow("LagrangianR:: copy constructor, inconsistent relation types for copy");

  isOutputPlugged = false;
  isInputPlugged  = false;

  const LagrangianR * lnlr = static_cast<const LagrangianR*>(&newLNLR);
  LagrangianRelationType = lnlr->getLagrangianRelationType();
  string plugin;
  // --- h ---
  string hPluginName = lnlr->getHFunctionName();
  setComputeHFunction(cShared.getPluginName(hPluginName), cShared.getPluginFunctionName(hPluginName));
  if (cShared.getPluginName(hPluginName) == "DefaultPlugin.so")
    isHPlugged = false;

  if (lnlr->getParametersListPtr() != NULL)
  {
    parametersList = new SimpleVector(lnlr->getParametersList());
    isParametersListAllocatedIn = true;
  }
  else
  {
    // Default value for parametersList: simpleVector of size 1
    parametersList = new SimpleVector(1);
    (*parametersList)(0) = 1.0;
    isParametersListAllocatedIn = true;
  }

  // --- G ---
  if (LagrangianRelationType == "holonom")
  {
    isGPlugged.resize(1, true);
    G.resize(1, NULL);
    isGAllocatedIn.resize(1, false);
    GFunctionName.resize(1, "none");
    if (lnlr->getGPtr() != NULL)
    {
      G[0] = new SiconosMatrix(lnlr->getG());
      isGAllocatedIn[0] = true;
    }
    plugin = lnlr->getGFunctionName();
    setComputeGFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin), 0);
    if (cShared.getPluginName(plugin) == "DefaultPlugin.so")
      isGPlugged[0] = false;
  }

  else if (LagrangianRelationType == "non-holonom" || LagrangianRelationType == "holonom+lambda")
  {
    isGPlugged.resize(2, true);
    G.resize(2, NULL);
    isGAllocatedIn.resize(2, false);
    GFunctionName.resize(2, "none");
    if (lnlr->getGPtr(0) != NULL)
    {
      G[0] = new SiconosMatrix(lnlr->getG(0));
      isGAllocatedIn[0] = true;
    }
    plugin = lnlr->getGFunctionName(0);
    setComputeGFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin), 0);
    if (cShared.getPluginName(plugin) == "DefaultPlugin.so")
      isGPlugged[0] = false;


    if (lnlr->getGPtr(1) != NULL)
    {
      G[1] = new SiconosMatrix(lnlr->getG(1));
      isGAllocatedIn[1] = true;
      isGPlugged[1] = false;
    }
    plugin = lnlr->getGFunctionName(1);
    setComputeGFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin), 1);
    if (cShared.getPluginName(plugin) == "DefaultPlugin.so")
      isGPlugged[0] = false;
  }
  else
    RuntimeException::selfThrow("LagrangianR:: copy constructor, not yet implemented for problem of type" + LagrangianRelationType);
}

LagrangianR::~LagrangianR()
{
  for (unsigned int i = 0; i < G.size(); i++)
  {
    if (isGAllocatedIn[i]) delete G[i];
    G[i] = NULL;
  }
  if (isParametersListAllocatedIn) delete parametersList;
  parametersList = NULL;
  h0Ptr = NULL;
  h1Ptr = NULL;
  G0Ptr = NULL;
  G10Ptr = NULL;
  G11Ptr = NULL;
  h2Ptr = NULL;
  G20Ptr = NULL;
  G21Ptr = NULL;
}

void LagrangianR::setLagrangianRelationType(const string & type)
{

  LagrangianRelationType = type;
  if (LagrangianRelationType == "holonom")
  {
    isGPlugged.resize(1, false);
    G.resize(1, NULL);
    isGAllocatedIn.resize(1, false);
    GFunctionName.resize(1, "none");
    setComputeGFunction("DefaultPlugin.so", "G0");
    isGPlugged[0] = false;
    setComputeHFunction("DefaultPlugin.so", "h0");
    isHPlugged = false;
  }
  else if (LagrangianRelationType == "non-holonom" || LagrangianRelationType == "holonom+lambda")
  {
    isGPlugged.resize(2, false);
    G.resize(2, NULL);
    isGAllocatedIn.resize(2, false);
    GFunctionName.resize(2, "none");
  }
  else
    RuntimeException::selfThrow("LagrangianR:: copy constructor, not yet implemented for problem of type" + LagrangianRelationType);


}

void LagrangianR::setGVector(const vector<SiconosMatrix*>& newVector)
{
  if (interaction == NULL)
    RuntimeException::selfThrow("LagrangianR:: setGVector, no interaction linked to the current relation");

  unsigned int sizeY = interaction->getNInteraction();
  unsigned int sizeQ = interaction->getSizeOfDS();

  for (unsigned int i = 0; i < G.size(); i++)
  {
    if (isGAllocatedIn[i]) delete G[i];
    G[i] = NULL;
  }
  G.resize(newVector.size(), NULL);
  isGAllocatedIn.resize(G.size(), false);
  if (LagrangianRelationType == "holonom+lambda")
  {
    if (G.size() != 2)
      RuntimeException::selfThrow("LagrangianR - setGVector: wrong size for G vector in holonom+lambda case");
    if (newVector[0]->size(0) != sizeQ || newVector[0]->size(1) != sizeY)
      RuntimeException::selfThrow("LagrangianR - setGVector: inconsistent input matrix size ");
    G[0] = new SiconosMatrix(*(newVector[0]));
    isGAllocatedIn[0] = true;
    if (newVector[1]->size(0) != sizeY || newVector[0]->size(1) != sizeY)
      RuntimeException::selfThrow("LagrangianR - setGVector: inconsistent input matrix size ");
    G[1] = new SiconosMatrix(*(newVector[1]));
    isGAllocatedIn[1] = true;
  }
  else
  {
    for (unsigned int i = 0; i < G.size(); i++)
    {
      if (newVector[i] != NULL)
      {
        if (newVector[i]->size(0) != sizeQ || newVector[i]->size(1) != sizeY)
          RuntimeException::selfThrow("LagrangianR - setGVector: inconsistent input matrix size ");
        G[i] = new SiconosMatrix(*(newVector[i]));
        isGAllocatedIn[i] = true;
      }
    }
  }
}

void LagrangianR::setG(const SiconosMatrix& newValue, const unsigned int & index)
{
  if (interaction == NULL)
    RuntimeException::selfThrow("LagrangianR:: setG, no interaction linked to the current relation");
  unsigned int sizeY = interaction->getNInteraction();
  unsigned int sizeQ = interaction->getSizeOfDS();
  if (LagrangianRelationType == "holonom+lambda" && index == 1)
  {
    if (newValue.size(1) != sizeY || newValue.size(0) != sizeY)
      RuntimeException::selfThrow("LagrangianR - setG: inconsistent input matrix size ");
  }
  else
  {
    if (newValue.size(1) != sizeQ || newValue.size(0) != sizeY)
      RuntimeException::selfThrow("LagrangianR - setG: inconsistent input matrix size ");
  }

  if (index >= G.size())
    RuntimeException::selfThrow("LagrangianR:: setG(mat,index), index out of range. Use rather setGVector?");

  if (G[index] == NULL)
  {
    G[index] = new SiconosMatrix(newValue);
    isGAllocatedIn[index] = true;
  }
  else
    *(G[index]) = newValue;
  isGPlugged[index] = false;
}

void LagrangianR::setGPtr(SiconosMatrix *newPtr, const unsigned int & index)
{
  if (interaction == NULL)
    RuntimeException::selfThrow("LagrangianR:: setGPtr, no interaction linked to the current relation");
  unsigned int sizeY = interaction->getNInteraction();
  unsigned int sizeQ = interaction->getSizeOfDS();
  if (LagrangianRelationType == "holonom+lambda" && index == 1)
  {
    if (newPtr->size(1) != sizeY || newPtr->size(0) != sizeY)
      RuntimeException::selfThrow("LagrangianR - setGPtr: inconsistent input matrix size ");
  }
  else
  {
    if (newPtr->size(1) != sizeQ || newPtr->size(0) != sizeY)
      RuntimeException::selfThrow("LagrangianR - setGPtr: inconsistent input matrix size ");
  }

  if (index >= G.size())
    RuntimeException::selfThrow("LagrangianR:: setGPtr(mat,index), index out of range. Use rather setGVector?");

  if (isGAllocatedIn[index]) delete G[index];
  G[index] = newPtr;
  isGAllocatedIn[index] = false;
  isGPlugged[index] = false;
}

void LagrangianR::setComputeHFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianR::setComputeHFunction\n");

  if (LagrangianRelationType == "holonom")
    cShared.setFunction(&h0Ptr, pluginPath, functionName);
  else if (LagrangianRelationType == "non-holonom")
    cShared.setFunction(&h1Ptr, pluginPath, functionName);
  else if (LagrangianRelationType == "holonom+lambda")
    cShared.setFunction(&h2Ptr, pluginPath, functionName);
  else
    RuntimeException::selfThrow("LagrangianR:: setComputeHFunction,  not yet implemented for this type of constraints");

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  hFunctionName = plugin + ":" + functionName;
  isHPlugged = true;
  OUT("LagrangianR::setComputeHFunction\n");
}

void LagrangianR::setComputeGFunction(const string& pluginPath, const string& functionName, const unsigned int & index)
{
  IN("LagrangianR::setComputeGFunction\n");

  if (index >= G.size())
    RuntimeException::selfThrow("LagrangianR:: setComputeGFunction, index out of range. Use rather setGVector?");

  if (interaction != NULL)
  {
    unsigned int sizeY = interaction->getNInteraction();
    unsigned int sizeQ = interaction->getSizeOfDS();

    if (G[index] == NULL)
    {
      if (LagrangianRelationType == "holonom+lambda" && index == 1)
        G[index] = new SiconosMatrix(sizeY, sizeY);
      else
        G[index] = new SiconosMatrix(sizeY, sizeQ);
      isGAllocatedIn[index] = true;
      G[index]->zero();
    }
  }
  else
    cout << "LagrangianR::setComputeGFunction: warning, no interaction linked with relation, can not set size of G matrix" << endl;

  if (LagrangianRelationType == "holonom")
    cShared.setFunction(&G0Ptr, pluginPath, functionName);
  else if (LagrangianRelationType == "non-holonom")
  {
    if (index == 0)
      cShared.setFunction(&G10Ptr, pluginPath, functionName);
    else if (index == 1)
      cShared.setFunction(&G11Ptr, pluginPath, functionName);
    else
      RuntimeException::selfThrow("LagrangianR:: setComputeGFunction, index out of range");
  }
  else if (LagrangianRelationType == "holonom+lambda")
  {
    if (index == 0)
      cShared.setFunction(&G20Ptr, pluginPath, functionName);
    else if (index == 1)
      cShared.setFunction(&G21Ptr, pluginPath, functionName);
    else
      RuntimeException::selfThrow("LagrangianR:: setComputeGFunction, index out of range");
  }
  else
    RuntimeException::selfThrow("LagrangianR:: setComputeGFunction,  not yet implemented for this type of constraints");

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  GFunctionName[index] = plugin + ":" + functionName;
  isGPlugged[index] = true;
  OUT("LagrangianR::setComputeGFunction\n");
}

void LagrangianR::setParametersList(const SimpleVector& newValue)
{
  if (isParametersListAllocatedIn) delete parametersList;
  parametersList = new SimpleVector(newValue);
  isParametersListAllocatedIn = true;
}

void LagrangianR::setParametersListPtr(SimpleVector *newPtr)
{
  if (isParametersListAllocatedIn) delete parametersList;
  parametersList = newPtr;
  isParametersListAllocatedIn = false;
}

void LagrangianR::computeH(const double& time)
{
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
  SimpleVector *lambda = interaction->getLambdaPtr(0);

  if (LagrangianRelationType == "holonom")
  {
    if (h0Ptr == NULL)
      RuntimeException::selfThrow("LagrangianR:computeH() is not linked to a plugin function");
    h0Ptr(&sizeQ, &(*qTmp)(0) , &sizeY, &(*y)(0), &(*parametersList)(0));
  }
  else if (LagrangianRelationType == "non-holonom")
  {
    if (h1Ptr == NULL)
      RuntimeException::selfThrow("LagrangianR:computeH() is not linked to a plugin function");
    h1Ptr(&sizeQ, &(*qTmp)(0), &time, &sizeY,  &(*y)(0), &(*parametersList)(0));
  }
  else if (LagrangianRelationType == "holonom+lambda")
  {
    if (h2Ptr == NULL)
      RuntimeException::selfThrow("LagrangianR:computeH() is not linked to a plugin function");
    h2Ptr(&sizeQ, &(*qTmp)(0), &(*lambda)(0), &sizeY,  &(*y)(0), &(*parametersList)(0));
  }
  else
    RuntimeException::selfThrow("LagrangianR::computeH,  not yet implemented for this type of constraints");
}

void LagrangianR::computeG(const double & time, const unsigned int & index)
{

  if (index >= G.size())
    RuntimeException::selfThrow("LagrangianR:: computeG(), index out of range. Use rather setGVector?");
  if (interaction == NULL)
    RuntimeException::selfThrow("LagrangianR:computeG(): no interaction linked with relation");

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
      RuntimeException::selfThrow("LagrangianLinearR::computeG not yet implemented for dynamical system of type: " + vDS[i]->getType());

    // convert vDS systems into LagrangianDS and put them in vLDS
    vLDS.push_back(static_cast<LagrangianDS*>(vDS[i]));

    // Put q of each DS into a composite
    // Warning: use copy constructors (add function), no link between pointers
    qTmp->add(vLDS[i]->getQ());
  }

  // get vector lambda of the current interaction
  SimpleVector *lambda = interaction->getLambdaPtr(0);

  if (LagrangianRelationType == "holonom")
  {
    SiconosMatrix * Gtmp = G[0];
    if (G0Ptr == NULL)
      RuntimeException::selfThrow("computeG() is not linked to a plugin function");
    G0Ptr(&sizeQ, &(*qTmp)(0), &sizeY, &(*Gtmp)(0, 0), &(*parametersList)(0));
  }
  else if (LagrangianRelationType == "non-holonom")
  {
    SiconosMatrix * Gtmp = G[index];
    if (index == 0)
    {
      if (G10Ptr == NULL)
        RuntimeException::selfThrow("computeG() is not linked to a plugin function");
      G10Ptr(&sizeQ, &(*qTmp)(0), &time, &sizeY, &(*Gtmp)(0, 0), &(*parametersList)(0));
    }
    else if (index == 1)
    {
      if (G11Ptr == NULL)
        RuntimeException::selfThrow("computeG() is not linked to a plugin function");
      G11Ptr(&sizeQ, &(*qTmp)(0), &time, &sizeY, &(*Gtmp)(0, 0), &(*parametersList)(0));
    }
    else
      RuntimeException::selfThrow("LagrangianR::computeG, index out of range");
  }
  else if (LagrangianRelationType == "holonom+lambda")
  {
    SiconosMatrix * Gtmp = G[index];
    if (index == 0)
    {
      if (G20Ptr == NULL)
        RuntimeException::selfThrow("computeG() is not linked to a plugin function");
      G20Ptr(&sizeQ, &(*qTmp)(0), &(*lambda)(0), &sizeY, &(*Gtmp)(0, 0), &(*parametersList)(0));
    }
    else if (index == 1)
    {
      if (G21Ptr == NULL)
        RuntimeException::selfThrow("computeG() is not linked to a plugin function");
      G21Ptr(&sizeQ, &(*qTmp)(0), &(*lambda)(0), &sizeY, &(*Gtmp)(0, 0), &(*parametersList)(0));
    }
    else
      RuntimeException::selfThrow("LagrangianR::computeG, index out of range");
  }
  else
    RuntimeException::selfThrow("LagrangianR::computeG,  not yet implemented for this type of constraints");
}

void LagrangianR::computeOutput(const double& time)
{
  if (interaction == NULL)
    RuntimeException::selfThrow("LagrangianR::computeOutput, no interaction linked with this relation");

  unsigned int sizeY = interaction->getNInteraction();
  unsigned int sizeQ = interaction->getSizeOfDS();

  // Get the DS concerned by the interaction of this relation
  vector<DynamicalSystem*> vDS = interaction->getDynamicalSystems();
  vector<LagrangianDS*> vLDS;

  unsigned int size = vDS.size(), i;
  CompositeVector *qTmp = new CompositeVector();
  CompositeVector *velocityTmp = new CompositeVector();
  for (i = 0; i < size; i++)
  {
    // check dynamical system type
    if (vDS[i]->getType() != LTIDS && vDS[i]->getType() != LNLDS)
      RuntimeException::selfThrow("LagrangianLinearR::computeOutput not yet implemented for dynamical system of type: " + vDS[i]->getType());

    // convert vDS systems into LagrangianDS and put them in vLDS
    vLDS.push_back(static_cast<LagrangianDS*>(vDS[i]));

    // Put q and velocity of each DS into a composite
    // Warning: use copy constructors (add function), no link between pointers
    qTmp->add(vLDS[i]->getQ());
    velocityTmp->add(vLDS[i]->getVelocity());
  }

  // get y and yDot of the interaction
  SimpleVector *y = interaction->getYPtr(0);
  SimpleVector *yDot = interaction->getYPtr(1);
  SimpleVector *lambda = interaction->getLambdaPtr(0);
  SimpleVector *lambdaDot = interaction->getLambdaPtr(1);

  // y = h(...) and yDot = GqDot
  if (LagrangianRelationType == "holonom")
  {
    if (h0Ptr == NULL)
      RuntimeException::selfThrow("LagrangianR:computeOutput() h0 is not linked to a plugin function");
    h0Ptr(&sizeQ, &(*qTmp)(0) , &sizeY, &(*y)(0), &(*parametersList)(0));
    SiconosMatrix * Gtmp = G[0];
    if (G0Ptr == NULL)
      RuntimeException::selfThrow("computeG() is not linked to a plugin function");
    G0Ptr(&sizeQ, &(*qTmp)(0), &sizeY, &(*Gtmp)(0, 0), &(*parametersList)(0));
    *yDot = *Gtmp * *velocityTmp;
  }
  // y = h(...) and yDot = G10qDot + G11
  else if (LagrangianRelationType == "non-holonom")
  {
    if (h1Ptr == NULL)
      RuntimeException::selfThrow("LagrangianR:computeOutput() h1 is not linked to a plugin function");
    h1Ptr(&sizeQ, &(*qTmp)(0), &time, &sizeY,  &(*y)(0), &(*parametersList)(0));
    SiconosMatrix * G0tmp = G[0];
    SiconosMatrix * G1tmp = G[1];

    if (G10Ptr == NULL)
      RuntimeException::selfThrow("computeG() is not linked to a plugin function");
    G10Ptr(&sizeQ, &(*qTmp)(0), &time, &sizeY, &(*G0tmp)(0, 0), &(*parametersList)(0));

    if (G11Ptr == NULL)
      RuntimeException::selfThrow("computeG() is not linked to a plugin function");
    G11Ptr(&sizeQ, &(*qTmp)(0), &time, &sizeY, &(*G1tmp)(0, 0), &(*parametersList)(0));

    // warning: G1 is a matrix
    SimpleVector * Id = new SimpleVector(1);
    (*Id)(0) = 1.0;
    *yDot = *G0tmp * *velocityTmp + *G1tmp * *Id;
    delete Id;
  }
  else if (LagrangianRelationType == "holonom+lambda")
  {
    if (h1Ptr == NULL)
      RuntimeException::selfThrow("LagrangianR:computeOutput() h2 is not linked to a plugin function");
    h2Ptr(&sizeQ, &(*qTmp)(0), &(*lambda)(0), &sizeY,  &(*y)(0), &(*parametersList)(0));
    SiconosMatrix * G0tmp = G[0];
    SiconosMatrix * G1tmp = G[1];

    if (G10Ptr == NULL)
      RuntimeException::selfThrow("computeG() is not linked to a plugin function");
    G20Ptr(&sizeQ, &(*qTmp)(0), &(*lambda)(0), &sizeY, &(*G0tmp)(0, 0), &(*parametersList)(0));

    if (G11Ptr == NULL)
      RuntimeException::selfThrow("computeG() is not linked to a plugin function");
    G21Ptr(&sizeQ, &(*qTmp)(0), &(*lambda)(0), &sizeY, &(*G1tmp)(0, 0), &(*parametersList)(0));

    *yDot = *G0tmp * *velocityTmp + *G1tmp * *lambdaDot;
  }
  else
    RuntimeException::selfThrow("LagrangianR::computeOutput(),  not yet implemented for this type of constraints");

  // free memory
  delete qTmp;
  delete velocityTmp;
}
void LagrangianR::computeFreeOutput(const double& time)
{
  // warning: the only difference with computeOutput is that we get qFREE and velocityFREE instead of q and velocity.
  if (interaction == NULL)
    RuntimeException::selfThrow("LagrangianR::computeOutput, no interaction linked with this relation");

  unsigned int sizeY = interaction->getNInteraction();
  unsigned int sizeQ = interaction->getSizeOfDS();

  // Get the DS concerned by the interaction of this relation
  vector<DynamicalSystem*> vDS = interaction->getDynamicalSystems();
  vector<LagrangianDS*> vLDS;

  unsigned int size = vDS.size(), i;
  CompositeVector *qTmp = new CompositeVector();
  CompositeVector *velocityTmp = new CompositeVector();
  for (i = 0; i < size; i++)
  {
    // check dynamical system type
    if (vDS[i]->getType() != LTIDS && vDS[i]->getType() != LNLDS)
      RuntimeException::selfThrow("LagrangianLinearR::computeOutput not yet implemented for dynamical system of type: " + vDS[i]->getType());

    // convert vDS systems into LagrangianDS and put them in vLDS
    vLDS.push_back(static_cast<LagrangianDS*>(vDS[i]));

    // Put q and velocity of each DS into a composite
    // Warning: use copy constructors (add function), no link between pointers
    qTmp->add(vLDS[i]->getQFree());
    velocityTmp->add(vLDS[i]->getVelocityFree());
  }

  // get y and yDot of the interaction
  SimpleVector *y = interaction->getYPtr(0);
  SimpleVector *yDot = interaction->getYPtr(1);
  SimpleVector *lambda = interaction->getLambdaPtr(0);
  SimpleVector *lambdaDot = interaction->getLambdaPtr(1);

  // y = h(...) and yDot = GqDot
  if (LagrangianRelationType == "holonom")
  {
    if (h0Ptr == NULL)
      RuntimeException::selfThrow("LagrangianR:computeOutput() h0 is not linked to a plugin function");
    h0Ptr(&sizeQ, &(*qTmp)(0) , &sizeY, &(*y)(0), &(*parametersList)(0));
    SiconosMatrix * Gtmp = G[0];
    if (G0Ptr == NULL)
      RuntimeException::selfThrow("computeG() is not linked to a plugin function");
    G0Ptr(&sizeQ, &(*qTmp)(0), &sizeY, &(*Gtmp)(0, 0), &(*parametersList)(0));

    *yDot = *Gtmp * *velocityTmp;
  }

  // y = h(...) and yDot = G10qDot + G11
  else if (LagrangianRelationType == "non-holonom")
  {
    if (h1Ptr == NULL)
      RuntimeException::selfThrow("LagrangianR:computeOutput() h1 is not linked to a plugin function");
    h1Ptr(&sizeQ, &(*qTmp)(0), &time, &sizeY,  &(*y)(0), &(*parametersList)(0));
    SiconosMatrix * G0tmp = G[0];
    SiconosMatrix * G1tmp = G[0];

    if (G10Ptr == NULL)
      RuntimeException::selfThrow("computeG() is not linked to a plugin function");
    G10Ptr(&sizeQ, &(*qTmp)(0), &time, &sizeY, &(*G0tmp)(0, 0), &(*parametersList)(0));

    if (G11Ptr == NULL)
      RuntimeException::selfThrow("computeG() is not linked to a plugin function");
    G11Ptr(&sizeQ, &(*qTmp)(0), &time, &sizeY, &(*G1tmp)(0, 0), &(*parametersList)(0));

    // warning: G1 is a matrix
    SimpleVector * Id = new SimpleVector(1);
    (*Id)(0) = 1.0;
    *yDot = *G0tmp * *velocityTmp + *G1tmp * *Id;
  }
  else if (LagrangianRelationType == "holonom+lambda")
  {
    if (h1Ptr == NULL)
      RuntimeException::selfThrow("LagrangianR:computeOutput() h2 is not linked to a plugin function");
    h2Ptr(&sizeQ, &(*qTmp)(0), &(*lambda)(0), &sizeY,  &(*y)(0), &(*parametersList)(0));
    SiconosMatrix * G0tmp = G[0];
    SiconosMatrix * G1tmp = G[1];

    if (G10Ptr == NULL)
      RuntimeException::selfThrow("computeG() is not linked to a plugin function");
    G20Ptr(&sizeQ, &(*qTmp)(0), &(*lambda)(0), &sizeY, &(*G0tmp)(0, 0), &(*parametersList)(0));

    if (G11Ptr == NULL)
      RuntimeException::selfThrow("computeG() is not linked to a plugin function");
    G21Ptr(&sizeQ, &(*qTmp)(0), &(*lambda)(0), &sizeY, &(*G1tmp)(0, 0), &(*parametersList)(0));

    *yDot = *G0tmp * *velocityTmp + *G1tmp * *lambdaDot;
  }
  else
    RuntimeException::selfThrow("LagrangianR::computeOutput(),  not yet implemented for this type of constraints");


  // free memory
  delete qTmp;
  delete velocityTmp;
}

void LagrangianR::computeInput(const double& time)
{
  IN("LagrangianLinearR::computeInput\n");
  if (interaction == NULL)
    RuntimeException::selfThrow("LagrangianLinearR::computeInput, no interaction linked with this relation");

  // Get the DS concerned by the interaction of this relation
  vector<DynamicalSystem*> vDS = interaction->getDynamicalSystems();
  vector<LagrangianDS*> vLDS;
  unsigned int numberDS = vDS.size(), i;
  vLDS.resize(numberDS);

  CompositeVector *p = new CompositeVector();
  string typeDS;

  for (i = 0; i < numberDS; i++)
  {
    // check dynamical system type
    typeDS = vDS[i] ->getType();
    if (typeDS != LTIDS && typeDS != LNLDS)
      RuntimeException::selfThrow("LagrangianLinearR::computeInput not yet implemented for this type of dynamical system " + typeDS);

    // convert vDS systems into LagrangianDS and put them in vLDS
    vLDS[i] = static_cast<LagrangianDS*>(vDS[i]);

    // Put p of each DS into a composite
    // Warning: use addPtr -> link between pointers
    p->addPtr(vLDS[i]->getPPtr());
  }

  // get lambda of the concerned interaction
  SimpleVector *lambda = interaction->getLambdaPtr(1);

  if (LagrangianRelationType == "holonom")
  {
    computeG(time, 0);
    *p += matTransVecMult(*(G[0]), *lambda);

  }
  else
    RuntimeException::selfThrow("LagrangianR::computeInput,  not yet implemented for constraints of type" + LagrangianRelationType);
}

void LagrangianR::getGBlockDS(DynamicalSystem * ds, SiconosMatrix& Block, const unsigned & index) const
{
  unsigned int k = 0;
  vector<DynamicalSystem*> vDS = interaction ->getDynamicalSystems();
  vector<DynamicalSystem*>::iterator itDS;
  itDS = vDS.begin();

  // look for ds
  while (*itDS != ds && itDS != vDS.end())
  {
    k += (*itDS)->getN() / 2;
    itDS++;
  }
  // check dimension
  if ((*itDS)->getN() / 2 != Block.size(1))
    RuntimeException::selfThrow("LagrangianR - getGBlockDS: inconsistent sizes between GBlock and DS");

  // get block
  unsigned int l = k + (*itDS)->getN() / 2 - 1;
  vector<unsigned int> index_list(4);
  index_list[0] = 0;
  index_list[1] = G[index]->size(0) - 1;
  index_list[2] = k;
  index_list[3] = l;
  G[index]->getBlock(index_list, Block);
}

void LagrangianR::getGBlockDS(const int& DSNumber, SiconosMatrix& Block, const unsigned & index) const
{
  unsigned int k = 0;

  vector<DynamicalSystem*> vDS = interaction ->getDynamicalSystems();

  vector<DynamicalSystem*>::iterator itDS;
  itDS = vDS.begin();

  // look for DS number DSNumber ...
  while ((*itDS)->getNumber() != DSNumber && itDS != vDS.end())
  {
    k += (*itDS)->getN() / 2;
    itDS++;
  }

  // check dimension
  if ((*itDS)->getN() / 2 != Block.size(1))
    RuntimeException::selfThrow("LagrangianR - getGlockDS: inconsistent sizes between GBlock and DS");

  // get block
  unsigned int l = k + (*itDS)->getN() / 2 - 1;
  vector<unsigned int> index_list(4);
  index_list[0] = 0;
  index_list[1] = G[index]->size(0) - 1;
  index_list[2] = k;
  index_list[3] = l;
  G[index]->getBlock(index_list, Block);
}

void LagrangianR::saveRelationToXML() const
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
  if (isHPlugged) cout << "h is plugged on function: " << hFunctionName << endl;
  for (unsigned int i = 0; i < G.size(); i++)
  {
    if (isGPlugged[i]) cout << "G" << i << " is plugged on function: " << GFunctionName[i] << endl;
    cout << " G" << i << ": " << endl;
    if (G[i] != NULL)
      G[i] -> display();
    else
      cout << " -> NULL " << endl;
  }
  cout << "===================================== " << endl;
  OUT("LagrangianR::display\n");
}
