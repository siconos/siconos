/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
#include "LagrangianDS.h"
using namespace std;

// Private function to set linked with members of Dynamical top class
void LagrangianDS::connectToDS()
{
  n = 2 * ndof;

  // x and related vectors
  x = new BlockVector(q, velocity);
  isXAllocatedIn[1] = true;
  x0 = new BlockVector(q0, velocity0);
  isXAllocatedIn[0] = true;
  xFree = new BlockVector(qFree, velocityFree);
  isXAllocatedIn[5] = true;

  // r
  r = new BlockVector(NULL, p);
  isRAllocatedIn[0] = true;
  r->zero();

  // xDot/VectorField
  xDot = new BlockVector(velocity, NULL);
  isXAllocatedIn[3] = true;

  // jacobianX
  workMatrix["zero-matrix"] = new SimpleMatrix(ndof, ndof);
  workMatrix["Id-matrix"] = new SimpleMatrix(ndof, ndof);
  workMatrix["Id-matrix"]->eye();
  jacobianX = new BlockMatrix(workMatrix["zero-matrix"], workMatrix["Id-matrix"], jacobianQFInt, jacobianVelocityFInt);

  // no setComputeVectorField or setJacobianXFunction, since computeVectorField and computeJacobianX are overloaded in the present class
}

// --- Constructor from an xml file ---
LagrangianDS::LagrangianDS(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* newNsds):
  DynamicalSystem(), ndof(0), q(NULL), q0(NULL), qFree(NULL), qMemory(NULL), velocity(NULL), velocity0(NULL),
  velocityFree(NULL), velocityMemory(NULL), p(NULL), mass(NULL), fInt(NULL), fExt(NULL),
  NNL(NULL), jacobianQFInt(NULL), jacobianVelocityFInt(NULL), jacobianQNNL(NULL), jacobianVelocityNNL(NULL),
  isPAllocatedIn(true), isMassAllocatedIn(true),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL), computeJacobianQFIntPtr(NULL),
  computeJacobianVelocityFIntPtr(NULL), computeJacobianQNNLPtr(NULL), computeJacobianVelocityNNLPtr(NULL), isDSup(false)
{
  if (dsXML == NULL)
    RuntimeException::selfThrow("LagrangianDS::LagrangianDS - DynamicalSystemXML paramater must not be NULL");

  if (newNsds != NULL) nsds = newNsds;

  // --- Dynamical System top-class members ---
  DSType = LNLDS;
  dsxml = dsXML;
  number = dsxml->getNumber();
  if (dsxml->hasId() == true) id = dsxml->getId();
  if (dsxml->hasStepsInMemory() == true) stepsInMemory = dsxml->getStepsInMemory();

  // -- Boundary conditions --
  fillBoundaryConditionsFromXml();
  // -- DS input-output --
  fillDsioFromXml();

  // --- Lagrangian class members ---
  // -- Lagrangian  xml object --
  LagrangianDSXML* lgptr = static_cast <LagrangianDSXML*>(dsxml);

  // -- Size of the system and number of degrees of freedom --
  ndof = lgptr->getNdof();

  // -- Vector and matrix members memory allocation --

  // q0 - velocity0
  q0 = new SimpleVector(lgptr->getQ0()); // required node in xml file
  if (q0->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS::xml constructor - size of input q0 differs from ndof");
  velocity0 = new SimpleVector(lgptr->getVelocity0()); // required node in xml file
  if (velocity0->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS::xml constructor - size of input velocity0 differs from ndof");

  // q and velocity
  if (lgptr->hasQ())
    q = new SimpleVector(lgptr->getQ()); // Means that first q is different from q0 ?? Strange case ...
  else
    q = new SimpleVector(*q0);           // initialize q with q0
  if (lgptr->hasVelocity())
    velocity = new SimpleVector(lgptr->getVelocity0()); // same remark as for q
  else
    velocity = new SimpleVector(*velocity0);

  qFree = new SimpleVector(ndof);
  velocityFree = new SimpleVector(ndof);

  for (unsigned int i = 0; i < 3; ++i)
  {
    isQAllocatedIn.push_back(true);
    isVelocityAllocatedIn.push_back(true);
  }

  // p
  p = new SimpleVector(ndof);

  // Memories
  isQAllocatedIn.push_back(false); // velocityMemory
  isVelocityAllocatedIn.push_back(false); // velocityMemory
  if (lgptr->hasQMemory())   // qMemory
  {
    qMemory = new SiconosMemory(lgptr->getQMemoryXML());
    isQAllocatedIn[3] = true;
  }
  if (lgptr->hasVelocityMemory())   // velocityMemory
  {
    velocityMemory = new SiconosMemory(lgptr->getVelocityMemoryXML());
    isVelocityAllocatedIn[3] = true;
  }

  // set default parameters list for plug-in and check-plug-in to false
  for (unsigned int i = 0; i < 8; ++i)
  {
    isLDSPlugin.push_back(false);
    // parametersList is set to default (vector of size 1 set to 0)
    parametersList.push_back(NULL);
    isParametersListAllocatedIn.push_back(true);
  }
  vector<SimpleVector*>::iterator iter;
  for (iter = parametersList.begin(); iter != parametersList.end(); iter++)
  {
    *iter = new SimpleVector(1);
    (*iter)->zero();
  }

  string plugin;

  // mass
  mass = new SimpleMatrix(ndof, ndof);
  if (lgptr->isMPlugin()) // if mass is plugged
  {
    plugin = lgptr->getMPlugin();
    setComputeMassFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    isLDSPlugin[0] = true;
  }
  else *mass = lgptr->getMMatrix();

  // fInt
  fInt = new SimpleVector(ndof);
  if (lgptr->hasFint())  // if fInt is given
  {
    if (lgptr->isFintPlugin()) // if fInt is plugged
    {
      plugin = lgptr->getFintPlugin();
      setComputeFIntFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      isLDSPlugin[1] = true;
    }
    else *fInt = lgptr->getFintVector();
  }

  // fExt
  fExt = new SimpleVector(ndof);
  if (lgptr->hasFext())  // if fExt is given
  {
    if (lgptr->isFextPlugin())// if fExt is plugged
    {
      plugin = lgptr->getFextPlugin();
      setComputeFExtFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      isLDSPlugin[2] = true;
    }
    else *fExt = lgptr->getFextVector();
  }

  // NNL
  NNL = new SimpleVector(ndof);
  if (lgptr ->hasNNL())// if NNL is given
  {
    if (lgptr->isNNLPlugin())// if NNL is plugged
    {
      plugin = lgptr->getNNLPlugin();
      setComputeNNLFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      isLDSPlugin[3] = true;
    }
    else *NNL = lgptr->getNNLVector();
  }

  // set "check" functions for variables above
  for (unsigned int i = 0; i < 3; ++i) areForcesAllocatedIn.push_back(true);

  // jacobian of fInt according to q
  jacobianQFInt = new SimpleMatrix(ndof, ndof);
  if (lgptr ->hasJacobianQFint())// if given
  {
    if (lgptr->isJacobianQFintPlugin())// if is plugged
    {
      plugin = lgptr->getJacobianQFintPlugin();
      setComputeJacobianQFIntFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      isLDSPlugin[4] = true;
    }
    else *jacobianQFInt = lgptr->getJacobianQFintMatrix();
  }

  // jacobian of fInt according to velocity
  jacobianVelocityFInt = new SimpleMatrix(ndof, ndof);
  if (lgptr ->hasJacobianVelocityFint())// if given
  {
    if (lgptr->isJacobianVelocityFintPlugin())// if is plugged
    {
      plugin = lgptr->getJacobianVelocityFintPlugin();
      setComputeJacobianVelocityFIntFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      isLDSPlugin[5] = true;
    }
    else *jacobianVelocityFInt = lgptr->getJacobianVelocityFintMatrix();
  }

  // jacobian of NNL according to q
  jacobianQNNL = new SimpleMatrix(ndof, ndof);
  if (lgptr -> hasJacobianQNNL()) // if given
  {
    if (lgptr->isJacobianQNNLPlugin())// if is plugged
    {
      plugin = lgptr->getJacobianQNNLPlugin();
      setComputeJacobianQNNLFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      isLDSPlugin[6] = true;
    }
    else *jacobianQNNL = lgptr->getJacobianQNNLMatrix();
  }

  // jacobian of NNL according to velocity
  jacobianVelocityNNL = new SimpleMatrix(ndof, ndof);
  if (lgptr ->  hasJacobianVelocityNNL())// if given
  {
    if (lgptr->isJacobianVelocityNNLPlugin())// if is plugged
    {
      plugin = lgptr->getJacobianVelocityNNLPlugin();
      setComputeJacobianVelocityNNLFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      isLDSPlugin[7] = true;
    }
    else *jacobianVelocityNNL = lgptr->getJacobianVelocityNNLMatrix();
  }

  // set "check" functions for jacobian matrices
  for (unsigned int i = 0; i < 4; ++i) isJacobianAllocatedIn.push_back(true);

  // set DynamicalClass members
  connectToDS();
}

// From a set of data; Mass filled-in directly from a siconosMatrix
LagrangianDS::LagrangianDS(const int& newNumber, const unsigned int& newNdof,
                           const SimpleVector& newQ0, const SimpleVector& newVelocity0,
                           const SiconosMatrix& newMass):
  DynamicalSystem(), ndof(newNdof), q(NULL), q0(NULL), qFree(NULL), qMemory(NULL), velocity(NULL), velocity0(NULL),
  velocityFree(NULL), velocityMemory(NULL),   p(NULL), mass(NULL), fInt(NULL), fExt(NULL),
  NNL(NULL), jacobianQFInt(NULL), jacobianVelocityFInt(NULL), jacobianQNNL(NULL), jacobianVelocityNNL(NULL),
  isPAllocatedIn(true), isMassAllocatedIn(true),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL), computeJacobianQFIntPtr(NULL),
  computeJacobianVelocityFIntPtr(NULL), computeJacobianQNNLPtr(NULL), computeJacobianVelocityNNLPtr(NULL), isDSup(false)
{
  // --- DS BASE-CLASS MEMBERS ---
  // --- Settings and xml load ---
  DSType = LNLDS;
  number = newNumber;
  n = 2 * ndof;

  // --- LAGRANGIAN INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --
  mass = new SimpleMatrix(newMass);
  q0 = new SimpleVector(newQ0);
  q = new SimpleVector(*q0);
  qFree = new SimpleVector(ndof);

  unsigned int i;
  for (i = 0; i < 3; ++i)
  {
    isQAllocatedIn.push_back(true);
    isVelocityAllocatedIn.push_back(true);
  }

  isQAllocatedIn.push_back(false); // qMemory
  isVelocityAllocatedIn.push_back(false); // velocityMemory
  velocity0 = new SimpleVector(newVelocity0);
  velocity = new SimpleVector(*velocity0);
  velocityFree = new SimpleVector(ndof);

  p = new SimpleVector(ndof);

  fInt = new SimpleVector(ndof);
  fExt = new SimpleVector(ndof);

  // parametersList is set to default (vector of size 1 set to 0)
  for (i = 0; i < 8; ++i)
  {
    isLDSPlugin.push_back(false);
    // parametersList is set to default (vector of size 1 set to 0)
    parametersList.push_back(NULL);
    isParametersListAllocatedIn.push_back(true);
  }
  vector<SimpleVector*>::iterator iter;
  for (iter = parametersList.begin(); iter != parametersList.end(); iter++)
  {
    *iter = new SimpleVector(1);
    (*iter)->zero();
  }

  NNL = new SimpleVector(ndof);
  areForcesAllocatedIn.resize(3, true);
  jacobianQFInt = new SimpleMatrix(ndof, ndof);
  jacobianVelocityFInt = new SimpleMatrix(ndof, ndof);
  jacobianQNNL = new SimpleMatrix(ndof, ndof);
  jacobianVelocityNNL = new SimpleMatrix(ndof, ndof);
  for (i = 0; i < 4; ++i) isJacobianAllocatedIn.push_back(true);

  // === Settings and links for master dynamical system variables ===
  connectToDS();

  setComputeMassFunction("DefaultPlugin.so", "computeMass");
  setComputeFIntFunction("DefaultPlugin.so", "computeFInt");
  setComputeFExtFunction("DefaultPlugin.so", "computeFExt");
  setComputeNNLFunction("DefaultPlugin.so", "computeNNL");
  setComputeJacobianQFIntFunction("DefaultPlugin.so", "computeJacobianQFInt");
  setComputeJacobianVelocityFIntFunction("DefaultPlugin.so", "computeJacobianVelocityFInt");
  setComputeJacobianQNNLFunction("DefaultPlugin.so", "computeJacobianQNNL");
  setComputeJacobianVelocityNNLFunction("DefaultPlugin.so", "computeJacobianVelocityNNL");
}

// From a set of data - Mass loaded from a plugin
LagrangianDS::LagrangianDS(const int& newNumber, const unsigned int& newNdof,
                           const SimpleVector& newQ0, const SimpleVector& newVelocity0, const string& massName):
  DynamicalSystem(), ndof(newNdof), q(NULL), q0(NULL), qFree(NULL), qMemory(NULL), velocity(NULL),
  velocity0(NULL), velocityFree(NULL), velocityMemory(NULL), p(NULL), mass(NULL),
  fInt(NULL), fExt(NULL), NNL(NULL), jacobianQFInt(NULL), jacobianVelocityFInt(NULL),
  jacobianQNNL(NULL), jacobianVelocityNNL(NULL), isPAllocatedIn(true), isMassAllocatedIn(true),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL), computeJacobianQFIntPtr(NULL),
  computeJacobianVelocityFIntPtr(NULL), computeJacobianQNNLPtr(NULL), computeJacobianVelocityNNLPtr(NULL), isDSup(false)
{
  // --- DS BASE-CLASS MEMBERS ---
  // --- Settings and xml load ---
  DSType = LNLDS;
  number = newNumber;

  // -- plugins --
  string plugin;
  // VectorField
  // JacobianX

  // --- LAGRANGIAN INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --
  q0 = new SimpleVector(newQ0);
  q = new SimpleVector(*q0);
  qFree = new SimpleVector(ndof);

  unsigned int i;
  for (i = 0; i < 3; ++i)
  {
    isQAllocatedIn.push_back(true);
    isVelocityAllocatedIn.push_back(true);
  }
  isQAllocatedIn.push_back(false); // qMemory
  isVelocityAllocatedIn.push_back(false); // velocityMemory
  velocity0 = new SimpleVector(newVelocity0);
  velocity = new SimpleVector(*velocity);
  velocityFree = new SimpleVector(ndof);

  p = new SimpleVector(ndof);

  mass = new SimpleMatrix(ndof, ndof);
  setComputeMassFunction(cShared.getPluginName(massName), cShared.getPluginFunctionName(massName));
  fInt = new SimpleVector(ndof);
  fExt = new SimpleVector(ndof);

  // parametersList is set to default (vector of size 1 set to 0)
  for (i = 0; i < 8; ++i)
  {
    isLDSPlugin.push_back(false);
    // parametersList is set to default (vector of size 1 set to 0)
    parametersList.push_back(NULL);
    isParametersListAllocatedIn.push_back(true);
  }
  vector<SimpleVector*>::iterator iter;
  for (iter = parametersList.begin(); iter != parametersList.end(); iter++)
  {
    *iter = new SimpleVector(1);
    (*iter)->zero();
  }

  NNL = new SimpleVector(ndof);
  areForcesAllocatedIn.resize(3, true);
  jacobianQFInt = new SimpleMatrix(ndof, ndof);
  jacobianVelocityFInt = new SimpleMatrix(ndof, ndof);
  jacobianQNNL = new SimpleMatrix(ndof, ndof);
  jacobianVelocityNNL = new SimpleMatrix(ndof, ndof);
  for (i = 0; i < 4; ++i) isJacobianAllocatedIn.push_back(true);


  // --- x, xDot and xFree update ---

  connectToDS();

  //   --- default plug-in ---
  setComputeFIntFunction("DefaultPlugin.so", "computeFInt");
  setComputeFExtFunction("DefaultPlugin.so", "computeFExt");
  setComputeNNLFunction("DefaultPlugin.so", "computeNNL");
  setComputeJacobianQFIntFunction("DefaultPlugin.so", "computeJacobianQFInt");
  setComputeJacobianVelocityFIntFunction("DefaultPlugin.so", "computeJacobianVelocityFInt");
  setComputeJacobianQNNLFunction("DefaultPlugin.so", "computeJacobianQNNL");
  setComputeJacobianVelocityNNLFunction("DefaultPlugin.so", "computeJacobianVelocityNNL");
}

// copy constructor
LagrangianDS::LagrangianDS(const DynamicalSystem & newDS):
  DynamicalSystem(), ndof(0), q(NULL), q0(NULL), qFree(NULL), qMemory(NULL),
  velocity(NULL), velocity0(NULL), velocityFree(NULL), velocityMemory(NULL), p(NULL), mass(NULL),
  fInt(NULL), fExt(NULL), NNL(NULL), jacobianQFInt(NULL), jacobianVelocityFInt(NULL),
  jacobianQNNL(NULL), jacobianVelocityNNL(NULL), isPAllocatedIn(true), isMassAllocatedIn(true),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL), computeJacobianQFIntPtr(NULL),
  computeJacobianVelocityFIntPtr(NULL), computeJacobianQNNLPtr(NULL), computeJacobianVelocityNNLPtr(NULL), isDSup(false)
{
  if (newDS.getType() != LNLDS && newDS.getType() != LTIDS)
    RuntimeException::selfThrow("LagrangianDS - copy constructor: try to copy into a LagrangianDS a DS of type: " + newDS.getType());

  DSType = LNLDS;

  // convert newDS to lagrangianDS by keeping const options
  const LagrangianDS * lnlds = static_cast<const LagrangianDS*>(&newDS);

  ndof = lnlds->getNdof();

  mass = new SimpleMatrix(lnlds->getMass());

  q0 = new SimpleVector(lnlds->getQ0());
  q = new SimpleVector(lnlds->getQ());
  qFree = new SimpleVector(lnlds->getQFree());

  unsigned int i;
  for (i = 0; i < 3; ++i)
  {
    isQAllocatedIn.push_back(true);
    isVelocityAllocatedIn.push_back(true);
  }
  isQAllocatedIn.push_back(false); // qMemory
  isVelocityAllocatedIn.push_back(false); // velocityMemory
  if (lnlds->getQMemoryPtr() != NULL)
  {
    qMemory = new SiconosMemory(lnlds->getQMemory());
    isQAllocatedIn[3] = true;
  }

  velocity0 = new SimpleVector(lnlds->getVelocity0());
  velocity  = new SimpleVector(lnlds->getVelocity0());
  velocityFree = new SimpleVector(lnlds->getVelocityFree());
  if (lnlds->getVelocityMemoryPtr() != NULL)
  {
    velocityMemory = new SiconosMemory(lnlds->getVelocityMemory());
    isVelocityAllocatedIn[3] = true;
  }

  p = new SimpleVector(lnlds->getP());

  areForcesAllocatedIn.resize(3, false);
  if (lnlds->getFIntPtr() != NULL)
  {
    fInt = new SimpleVector(lnlds->getFInt());
    areForcesAllocatedIn[0] = true;
  }

  setParametersListVector(lnlds->getParametersListVector());

  if (lnlds->getFExtPtr() != NULL)
  {
    fExt = new SimpleVector(lnlds->getFExt());
    areForcesAllocatedIn[1] = true;
  }

  if (lnlds->getNNLPtr() != NULL)
  {
    NNL = new SimpleVector(lnlds->getNNL());
    areForcesAllocatedIn[2] = true;
  }

  isJacobianAllocatedIn.resize(4, false);

  if (lnlds->getJacobianQFIntPtr() != NULL)
  {
    jacobianQFInt = new SimpleMatrix(lnlds->getJacobianQFInt());
    isJacobianAllocatedIn[0] = true;
  }
  if (lnlds->getJacobianVelocityFIntPtr() != NULL)
  {
    jacobianVelocityFInt = new SimpleMatrix(lnlds->getJacobianVelocityFInt());
    isJacobianAllocatedIn[1] = true;
  }
  if (lnlds->getJacobianQNNLPtr() != NULL)
  {
    jacobianQNNL = new SimpleMatrix(lnlds->getJacobianQNNL());
    isJacobianAllocatedIn[2] = true;
  }
  if (lnlds->getJacobianVelocityNNLPtr() != NULL)
  {
    jacobianVelocityNNL = new SimpleMatrix(lnlds->getJacobianVelocityNNL());
    isJacobianAllocatedIn[3] = true;
  }

  isLDSPlugin = lnlds->getIsLDSPlugin();
  string pluginPath, functionName;
  if (isLDSPlugin[0])
  {
    massFunctionName = lnlds -> getMassFunctionName();
    functionName = cShared.getPluginFunctionName(massFunctionName);
    pluginPath  = cShared.getPluginName(massFunctionName);
    setComputeMassFunction(pluginPath, functionName);
  }

  if (isLDSPlugin[1])
  {
    fIntFunctionName = lnlds -> getFIntFunctionName();
    functionName = cShared.getPluginFunctionName(fIntFunctionName);
    pluginPath  = cShared.getPluginName(fIntFunctionName);
    setComputeFIntFunction(pluginPath, functionName);
  }

  if (isLDSPlugin[2])
  {
    fExtFunctionName = lnlds -> getFExtFunctionName();
    functionName = cShared.getPluginFunctionName(fExtFunctionName);
    pluginPath  = cShared.getPluginName(fExtFunctionName);
    setComputeFExtFunction(pluginPath, functionName);
  }

  if (isLDSPlugin[3])
  {
    NNLFunctionName = lnlds -> getNNLFunctionName();
    functionName = cShared.getPluginFunctionName(NNLFunctionName);
    pluginPath  = cShared.getPluginName(NNLFunctionName);
    setComputeNNLFunction(pluginPath, functionName);
  }

  if (isLDSPlugin[4])
  {
    jacobianQFIntFunctionName = lnlds -> getJacobianQFIntFunctionName();
    functionName = cShared.getPluginFunctionName(jacobianQFIntFunctionName);
    pluginPath  = cShared.getPluginName(jacobianQFIntFunctionName);
    setComputeJacobianQFIntFunction(pluginPath, functionName);
  }

  if (isLDSPlugin[5])
  {
    jacobianVelocityFIntFunctionName = lnlds -> getJacobianVelocityFIntFunctionName();
    functionName = cShared.getPluginFunctionName(jacobianVelocityFIntFunctionName);
    pluginPath  = cShared.getPluginName(jacobianVelocityFIntFunctionName);
    setComputeJacobianVelocityFIntFunction(pluginPath, functionName);
  }

  if (isLDSPlugin[6])
  {
    jacobianQNNLFunctionName = lnlds -> getJacobianQNNLFunctionName();
    functionName = cShared.getPluginFunctionName(jacobianQNNLFunctionName);
    pluginPath  = cShared.getPluginName(jacobianQNNLFunctionName);
    setComputeJacobianQNNLFunction(pluginPath, functionName);
  }

  if (isLDSPlugin[7])
  {
    jacobianVelocityNNLFunctionName = lnlds -> getJacobianVelocityNNLFunctionName();
    functionName = cShared.getPluginFunctionName(jacobianVelocityNNLFunctionName);
    pluginPath  = cShared.getPluginName(jacobianVelocityNNLFunctionName);
    setComputeJacobianVelocityNNLFunction(pluginPath, functionName);
  }

  // Set variables of top-class DynamicalSystem
  connectToDS();
}

// Destructor
LagrangianDS::~LagrangianDS()
{
  if (isQAllocatedIn[0]) delete q;
  q = NULL;
  if (isQAllocatedIn[1]) delete q0;
  q0 = NULL;
  if (isQAllocatedIn[2]) delete qFree;
  qFree = NULL;
  if (isQAllocatedIn[3]) delete qMemory;
  qMemory = NULL;
  if (isVelocityAllocatedIn[0]) delete velocity ;
  velocity = NULL;
  if (isVelocityAllocatedIn[1])delete velocity0 ;
  velocity0 = NULL;
  if (isVelocityAllocatedIn[2])delete velocityFree ;
  velocityFree = NULL;
  if (isVelocityAllocatedIn[3])delete velocityMemory;
  velocityMemory = NULL;
  if (isPAllocatedIn) delete p ;
  p = NULL;
  if (isMassAllocatedIn) delete mass;
  mass = NULL;
  if (areForcesAllocatedIn[0])delete fInt ;
  fInt = NULL;
  if (areForcesAllocatedIn[1])delete fExt ;
  fExt = NULL;
  if (areForcesAllocatedIn[2])delete NNL ;
  NNL = NULL;
  if (isJacobianAllocatedIn[0])delete jacobianQFInt  ;
  jacobianQFInt = NULL;
  if (isJacobianAllocatedIn[1])delete  jacobianVelocityFInt ;
  jacobianVelocityFInt = NULL;
  if (isJacobianAllocatedIn[2])delete jacobianQNNL ;
  jacobianQNNL = NULL;
  if (isJacobianAllocatedIn[3])delete jacobianVelocityNNL ;
  jacobianVelocityNNL = NULL;

  for (unsigned int i = 0; i < parametersList.size(); i++)
  {
    if (isParametersListAllocatedIn[i]) delete parametersList[i];
    parametersList[i] = NULL;
  }
  isParametersListAllocatedIn.resize(8, false);

  // clean workMatrix map. Warning: if set/get functions are added for this object,
  // add a isAlloc. deque to check in-class allocation.
  map<string, SiconosMatrix*>::iterator it;
  for (it = workMatrix.begin(); it != workMatrix.end(); ++it)
    if (it->second != NULL) delete it->second;
}

void LagrangianDS::initialize(const double& time, const unsigned int& sizeOfMemory)
{
  // set q and velocity to q0 and velocity0
  *q = *q0;
  *velocity = *velocity0;

  // reset r and free vectors
  qFree->zero();
  velocityFree->zero();
  p->zero();

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  // compute plug-in values for mass, FInt etc ...
  isConstant["mass"] = false;
  SimpleVector* param;
  if (isLDSPlugin[0]) // mass
  {
    if (computeMassPtr == NULL)
      RuntimeException::selfThrow("LagrangianDS initialize, computeMass() is not linked to a plugin function");
    param = parametersList[0];
    computeMassPtr(ndof, &(*q)(0), &(*mass)(0, 0), &(*param)(0));
  }
  if (isLDSPlugin[1]) // FInt
  {
    if (computeFIntPtr == NULL)
      RuntimeException::selfThrow("LagrangianDS initialize, computeFInt() is not linked to a plugin function");
    param = parametersList[1];
    computeFIntPtr(ndof, &time, &(*q)(0), &(*velocity)(0), &(*fInt)(0), &(*param)(0));
  }
  if (isLDSPlugin[2]) // FExt
  {
    if (computeFExtPtr == NULL)
      RuntimeException::selfThrow("LagrangianDS initialize, computeFExt() is not linked to a plugin function");
    param = parametersList[2];
    computeFExtPtr(ndof, &time, &(*fExt)(0), &(*param)(0));
  }
  if (isLDSPlugin[3]) // NNL
  {
    if (computeNNLPtr == NULL)
      RuntimeException::selfThrow("LagrangianDS initialize, computeNNL() is not linked to a plugin function");
    param = parametersList[3];
    computeNNLPtr(ndof, &(*q)(0), &(*velocity)(0), &(*NNL)(0), &(*param)(0));
  }
  if (isLDSPlugin[4]) // nablaQFInt
  {
    if (computeJacobianQFIntPtr == NULL)
      RuntimeException::selfThrow("LagrangianDS initialize, computeJacobianQFInt() is not linked to a plugin function");
    param = parametersList[4];
    computeJacobianQFIntPtr(ndof, &time, &(*q)(0), &(*velocity)(0), &(*jacobianQFInt)(0, 0), &(*param)(0));
  }
  if (isLDSPlugin[5]) // nablaVFInt
  {
    if (computeJacobianVelocityFIntPtr == NULL)
      RuntimeException::selfThrow("LagrangianDS initialize, computeJacobianVelocityFInt() is not linked to a plugin function");
    param = parametersList[5];
    computeJacobianVelocityFIntPtr(ndof, &time, &(*q)(0), &(*velocity)(0), &(*jacobianVelocityFInt)(0, 0), &(*param)(0));
  }
  if (isLDSPlugin[6]) // nablaQNNL
  {
    if (computeJacobianQNNLPtr == NULL)
      RuntimeException::selfThrow("LagrangianDS initialize, computeJacobianQNNL() is not linked to a plugin function");
    param = parametersList[6];
    computeJacobianQNNLPtr(ndof, &(*q)(0), &(*velocity)(0), &(*jacobianQNNL)(0, 0), &(*param)(0));
  }
  if (isLDSPlugin[7]) // nablaVNNL
  {
    if (computeJacobianVelocityNNLPtr == NULL)
      RuntimeException::selfThrow("LagrangianDS initialize, computeJacobianVelocityNNL() is not linked to a plugin function");
    param = parametersList[7];
    computeJacobianVelocityNNLPtr(ndof, &(*q)(0), &(*velocity)(0), &(*jacobianVelocityNNL)(0, 0), &(*param)(0));
  }

  // vectorField
  // note that xDot(0) = velocity, with pointer link, must already be set.
  SiconosVector* tmp = xDot->getVectorPtr(1); // Pointer link!
  // Compute M-1
  workMatrix["inverseOfMass"] = new SimpleMatrix(*mass);
  workMatrix["inverseOfMass"]->PLUFactorizationInPlace();
  workMatrix["inverseOfMass"]->PLUInverseInPlace();

  *tmp = *(workMatrix["inverseOfMass"]) * (*fExt - *fInt - *NNL);
  // todo: use linearSolve to avoid inversion ? Or save M-1 to avoid re-computation. See this when "M" will be added in DS or LDS.

  // jacobianX

  SiconosMatrix * tmp2 = jacobianX->getBlockPtr(1, 0); // Pointer link!
  // !!! jacobian of M according to q is not take into account at the time !!!
  *tmp2 = -1 * *workMatrix["inverseOfMass"] * (*jacobianQFInt + *jacobianQNNL);
  tmp2 = jacobianX->getBlockPtr(1, 1); // Pointer link!
  *tmp2 = -1 * *workMatrix["inverseOfMass"] * (*jacobianVelocityFInt + *jacobianVelocityNNL);

  isDSup = true;

  // \todo: control terms handling
}

// --- GETTERS/SETTERS ---

void LagrangianDS::setQ(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ: inconsistent input vector size ");

  if (q == NULL)
  {
    q = new SimpleVector(newValue);
    isQAllocatedIn[0] = true;
  }
  else
    *q = newValue;
  isDSup = false;
}

void LagrangianDS::setQPtr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQPtr: inconsistent input vector size ");

  if (isQAllocatedIn[0]) delete q;
  q = newPtr;
  isQAllocatedIn[0] = false;
  isDSup = false;
}

void LagrangianDS::setQ0(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ0: inconsistent input vector size ");

  if (q0 == NULL)
  {
    q0 = new SimpleVector(newValue);
    isQAllocatedIn[1] = true;
  }
  else
    *q0 = newValue;
}

void LagrangianDS::setQ0Ptr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ0Ptr: inconsistent input vector size ");

  if (isQAllocatedIn[1]) delete q0;
  q0 = newPtr;
  isQAllocatedIn[1] = false;
}

void LagrangianDS::setQFree(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQFree: inconsistent input vector size ");

  if (qFree == NULL)
  {
    qFree = new SimpleVector(newValue);
    isQAllocatedIn[2] = true;
  }
  else
    *qFree = newValue;
}

void LagrangianDS::setQFreePtr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQFreePtr: inconsistent input vector size ");

  if (isQAllocatedIn[2]) delete qFree;
  qFree = newPtr;
  isQAllocatedIn[2] = false;
}

void LagrangianDS::setQMemory(const SiconosMemory& newValue)
{
  if (qMemory == NULL)
  {
    qMemory = new SiconosMemory(newValue);
    isQAllocatedIn[3] = true;
  }
  else
    *qMemory = newValue;
}

void LagrangianDS::setQMemoryPtr(SiconosMemory * newPtr)
{
  if (isQAllocatedIn[3]) delete qMemory;
  qMemory = newPtr;
  isQAllocatedIn[3] = false;
}

void LagrangianDS::setVelocity(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity: inconsistent input vector size ");

  if (velocity == NULL)
  {
    velocity = new SimpleVector(newValue);
    isVelocityAllocatedIn[0] = true;
  }
  else
    *velocity = newValue;
  isDSup = false;
}

void LagrangianDS::setVelocityPtr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocityPtr: inconsistent input vector size ");

  if (isVelocityAllocatedIn[0]) delete velocity;
  velocity = newPtr;
  isVelocityAllocatedIn[0] = false;
  isDSup = false;
}

void LagrangianDS::setVelocity0(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity0: inconsistent input vector size ");

  if (velocity0 == NULL)
  {
    velocity0 = new SimpleVector(newValue);
    isVelocityAllocatedIn[1] = true;
  }
  else
    *velocity0 = newValue;
}

void LagrangianDS::setVelocity0Ptr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity0Ptr: inconsistent input vector size ");

  if (isVelocityAllocatedIn[1]) delete velocity0;
  velocity0 = newPtr;
  isVelocityAllocatedIn[1] = false;
}

void LagrangianDS::setVelocityFree(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocityFree: inconsistent input vector size ");

  if (velocityFree == NULL)
  {
    velocityFree = new SimpleVector(newValue);
    isVelocityAllocatedIn[2] = true;
  }
  else
    *velocityFree = newValue;
}

void LagrangianDS::setVelocityFreePtr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocityFreePtr: inconsistent input vector size ");

  if (isVelocityAllocatedIn[2]) delete velocityFree;
  velocityFree = newPtr;
  isVelocityAllocatedIn[2] = false;
}

void LagrangianDS::setVelocityMemory(const SiconosMemory& newValue)
{
  if (velocityMemory == NULL)
  {
    velocityMemory = new SiconosMemory(newValue);
    isVelocityAllocatedIn[3] = true;
  }
  else
    *velocityMemory = newValue;
}

void LagrangianDS::setVelocityMemoryPtr(SiconosMemory * newPtr)
{
  if (isVelocityAllocatedIn[3]) delete velocityMemory;
  velocityMemory = newPtr;
  isVelocityAllocatedIn[3] = false;
}

void LagrangianDS::setP(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setP: inconsistent input vector size ");

  if (p == NULL)
  {
    p = new SimpleVector(newValue);
    isPAllocatedIn = true;
  }
  else
    *p = newValue;
}

void LagrangianDS::setPPtr(SimpleVector *newPtr)
{

  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setPPtr: inconsistent input vector size ");

  if (isPAllocatedIn) delete p;
  p = newPtr;
  isPAllocatedIn = false;
}

void LagrangianDS::setMass(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setMass: inconsistent input matrix size ");

  if (mass == NULL)
  {
    mass = new SimpleMatrix(newValue);
    isMassAllocatedIn = true;
  }
  else
    *mass = newValue;
  isLDSPlugin[0] = false;
  isDSup = false;
}

void LagrangianDS::setMassPtr(SiconosMatrix *newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setMassPtr: inconsistent input matrix size ");

  if (isMassAllocatedIn) delete mass;
  mass = newPtr;
  isMassAllocatedIn = false;
  isLDSPlugin[0] = false;
  isDSup = false;
}

void LagrangianDS::setFInt(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setFInt: inconsistent dimensions with problem size for input vector FInt");

  if (fInt == NULL)
  {
    fInt = new SimpleVector(newValue);
    areForcesAllocatedIn[0] = true;
  }
  else
    *fInt = newValue;
  isLDSPlugin[1] = false;
  isDSup = false;
}

void LagrangianDS::setFIntPtr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setFIntPtr: inconsistent input matrix size ");

  if (areForcesAllocatedIn[0]) delete fInt;
  fInt = newPtr;
  areForcesAllocatedIn[0] = false;
  isLDSPlugin[1] = false;
  isDSup = false;
}

void LagrangianDS::setFExt(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setFExt: inconsistent dimensions with problem size for input vector FExt");

  if (fExt == NULL)
  {
    fExt = new SimpleVector(newValue);
    areForcesAllocatedIn[1] = true;
  }
  else
    *fExt = newValue;
  isLDSPlugin[2] = false;
  isDSup = false;
}

void LagrangianDS::setFExtPtr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setFIntPtr: inconsistent input matrix size ");
  if (areForcesAllocatedIn[1]) delete fExt;
  fExt = newPtr;
  areForcesAllocatedIn[1] = false;
  isLDSPlugin[2] = false;
  isDSup = false;
}

void LagrangianDS::setNNL(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setNNL: inconsistent dimensions with problem size for input vector NNL");

  if (NNL == NULL)
  {
    NNL = new SimpleVector(ndof);
    areForcesAllocatedIn[2] = true;
  }
  else
    *NNL = newValue;
  isLDSPlugin[3] = false;
  isDSup = false;
}

void LagrangianDS::setNNLPtr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setNNLPtr: inconsistent input matrix size ");

  if (areForcesAllocatedIn[2]) delete NNL;
  NNL = newPtr;
  areForcesAllocatedIn[2] = false;
  isLDSPlugin[3] = false;
  isDSup = false;
}

void LagrangianDS::setJacobianQFInt(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianQFInt: inconsistent dimensions with problem size for input matrix JacobianQFInt");

  if (jacobianQFInt == NULL)
  {
    jacobianQFInt = new SimpleMatrix(newValue);
    isJacobianAllocatedIn[0] = true;
  }
  else
    *jacobianQFInt = newValue;
  isLDSPlugin[4] = false;
  isDSup = false;
}

void LagrangianDS::setJacobianQFIntPtr(SiconosMatrix *newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianQFIntPtr: inconsistent input matrix size ");

  if (isJacobianAllocatedIn[0]) delete jacobianQFInt;
  jacobianQFInt = newPtr;
  isJacobianAllocatedIn[0] = false;
  isLDSPlugin[4] = false;
  isDSup = false;
}

void LagrangianDS::setJacobianVelocityFInt(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianVelocityFInt: inconsistent dimensions with problem size for input matrix JacobianVelocityFInt");

  if (jacobianVelocityFInt == NULL)
  {
    jacobianVelocityFInt = new SimpleMatrix(newValue);
    isJacobianAllocatedIn[1] = true;
  }
  else
    *jacobianVelocityFInt = newValue;
  isLDSPlugin[5] = false;
  isDSup = false;
}

void LagrangianDS::setJacobianVelocityFIntPtr(SiconosMatrix *newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianVelocityFIntPtr: inconsistent input matrix size ");

  if (isJacobianAllocatedIn[1]) delete jacobianVelocityFInt;
  jacobianVelocityFInt = newPtr;
  isJacobianAllocatedIn[1] = false;
  isLDSPlugin[5] = false;
  isDSup = false;
}

void LagrangianDS::setJacobianQNNL(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianQNNL: inconsistent dimensions with problem size for input matrix JacobianQNNL");

  if (jacobianQNNL == NULL)
  {
    jacobianQNNL = new SimpleMatrix(newValue);
    isJacobianAllocatedIn[2] = true;
  }
  else
    *jacobianQNNL = newValue;
  isLDSPlugin[6] = false;
  isDSup = false;
}

void LagrangianDS::setJacobianQNNLPtr(SiconosMatrix *newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianQNNLPtr: inconsistent input matrix size ");

  if (isJacobianAllocatedIn[2]) delete jacobianQNNL;
  jacobianQNNL = newPtr;
  isJacobianAllocatedIn[2] = false;
  isLDSPlugin[6] = false;
  isDSup = false;
}

void  LagrangianDS::setJacobianVelocityNNL(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianVelocityNNL: inconsistent dimensions with problem size for input matrix JacobianVelocityNNL");

  if (jacobianVelocityNNL == NULL)
  {
    jacobianVelocityNNL = new SimpleMatrix(newValue);
    isJacobianAllocatedIn[3] = true;
  }
  else
    *jacobianVelocityNNL = newValue;
  isLDSPlugin[7] = false;
  isDSup = false;
}

void LagrangianDS::setJacobianVelocityNNLPtr(SiconosMatrix *newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianVelocityNNLPtr: inconsistent input matrix size ");

  if (isJacobianAllocatedIn[3]) delete jacobianVelocityNNL;
  jacobianVelocityNNL = newPtr;
  isJacobianAllocatedIn[3] = false;
  isLDSPlugin[7] = false;
  isDSup = false;
}


// --- Plugins related functions ---
void LagrangianDS::computeMass()
{
  if (!isConstant["mass"])
  {
    if (computeMassPtr == NULL)
      RuntimeException::selfThrow("computeMass() is not linked to a plugin function");
    SimpleVector* param = parametersList[0];
    computeMassPtr(ndof, &(*q)(0), &(*mass)(0, 0), &(*param)(0));
  }
}

void LagrangianDS::computeMass(SimpleVector *q2)
{
  if (!isConstant["mass"])
  {
    if (computeMassPtr == NULL)
      RuntimeException::selfThrow("computeMass() is not linked to a plugin function");

    SimpleVector* param = parametersList[0];
    computeMassPtr(ndof, &(*q2)(0), &(*mass)(0, 0), &(*param)(0));
  }
}

void LagrangianDS::computeFInt(const double& time)
{
  if (computeFIntPtr == NULL)
    RuntimeException::selfThrow("computeFInt() is not linked to a plugin function");

  SimpleVector* param = parametersList[1];
  computeFIntPtr(ndof, &time, &(*q)(0), &(*velocity)(0), &(*fInt)(0), &(*param)(0));
}
void LagrangianDS::computeFInt(const double& time, SimpleVector *q2, SimpleVector *velocity2)
{
  if (computeFIntPtr == NULL)
    RuntimeException::selfThrow("computeFInt() is not linked to a plugin function");

  SimpleVector* param = parametersList[1];
  computeFIntPtr(ndof, &time, &(*q2)(0), &(*velocity2)(0), &(*fInt)(0), &(*param)(0));
}

void LagrangianDS::computeFExt(const double& time)
{
  if (computeFExtPtr == NULL)
    RuntimeException::selfThrow("computeFExt() is not linked to a plugin function");

  SimpleVector* param = parametersList[2];
  computeFExtPtr(ndof, &time, &(*fExt)(0), &(*param)(0));
}

void LagrangianDS::computeNNL()
{
  if (computeNNLPtr == NULL)
    RuntimeException::selfThrow("computeQ() is not linked to a plugin function");

  SimpleVector* param = parametersList[3];
  computeNNLPtr(ndof, &(*q)(0), &(*velocity)(0), &(*NNL)(0), &(*param)(0));
}

void LagrangianDS::computeNNL(SimpleVector *q2, SimpleVector *velocity2)
{
  if (computeNNLPtr == NULL)
    RuntimeException::selfThrow("computeQ() is not linked to a plugin function");

  SimpleVector* param = parametersList[3];
  computeNNLPtr(ndof, &(*q2)(0), &(*velocity2)(0), &(*NNL)(0), &(*param)(0));
}

void LagrangianDS::computeJacobianQFInt(const double& time)
{
  if (computeJacobianQFIntPtr == NULL)
    RuntimeException::selfThrow("computeJacobianQFInt() is not linked to a plugin function");

  SimpleVector* param = parametersList[4];
  computeJacobianQFIntPtr(ndof, &time, &(*q)(0), &(*velocity)(0), &(*jacobianQFInt)(0, 0), &(*param)(0));
}

void LagrangianDS::computeJacobianQFInt(const double& time, SimpleVector *q2, SimpleVector *velocity2)
{
  if (computeJacobianQFIntPtr == NULL)
    RuntimeException::selfThrow("computeJacobianQFInt() is not linked to a plugin function");

  SimpleVector* param = parametersList[4];
  computeJacobianQFIntPtr(ndof, &time, &(*q2)(0), &(*velocity2)(0), &(*jacobianQFInt)(0, 0), &(*param)(0));
}

void LagrangianDS::computeJacobianVelocityFInt(const double & time)
{
  if (computeJacobianVelocityFIntPtr == NULL)
    RuntimeException::selfThrow("computeJacobianVelocityFInt() is not linked to a plugin function");

  SimpleVector* param = parametersList[5];
  computeJacobianVelocityFIntPtr(ndof, &time, &(*q)(0), &(*velocity)(0), &(*jacobianVelocityFInt)(0, 0), &(*param)(0));
}
void LagrangianDS::computeJacobianVelocityFInt(const double & time, SimpleVector *q2, SimpleVector *velocity2)
{
  if (computeJacobianVelocityFIntPtr == NULL)
    RuntimeException::selfThrow("computeJacobianVelocityFInt() is not linked to a plugin function");

  SimpleVector* param = parametersList[5];
  computeJacobianVelocityFIntPtr(ndof, &time, &(*q2)(0), &(*velocity2)(0), &(*jacobianVelocityFInt)(0, 0), &(*param)(0));
}

void LagrangianDS::computeJacobianQNNL()
{
  if (computeJacobianQNNLPtr == NULL)
    RuntimeException::selfThrow("computeJacobianQNNL() is not linked to a plugin function");

  SimpleVector* param = parametersList[6];
  computeJacobianQNNLPtr(ndof, &(*q)(0), &(*velocity)(0), &(*jacobianQNNL)(0, 0), &(*param)(0));
}

void LagrangianDS::computeJacobianQNNL(SimpleVector *q2, SimpleVector *velocity2)
{
  if (computeJacobianQNNLPtr == NULL)
    RuntimeException::selfThrow("computeJacobianQNNL() is not linked to a plugin function");

  SimpleVector* param = parametersList[6];
  computeJacobianQNNLPtr(ndof, &(*q2)(0), &(*velocity2)(0), &(*jacobianQNNL)(0, 0), &(*param)(0));
}

void LagrangianDS::computeJacobianVelocityNNL()
{
  if (computeJacobianVelocityNNLPtr == NULL)
    RuntimeException::selfThrow("computeJacobianVelocityNNL() is not linked to a plugin function");

  SimpleVector* param = parametersList[7];
  computeJacobianVelocityNNLPtr(ndof, &(*q)(0), &(*velocity)(0), &(*jacobianVelocityNNL)(0, 0), &(*param)(0));
}
void LagrangianDS::computeJacobianVelocityNNL(SimpleVector *q2, SimpleVector *velocity2)
{
  if (computeJacobianVelocityNNLPtr == NULL)
    RuntimeException::selfThrow("computeJacobianVelocityNNL() is not linked to a plugin function");

  SimpleVector* param = parametersList[7];
  computeJacobianVelocityNNLPtr(ndof, &(*q2)(0), &(*velocity2)(0), &(*jacobianVelocityNNL)(0, 0), &(*param)(0));
}

void LagrangianDS::computeVectorField(const double& time)
{
  SiconosVector* tmp = xDot->getVectorPtr(1); // Pointer link!
  if (!isDSup) // if it is necessary to re-compute mass, FInt ..., ie if they have not been compiled during the present time step
  {
    computeMass();
    computeFInt(time);
    computeFExt(time);
    computeNNL();
    // note that xDot(0) = velocity with pointer link must already be set.

    // Compute M-1
    map<string, SiconosMatrix*>::iterator it = workMatrix.find("inverseOfMass");
    if (it == workMatrix.end()) // if it is the first call to computeVectorField
    {
      workMatrix["inverseOfMass"] = new SimpleMatrix(*mass);
      workMatrix["inverseOfMass"]->PLUFactorizationInPlace();
      workMatrix["inverseOfMass"]->PLUInverseInPlace();
    }
    else if (!isConstant["mass"])
    {
      *(workMatrix["inverseOfMass"]) = *mass;
      if (!((workMatrix["inverseOfMass"])->isFactorized()))
        workMatrix["inverseOfMass"]->PLUFactorizationInPlace();
      workMatrix["inverseOfMass"]->PLUInverseInPlace();
    }
    isDSup = true;
  }
  *tmp = *(workMatrix["inverseOfMass"]) * (*fExt - *fInt - *NNL);
  // todo: use linearSolve to avoid inversion ? Or save M-1 to avoid re-computation. See this when "M" will be added in DS or LDS.
}

void LagrangianDS::computeJacobianX(const double& time)
{
  SiconosMatrix* tmp = jacobianX->getBlockPtr(1, 0); // Pointer link!
  if (!isDSup) // if it is necessary to re-compute mass, FInt ..., ie if they have not been compiled during the present time step
  {
    computeMass();
    computeFInt(time);
    computeFExt(time);
    computeNNL();
    // note that xDot(0) = velocity with pointer link must already be set.

    // Compute M-1 if required
    map<string, SiconosMatrix*>::iterator it = workMatrix.find("inverseOfMass");
    if (it == workMatrix.end()) // if it is the first call to computeVectorField or computeJacobianX
    {
      workMatrix["inverseOfMass"] = new SimpleMatrix(*mass);
      workMatrix["inverseOfMass"]->PLUFactorizationInPlace();
      workMatrix["inverseOfMass"]->PLUInverseInPlace();
    }
    else if (!isConstant["mass"]) // add a flag to avoid re-computation of M when vectorField has been called before  ...
    {
      *(workMatrix["inverseOfMass"]) = *mass;
      if (!workMatrix["inverseOfMass"]->isFactorized())
        workMatrix["inverseOfMass"]->PLUFactorizationInPlace();
      workMatrix["inverseOfMass"]->PLUInverseInPlace();
    }

    // !!! jacobian of M according to q is not take into account at the time !!!
    isDSup = true;
  }

  *tmp = -1 * *workMatrix["inverseOfMass"] * (*jacobianQFInt + *jacobianQNNL);
  tmp = jacobianX->getBlockPtr(1, 1); // Pointer link!
  *tmp = -1 * *workMatrix["inverseOfMass"] * (*jacobianVelocityFInt + *jacobianVelocityNNL);
}

void LagrangianDS::setVectorFieldFunction(const std::string & pluginPath, const std::string& functionName)
{
  cout << " /!\\ LagrangianDS setVectorFieldFunction: useless function call. Set FInt, NNL ...plug-in /!\\ ." << endl;
}

void LagrangianDS::setComputeJacobianXFunction(const std::string & pluginPath, const std::string & functionName)
{
  cout << " /!\\ LagrangianDS setVectorFieldFunction: useless function call. Set FInt, NNL ... plug-in /!\\ ." << endl;
}

void LagrangianDS::setComputeMassFunction(const string& pluginPath, const string& functionName)
{
  if (mass == NULL)
  {
    mass = new SimpleMatrix(ndof, ndof);
    isMassAllocatedIn = true;
  }

  cShared.setFunction(&computeMassPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  massFunctionName = plugin + ":" + functionName;
  isLDSPlugin[0] = true;
}

void LagrangianDS::setComputeFIntFunction(const string& pluginPath, const string& functionName)
{
  if (fInt == NULL)
  {
    fInt = new SimpleVector(ndof);
    areForcesAllocatedIn[0] = true;
  }

  cShared.setFunction(&computeFIntPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  fIntFunctionName = plugin + ":" + functionName;
  isLDSPlugin[1] = true;
}

void LagrangianDS::setComputeFExtFunction(const string& pluginPath, const string& functionName)
{
  if (fExt == NULL)
  {
    fExt = new SimpleVector(ndof);
    areForcesAllocatedIn[1] = true;
  }

  cShared.setFunction(&computeFExtPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  fExtFunctionName = plugin + ":" + functionName;
  isLDSPlugin[2] = true;
}

void LagrangianDS::setComputeNNLFunction(const string& pluginPath, const string& functionName)
{
  if (NNL == NULL)
  {
    NNL = new SimpleVector(ndof);
    areForcesAllocatedIn[2] = true;
  }

  cShared.setFunction(&computeNNLPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  NNLFunctionName = plugin + ":" + functionName;
  isLDSPlugin[3] = true;
}

void LagrangianDS::setComputeJacobianQFIntFunction(const string& pluginPath, const string& functionName)
{
  if (jacobianQFInt == NULL)
  {
    jacobianQFInt = new SimpleMatrix(ndof, ndof);
    isJacobianAllocatedIn[0] = true;
  }

  cShared.setFunction(&computeJacobianQFIntPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  jacobianQFIntFunctionName = plugin + ":" + functionName;
  isLDSPlugin[4] = true;
}

void LagrangianDS::setComputeJacobianVelocityFIntFunction(const string& pluginPath, const string& functionName)
{
  if (jacobianVelocityFInt == NULL)
  {
    jacobianVelocityFInt = new SimpleMatrix(ndof, ndof);
    isJacobianAllocatedIn[1] = true;
  }

  cShared.setFunction(&computeJacobianVelocityFIntPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  jacobianVelocityFIntFunctionName = plugin + ":" + functionName;
  isLDSPlugin[5] = true;
}

void LagrangianDS::setComputeJacobianQNNLFunction(const string& pluginPath, const string& functionName)
{
  if (jacobianQNNL == NULL)
  {
    jacobianQNNL = new SimpleMatrix(ndof, ndof);
    isJacobianAllocatedIn[2] = true;
  }

  cShared.setFunction(&computeJacobianQNNLPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  jacobianQNNLFunctionName = plugin + ":" + functionName;
  isLDSPlugin[6] = true;
}

void LagrangianDS::setComputeJacobianVelocityNNLFunction(const string& pluginPath, const string& functionName)
{
  if (jacobianVelocityNNL == NULL)
  {
    jacobianVelocityNNL = new SimpleMatrix(ndof, ndof);
    isJacobianAllocatedIn[3] = true;
  }

  cShared.setFunction(&computeJacobianVelocityNNLPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  jacobianVelocityNNLFunctionName = plugin + ":" + functionName;
  isLDSPlugin[7] = true;
}

void LagrangianDS::setParametersListVector(const std::vector<SimpleVector*>& newVector)
{
  // copy!!
  for (unsigned int i = 0; i < parametersList.size(); ++i)
  {
    if (isParametersListAllocatedIn[i]) delete parametersList[i];
    *(parametersList[i]) = *(newVector[i]);
    isParametersListAllocatedIn[i] = true;
  }
}

void LagrangianDS::setParametersList(const SimpleVector& newValue, const unsigned int & index)
{
  if (isParametersListAllocatedIn[index]) delete parametersList[index];
  parametersList[index] = new SimpleVector(newValue);
  isParametersListAllocatedIn[index] = true;
}

void LagrangianDS::setParametersListPtr(SimpleVector *newPtr, const unsigned int & index)
{
  if (isParametersListAllocatedIn[index]) delete parametersList[index];
  parametersList[index] = newPtr;
  isParametersListAllocatedIn[index] = false;
}

void LagrangianDS::saveDSToXML()
{
  //--- general DS data---
  saveDSDataToXML();
  // --- other data ---
  if (dsxml != NULL)
  {
    LagrangianDSXML* lgptr = static_cast <LagrangianDSXML*>(dsxml);
    lgptr->setNdof(ndof);
    lgptr->setMPlugin(massFunctionName);
    lgptr->setQ(q);
    lgptr->setQ0(q0);
    lgptr->setQMemory(qMemory);
    lgptr->setVelocity(velocity);
    lgptr->setVelocity0(velocity0);
    lgptr->setVelocityMemory(velocityMemory);

    // FExt
    if (lgptr->hasFext())
    {
      if (!lgptr->isFextPlugin())
        lgptr->setFextVector(fExt);
    }
    else
      lgptr->setFextPlugin(fExtFunctionName);

    // FInt
    if (lgptr->hasFint())
    {
      if (!lgptr->isFintPlugin())
        if (fInt->size() > 0)
          lgptr->setFintVector(fInt);
        else cout << "Warning : Fint can't be saved, the Fint vector is not defined." << endl;
    }
    else
      lgptr->setFintPlugin(fIntFunctionName);

    // JacobianQFInt
    if (lgptr->hasJacobianQFint())
    {
      if (!lgptr->isJacobianQFintPlugin())
        lgptr->setJacobianQFintMatrix(jacobianQFInt);
    }
    else
      lgptr->setJacobianQFintPlugin(jacobianQFIntFunctionName);

    // JacobianVelocityFInt
    if (lgptr->hasJacobianVelocityFint())
    {
      if (!lgptr->isJacobianVelocityFintPlugin())
        lgptr->setJacobianVelocityFintMatrix(jacobianVelocityFInt);
    }
    else
      lgptr->setJacobianVelocityFintPlugin(jacobianVelocityFIntFunctionName);

    // JacobianQNNL
    if (lgptr->hasJacobianQNNL())
    {
      if (!lgptr->isJacobianQNNLPlugin())
        lgptr->setJacobianQNNLMatrix(jacobianQNNL);
    }
    else
      lgptr->setJacobianQNNLPlugin(jacobianQNNLFunctionName);

    // JacobianVelocityNNLFunction
    if (lgptr->hasJacobianVelocityNNL())
    {
      if (!lgptr->isJacobianVelocityNNLPlugin())
        lgptr->setJacobianVelocityNNLMatrix(jacobianVelocityNNL);
    }
    else
      lgptr->setJacobianVelocityNNLPlugin(jacobianVelocityNNLFunctionName);

    // NNL
    if (lgptr->hasNNL())
    {
      if (!lgptr->isNNLPlugin())
        lgptr->setNNLVector(NNL);
    }
    else
      lgptr->setNNLPlugin(NNLFunctionName);
  }
  else RuntimeException::selfThrow("LagrangianDS::saveDSToXML - object DynamicalSystemXML does not exist");
}

void LagrangianDS::display() const
{
  cout << "===== Lagrangian System display ===== " << endl;
  DynamicalSystem::display();
  cout << "____ data of the LagrangianDS " << endl;
  cout << "| ndof : " << ndof << endl;
  cout << "| q " << endl;
  if (q != NULL) q->display();
  else cout << "-> NULL" << endl;
  cout << "| q0 " << endl;
  if (q0 != NULL) q0->display();
  else cout << "-> NULL" << endl;
  cout << "| qFree " << endl;
  if (qFree != NULL) qFree->display();
  else cout << "-> NULL" << endl;
  cout << "| velocity " << endl;
  if (velocity != NULL) velocity->display();
  else cout << "-> NULL" << endl;
  cout << "| velocity0 " << endl;
  if (velocity0 != NULL) velocity0->display();
  else cout << "-> NULL" << endl;
  cout << "| velocityFree " << endl;
  if (velocityFree != NULL) velocityFree->display();
  else cout << "-> NULL" << endl;
  cout << "| p " << endl;
  if (p != NULL) p->display();
  else cout << "-> NULL" << endl;
  cout << "===================================== " << endl;
}

// --- Functions for memory handling ---
void LagrangianDS::initMemory(const unsigned int& steps)
{
  DynamicalSystem::initMemory(steps);
  // warning : depends on xml loading ??? To be reviewed
  //qMemory = new SiconosMemory( steps, qMemory.getSiconosMemoryXML() );
  //velocityMemory = new SiconosMemory( steps, velocityMemory.getSiconosMemoryXML() );
  if (isQAllocatedIn[3]) delete qMemory;
  qMemory = new SiconosMemory(steps);
  isQAllocatedIn[3] = true;
  if (isVelocityAllocatedIn[3]) delete velocityMemory;
  velocityMemory = new SiconosMemory(steps);
  isVelocityAllocatedIn[3] = true;
  swapInMemory();
}

void LagrangianDS::swapInMemory()
{
  DynamicalSystem::swapInMemory();
  qMemory->swap(*q);
  velocityMemory->swap(*velocity);
  // initialization of the reaction force due to the non smooth law
  p->zero();
}

LagrangianDS* LagrangianDS::convert(DynamicalSystem* ds)
{
  LagrangianDS* lnlds = dynamic_cast<LagrangianDS*>(ds);
  return lnlds;
}

double LagrangianDS::dsConvergenceIndicator()
{
  double dsCvgIndic;
  // Velocity is used to calculate the indicator.
  SimpleVector *diff = new SimpleVector(velocity->size());
  // Compute difference between present and previous Newton steps
  SimpleVector * valRef = tmpWorkVector["LagNLDSMoreau"];
  *diff =  *(getVelocityPtr()) - *valRef;
  if (valRef->norm() != 0)
    dsCvgIndic = diff->norm() / (valRef->norm());
  else
    dsCvgIndic = diff->norm();
  delete diff;
  return (dsCvgIndic);
}

// -- Default constructor --
LagrangianDS::LagrangianDS():
  DynamicalSystem(), ndof(0), q(NULL), q0(NULL), qFree(NULL), qMemory(NULL), velocity(NULL), velocity0(NULL),
  velocityFree(NULL), velocityMemory(NULL), p(NULL), mass(NULL), fInt(NULL), fExt(NULL), NNL(NULL), jacobianQFInt(NULL),
  jacobianVelocityFInt(NULL), jacobianQNNL(NULL), jacobianVelocityNNL(NULL),
  isPAllocatedIn(false), isMassAllocatedIn(false),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL), computeJacobianQFIntPtr(NULL),
  computeJacobianVelocityFIntPtr(NULL), computeJacobianQNNLPtr(NULL), computeJacobianVelocityNNLPtr(NULL)
{
  DSType = LNLDS;
  isQAllocatedIn.resize(4, false);
  isVelocityAllocatedIn.resize(4, false);
  areForcesAllocatedIn.resize(3, false);
  isJacobianAllocatedIn.resize(4, false);
  isLDSPlugin.resize(8, false);
  // parametersList is set to default (vector of size 1 set to 0)
  parametersList.resize(8, NULL);
  vector<SimpleVector*>::iterator iter;
  for (iter = parametersList.begin(); iter != parametersList.end(); iter++)
  {
    *iter = new SimpleVector(1);
    (*iter)->zero();
  }
  isParametersListAllocatedIn.resize(8, true);

  // --- plugins connected to DefaultPlugin ---
  setComputeMassFunction("DefaultPlugin.so", "computeMass");
  setComputeFIntFunction("DefaultPlugin.so", "computeFInt");
  setComputeFExtFunction("DefaultPlugin.so", "computeFExt");
  setComputeNNLFunction("DefaultPlugin.so", "computeNNL");
  setComputeJacobianQFIntFunction("DefaultPlugin.so", "computeJacobianQFInt");
  setComputeJacobianVelocityFIntFunction("DefaultPlugin.so", "computeJacobianVelocityFInt");
  setComputeJacobianQNNLFunction("DefaultPlugin.so", "computeJacobianQNNL");
  setComputeJacobianVelocityNNLFunction("DefaultPlugin.so", "computeJacobianVelocityNNL");
}
