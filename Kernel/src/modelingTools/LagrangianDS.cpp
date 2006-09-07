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
#include "LagrangianDS.h"
using namespace std;

// Private function to set linked with members of Dynamical top class
void LagrangianDS::connectToDS()
{
  n = 2 * ndof;

  // x and related vectors
  x = new BlockVector(q, velocity);
  isAllocatedIn["x"] = true;

  x0 = new BlockVector(q0, velocity0);
  isAllocatedIn["x0"] = true;
  xFree = new BlockVector(qFree, velocityFree);
  isAllocatedIn["xFree"] = true;

  // r
  r = new BlockVector(NULL, p[2]);
  isAllocatedIn["r"] = true;
  r->zero();

  // rhs // add conditional allocation ??
  rhs = new BlockVector(velocity, NULL);
  isAllocatedIn["rhs"] = true;

  // jacobianXF

  if (jacobianQFInt != NULL || jacobianQNNL != NULL || jacobianVelocityFInt != NULL || jacobianVelocityNNL != NULL)
  {
    workMatrix["zero-matrix"] = new SimpleMatrix(ndof, ndof);
    workMatrix["Id-matrix"] = new SimpleMatrix(ndof, ndof);
    workMatrix["jacob-block10"] = new SimpleMatrix(ndof, ndof);
    workMatrix["jacob-block11"] = new SimpleMatrix(ndof, ndof);
    workMatrix["Id-matrix"]->eye();
    jacobianXF = new BlockMatrix(workMatrix["zero-matrix"], workMatrix["Id-matrix"], workMatrix["jacob-block10"], workMatrix["jacob-block11"]);
    isAllocatedIn["jacobianXF"] = true;
  }
  else isAllocatedIn["jacobianXF"] = false;

  isPlugin["f"] = false;
  isPlugin["jacobianXF"] = false;
}

void LagrangianDS::initAllocationFlags(const bool in)
{
  if (in) // initialize flag to true for required input data and to false for optional ones
  {
    isAllocatedIn["q0"] = true;
    isAllocatedIn["q"] = true;
    isAllocatedIn["qFree"] = false;
    isAllocatedIn["qMemory"] = false;
    isAllocatedIn["velocity0"] = true;
    isAllocatedIn["velocity"] = true;
    isAllocatedIn["velocityFree"] = false;
    isAllocatedIn["velocityMemory"] = false;
    isAllocatedIn["p0"] = false;
    isAllocatedIn["p1"] = false;
    isAllocatedIn["p2"] = false;
    isAllocatedIn["mass"] = true;
    isAllocatedIn["fInt"] = false;
    isAllocatedIn["fExt"] = false;
    isAllocatedIn["NNL"] = false;
    isAllocatedIn["jacobianQFInt"] = false;
    isAllocatedIn["jacobianVelocityFInt"] = false;
    isAllocatedIn["jacobianQNNL"] = false;
    isAllocatedIn["jacobianVelocityNNL"] = false;
  }
  else
  {
    isAllocatedIn["q0"] = false;
    isAllocatedIn["q"] = false;
    isAllocatedIn["qFree"] = false;
    isAllocatedIn["qMemory"] = false;
    isAllocatedIn["velocity0"] = false;
    isAllocatedIn["velocity"] = false;
    isAllocatedIn["velocityFree"] = false;
    isAllocatedIn["velocityMemory"] = false;
    isAllocatedIn["p0"] = false;
    isAllocatedIn["p1"] = false;
    isAllocatedIn["p2"] = false;
    isAllocatedIn["mass"] = false;
    isAllocatedIn["fInt"] = false;
    isAllocatedIn["fExt"] = false;
    isAllocatedIn["NNL"] = false;
    isAllocatedIn["jacobianQFInt"] = false;
    isAllocatedIn["jacobianVelocityFInt"] = false;
    isAllocatedIn["jacobianQNNL"] = false;
    isAllocatedIn["jacobianVelocityNNL"] = false;
  }
}

void LagrangianDS::initPluginFlags(const bool val)
{
  isPlugin["mass"] = val;
  isPlugin["fInt"] = val;
  isPlugin["fExt"] = val;
  isPlugin["NNL"] = val;
  isPlugin["jacobianQFInt"] = val;
  isPlugin["jacobianVelocityFInt"] = val;
  isPlugin["jacobianQNNL"] = val;
  isPlugin["jacobianVelocityNNL"] = val;
}

// -- Default constructor --
LagrangianDS::LagrangianDS():
  DynamicalSystem(), ndof(0), q(NULL), q0(NULL), qFree(NULL), qMemory(NULL), velocity(NULL), velocity0(NULL),
  velocityFree(NULL), velocityMemory(NULL), mass(NULL), fInt(NULL), fExt(NULL), NNL(NULL), jacobianQFInt(NULL),
  jacobianVelocityFInt(NULL), jacobianQNNL(NULL), jacobianVelocityNNL(NULL),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL), computeJacobianQFIntPtr(NULL),
  computeJacobianVelocityFIntPtr(NULL), computeJacobianQNNLPtr(NULL), computeJacobianVelocityNNLPtr(NULL)
{
  DSType = LNLDS;
  initAllocationFlags(false);
  initPluginFlags(false);
  p.resize(3, NULL);
  // !!! No plug-in connection !!!
}

// --- Constructor from an xml file ---
LagrangianDS::LagrangianDS(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* newNsds):
  DynamicalSystem(), ndof(0), q(NULL), q0(NULL), qFree(NULL), qMemory(NULL), velocity(NULL), velocity0(NULL),
  velocityFree(NULL), velocityMemory(NULL), mass(NULL), fInt(NULL), fExt(NULL),
  NNL(NULL), jacobianQFInt(NULL), jacobianVelocityFInt(NULL), jacobianQNNL(NULL), jacobianVelocityNNL(NULL),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL), computeJacobianQFIntPtr(NULL),
  computeJacobianVelocityFIntPtr(NULL), computeJacobianQNNLPtr(NULL), computeJacobianVelocityNNLPtr(NULL)
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

  initAllocationFlags(); // Default configuration

  // Memories
  if (lgptr->hasQMemory())   // qMemory
  {
    qMemory = new SiconosMemory(lgptr->getQMemoryXML());
    isAllocatedIn["qMemory"] = true;
  }

  if (lgptr->hasVelocityMemory())   // velocityMemory
  {
    velocityMemory = new SiconosMemory(lgptr->getVelocityMemoryXML());
    isAllocatedIn["velocityMemory"] = true;
  }

  string plugin;
  initPluginFlags(false);

  // mass
  mass = new SimpleMatrix(ndof, ndof);
  if (lgptr->isMPlugin()) // if mass is plugged
  {
    plugin = lgptr->getMPlugin();
    setComputeMassFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
  }
  else *mass = lgptr->getMMatrix();

  // === Optional inputs ===

  // fInt
  if (lgptr->hasFint())  // if fInt is given
  {
    fInt = new SimpleVector(ndof);
    isAllocatedIn["fInt"] = true;
    if (lgptr->isFintPlugin()) // if fInt is plugged
    {
      plugin = lgptr->getFintPlugin();
      setComputeFIntFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else *fInt = lgptr->getFintVector();
  }

  // fExt
  if (lgptr->hasFext())  // if fExt is given
  {
    fExt = new SimpleVector(ndof);
    isAllocatedIn["fExt"] = true;
    if (lgptr->isFextPlugin())// if fExt is plugged
    {
      plugin = lgptr->getFextPlugin();
      setComputeFExtFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else *fExt = lgptr->getFextVector();
  }

  // NNL
  if (lgptr ->hasNNL())// if NNL is given
  {
    NNL = new SimpleVector(ndof);
    isAllocatedIn["NNL"] = true;
    if (lgptr->isNNLPlugin())// if NNL is plugged
    {
      plugin = lgptr->getNNLPlugin();
      setComputeNNLFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else *NNL = lgptr->getNNLVector();
  }

  // jacobian of fInt according to q
  if (lgptr ->hasJacobianQFint())// if given
  {
    jacobianQFInt = new SimpleMatrix(ndof, ndof);
    isAllocatedIn["jacobianQFInt"] = true;
    if (lgptr->isJacobianQFintPlugin())// if is plugged
    {
      plugin = lgptr->getJacobianQFintPlugin();
      setComputeJacobianQFIntFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else *jacobianQFInt = lgptr->getJacobianQFintMatrix();
  }

  // jacobian of fInt according to velocity
  if (lgptr ->hasJacobianVelocityFint())// if given
  {
    jacobianVelocityFInt = new SimpleMatrix(ndof, ndof);
    isAllocatedIn["jacobianVelocityFInt"] = true;
    if (lgptr->isJacobianVelocityFintPlugin())// if is plugged
    {
      plugin = lgptr->getJacobianVelocityFintPlugin();
      setComputeJacobianVelocityFIntFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else *jacobianVelocityFInt = lgptr->getJacobianVelocityFintMatrix();
  }

  // jacobian of NNL according to q
  if (lgptr -> hasJacobianQNNL()) // if given
  {
    jacobianQNNL = new SimpleMatrix(ndof, ndof);
    isAllocatedIn["jacobianQNNL"] = true;
    if (lgptr->isJacobianQNNLPlugin())// if is plugged
    {
      plugin = lgptr->getJacobianQNNLPlugin();
      setComputeJacobianQNNLFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else *jacobianQNNL = lgptr->getJacobianQNNLMatrix();
  }

  // jacobian of NNL according to velocity
  if (lgptr ->  hasJacobianVelocityNNL())// if given
  {
    jacobianVelocityNNL = new SimpleMatrix(ndof, ndof);
    isAllocatedIn["jacobianVelocityNNL"] = true;
    if (lgptr->isJacobianVelocityNNLPlugin())// if is plugged
    {
      plugin = lgptr->getJacobianVelocityNNLPlugin();
      setComputeJacobianVelocityNNLFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else *jacobianVelocityNNL = lgptr->getJacobianVelocityNNLMatrix();
  }
  p.resize(3, NULL);
}

// From a set of data; Mass filled-in directly from a siconosMatrix
LagrangianDS::LagrangianDS(const int newNumber, const unsigned int newNdof,
                           const SimpleVector& newQ0, const SimpleVector& newVelocity0,
                           const SiconosMatrix& newMass):
  DynamicalSystem(), ndof(newNdof), q(NULL), q0(NULL), qFree(NULL), qMemory(NULL), velocity(NULL), velocity0(NULL),
  velocityFree(NULL), velocityMemory(NULL), mass(NULL), fInt(NULL), fExt(NULL),
  NNL(NULL), jacobianQFInt(NULL), jacobianVelocityFInt(NULL), jacobianQNNL(NULL), jacobianVelocityNNL(NULL),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL), computeJacobianQFIntPtr(NULL),
  computeJacobianVelocityFIntPtr(NULL), computeJacobianQNNLPtr(NULL), computeJacobianVelocityNNLPtr(NULL)
{
  // --- DS BASE-CLASS MEMBERS ---
  // --- Settings and xml load ---
  DSType = LNLDS;
  number = newNumber;

  // --- LAGRANGIAN INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --
  mass = new SimpleMatrix(newMass);
  q0 = new SimpleVector(newQ0);
  q = new SimpleVector(*q0);

  velocity0 = new SimpleVector(newVelocity0);
  velocity = new SimpleVector(*velocity0);

  // set allocation flags: true for required input, false for others
  initAllocationFlags(); // Default

  //  fInt = new SimpleVector(ndof);
  //  fExt = new SimpleVector(ndof);
  //   NNL = new SimpleVector(ndof);
  //   jacobianQFInt = new SimpleMatrix(ndof,ndof);
  //   jacobianVelocityFInt = new SimpleMatrix(ndof,ndof);
  //   jacobianQNNL = new SimpleMatrix(ndof,ndof);
  //   jacobianVelocityNNL = new SimpleMatrix(ndof,ndof);



  setComputeMassFunction("DefaultPlugin.so", "computeMass");
  //   setComputeFIntFunction("DefaultPlugin.so", "computeFInt");
  //   setComputeFExtFunction("DefaultPlugin.so", "computeFExt");
  //   setComputeNNLFunction("DefaultPlugin.so", "computeNNL");
  //   setComputeJacobianQFIntFunction("DefaultPlugin.so", "computeJacobianQFInt");
  //   setComputeJacobianVelocityFIntFunction("DefaultPlugin.so", "computeJacobianVelocityFInt");
  //   setComputeJacobianQNNLFunction("DefaultPlugin.so", "computeJacobianQNNL");
  //   setComputeJacobianVelocityNNLFunction("DefaultPlugin.so", "computeJacobianVelocityNNL");
  p.resize(3, NULL);
}

// From a set of data - Mass loaded from a plugin
LagrangianDS::LagrangianDS(const int newNumber, const unsigned int newNdof,
                           const SimpleVector& newQ0, const SimpleVector& newVelocity0, const string massName):
  DynamicalSystem(), ndof(newNdof), q(NULL), q0(NULL), qFree(NULL), qMemory(NULL), velocity(NULL),
  velocity0(NULL), velocityFree(NULL), velocityMemory(NULL), mass(NULL),
  fInt(NULL), fExt(NULL), NNL(NULL), jacobianQFInt(NULL), jacobianVelocityFInt(NULL),
  jacobianQNNL(NULL), jacobianVelocityNNL(NULL),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL), computeJacobianQFIntPtr(NULL),
  computeJacobianVelocityFIntPtr(NULL), computeJacobianQNNLPtr(NULL), computeJacobianVelocityNNLPtr(NULL)
{
  // --- DS BASE-CLASS MEMBERS ---
  // --- Settings and xml load ---
  DSType = LNLDS;
  number = newNumber;

  // --- LAGRANGIAN INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --
  q0 = new SimpleVector(newQ0);
  q = new SimpleVector(*q0);

  velocity0 = new SimpleVector(newVelocity0);
  velocity = new SimpleVector(*velocity0);

  mass = new SimpleMatrix(ndof, ndof);
  setComputeMassFunction(cShared.getPluginName(massName), cShared.getPluginFunctionName(massName));
  //   fInt = new SimpleVector(ndof);
  //   fExt = new SimpleVector(ndof);

  //   NNL = new SimpleVector(ndof);
  //   jacobianQFInt = new SimpleMatrix(ndof,ndof);
  //   jacobianVelocityFInt = new SimpleMatrix(ndof,ndof);
  //   jacobianQNNL = new SimpleMatrix(ndof,ndof);
  //   jacobianVelocityNNL = new SimpleMatrix(ndof,ndof);

  initAllocationFlags();// default

  //   --- default plug-in ---
  //   setComputeFIntFunction("DefaultPlugin.so", "computeFInt");
  //   setComputeFExtFunction("DefaultPlugin.so", "computeFExt");
  //   setComputeNNLFunction("DefaultPlugin.so", "computeNNL");
  //   setComputeJacobianQFIntFunction("DefaultPlugin.so", "computeJacobianQFInt");
  //   setComputeJacobianVelocityFIntFunction("DefaultPlugin.so", "computeJacobianVelocityFInt");
  //   setComputeJacobianQNNLFunction("DefaultPlugin.so", "computeJacobianQNNL");
  //   setComputeJacobianVelocityNNLFunction("DefaultPlugin.so", "computeJacobianVelocityNNL");
  p.resize(3, NULL);
}

// copy constructor
LagrangianDS::LagrangianDS(const DynamicalSystem & newDS):
  DynamicalSystem(), ndof(0), q(NULL), q0(NULL), qFree(NULL), qMemory(NULL),
  velocity(NULL), velocity0(NULL), velocityFree(NULL), velocityMemory(NULL), mass(NULL),
  fInt(NULL), fExt(NULL), NNL(NULL), jacobianQFInt(NULL), jacobianVelocityFInt(NULL),
  jacobianQNNL(NULL), jacobianVelocityNNL(NULL),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL), computeJacobianQFIntPtr(NULL),
  computeJacobianVelocityFIntPtr(NULL), computeJacobianQNNLPtr(NULL), computeJacobianVelocityNNLPtr(NULL)
{
  if (newDS.getType() != LNLDS)
    RuntimeException::selfThrow("LagrangianDS - copy constructor: try to copy into a LagrangianDS a DS of type: " + newDS.getType());

  DSType = LNLDS;
  id = "copy";
  stepsInMemory = newDS.getStepsInMemory();

  // convert newDS to lagrangianDS and keeps const options
  const LagrangianDS * lnlds = static_cast<const LagrangianDS*>(&newDS);

  ndof = lnlds->getNdof();

  initAllocationFlags();

  mass = new SimpleMatrix(lnlds->getMass());

  q0 = new SimpleVector(lnlds->getQ0());
  q = new SimpleVector(lnlds->getQ());

  if (lnlds->getQMemoryPtr() != NULL)
    qMemory = new SiconosMemory(lnlds->getQMemory());
  else isAllocatedIn["qMemory"] = false;

  velocity0 = new SimpleVector(lnlds->getVelocity0());
  velocity  = new SimpleVector(lnlds->getVelocity0());

  if (lnlds->getVelocityMemoryPtr() != NULL)
    velocityMemory = new SiconosMemory(lnlds->getVelocityMemory());
  else isAllocatedIn["velocityMemory"] = false;

  // Warning: p is not copied !!!

  if (lnlds->getFIntPtr() != NULL)
  {
    fInt = new SimpleVector(lnlds->getFInt());
    isAllocatedIn["fInt"] = true;
  }

  if (lnlds->getFExtPtr() != NULL)
  {
    fExt = new SimpleVector(lnlds->getFExt());
    isAllocatedIn["fExt"] = true;
  }

  if (lnlds->getNNLPtr() != NULL)
  {
    NNL = new SimpleVector(lnlds->getNNL());
    isAllocatedIn["NNL"] = true;
  }

  if (lnlds->getJacobianQFIntPtr() != NULL)
  {
    jacobianQFInt = new SimpleMatrix(lnlds->getJacobianQFInt());
    isAllocatedIn["jacobianQFInt"] = true;
  }

  if (lnlds->getJacobianVelocityFIntPtr() != NULL)
  {
    jacobianVelocityFInt = new SimpleMatrix(lnlds->getJacobianVelocityFInt());
    isAllocatedIn["jacobianVelocityFInt"] = true;
  }

  if (lnlds->getJacobianQNNLPtr() != NULL)
  {
    jacobianQNNL = new SimpleMatrix(lnlds->getJacobianQNNL());
    isAllocatedIn["jacobianQNNL"] = true;
  }

  if (lnlds->getJacobianVelocityNNLPtr() != NULL)
  {
    jacobianVelocityNNL = new SimpleMatrix(lnlds->getJacobianVelocityNNL());
    isAllocatedIn["jacobianVelocityNNL"] = true;
  }

  string pluginPath, functionName;
  isPlugin = lnlds->getIsPlugin();
  setParameters(newDS.getParameters());   // Copy !!

  if (isPlugin["mass"])
  {
    massFunctionName = lnlds -> getMassFunctionName();
    functionName = cShared.getPluginFunctionName(massFunctionName);
    pluginPath  = cShared.getPluginName(massFunctionName);
    setComputeMassFunction(pluginPath, functionName);
  }

  if (isPlugin["fInt"])
  {
    fIntFunctionName = lnlds -> getFIntFunctionName();
    functionName = cShared.getPluginFunctionName(fIntFunctionName);
    pluginPath  = cShared.getPluginName(fIntFunctionName);
    setComputeFIntFunction(pluginPath, functionName);
  }

  if (isPlugin["fExt"])
  {
    fExtFunctionName = lnlds -> getFExtFunctionName();
    functionName = cShared.getPluginFunctionName(fExtFunctionName);
    pluginPath  = cShared.getPluginName(fExtFunctionName);
    setComputeFExtFunction(pluginPath, functionName);
  }

  if (isPlugin["NNL"])
  {
    NNLFunctionName = lnlds -> getNNLFunctionName();
    functionName = cShared.getPluginFunctionName(NNLFunctionName);
    pluginPath  = cShared.getPluginName(NNLFunctionName);
    setComputeNNLFunction(pluginPath, functionName);
  }

  if (isPlugin["jacobianQFInt"])
  {
    jacobianQFIntFunctionName = lnlds -> getJacobianQFIntFunctionName();
    functionName = cShared.getPluginFunctionName(jacobianQFIntFunctionName);
    pluginPath  = cShared.getPluginName(jacobianQFIntFunctionName);
    setComputeJacobianQFIntFunction(pluginPath, functionName);
  }

  if (isPlugin["jacobianVelocityFInt"])
  {
    jacobianVelocityFIntFunctionName = lnlds -> getJacobianVelocityFIntFunctionName();
    functionName = cShared.getPluginFunctionName(jacobianVelocityFIntFunctionName);
    pluginPath  = cShared.getPluginName(jacobianVelocityFIntFunctionName);
    setComputeJacobianVelocityFIntFunction(pluginPath, functionName);
  }

  if (isPlugin["jacobianQNNL"])
  {
    jacobianQNNLFunctionName = lnlds -> getJacobianQNNLFunctionName();
    functionName = cShared.getPluginFunctionName(jacobianQNNLFunctionName);
    pluginPath  = cShared.getPluginName(jacobianQNNLFunctionName);
    setComputeJacobianQNNLFunction(pluginPath, functionName);
  }

  if (isPlugin["jacobianVelocityNNL"])
  {
    jacobianVelocityNNLFunctionName = lnlds -> getJacobianVelocityNNLFunctionName();
    functionName = cShared.getPluginFunctionName(jacobianVelocityNNLFunctionName);
    pluginPath  = cShared.getPluginName(jacobianVelocityNNLFunctionName);
    setComputeJacobianVelocityNNLFunction(pluginPath, functionName);
  }
  p.resize(3, NULL);
}

// Destructor
LagrangianDS::~LagrangianDS()
{
  if (isAllocatedIn["q"]) delete q;
  q = NULL;
  if (isAllocatedIn["q0"]) delete q0;
  q0 = NULL;
  if (isAllocatedIn["qFree"]) delete qFree;
  qFree = NULL;
  if (isAllocatedIn["qMemory"]) delete qMemory;
  qMemory = NULL;
  if (isAllocatedIn["velocity"]) delete velocity ;
  velocity = NULL;
  if (isAllocatedIn["velocity0"])delete velocity0 ;
  velocity0 = NULL;
  if (isAllocatedIn["velocityFree"])delete velocityFree ;
  velocityFree = NULL;
  if (isAllocatedIn["velocityMemory"])delete velocityMemory;
  velocityMemory = NULL;
  if (isAllocatedIn["p0"]) delete p[0] ;
  p[0] = NULL;
  if (isAllocatedIn["p1"]) delete p[1] ;
  p[1] = NULL;
  if (isAllocatedIn["p2"]) delete p[2] ;
  p[2] = NULL;
  if (isAllocatedIn["mass"]) delete mass;
  mass = NULL;
  if (isAllocatedIn["fInt"])delete fInt ;
  fInt = NULL;
  if (isAllocatedIn["fExt"])delete fExt ;
  fExt = NULL;
  if (isAllocatedIn["NNL"])delete NNL ;
  NNL = NULL;
  if (isAllocatedIn["jacobianQFInt"])delete jacobianQFInt  ;
  jacobianQFInt = NULL;
  if (isAllocatedIn["jacobianVelocityFInt"])delete  jacobianVelocityFInt ;
  jacobianVelocityFInt = NULL;
  if (isAllocatedIn["jacobianQNNL"])delete jacobianQNNL ;
  jacobianQNNL = NULL;
  if (isAllocatedIn["jacobianVelocityNNL"])delete jacobianVelocityNNL ;
  jacobianVelocityNNL = NULL;

  // clean workMatrix map. Warning: if set/get functions are added for this object,
  // add a isAlloc. deque to check in-class allocation.
  map<string, SiconosMatrix*>::iterator it;
  for (it = workMatrix.begin(); it != workMatrix.end(); ++it)
    if (it->second != NULL) delete it->second;
}

bool LagrangianDS::checkDynamicalSystem()
{
  bool output = true;
  // ndof
  if (ndof == 0)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - number of degrees of freedom is equal to 0.");
    output = false;
  }

  // q0 and velocity0 != NULL
  if (q0 == NULL || velocity0 == NULL)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - initial conditions are badly set.");
    output = false;
  }

  // Mass
  if (mass == NULL)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - Mass not set.");
    output = false;
  }

  // fInt
  if ((fInt != NULL  && isPlugin["fInt"]) && (jacobianQFInt == NULL || jacobianVelocityFInt == NULL))
    // ie if fInt is defined and not constant => its Jacobian must be defined (but not necessarily plugged)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - You defined fInt but not its Jacobian (according to q and velocity).");
    output = false;
  }

  // NNL
  if ((NNL != NULL  && isPlugin["NNL"]) && (jacobianQNNL == NULL || jacobianVelocityNNL == NULL))
    // ie if NNL is defined and not constant => its Jacobian must be defined (but not necessarily plugged)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - You defined NNL but not its Jacobian (according to q and velocity).");
    output = false;
  }

  if (isPlugin["f"] || isPlugin["jacobianXF"])
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - f and/or its Jacobian can not be plugged for a Lagrangian system.");
    output = false;
  }
  return output;
}

void LagrangianDS::initFreeVectors(const string type)
{
  if (type == "TimeStepping")
  {
    qFree = new SimpleVector(ndof);
    velocityFree = new SimpleVector(ndof);
    isAllocatedIn["qFree"] = true;
    isAllocatedIn["velocityFree"] = true;
  }
  else if (type == "EventDriven")
  {
    qFree = q;
    velocityFree = velocity;
    isAllocatedIn["qFree"] = false;
    isAllocatedIn["velocityFree"] = false;
  }
  else
    RuntimeException::selfThrow("LagrangianDS::initFreeVectors(simulationType) - Unknown simulation type.");

  qFree->zero();
  velocityFree->zero();
}

// TEMPORARY FUNCTION: Must be called before this->initialize
void LagrangianDS::initP(const string simulationType)
{
  if (simulationType == "TimeStepping")
  {
    p[1] = new SimpleVector(ndof);
    isAllocatedIn["p1"] = true;
    p[2] = p[1];
    isAllocatedIn["p2"] = false;
  }
  else if (simulationType == "EventDriven")
  {
    p[1] = new SimpleVector(ndof);
    isAllocatedIn["p1"] = true;
    p[2] = new SimpleVector(ndof);
    isAllocatedIn["p2"] = true;
  }
}

void LagrangianDS::initialize(const string simulationType, const double time, const unsigned int sizeOfMemory)
{
  initFreeVectors(simulationType);

  initP(simulationType);

  // Set variables of top-class DynamicalSystem
  connectToDS(); // note that connection can not be done during constructor call, since user can complete the ds after (add plugin or anything else).
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;

  // set q and velocity to q0 and velocity0
  *q = *q0;
  *velocity = *velocity0;

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  // compute plug-in values for mass, FInt etc ...
  computeMass();

  // rhs
  // note that rhs(0) = velocity, with pointer link, must already be set.
  SiconosVector* vField = rhs->getVectorPtr(1); // Pointer link!
  vField->zero();

  // Compute M-1
  workMatrix["inverseOfMass"] = new SimpleMatrix(*mass);
  workMatrix["inverseOfMass"]->PLUFactorizationInPlace();
  workMatrix["inverseOfMass"]->PLUInverseInPlace();

  bool flag = false;
  if (fExt != NULL)
  {
    computeFExt(time);
    *vField += *fExt;
    flag = true;
  }
  if (fInt != NULL)
  {
    computeFInt(time);
    *vField -= *fInt;
    flag = true;
  }
  if (NNL != NULL)
  {
    computeNNL();
    *vField -= *NNL;
    flag = true;
  }
  if (flag) // else vField = 0
    *vField = *(workMatrix["inverseOfMass"])**vField;
  // todo: use linearSolve to avoid inversion ? Or save M-1 to avoid re-computation. See this when "M" will be added in DS or LDS.

  // jacobianXF
  if (jacobianQFInt != NULL || jacobianQNNL != NULL || jacobianVelocityFInt != NULL || jacobianVelocityNNL != NULL)
  {
    workMatrix["jacob-block10"]->zero();
    workMatrix["jacob-block11"]->zero();
  }

  flag = false;
  if (jacobianQFInt != NULL)
  {
    computeJacobianQFInt(time);
    *workMatrix["jacob-block10"] -= *jacobianQFInt;
    flag = true;
  }
  if (jacobianQNNL != NULL)
  {
    computeJacobianQNNL();
    *workMatrix["jacob-block10"] -= *jacobianQNNL;
    flag = true;
  }
  if (flag)
    *workMatrix["jacob-block10"] = *workMatrix["inverseOfMass"] **workMatrix["jacob-block10"];

  // !!! jacobian of M according to q is not take into account at the time !!!
  flag = false;
  if (jacobianVelocityFInt != NULL)
  {
    computeJacobianVelocityFInt(time);
    *workMatrix["jacob-block11"] -= *jacobianVelocityFInt;
    flag = true;
  }
  if (jacobianVelocityNNL != NULL)
  {
    computeJacobianVelocityNNL();
    *workMatrix["jacob-block11"] -= *jacobianVelocityNNL;
    flag = true;
  }
  if (flag)
    *workMatrix["jacob-block11"] = *workMatrix["inverseOfMass"] * *workMatrix["jacob-block11"];

  // \todo: control terms handling
}

void LagrangianDS::update(const double time)
{

  computeRhs(time, 0);

  // jacobianXF
  if (jacobianQFInt != NULL || jacobianQNNL != NULL || jacobianVelocityFInt != NULL || jacobianVelocityNNL != NULL)
  {
    workMatrix["jacob-block10"]->zero();
    workMatrix["jacob-block11"]->zero();
  }

  bool flag = false;
  if (jacobianQFInt != NULL)
  {
    computeJacobianQFInt(time);
    *workMatrix["jacob-block10"] -= *jacobianQFInt;
    flag = true;
  }
  if (jacobianQNNL != NULL)
  {
    computeJacobianQNNL();
    *workMatrix["jacob-block10"] -= *jacobianQNNL;
    flag = true;
  }
  if (flag)
    *workMatrix["jacob-block10"] = *workMatrix["inverseOfMass"] **workMatrix["jacob-block10"];

  // !!! jacobian of M according to q is not take into account at the time !!!
  flag = false;
  if (jacobianVelocityFInt != NULL)
  {
    computeJacobianVelocityFInt(time);
    *workMatrix["jacob-block11"] -= *jacobianVelocityFInt;
    flag = true;
  }
  if (jacobianVelocityNNL != NULL)
  {
    computeJacobianVelocityNNL();
    *workMatrix["jacob-block11"] -= *jacobianVelocityNNL;
    flag = true;
  }
  if (flag)
    *workMatrix["jacob-block11"] = *workMatrix["inverseOfMass"] * *workMatrix["jacob-block11"];

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
    isAllocatedIn["q"] = true;
  }
  else
    *q = newValue;
}

void LagrangianDS::setQPtr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQPtr: inconsistent input vector size ");

  if (isAllocatedIn["q"]) delete q;
  q = newPtr;
  isAllocatedIn["q"] = false;
}

void LagrangianDS::setQ0(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ0: inconsistent input vector size ");

  if (q0 == NULL)
  {
    q0 = new SimpleVector(newValue);
    isAllocatedIn["q0"] = true;
  }
  else
    *q0 = newValue;
}

void LagrangianDS::setQ0Ptr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ0Ptr: inconsistent input vector size ");

  if (isAllocatedIn["q0"]) delete q0;
  q0 = newPtr;
  isAllocatedIn["q0"] = false;
}

void LagrangianDS::setQFree(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQFree: inconsistent input vector size ");

  if (qFree == NULL)
  {
    qFree = new SimpleVector(newValue);
    isAllocatedIn["qFree"] = true;
  }
  else
    *qFree = newValue;
}

void LagrangianDS::setQFreePtr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQFreePtr: inconsistent input vector size ");

  if (isAllocatedIn["qFree"]) delete qFree;
  qFree = newPtr;
  isAllocatedIn["qFree"] = false;
}

void LagrangianDS::setQMemory(const SiconosMemory& newValue)
{
  if (qMemory == NULL)
  {
    qMemory = new SiconosMemory(newValue);
    isAllocatedIn["qMemory"] = true;
  }
  else
    *qMemory = newValue;
}

void LagrangianDS::setQMemoryPtr(SiconosMemory * newPtr)
{
  if (isAllocatedIn["qMemory"]) delete qMemory;
  qMemory = newPtr;
  isAllocatedIn["qMemory"] = false;
}

void LagrangianDS::setVelocity(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity: inconsistent input vector size ");

  if (velocity == NULL)
  {
    velocity = new SimpleVector(newValue);
    isAllocatedIn["velocity"] = true;
  }
  else
    *velocity = newValue;
}

void LagrangianDS::setVelocityPtr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocityPtr: inconsistent input vector size ");

  if (isAllocatedIn["velocity"]) delete velocity;
  velocity = newPtr;
  isAllocatedIn["velocity"] = false;
}

void LagrangianDS::setVelocity0(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity0: inconsistent input vector size ");

  if (velocity0 == NULL)
  {
    velocity0 = new SimpleVector(newValue);
    isAllocatedIn["velocity0"] = true;
  }
  else
    *velocity0 = newValue;
}

void LagrangianDS::setVelocity0Ptr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity0Ptr: inconsistent input vector size ");

  if (isAllocatedIn["velocity0"]) delete velocity0;
  velocity0 = newPtr;
  isAllocatedIn["velocity0"] = false;
}

void LagrangianDS::setVelocityFree(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocityFree: inconsistent input vector size ");

  if (velocityFree == NULL)
  {
    velocityFree = new SimpleVector(newValue);
    isAllocatedIn["velocityFree"] = true;
  }
  else
    *velocityFree = newValue;
}

void LagrangianDS::setVelocityFreePtr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocityFreePtr: inconsistent input vector size ");

  if (isAllocatedIn["velocityFree"]) delete velocityFree;
  velocityFree = newPtr;
  isAllocatedIn["velocityFree"] = false;
}

SimpleVector* LagrangianDS::getAccelerationPtr() const
{
  return static_cast<SimpleVector*>(rhs->getVectorPtr(1));
}

void LagrangianDS::setVelocityMemory(const SiconosMemory& newValue)
{
  if (velocityMemory == NULL)
  {
    velocityMemory = new SiconosMemory(newValue);
    isAllocatedIn["velocityMemory"] = true;
  }
  else
    *velocityMemory = newValue;
}

void LagrangianDS::setVelocityMemoryPtr(SiconosMemory * newPtr)
{
  if (isAllocatedIn["velocityMemory"]) delete velocityMemory;
  velocityMemory = newPtr;
  isAllocatedIn["velocityMemory"] = false;
}

void LagrangianDS::setP(const SimpleVector& newValue, const unsigned int level)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setP: inconsistent input vector size ");

  if (p[level] == NULL)
  {
    p[level] = new SimpleVector(newValue);
    string stringValue;
    stringstream sstr;
    sstr << level;
    sstr >> stringValue;
    stringValue = "p" + stringValue;
    isAllocatedIn[stringValue] = true;
  }
  else
    *(p[level]) = newValue;
}

void LagrangianDS::setPPtr(SimpleVector *newPtr, const unsigned int level)
{

  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setPPtr: inconsistent input vector size ");

  string stringValue;
  stringstream sstr;
  sstr << level;
  sstr >> stringValue;
  stringValue = "p" + stringValue;

  if (isAllocatedIn[stringValue]) delete p[level];
  p[level] = newPtr;
  isAllocatedIn[stringValue] = false;
}

SiconosMatrix * LagrangianDS::getInverseOfMassPtr()
{
  if (workMatrix.find("inverseOfMass") == workMatrix.end()) // ie if inverse of Mass has not been computed
    RuntimeException::selfThrow("LagrangianDS - getInverseOfMassPtr: you need to compute this matrix first.");

  // Warning: in this function we do not checked that tha matrix is up to date. If M depends on q, it may require a recomputation before the get.
  return workMatrix["inverseOfMass"];
}

void LagrangianDS::setMass(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setMass: inconsistent input matrix size ");

  if (mass == NULL)
  {
    mass = new SimpleMatrix(newValue);
    isAllocatedIn["mass"] = true;
  }
  else
    *mass = newValue;
  isPlugin["mass"] = false;
}

void LagrangianDS::setMassPtr(SiconosMatrix *newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setMassPtr: inconsistent input matrix size ");

  if (isAllocatedIn["mass"]) delete mass;
  mass = newPtr;
  isAllocatedIn["mass"] = false;
  isPlugin["mass"] = false;
}

void LagrangianDS::setFInt(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setFInt: inconsistent dimensions with problem size for input vector FInt");

  if (fInt == NULL)
  {
    fInt = new SimpleVector(newValue);
    isAllocatedIn["fInt"] = true;
  }
  else
    *fInt = newValue;
  isPlugin["fInt"] = false;
}

void LagrangianDS::setFIntPtr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setFIntPtr: inconsistent input matrix size ");

  if (isAllocatedIn["fInt"]) delete fInt;
  fInt = newPtr;
  isAllocatedIn["fInt"] = false;
  isPlugin["fInt"] = false;
}

void LagrangianDS::setFExt(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setFExt: inconsistent dimensions with problem size for input vector FExt");

  if (fExt == NULL)
  {
    fExt = new SimpleVector(newValue);
    isAllocatedIn["fExt"] = true;
  }
  else
    *fExt = newValue;
  isPlugin["fExt"] = false;
}

void LagrangianDS::setFExtPtr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setFIntPtr: inconsistent input matrix size ");
  if (isAllocatedIn["fExt"]) delete fExt;
  fExt = newPtr;
  isAllocatedIn["fExt"] = false;
  isPlugin["fExt"] = false;
}

void LagrangianDS::setNNL(const SimpleVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setNNL: inconsistent dimensions with problem size for input vector NNL");

  if (NNL == NULL)
  {
    NNL = new SimpleVector(ndof);
    isAllocatedIn["NNL"] = true;
  }
  else
    *NNL = newValue;
  isPlugin["NNL"] = false;
}

void LagrangianDS::setNNLPtr(SimpleVector *newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setNNLPtr: inconsistent input matrix size ");

  if (isAllocatedIn["NNL"]) delete NNL;
  NNL = newPtr;
  isAllocatedIn["NNL"] = false;
  isPlugin["NNL"] = false;
}

void LagrangianDS::setJacobianQFInt(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianQFInt: inconsistent dimensions with problem size for input matrix JacobianQFInt");

  if (jacobianQFInt == NULL)
  {
    jacobianQFInt = new SimpleMatrix(newValue);
    isAllocatedIn["jacobianQFInt"] = true;
  }
  else
    *jacobianQFInt = newValue;
  isPlugin["jacobianQFInt"] = false;
}

void LagrangianDS::setJacobianQFIntPtr(SiconosMatrix *newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianQFIntPtr: inconsistent input matrix size ");

  if (isAllocatedIn["jacobianQFInt"]) delete jacobianQFInt;
  jacobianQFInt = newPtr;
  isAllocatedIn["jacobianQFInt"] = false;
  isPlugin["jacobianQFInt"] = false;
}

void LagrangianDS::setJacobianVelocityFInt(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianVelocityFInt: inconsistent dimensions with problem size for input matrix JacobianVelocityFInt");

  if (jacobianVelocityFInt == NULL)
  {
    jacobianVelocityFInt = new SimpleMatrix(newValue);
    isAllocatedIn["jacobianVelocityFInt"] = true;
  }
  else
    *jacobianVelocityFInt = newValue;
  isPlugin["jacobianVelocityFInt"] = false;
}

void LagrangianDS::setJacobianVelocityFIntPtr(SiconosMatrix *newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianVelocityFIntPtr: inconsistent input matrix size ");

  if (isAllocatedIn["jacobianVelocityFInt"]) delete jacobianVelocityFInt;
  jacobianVelocityFInt = newPtr;
  isAllocatedIn["jacobianVelocityFInt"] = false;
  isPlugin["jacobianVelocityFInt"] = false;
}

void LagrangianDS::setJacobianQNNL(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianQNNL: inconsistent dimensions with problem size for input matrix JacobianQNNL");

  if (jacobianQNNL == NULL)
  {
    jacobianQNNL = new SimpleMatrix(newValue);
    isAllocatedIn["jacobianQNNL"] = true;
  }
  else
    *jacobianQNNL = newValue;
  isPlugin["jacobianQNNL"] = false;
}

void LagrangianDS::setJacobianQNNLPtr(SiconosMatrix *newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianQNNLPtr: inconsistent input matrix size ");

  if (isAllocatedIn["jacobianQNNL"]) delete jacobianQNNL;
  jacobianQNNL = newPtr;
  isAllocatedIn["jacobianQNNL"] = false;
  isPlugin["jacobianQNNL"] = false;
}

void  LagrangianDS::setJacobianVelocityNNL(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianVelocityNNL: inconsistent dimensions with problem size for input matrix JacobianVelocityNNL");

  if (jacobianVelocityNNL == NULL)
  {
    jacobianVelocityNNL = new SimpleMatrix(newValue);
    isAllocatedIn["jacobianVelocityNNL"] = true;
  }
  else
    *jacobianVelocityNNL = newValue;
  isPlugin["jacobianVelocityNNL"] = false;
}

void LagrangianDS::setJacobianVelocityNNLPtr(SiconosMatrix *newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianVelocityNNLPtr: inconsistent input matrix size ");

  if (isAllocatedIn["jacobianVelocityNNL"]) delete jacobianVelocityNNL;
  jacobianVelocityNNL = newPtr;
  isAllocatedIn["jacobianVelocityNNL"] = false;
  isPlugin["jacobianVelocityNNL"] = false;
}


// --- Plugins related functions ---
void LagrangianDS::setComputeMassFunction(const string pluginPath, const string functionName)
{
  if (mass == NULL)
  {
    mass = new SimpleMatrix(ndof, ndof);
    isAllocatedIn["mass"] = true;
  }

  cShared.setFunction(&computeMassPtr, pluginPath, functionName);

  initParameter("mass");

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  massFunctionName = plugin + ":" + functionName;
  isPlugin["mass"] = true;
}

void LagrangianDS::setComputeFIntFunction(const string pluginPath, const string functionName)
{
  if (fInt == NULL)
  {
    fInt = new SimpleVector(ndof);
    isAllocatedIn["fInt"] = true;
  }

  cShared.setFunction(&computeFIntPtr, pluginPath, functionName);

  initParameter("fInt");

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  fIntFunctionName = plugin + ":" + functionName;
  isPlugin["fInt"] = true;
}

void LagrangianDS::setComputeFExtFunction(const string pluginPath, const string functionName)
{
  if (fExt == NULL)
  {
    fExt = new SimpleVector(ndof);
    isAllocatedIn["fExt"] = true;
  }

  cShared.setFunction(&computeFExtPtr, pluginPath, functionName);

  initParameter("fExt");

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  fExtFunctionName = plugin + ":" + functionName;
  isPlugin["fExt"] = true;
}

void LagrangianDS::setComputeNNLFunction(const string pluginPath, const string functionName)
{
  if (NNL == NULL)
  {
    NNL = new SimpleVector(ndof);
    isAllocatedIn["NNL"] = true;
  }

  cShared.setFunction(&computeNNLPtr, pluginPath, functionName);

  initParameter("NNL");

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  NNLFunctionName = plugin + ":" + functionName;
  isPlugin["NNL"] = true;
}

void LagrangianDS::setComputeJacobianQFIntFunction(const string pluginPath, const string functionName)
{
  if (jacobianQFInt == NULL)
  {
    jacobianQFInt = new SimpleMatrix(ndof, ndof);
    isAllocatedIn["jacobianQFInt"] = true;
  }

  cShared.setFunction(&computeJacobianQFIntPtr, pluginPath, functionName);

  initParameter("jacobianQFInt");

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  jacobianQFIntFunctionName = plugin + ":" + functionName;
  isPlugin["jacobianQFInt"] = true;
}

void LagrangianDS::setComputeJacobianVelocityFIntFunction(const string pluginPath, const string functionName)
{
  if (jacobianVelocityFInt == NULL)
  {
    jacobianVelocityFInt = new SimpleMatrix(ndof, ndof);
    isAllocatedIn["jacobianVelocityFInt"] = true;
  }

  cShared.setFunction(&computeJacobianVelocityFIntPtr, pluginPath, functionName);

  initParameter("jacobianVelocityFInt");

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  jacobianVelocityFIntFunctionName = plugin + ":" + functionName;
  isPlugin["jacobianVelocityFInt"] = true;
}

void LagrangianDS::setComputeJacobianQNNLFunction(const string pluginPath, const string functionName)
{
  if (jacobianQNNL == NULL)
  {
    jacobianQNNL = new SimpleMatrix(ndof, ndof);
    isAllocatedIn["jacobianQNNL"] = true;
  }

  cShared.setFunction(&computeJacobianQNNLPtr, pluginPath, functionName);

  initParameter("jacobianQNNL");

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  jacobianQNNLFunctionName = plugin + ":" + functionName;
  isPlugin["jacobianQNNL"] = true;
}

void LagrangianDS::setComputeJacobianVelocityNNLFunction(const string pluginPath, const string functionName)
{
  if (jacobianVelocityNNL == NULL)
  {
    jacobianVelocityNNL = new SimpleMatrix(ndof, ndof);
    isAllocatedIn["jacobianVelocityNNL"] = true;
  }

  cShared.setFunction(&computeJacobianVelocityNNLPtr, pluginPath, functionName);

  initParameter("jacobianVelocityNNL");

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  jacobianVelocityNNLFunctionName = plugin + ":" + functionName;
  isPlugin["jacobianVelocityNNL"] = true;
}

void LagrangianDS::computeMass()
{
  if (isPlugin["mass"])
  {
    if (computeMassPtr == NULL)
      RuntimeException::selfThrow("computeMass() is not linked to a plugin function");
    SimpleVector* param = parametersList["mass"];
    computeMassPtr(ndof, &(*q)(0), &(*mass)(0, 0), &(*param)(0));
  }
}

void LagrangianDS::computeMass(SimpleVector *q2)
{
  if (isPlugin["mass"])
  {
    if (computeMassPtr == NULL)
      RuntimeException::selfThrow("computeMass() is not linked to a plugin function");

    SimpleVector* param = parametersList["mass"];
    computeMassPtr(ndof, &(*q2)(0), &(*mass)(0, 0), &(*param)(0));
  }
}

void LagrangianDS::computeFInt(const double time)
{
  if (isPlugin["fInt"])
  {
    if (computeFIntPtr == NULL)
      RuntimeException::selfThrow("computeFInt() is not linked to a plugin function");

    SimpleVector* param = parametersList["fInt"];
    computeFIntPtr(ndof, &time, &(*q)(0), &(*velocity)(0), &(*fInt)(0), &(*param)(0));
  }
}
void LagrangianDS::computeFInt(const double time, SimpleVector *q2, SimpleVector *velocity2)
{
  if (isPlugin["fInt"])
  {
    if (computeFIntPtr == NULL)
      RuntimeException::selfThrow("computeFInt() is not linked to a plugin function");

    SimpleVector* param = parametersList["fInt"];
    computeFIntPtr(ndof, &time, &(*q2)(0), &(*velocity2)(0), &(*fInt)(0), &(*param)(0));
  }
}

void LagrangianDS::computeFExt(const double time)
{
  if (isPlugin["fExt"])
  {
    if (computeFExtPtr == NULL)
      RuntimeException::selfThrow("computeFExt() is not linked to a plugin function");

    SimpleVector* param = parametersList["fExt"];
    computeFExtPtr(ndof, &time, &(*fExt)(0), &(*param)(0));
  }
}

void LagrangianDS::computeNNL()
{
  if (isPlugin["NNL"])
  {
    if (computeNNLPtr == NULL)
      RuntimeException::selfThrow("computeQ() is not linked to a plugin function");

    SimpleVector* param = parametersList["NNL"];
    computeNNLPtr(ndof, &(*q)(0), &(*velocity)(0), &(*NNL)(0), &(*param)(0));
  }
}

void LagrangianDS::computeNNL(SimpleVector *q2, SimpleVector *velocity2)
{
  if (isPlugin["NNL"])
  {
    if (computeNNLPtr == NULL)
      RuntimeException::selfThrow("computeQ() is not linked to a plugin function");

    SimpleVector* param = parametersList["NNL"];
    computeNNLPtr(ndof, &(*q2)(0), &(*velocity2)(0), &(*NNL)(0), &(*param)(0));
  }
}

void LagrangianDS::computeJacobianQFInt(const double time)
{
  if (isPlugin["jacobianQFInt"])
  {
    if (computeJacobianQFIntPtr == NULL)
      RuntimeException::selfThrow("computeJacobianQFInt() is not linked to a plugin function");

    SimpleVector* param = parametersList["jacobianQFInt"];
    computeJacobianQFIntPtr(ndof, &time, &(*q)(0), &(*velocity)(0), &(*jacobianQFInt)(0, 0), &(*param)(0));
  }
}

void LagrangianDS::computeJacobianQFInt(const double time, SimpleVector *q2, SimpleVector *velocity2)
{
  if (isPlugin["jacobianQFInt"])
  {
    if (computeJacobianQFIntPtr == NULL)
      RuntimeException::selfThrow("computeJacobianQFInt() is not linked to a plugin function");

    SimpleVector* param = parametersList["jacobianQFInt"];
    computeJacobianQFIntPtr(ndof, &time, &(*q2)(0), &(*velocity2)(0), &(*jacobianQFInt)(0, 0), &(*param)(0));
  }
}

void LagrangianDS::computeJacobianVelocityFInt(const double  time)
{
  if (isPlugin["jacobianVelocityFInt"])
  {
    if (computeJacobianVelocityFIntPtr == NULL)
      RuntimeException::selfThrow("computeJacobianVelocityFInt() is not linked to a plugin function");

    SimpleVector* param = parametersList["jacobianVelocityFInt"];
    computeJacobianVelocityFIntPtr(ndof, &time, &(*q)(0), &(*velocity)(0), &(*jacobianVelocityFInt)(0, 0), &(*param)(0));
  }
}
void LagrangianDS::computeJacobianVelocityFInt(const double  time, SimpleVector *q2, SimpleVector *velocity2)
{
  if (isPlugin["jacobianVelocityFInt"])
  {
    if (computeJacobianVelocityFIntPtr == NULL)
      RuntimeException::selfThrow("computeJacobianVelocityFInt() is not linked to a plugin function");

    SimpleVector* param = parametersList["jacobianVelocityFInt"];
    computeJacobianVelocityFIntPtr(ndof, &time, &(*q2)(0), &(*velocity2)(0), &(*jacobianVelocityFInt)(0, 0), &(*param)(0));
  }
}

void LagrangianDS::computeJacobianQNNL()
{
  if (isPlugin["jacobianQNNL"])
  {
    if (computeJacobianQNNLPtr == NULL)
      RuntimeException::selfThrow("computeJacobianQNNL() is not linked to a plugin function");

    SimpleVector* param = parametersList["jacobianQNNL"];
    computeJacobianQNNLPtr(ndof, &(*q)(0), &(*velocity)(0), &(*jacobianQNNL)(0, 0), &(*param)(0));
  }
}

void LagrangianDS::computeJacobianQNNL(SimpleVector *q2, SimpleVector *velocity2)
{
  if (isPlugin["jacobianQNNL"])
  {
    if (computeJacobianQNNLPtr == NULL)
      RuntimeException::selfThrow("computeJacobianQNNL() is not linked to a plugin function");

    SimpleVector* param = parametersList["jacobianQNNL"];
    computeJacobianQNNLPtr(ndof, &(*q2)(0), &(*velocity2)(0), &(*jacobianQNNL)(0, 0), &(*param)(0));
  }
}

void LagrangianDS::computeJacobianVelocityNNL()
{
  if (isPlugin["jacobianVelocityNNL"])
  {
    if (computeJacobianVelocityNNLPtr == NULL)
      RuntimeException::selfThrow("computeJacobianVelocityNNL() is not linked to a plugin function");

    SimpleVector* param = parametersList["jacobianVelocityNNL"];
    computeJacobianVelocityNNLPtr(ndof, &(*q)(0), &(*velocity)(0), &(*jacobianVelocityNNL)(0, 0), &(*param)(0));
  }
}
void LagrangianDS::computeJacobianVelocityNNL(SimpleVector *q2, SimpleVector *velocity2)
{
  if (isPlugin["jacobianVelocityNNL"])
  {
    if (computeJacobianVelocityNNLPtr == NULL)
      RuntimeException::selfThrow("computeJacobianVelocityNNL() is not linked to a plugin function");

    SimpleVector* param = parametersList["jacobianVelocityNNL"];
    computeJacobianVelocityNNLPtr(ndof, &(*q2)(0), &(*velocity2)(0), &(*jacobianVelocityNNL)(0, 0), &(*param)(0));
  }
}

void LagrangianDS::computeInverseOfMass()
{
  // if this is the first call to inverse of mass ...
  map<string, SiconosMatrix*>::iterator it = workMatrix.find("inverseOfMass");
  if (it == workMatrix.end())
  {
    workMatrix["inverseOfMass"] = new SimpleMatrix(*mass);
    workMatrix["inverseOfMass"]->PLUFactorizationInPlace();
    workMatrix["inverseOfMass"]->PLUInverseInPlace();
  }
  else if (isPlugin["mass"]) // if inverse has already been compute but need recomputation since M is not constant
  {
    *(workMatrix["inverseOfMass"]) = *mass;
    if (!((workMatrix["inverseOfMass"])->isFactorized()))
      workMatrix["inverseOfMass"]->PLUFactorizationInPlace();
    workMatrix["inverseOfMass"]->PLUInverseInPlace();
  }
}

void LagrangianDS::computeRhs(const double time, const bool isDSup)
{
  // if isDSup == true, this means that there is no need to re-compute mass ...
  // note that rhs(0) = velocity, with pointer link, must already be set.
  SiconosVector* vField = rhs->getVectorPtr(1); // Pointer link!

  *vField = *(p[2]); // Warning: r update is done in Interactions/Relations

  if (!isDSup) // if it is necessary to re-compute mass, FInt ..., ie if they have not been compiled during the present time step
  {
    // compute M and M-1
    computeMass();
    computeInverseOfMass();

    if (fExt != NULL)
    {
      computeFExt(time);
      *vField += *fExt;
    }
    if (fInt != NULL)
    {
      computeFInt(time);
      *vField -= *fInt;
    }
    if (NNL != NULL)
    {
      computeNNL();
      *vField -= *NNL;
    }
  }
  else
  {
    if (fExt != NULL)
      *vField += *fExt;

    if (fInt != NULL)
      *vField -= *fInt;

    if (NNL != NULL)
      *vField -= *NNL;
  }

  *vField = *(workMatrix["inverseOfMass"])**vField ;
  // todo: use linearSolve to avoid inversion ? Or save M-1 to avoid re-computation. See this when "M" will be added in DS or LDS.

}

void LagrangianDS::computeJacobianXRhs(const double time, const bool isDSup)
{

  // jacobian_x of R ???? To do.


  // if isDSup == true, this means that there is no need to re-compute mass ...
  workMatrix["jacob-block10"]->zero();
  workMatrix["jacob-block11"]->zero();

  bool flag = false, flag2 = false;
  if (!isDSup)
  {
    // compute M and M-1
    computeMass();
    computeInverseOfMass();
    // !!! jacobian of M according to q is not taken into account at the time !!!
    if (jacobianQFInt != NULL)
    {
      computeJacobianQFInt(time);
      *workMatrix["jacob-block10"] -= *jacobianQFInt;
      flag = true;
    }
    if (jacobianQNNL != NULL)
    {
      computeJacobianQNNL();
      *workMatrix["jacob-block10"] -= *jacobianQNNL;
      flag = true;
    }
    if (jacobianVelocityFInt != NULL)
    {
      computeJacobianVelocityFInt(time);
      *workMatrix["jacob-block11"] -= *jacobianVelocityFInt;
      flag2 = true;
    }
    if (jacobianVelocityNNL != NULL)
    {
      computeJacobianVelocityNNL();
      *workMatrix["jacob-block11"] -= *jacobianVelocityNNL;
      flag2 = true;
    }
  }
  else
  {
    if (jacobianQFInt != NULL)
    {
      *workMatrix["jacob-block10"] -= *jacobianQFInt;
      flag = true;
    }
    if (jacobianQNNL != NULL)
    {
      *workMatrix["jacob-block10"] -= *jacobianQNNL;
      flag = true;
    }
    if (jacobianVelocityFInt != NULL)
    {
      *workMatrix["jacob-block11"] -= *jacobianVelocityFInt;
      flag2 = true;
    }
    if (jacobianVelocityNNL != NULL)
    {
      *workMatrix["jacob-block11"] -= *jacobianVelocityNNL;
      flag2 = true;
    }
  }

  if (flag) // else workMatrix["jacob-block10"] = 0
    *workMatrix["jacob-block10"] = *workMatrix["inverseOfMass"] * *workMatrix["jacob-block10"];
  if (flag2) // else workMatrix["jacob-block11"] = 0
    *workMatrix["jacob-block11"] = *workMatrix["inverseOfMass"] * *workMatrix["jacob-block11"];
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
    lgptr->setQ(*q);
    lgptr->setQ0(*q0);
    lgptr->setQMemory(*qMemory);
    lgptr->setVelocity(*velocity);
    lgptr->setVelocity0(*velocity0);
    lgptr->setVelocityMemory(*velocityMemory);

    // FExt
    if (lgptr->hasFext())
    {
      if (!lgptr->isFextPlugin())
        lgptr->setFextVector(*fExt);
    }
    else
      lgptr->setFextPlugin(fExtFunctionName);

    // FInt
    if (lgptr->hasFint())
    {
      if (!lgptr->isFintPlugin())
        if (fInt->size() > 0)
          lgptr->setFintVector(*fInt);
        else cout << "Warning : Fint can't be saved, the Fint vector is not defined." << endl;
    }
    else
      lgptr->setFintPlugin(fIntFunctionName);

    // JacobianQFInt
    if (lgptr->hasJacobianQFint())
    {
      if (!lgptr->isJacobianQFintPlugin())
        lgptr->setJacobianQFintMatrix(*jacobianQFInt);
    }
    else
      lgptr->setJacobianQFintPlugin(jacobianQFIntFunctionName);

    // JacobianVelocityFInt
    if (lgptr->hasJacobianVelocityFint())
    {
      if (!lgptr->isJacobianVelocityFintPlugin())
        lgptr->setJacobianVelocityFintMatrix(*jacobianVelocityFInt);
    }
    else
      lgptr->setJacobianVelocityFintPlugin(jacobianVelocityFIntFunctionName);

    // JacobianQNNL
    if (lgptr->hasJacobianQNNL())
    {
      if (!lgptr->isJacobianQNNLPlugin())
        lgptr->setJacobianQNNLMatrix(*jacobianQNNL);
    }
    else
      lgptr->setJacobianQNNLPlugin(jacobianQNNLFunctionName);

    // JacobianVelocityNNLFunction
    if (lgptr->hasJacobianVelocityNNL())
    {
      if (!lgptr->isJacobianVelocityNNLPlugin())
        lgptr->setJacobianVelocityNNLMatrix(*jacobianVelocityNNL);
    }
    else
      lgptr->setJacobianVelocityNNLPlugin(jacobianVelocityNNLFunctionName);

    // NNL
    if (lgptr->hasNNL())
    {
      if (!lgptr->isNNLPlugin())
        lgptr->setNNLVector(*NNL);
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
  if (p[2] != NULL) p[2]->display();
  else cout << "-> NULL" << endl;
  cout << "===================================== " << endl;
}

// --- Functions for memory handling ---
void LagrangianDS::initMemory(const unsigned int steps)
{
  DynamicalSystem::initMemory(steps);

  // warning : depends on xml loading ??? To be reviewed
  //qMemory = new SiconosMemory( steps, qMemory.getSiconosMemoryXML() );
  //velocityMemory = new SiconosMemory( steps, velocityMemory.getSiconosMemoryXML() );
  if (isAllocatedIn["qMemory"]) delete qMemory;
  qMemory = new SiconosMemory(steps);
  isAllocatedIn["qMemory"] = true;
  if (isAllocatedIn["velocityMemory"]) delete velocityMemory;
  velocityMemory = new SiconosMemory(steps);
  isAllocatedIn["velocityMemory"] = true;
  swapInMemory();
}

void LagrangianDS::swapInMemory()
{
  DynamicalSystem::swapInMemory();
  qMemory->swap(q);
  velocityMemory->swap(velocity);
  // initialization of the reaction force due to the non smooth law
  p[1]->zero();
  if (isAllocatedIn["p2"]) p[2]->zero();
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

void LagrangianDS::computeQFree(const double time, const unsigned int level, SiconosVector* qFreeOut)
{
  // to compute qFree, derivative number level. Mainly used in EventDriven to compute second derivative
  // of q for Output y computation.

  if (qFreeOut->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS::computeQFree - Wrong size for output (different from ndof)");


  if (level != 2)
    RuntimeException::selfThrow("LagrangianDS::computeQFree - Only implemented for second derivative.");

  // Warning: we suppose that all other operators are up to date (FInt, FExt ...)

  qFreeOut->zero();
  if (fInt != NULL)
    *qFreeOut += *fInt;
  if (fExt != NULL)
    *qFreeOut -= *fExt;
  if (NNL != NULL)
    *qFreeOut += *NNL;

  *qFreeOut = *workMatrix["inverseOfMass"]**qFreeOut;

}

void LagrangianDS::resetNonSmoothPart()
{
  p[1]->zero();
  if (isAllocatedIn["p2"]) p[2]->zero();
}
