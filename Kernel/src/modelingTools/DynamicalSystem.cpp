/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
#include "DynamicalSystem.h"
#include "DynamicalSystemXML.h"
#include "BlockVector.h"

using namespace std;

void DynamicalSystem::initAllocationFlags(bool in) // default in = true.
{
  isAllocatedIn["xMemory"] = in;
  isAllocatedIn["z"] = in;
  isAllocatedIn["jacobianXRhs"] = in;
  isAllocatedIn["x0"] = in;
  isAllocatedIn["x"] = in;
  isAllocatedIn["rhs"] = in;
  isAllocatedIn["g"] = in;
  isAllocatedIn["jacobianG0"] = in;
  isAllocatedIn["jacobianG1"] = in;
}

void DynamicalSystem::initPluginFlags(bool val)
{
  isPlugin["g"] = val;
  isPlugin["jacobianG0"] = val;
  isPlugin["jacobianG1"] = val;
}

// ===== CONSTRUCTORS =====

// Default constructor (protected)
DynamicalSystem::DynamicalSystem(const string& type):
  DSType(type), number(0), id("none"), nsds(NULL), n(0), x0(NULL), jacobianXRhs(NULL), z(NULL), g(NULL), computeGPtr(NULL),
  xMemory(NULL), stepsInMemory(1), dsxml(NULL)
{
  initAllocationFlags(false);
  x.resize(2, NULL);
}

// From XML file (warning: newNsds is optional, default = NULL)
DynamicalSystem::DynamicalSystem(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* newNsds):
  DSType(dsXML->getType()), number(dsXML->getNumber()), id("none"), nsds(newNsds), n(0), x0(NULL), jacobianXRhs(NULL), z(NULL), g(NULL), computeGPtr(NULL),
  xMemory(NULL), stepsInMemory(1), dsxml(dsXML)
{
  if (dsXML == NULL)  // Not really useful: if NULL, error in init. list of the present constructor ... It is also tested in NSDS xml constructor.
    RuntimeException::selfThrow("DynamicalSystem::DynamicalSystem - DynamicalSystemXML paramater must not be NULL");

  // Only the following data are set in this general constructor:
  //  - DSType
  //  - number and Id
  //  - nsds
  //  - z
  //  - stepsInMemory
  // All others are dependent of the derived class type.

  // === Id ===
  if (dsxml->hasId()) id = dsxml->getId();

  // === Initial conditions ===
  // Warning: n is set thanks to vector of initial conditions size. That will be set in derived classes constructor.

  initAllocationFlags(false); // default values
  x.resize(2, NULL);

  // z - Optional parameter.
  if (dsxml->hasZ())
  {
    z = new SimpleVector(dsxml->getZ());
    isAllocatedIn["z"] = true;
  }

  if (dsxml->hasStepsInMemory()) stepsInMemory = dsxml->getStepsInMemory();
}

// From a minimum set of data
DynamicalSystem::DynamicalSystem(const string& type, int newNumber, unsigned int newN):
  DSType(type), number(newNumber), id("none"), nsds(NULL), n(newN), x0(NULL), jacobianXRhs(NULL), z(NULL), g(NULL), computeGPtr(NULL),
  xMemory(NULL), stepsInMemory(1), dsxml(NULL)
{
  initAllocationFlags(false);
  x.resize(2, NULL);
}

// --- Destructor ---
DynamicalSystem::~DynamicalSystem()
{
  if (isAllocatedIn["x0"])delete x0;
  x0 = NULL ;
  if (isAllocatedIn["x"]) delete x[0];
  x[0] = NULL;
  if (isAllocatedIn["rhs"]) delete x[1];
  x[1] = NULL;
  if (isAllocatedIn["xMemory"]) delete xMemory;
  xMemory = NULL;
  if (isAllocatedIn["jacobianXRhs"]) delete jacobianXRhs;
  jacobianXRhs = NULL;
  if (isAllocatedIn["z"]) delete z;
  z = NULL;
  if (isAllocatedIn["g"]) delete g;
  g = NULL;
  if (isAllocatedIn["jacobianG0"]) delete jacobianG[0];
  if (isAllocatedIn["jacobianG1"]) delete jacobianG[1];
  jacobianG.clear();
  computeGPtr = NULL;
  // clean workMatrix map. Warning: if set/get functions are added for this object,
  // add a isAlloc. deque to check in-class allocation.
  map<string, SiconosMatrix*>::iterator it;
  for (it = workMatrix.begin(); it != workMatrix.end(); ++it)
    if (it->second != NULL)
      delete it->second;
}

bool DynamicalSystem::checkDynamicalSystem()
{
  bool output = true;
  // n
  if (n == 0)
  {
    RuntimeException::selfThrow("DynamicalSystem::checkDynamicalSystem - number of degrees of freedom is equal to 0.");
    output = false;
  }
  // x0 != NULL
  if (x0 == NULL)
  {
    RuntimeException::selfThrow("DynamicalSystem::checkDynamicalSystem - x0 not set.");
    output = false;
  }
  return output;
}

// Setters

void DynamicalSystem::setX0(const SiconosVector& newValue)
{
  // check dimensions ...
  if (newValue.size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setX0 - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (x0 != NULL)
    *x0 = newValue;

  else
  {
    if (newValue.isBlock())
      x0 = new BlockVector(newValue);
    else
      x0 = new SimpleVector(newValue);
    isAllocatedIn["x0"] = true;
  }
}

void DynamicalSystem::setX0Ptr(SiconosVector* newPtr)
{
  // check dimensions ...
  if (newPtr->size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setX0Ptr - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (isAllocatedIn["x0"]) delete x0;
  x0 = newPtr;
  isAllocatedIn["x0"] = false;
}

void DynamicalSystem::setX(const SiconosVector& newValue)
{
  // Warning: this only sets the value of x[0]
  // We suppose that both x and (*x)[0] are properly allocated.

  // check dimensions ...
  if (newValue.size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setX - inconsistent sizes between x input and n - Maybe you forget to set n?");

  if (x[0] == NULL)
  {
    x[0] = new SimpleVector(newValue);
    isAllocatedIn["x"] = true;
  }
  else
    *(x[0]) = newValue;
}

void DynamicalSystem::setXPtr(SiconosVector* newPtr)
{
  // Warning: this only sets the pointer x[0]

  // check dimensions ...
  if (newPtr->size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setXPtr - inconsistent sizes between x input and n - Maybe you forget to set n?");

  if (isAllocatedIn["x"])
    delete x[0];

  x[0] = newPtr;
  isAllocatedIn["x"] = false;
}

void DynamicalSystem::setRhs(const SiconosVector& newValue)
{
  // Warning: this only sets the value of x[1]

  // check dimensions ...
  if (newValue.size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setRhs - inconsistent sizes between x input and n - Maybe you forget to set n?");

  if (x[1] == NULL)
  {
    x[1] = new SimpleVector(newValue);
    isAllocatedIn["rhs"] = true;
  }
  else
    *(x[1]) = newValue;
}

void DynamicalSystem::setRhsPtr(SiconosVector* newPtr)
{
  // Warning: this only sets the pointer (*x)[1]

  // check dimensions ...
  if (newPtr->size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setRhsPtr - inconsistent sizes between x input and n - Maybe you forget to set n?");

  if (isAllocatedIn["rhs"])
    delete x[1];

  x[1] = newPtr;
  isAllocatedIn["rhs"] = false;
}

void DynamicalSystem::setJacobianXRhs(const SiconosMatrix& newValue)
{
  // check dimensions ...
  if (newValue.size(0) != n || newValue.size(1) != n)
    RuntimeException::selfThrow("DynamicalSystem::setJacobianXRhs - inconsistent sizes between jacobianXRhs input and n - Maybe you forget to set n?");

  if (jacobianXRhs != NULL)
    *jacobianXRhs = newValue;

  else
  {
    jacobianXRhs = new SimpleMatrix(newValue);
    isAllocatedIn["jacobianXRhs"] = true;
  }
  isPlugin["jacobianXRhs"] = false;
}

void DynamicalSystem::setJacobianXRhsPtr(SiconosMatrix *newPtr)
{
  // check dimensions ...
  if (newPtr->size(0) != n || newPtr->size(1) != n)
    RuntimeException::selfThrow("DynamicalSystem::setJacobianXRhsPtr - inconsistent sizes between jacobianXRhs input and n - Maybe you forget to set n?");

  if (isAllocatedIn["jacobianXRhs"]) delete jacobianXRhs;
  jacobianXRhs = newPtr;
  isAllocatedIn["jacobianXRhs"] = false;
  isPlugin["jacobianXRhs"] = false;
}

void DynamicalSystem::setZ(const SiconosVector& newValue)
{
  if (z != NULL)
  {
    if (newValue.size() != z->size())
      RuntimeException::selfThrow("DynamicalSystem::setZ - inconsistent sizes between input and existing z - To change z size use setZPtr.");
    *z = newValue;
  }
  else
  {
    if (newValue.isBlock())
      z = new BlockVector(newValue);
    else
      z = new SimpleVector(newValue);
    isAllocatedIn["z"] = true;
  }
}

void DynamicalSystem::setZPtr(SiconosVector* newPtr)
{
  if (isAllocatedIn["z"]) delete z;
  z = newPtr;
  isAllocatedIn["z"] = false;
}

void DynamicalSystem::setG(const SiconosVector& newValue)
{
  // check dimensions ...
  if (newValue.size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setG - inconsistent sizes between g input and n - Maybe you forget to set n?");

  if (g == NULL)
  {
    g = new SimpleVector(newValue);
    isAllocatedIn["g"] = true;
  }
  else
    *g = newValue;
}

void DynamicalSystem::setGPtr(SiconosVector* newPtr)
{
  // check dimensions ...
  if (newPtr->size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setGPtr - inconsistent sizes between g input and n - Maybe you forget to set n?");

  if (isAllocatedIn["g"])
    delete g;

  g = newPtr;
  isAllocatedIn["g"] = false;
}

void DynamicalSystem::setJacobianG(unsigned int i, const SiconosMatrix& newValue)
{
  if (newValue.size(0) != n || newValue.size(1) != n)
    RuntimeException::selfThrow("DynamicalSystem - setJacobianG: inconsistent dimensions with problem size for input matrix JacobianG");

  string name = "jacobianG" + toString<unsigned int>(i);
  if (jacobianG[i] == NULL)
  {
    jacobianG[i] = new SimpleMatrix(newValue);
    isAllocatedIn[name] = true;
  }
  else
    *jacobianG[i] = newValue;
  isPlugin[name] = false;
}

void DynamicalSystem::setJacobianGPtr(unsigned int i, SiconosMatrix *newPtr)
{
  if (newPtr->size(0) != n || newPtr->size(1) != n)
    RuntimeException::selfThrow("DynamicalSystem - setJacobianGPtr: inconsistent input matrix size ");

  string name = "jacobianG" + toString<unsigned int>(i);
  if (isAllocatedIn[name]) delete jacobianG[i];
  jacobianG[i] = newPtr;
  isAllocatedIn[name] = false;
  isPlugin[name] = false;
}

void DynamicalSystem::setXMemory(const SiconosMemory& newValue)
{
  if (xMemory != NULL)
  {
    if (newValue.getMemorySize() != xMemory->getMemorySize())
      RuntimeException::selfThrow("DynamicalSystem::setXMemory - inconsistent sizes between xMemory input and existing memorySize");
    else
      *xMemory = newValue;
  }
  else
  {
    xMemory = new SiconosMemory(newValue);
    isAllocatedIn["xMemory"] = true;
  }
}

void DynamicalSystem::setXMemoryPtr(SiconosMemory * newPtr)
{
  if (isAllocatedIn["xMemory"]) delete xMemory;
  xMemory = newPtr;
  isAllocatedIn["xMemory"] = false;
}

void DynamicalSystem::update(double time)
{
  computeRhs(time);
  computeJacobianXRhs(time);
}

// ===== MEMORY MANAGEMENT FUNCTIONS =====

void DynamicalSystem::initMemory(unsigned int steps)
{
  if (steps == 0)
    cout << "Warning : FirstOrderNonLinearDS::initMemory with size equal to zero" << endl;
  else
  {
    stepsInMemory = steps;
    if (isAllocatedIn["xMemory"]) delete xMemory;
    xMemory = new SiconosMemory(steps);
    isAllocatedIn["xMemory"] = true;
  }
}

void DynamicalSystem::setComputeGFunction(const string& pluginPath, const string& functionName)
{
  if (g == NULL)
  {
    g = new SimpleVector(n);
    isAllocatedIn["g"] = true;
  }

  cShared.setFunction(&computeGPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["g"] = plugin + ":" + functionName;
  isPlugin["g"] = true;
}

void DynamicalSystem::setComputeJacobianGFunction(unsigned int i, const string& pluginPath, const string& functionName)
{
  string name = "jacobianG" + toString<unsigned int>(i);

  if (jacobianG[i] == NULL)
  {
    jacobianG[i] = new SimpleMatrix(n, n);
    isAllocatedIn[name] = true;
  }

  cShared.setFunction(&computeJacobianGPtr[i], pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames[name] = plugin + ":" + functionName;
  isPlugin[name] = true;
}

void DynamicalSystem::computeG(double time)
{
  if (isPlugin["g"])
  {
    if (computeGPtr == NULL)
      RuntimeException::selfThrow("DynamicalSystem::computeG() is not linked to a plugin function");

    computeGPtr(time, n, &(*x[0])(0), &(*x[1])(0), &(*g)(0), z->size(), &(*z)(0));
  }
  else RuntimeException::selfThrow("DynamicalSystem::computeG - Not yet implemented for DS of type " + DSType);
}

void DynamicalSystem::computeJacobianG(unsigned int i, double time)
{
  string name = "jacobianG" + toString<unsigned int>(i);

  if (isPlugin[name])
  {
    if (computeJacobianGPtr[i] == NULL)
      RuntimeException::selfThrow("computeJacobianG(i,time) is not linked to a plugin function. i=" + i);

    (computeJacobianGPtr[i])(time, n, &(*x[0])(0), &(*x[1])(0), &(*jacobianG[i])(0, 0), z->size(), &(*z)(0));
  }
  else RuntimeException::selfThrow("DynamicalSystem::computeJacobianG - Not yet implemented for DS of type " + DSType);
}

// ===== XML MANAGEMENT FUNCTIONS =====

void DynamicalSystem::saveDSToXML()
{
  RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - Not yet tested for DS: do not use it.");

  if (dsxml == NULL)
    RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - The DynamicalSystemXML object doesn't exists");

  // --- Specific (depending on derived class) DS data ---
  saveSpecificDataToXML();

  // --- other data ---
}

// ===== MISCELLANEOUS ====

double DynamicalSystem::dsConvergenceIndicator()
{
  RuntimeException::selfThrow("DynamicalSystem:dsConvergenceIndicator - not yet implemented for Dynamical system type :" + DSType);
  return 0;
}

