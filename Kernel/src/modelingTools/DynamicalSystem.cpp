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

using namespace std;

void DynamicalSystem::initAllocationFlags(const bool in) // default in = true.
{
  isAllocatedIn["xMemory"] = false;
  isAllocatedIn["rMemory"] = false;
  isAllocatedIn["xFree"] = false;
  isAllocatedIn["u"] = false;
  isAllocatedIn["T"] = false;
  isAllocatedIn["M"] = false;
  isAllocatedIn["f"] = false;
  isAllocatedIn["jacobianXF"] = false;
  if (in) // set default, minimum required, configuration
  {
    isAllocatedIn["x0"] = true;
    isAllocatedIn["x"] = true;
    isAllocatedIn["rhs"] = true; // the right-hand side, \f$ \dot x\f$
    isAllocatedIn["jacobianXRhs"] = true;
    isAllocatedIn["r"] = true;
  }
  else // if not default, set all to false (default, no argument, constructor case)
  {
    isAllocatedIn["x0"] = false;
    isAllocatedIn["x"] = false;
    isAllocatedIn["rhs"] = true;
    isAllocatedIn["jacobianXRhs"] = false;
    isAllocatedIn["r"] = false;
  }
}

void DynamicalSystem::initPluginFlags(const bool val)
{
  isPlugin["f"] = val;
  isPlugin["jacobianXF"] = val;
  isPlugin["u"] = val;
  isPlugin["T"] = val;
}

void DynamicalSystem::initParameter(const std::string id)
{
  if (parametersList[id] == NULL)
  {
    parametersList[id] = new SimpleVector(1);
    string alloc = "parameter_for_" + id;
    isAllocatedIn[alloc] = true;
  }
}

// ===== CONSTRUCTORS =====

// Default constructor (protected)
DynamicalSystem::DynamicalSystem(const string type):
  DSType(type), nsds(NULL), number(0), id("none"), n(0), x0(NULL), xMemory(NULL),
  jacobianXRhs(NULL), xFree(NULL), r(NULL), rMemory(NULL), M(NULL), invM(NULL), f(NULL), jacobianXF(NULL),
  uSize(0), u(NULL), T(NULL), stepsInMemory(1), dsxml(NULL),
  computeFFunctionName("none"), computeJacobianXFFunctionName("none"), computeUFunctionName("none"),
  computeTFunctionName("none"), computeFPtr(NULL), computeJacobianXFPtr(NULL),
  computeUPtr(NULL), computeTPtr(NULL)
{
  initAllocationFlags(false);
  initPluginFlags(false);
  x.resize(2, NULL);
}

// From XML file (warning: newNsds is optional, default = NULL)
DynamicalSystem::DynamicalSystem(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* newNsds):
  DSType(NLDS), nsds(newNsds), number(0), id("none"), n(0), x0(NULL), xMemory(NULL),
  jacobianXRhs(NULL), xFree(NULL), r(NULL), rMemory(NULL), M(NULL), invM(NULL), f(NULL), jacobianXF(NULL),
  uSize(0), u(NULL), T(NULL), stepsInMemory(1), dsxml(dsXML),
  computeFFunctionName("none"), computeJacobianXFFunctionName("none"), computeUFunctionName("none"),
  computeTFunctionName("none"), computeFPtr(NULL), computeJacobianXFPtr(NULL),
  computeUPtr(NULL), computeTPtr(NULL)
{
  if (dsXML == NULL)
    RuntimeException::selfThrow("DynamicalSystem::DynamicalSystem - DynamicalSystemXML paramater must not be NULL");

  // === Id, number ... ===
  number = dsxml->getNumber();
  if (dsxml->hasId()) id = dsxml->getId();

  // === Initial conditions ===
  // Warning: n is set thanks to x0 size
  if (! dsxml->hasX0())
    RuntimeException::selfThrow("DynamicalSystem:: xml constructor, x0 is a required input");

  x0 = new SimpleVector(dsxml->getX0());
  n = x0->size();

  // === Current state (optional input) ===
  x.resize(2, NULL);
  // x is composed of two blocks of size n, (*x)[0] = \f$ x \f$ and (*x)[1]=\f$ \dot x \f$.
  if (dsxml->hasX())
    x[0] = new SimpleVector(dsxml->getX());
  else // (*x)[0] initialize with x0.
    x[0] = new SimpleVector(*x0);

  // build and initialize right-hand side
  x[1] = new SimpleVector(n);

  // (*x)[1] is set to zero during x allocation
  // Memory for its jacobian is allocated below:
  jacobianXRhs = new SimpleMatrix(n, n);

  // r
  r = new SimpleVector(n);

  initAllocationFlags(); // default values

  string plugin;

  // M - Optional parameter supposed to be a SimpleMatrix with xml
  if (dsxml->hasM())
  {
    M = new SimpleMatrix(dsxml->getM());
    isAllocatedIn["M"] = true;
  }

  // f and jacobianXF are required for DynamicalSystem but not for derived class.
  // Then we can not set exception if they are not given. By default, they are connected to defaultPlugin
  // and isPlugin corresponding flag has to be set in the derived class.

  if (dsxml->hasF())
  {
    f = new SimpleVector(n);
    isAllocatedIn["f"] = true;
    if (dsxml->isFPlugin())
    {
      plugin = dsxml->getFPlugin();
      setComputeFFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));     // this set isPlugin["f"] to true.
    }
    else
    {
      if (dsxml->getFVector().size() != n)
        RuntimeException::selfThrow("DynamicalSystem:: xml constructor, f size differs from n!");
      *f = dsxml->getFVector();
      isPlugin["f"] = false;
    }
  }
  else // set f plug-in to default.
    setComputeFFunction("DefaultPlugin.so", "computeF");
  // Warning: previous line sets isPlugin["f"] AND isAllocatedIn["f"] to true.

  if (dsxml->hasJacobianXF())
  {
    if (dsxml->isJacobianXFPlugin())
    {
      jacobianXF = new SimpleMatrix(n, n);
      plugin = dsxml->getJacobianXFPlugin();
      setComputeJacobianXFFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));     // this set isPlugin["jacobianXF"] to true
    }
    else // This means that jacobianXF is constant
    {
      jacobianXF = new SimpleMatrix(dsxml->getJacobianXFMatrix());
      if (jacobianXF->size(0) != n || jacobianXF->size(1) != n)
        RuntimeException::selfThrow("DynamicalSystem:: xml constructor, jacobianXF size differs from n!");
      isPlugin["jacobianXF"] = false;
    }
    isAllocatedIn["jacobianXF"] = true;
  }
  else // set jacobian f plug-in to default.
    setComputeJacobianXFFunction("DefaultPlugin.so", "computeJacobianXF");
  // Warning: previous line sets isPlugin["jacobianXF"] AND isAllocatedIn["jacobianXF"] to true.

  // Memory
  if (dsxml->hasStepsInMemory()) stepsInMemory = dsxml->getStepsInMemory();

  if (dsxml->hasXMemory())
  {
    xMemory = new SiconosMemory(dsxml->getXMemoryXML());
    isAllocatedIn["xMemory"] = true;
  }

  // --- u and T xml loading (optional) ---
  isPlugin["u"] = false;
  isPlugin["T"] = false;

  // if T then u required but u does not need T.

  if (dsxml->hasU())
  {
    // check if uSize has already been loaded
    if (dsxml->hasUSize())
      uSize = dsxml -> getUSize();
    else
      RuntimeException::selfThrow("DynamicalSystem:: xml constructor, uSize is a required input");

    u = new SimpleVector(uSize);
    isAllocatedIn["u"] = true;

    if (dsxml->isUPlugin())
    {
      plugin = dsxml->getUPlugin();
      setComputeUFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      isPlugin["u"] = true;
    }
    else
    {
      if (dsxml->getUVector().size() != uSize)
        RuntimeException::selfThrow("DynamicalSystem:: xml constructor, inconsistent values for uSize and u.");
      *u = dsxml->getUVector();
    }
  }

  if (dsxml->hasT())
  {
    if (u == NULL)
      RuntimeException::selfThrow("DynamicalSystem:: xml constructor, definition of T without u.");

    T = new SimpleMatrix(n, uSize);
    isAllocatedIn["T"] = true;

    // if T plugin, set compute function, else load matrix from xml
    if (dsxml->isTPlugin())
    {
      plugin = dsxml->getTPlugin();
      setComputeTFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      isPlugin["T"] = true;
    }
    else
    {
      if (dsxml->getTMatrix().size(1) != uSize || dsxml->getTMatrix().size(0) != n)
        RuntimeException::selfThrow("DynamicalSystem:: xml constructor, inconsistent values for uSize, n and T.");
      *T = dsxml->getTMatrix();
    }
  }

  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// From a minimum set of data
DynamicalSystem::DynamicalSystem(const int newNumber, const SiconosVector& newX0, const string fPlugin, const string jacobianXFPlugin):
  DSType(NLDS), nsds(NULL), number(newNumber), id("none"), n(newX0.size()), x0(NULL), xMemory(NULL),
  jacobianXRhs(NULL), xFree(NULL), r(NULL), rMemory(NULL), M(NULL), invM(NULL), f(NULL), jacobianXF(NULL), uSize(0), u(NULL),
  T(NULL), stepsInMemory(1), dsxml(NULL),
  computeFFunctionName("none"), computeJacobianXFFunctionName("none"), computeUFunctionName("none"),
  computeTFunctionName("none"), computeFPtr(NULL), computeJacobianXFPtr(NULL),
  computeUPtr(NULL), computeTPtr(NULL)
{

  // === The dimension of the problem is given newX0.size() ===

  // == Initial conditions ==
  x0 = new SimpleVector(newX0);

  // == Current state ==
  // x is composed of two blocks of size n, x[0] = \f$ x \f$ and x[1]=\f$ \dot x \f$.
  // x[0] initialized with x0.
  x.push_back(new SimpleVector(*x0));
  // x[1]
  x.push_back(new SimpleVector(n));

  // Memory for its jacobian is allocated below:
  jacobianXRhs = new SimpleMatrix(n, n);

  // == r ==
  r = new SimpleVector(n);

  initAllocationFlags();

  // == f and its jacobian ==
  // Allocation and link with the plug-in
  f = new SimpleVector(n);
  setComputeFFunction(cShared.getPluginName(fPlugin), cShared.getPluginFunctionName(fPlugin));
  jacobianXF = new SimpleMatrix(n, n);
  setComputeJacobianXFFunction(cShared.getPluginName(jacobianXFPlugin), cShared.getPluginFunctionName(jacobianXFPlugin));
  isAllocatedIn["f"] = true;
  isAllocatedIn["jacobianXF"] = true;
  // u and T are optional
  isPlugin["u"] = false;
  isPlugin["T"] = false;
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// copy constructor
DynamicalSystem::DynamicalSystem(const DynamicalSystem& newDS):
  DSType(NLDS), nsds(newDS.getNonSmoothDynamicalSystemPtr()), number(-2), id("copy"), n(newDS.getN()),
  x0(NULL), xMemory(NULL), jacobianXRhs(NULL),
  xFree(NULL), r(NULL), rMemory(NULL), M(NULL), invM(NULL), f(NULL), jacobianXF(NULL), uSize(newDS.getUSize()), u(NULL), T(NULL),
  stepsInMemory(newDS.getStepsInMemory()), dsxml(NULL),
  computeFFunctionName(newDS.getComputeFFunctionName()), computeJacobianXFFunctionName(newDS. getComputeJacobianXFFunctionName()),
  computeUFunctionName(newDS.getComputeUFunctionName()), computeTFunctionName(newDS.getComputeTFunctionName()),
  computeFPtr(NULL), computeJacobianXFPtr(NULL),
  computeUPtr(NULL), computeTPtr(NULL)
{

  cout << "!!! Warning: Dynamical System copy, do not forget to set id and number for the new system !!! " << endl;

  x0 = new SimpleVector(newDS.getX0());
  x.push_back(new SimpleVector(newDS.getX()));
  x.push_back(new SimpleVector(newDS.getRhs()));

  if (newDS.getJacobianXRhsPtr() != NULL)
    jacobianXRhs = new SimpleMatrix(newDS.getJacobianXRhs());
  else
    jacobianXRhs = new SimpleMatrix(n, n);

  if (newDS.getRPtr() != NULL)
    r = new SimpleVector(newDS.getR());
  else
    r = new SimpleVector(n);

  initAllocationFlags();

  f = new SimpleVector(newDS.getF());
  jacobianXF = new SimpleMatrix(newDS.getJacobianXF());
  isAllocatedIn["f"] = true;
  isAllocatedIn["jacobianXF"] = true;

  if (newDS.getXMemoryPtr() != NULL)
  {
    xMemory = new SiconosMemory(newDS.getXMemory());
    isAllocatedIn["xMemory"] = true;
  }

  if (newDS.getRMemoryPtr() != NULL)
  {
    rMemory = new SiconosMemory(newDS.getRMemory());
    isAllocatedIn["rMemory"] = true;
  }

  if (newDS.getUPtr() != NULL)
  {
    if (newDS.getUPtr()->isBlock())
      u = new BlockVector(newDS.getU());
    else
      u = new SimpleVector(newDS.getU());
    isAllocatedIn["u"] = true;
  }


  if (newDS.getTPtr() != NULL)
  {
    T = new SimpleMatrix(newDS.getT());
    isAllocatedIn["T"] = true ;
  }

  setParameters(newDS.getParameters());   // Copy !!

  string pluginPath, functionName;

  isPlugin = newDS.getIsPlugin();

  if (isPlugin["f"])
  {
    functionName = cShared.getPluginFunctionName(computeFFunctionName);
    pluginPath  = cShared.getPluginName(computeFFunctionName);
    setComputeFFunction(pluginPath, functionName);
  }

  if (isPlugin["jacobianXF"])
  {
    functionName = cShared.getPluginFunctionName(computeJacobianXFFunctionName);
    pluginPath  = cShared.getPluginName(computeJacobianXFFunctionName);
    setComputeJacobianXFFunction(pluginPath, functionName);
  }

  // plugin for u
  if (isPlugin["u"])
  {
    functionName = cShared.getPluginFunctionName(computeUFunctionName);
    pluginPath  = cShared.getPluginName(computeUFunctionName);
    setComputeUFunction(pluginPath, functionName);
  }
  // plugin for T
  if (isPlugin["T"])
  {
    functionName = cShared.getPluginFunctionName(computeTFunctionName);
    pluginPath  = cShared.getPluginName(computeTFunctionName);
    setComputeTFunction(pluginPath, functionName);
  }

  // xml link is not copied.
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
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
  if (isAllocatedIn["xFree"]) delete xFree;
  xFree = NULL;
  if (isAllocatedIn["r"]) delete r;
  r = NULL;
  if (isAllocatedIn["rMemory"]) delete rMemory;
  rMemory = NULL;
  if (isAllocatedIn["M"]) delete M;
  M = NULL;
  if (isAllocatedIn["f"]) delete f;
  f = NULL;
  if (isAllocatedIn["jacobianXF"]) delete jacobianXF;
  jacobianXF = NULL;
  if (isAllocatedIn["u"]) delete u;
  u = NULL;
  if (isAllocatedIn["T"]) delete T;
  T = NULL;

  map<string, SiconosVector*>::iterator it;
  for (it = parametersList.begin(); it != parametersList.end(); ++it)
  {
    string alloc = "parameter_for_" + it->first;
    if (isAllocatedIn[alloc]) delete it->second;
  }
  parametersList.clear();
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

  // f
  if (f == NULL)
  {
    RuntimeException::selfThrow("DynamicalSystem::checkDynamicalSystem - f not set.");
    output = false;
  }

  // jacobianXF
  if (jacobianXF == NULL)
  {
    RuntimeException::selfThrow("DynamicalSystem::checkDynamicalSystem - Jacobian of f according to x is a required input and is not set.");
    output = false;
  }
  return output;
}

void DynamicalSystem::initFreeVectors(const string type)
{
  if (type == "TimeStepping")
  {
    xFree = new SimpleVector(n);
    isAllocatedIn["xFree"] = true;
  }
  else if (type == "EventDriven")
  {
    xFree = x[0];
    isAllocatedIn["xFree"] = false;
  }
  else
    RuntimeException::selfThrow("DynamicalSystem::initFreeVectors(simulationType) - Unknown simulation type.");
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

void DynamicalSystem::setXFree(const SiconosVector& newValue)
{
  // check dimensions ...
  if (newValue.size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setXFree - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (xFree != NULL)
    *xFree = newValue;

  else
  {
    if (newValue.isBlock())
      xFree = new BlockVector(newValue);
    else
      xFree = new SimpleVector(newValue);
    isAllocatedIn["xFree"] = true;
  }
}

void DynamicalSystem::setXFreePtr(SiconosVector* newPtr)
{
  // check dimensions ...
  if (newPtr->size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setXFreePtr - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (isAllocatedIn["xFree"]) delete xFree;
  xFree = newPtr;
  isAllocatedIn["xFree"] = false;
}

void DynamicalSystem::setR(const SiconosVector& newValue)
{
  // check dimensions ...
  if (newValue.size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setR - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (r != NULL)
    *r = newValue;

  else
  {
    r = new SimpleVector(newValue);
    isAllocatedIn["r"] = true;
  }
}

void DynamicalSystem::setRPtr(SiconosVector *newPtr)
{
  // check dimensions ...
  if (newPtr->size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setRPtr - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (isAllocatedIn["r"]) delete r;
  r = newPtr;
  isAllocatedIn["r"] = false;
}

void DynamicalSystem::setRMemory(const SiconosMemory& newValue)
{
  if (rMemory != NULL)
  {
    if (newValue.getMemorySize() != rMemory->getMemorySize())
      RuntimeException::selfThrow("DynamicalSystem::setRMemory - inconsistent sizes between rMemory input and existing memorySize");
    else
      *rMemory = newValue;
  }
  else
  {
    rMemory = new SiconosMemory(newValue);
    isAllocatedIn["rMemory"] = true;
  }
}

void DynamicalSystem::setRMemoryPtr(SiconosMemory * newPtr)
{
  if (isAllocatedIn["rMemory"]) delete rMemory;
  rMemory = newPtr;
  isAllocatedIn["rMemory"] = false;
}

void DynamicalSystem::setM(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != n || newValue.size(1) != n)
    RuntimeException::selfThrow("DynamicalSystem::setM: inconsistent dimensions with problem size for input matrix M.");

  if (M == NULL)
  {
    M = new SimpleMatrix(n, n);
    isAllocatedIn["M"] = true;
  }
  *M = newValue;
}

void DynamicalSystem::setMPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn["M"]) delete M;
  M = newPtr;
  isAllocatedIn["M"] = false;
}

void DynamicalSystem::setInvM(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != n || newValue.size(1) != n)
    RuntimeException::selfThrow("DynamicalSystem::setInvM: inconsistent dimensions with problem size for input matrix.");

  if (invM == NULL)
  {
    invM = new SimpleMatrix(n, n);
    isAllocatedIn["invM"] = true;
  }
  *invM = newValue;
}

void DynamicalSystem::setInvMPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn["invM"]) delete invM;
  invM = newPtr;
  isAllocatedIn["invM"] = false;
}

void DynamicalSystem::setF(const SiconosVector& newValue)
{
  // check dimensions ...
  if (newValue.size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setF - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (f != NULL)
    *f = newValue;

  else
  {
    if (newValue.isBlock())
      f = new BlockVector(newValue);
    else
      f = new SimpleVector(newValue);
    isAllocatedIn["f"] = true;
  }
}

void DynamicalSystem::setFPtr(SiconosVector* newPtr)
{
  // check dimensions ...
  if (newPtr->size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setFPtr - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (isAllocatedIn["f"]) delete f;
  f = newPtr;
  isAllocatedIn["f"] = false;
}

void DynamicalSystem::setJacobianXF(const SiconosMatrix& newValue)
{
  // check dimensions ...
  if (newValue.size(0) != n || newValue.size(1) != n)
    RuntimeException::selfThrow("DynamicalSystem::setJacobianXF - inconsistent sizes between jacobianXF input and n - Maybe you forget to set n?");

  if (jacobianXF != NULL)
    *jacobianXF = newValue;

  else
  {
    jacobianXF = new SimpleMatrix(newValue);
    isAllocatedIn["jacobianXF"] = true;
  }
  isPlugin["jacobianXF"] = false;
}

void DynamicalSystem::setJacobianXFPtr(SiconosMatrix *newPtr)
{
  // check dimensions ...
  if (newPtr->size(0) != n || newPtr->size(1) != n)
    RuntimeException::selfThrow("DynamicalSystem::setJacobianXFPtr - inconsistent sizes between jacobianXF input and n - Maybe you forget to set n?");

  if (isAllocatedIn["jacobianXF"]) delete jacobianXF;
  jacobianXF = newPtr;
  isAllocatedIn["jacobianXF"] = false;
  isPlugin["jacobianXF"] = false;
}

void  DynamicalSystem::setUSize(const unsigned int newUSize)
{
  if (isAllocatedIn["u"]) delete u;
  uSize = newUSize;
  u = new SimpleVector(uSize);
  isAllocatedIn["u"] = true;
  u->zero();
}

// Three steps to set u:
//  - Check if uSize has been given (default value=0 in all constructors)
//  - Allocate memory for u, if necessary
//  - Set value for u
void DynamicalSystem::setU(const SiconosVector& newValue)
{
  if (uSize == 0 || newValue.size() != uSize)
    RuntimeException::selfThrow("DynamicalSystem::setU - inconsistent sizes between u input and uSize - Maybe you forget to set uSize?");
  if (u != NULL)
    *u = newValue;

  else
  {
    if (newValue.isBlock())
      u = new BlockVector(newValue);
    else
      u = new SimpleVector(newValue);
    isAllocatedIn["u"] = true;
  }
  isPlugin["u"] = false;
}

void DynamicalSystem::setUPtr(SiconosVector* newPtr)
{
  if (uSize == 0 || newPtr->size() != uSize)
    RuntimeException::selfThrow("DynamicalSystem::setUPtr - inconsistent sizes between u input and uSize - Maybe you forget to set uSize?");
  // check dimensions ...

  if (isAllocatedIn["u"]) delete u;
  u = newPtr;
  isAllocatedIn["u"] = false;
  isPlugin["u"] = false;
}

void DynamicalSystem::setT(const SiconosMatrix& newValue)
{
  // check dimensions ...
  if (uSize == 0 || newValue.size(1) != uSize || newValue.size(0) != n)
    RuntimeException::selfThrow("DynamicalSystem::setT - inconsistent sizes between T input, uSize and/or n - Maybe you forget to set n or uSize?");

  if (T != NULL)
    *T = newValue;
  else
  {
    T = new SimpleMatrix(newValue);
    isAllocatedIn["T"] = true;
  }
  isPlugin["T"] = false;
}

void DynamicalSystem::setTPtr(SiconosMatrix *newPtr)
{
  // check dimensions ...
  if (uSize == 0 || newPtr->size(1) != uSize || newPtr->size(0) != n)
    RuntimeException::selfThrow("DynamicalSystem::setTPtr - inconsistent sizes between T input, uSize and/or n - Maybe you forget to set n or uSize?");

  if (isAllocatedIn["T"]) delete T;
  T = newPtr;
  isAllocatedIn["T"] = false;
  isPlugin["T"] = false;
}

void DynamicalSystem::initialize(const string& simulationType, double time, unsigned int sizeOfMemory)
{
  initFreeVectors(simulationType);

  // reset x to x0, xFree and r to zero.
  xFree->zero();
  r->zero();
  *(x[0]) = *x0;

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  // compute initial values for f and jacobianXF, initialize right-hand side.
  x[1]->zero();
  computeRhs(time); // this will compute, if required, f, u and T.

  jacobianXRhs->zero();
  computeJacobianXRhs(time);

}

void DynamicalSystem::update(const double time)
{

  // compute initial values for f and jacobianXF, initialize right-hand side.
  x[1]->zero();
  computeRhs(time); // this will compute, if required, f, u and T.

  jacobianXRhs->zero();
  computeJacobianXRhs(time);

}

// ===== MEMORY MANAGEMENT FUNCTIONS =====

void DynamicalSystem::initMemory(const unsigned int steps)
{
  if (steps == 0)
    cout << "Warning : DynamicalSystem::initMemory with size equal to zero" << endl;
  else
  {
    stepsInMemory = steps;
    if (isAllocatedIn["xMemory"]) delete xMemory;
    xMemory = new SiconosMemory(steps);
    isAllocatedIn["xMemory"] = true;

    if (isAllocatedIn["rMemory"]) delete rMemory;
    rMemory = new SiconosMemory(steps);
    isAllocatedIn["rMemory"] = true;
  }
}

void DynamicalSystem::swapInMemory()
{
  xMemory->swap(x[0]);
  rMemory->swap(r);
}

// ===== COMPUTE PLUGINS FUNCTIONS =====

void DynamicalSystem::setComputeFFunction(const string pluginPath, const string functionName)
{
  if (f == NULL)
  {
    f = new SimpleVector(n);
    isAllocatedIn["f"] = true ;
  }

  initParameter("f");

  computeFPtr = NULL;
  cShared.setFunction(&computeFPtr, pluginPath, functionName);
  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  computeFFunctionName = plugin + ":" + functionName;
  isPlugin["f"] = true;
}

void DynamicalSystem::setComputeJacobianXFFunction(const string pluginPath, const string functionName)
{
  if (jacobianXF == NULL)
  {
    jacobianXF = new SimpleMatrix(n, n);
    isAllocatedIn["jacobianXF"] = true ;
  }

  initParameter("jacobianXF");

  computeJacobianXFPtr = NULL;
  cShared.setFunction(&computeJacobianXFPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  computeJacobianXFFunctionName = plugin + ":" + functionName;
  isPlugin["jacobianXF"] = true;
}

void DynamicalSystem::setComputeUFunction(const string pluginPath, const string functionName)
{
  // since u is not allocated by default, memory must be reserved for it
  if (uSize == 0)
    RuntimeException::selfThrow("DynamicalSystem::setComputeUFunction - uSize is equal to 0 - Maybe you forget to set it?");

  if (u == NULL)
  {
    u = new SimpleVector(uSize);
    isAllocatedIn["u"] = true;
    u->zero();
  }

  initParameter("u");

  computeUPtr = NULL;
  cShared.setFunction(&computeUPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  computeUFunctionName = plugin + ":" + functionName;
  isPlugin["u"] = true;
}

void DynamicalSystem::setComputeTFunction(const string pluginPath, const string functionName)
{
  // since T is not allocated by default, memory must be reserved for it
  if (uSize == 0)
    RuntimeException::selfThrow("DynamicalSystem::setComputeUFunction - uSize is equal to 0 - Maybe you forget to set it?");

  if (T == NULL)
  {
    T = new SimpleMatrix(n, uSize);
    isAllocatedIn["T"] = true;
    T->zero();
  }

  initParameter("T");

  computeTPtr = NULL;
  cShared.setFunction(&computeTPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  computeTFunctionName = plugin + ":" + functionName;
  isPlugin["T"] = true;
}

void DynamicalSystem::setParameters(const std::map<string, SiconosVector*>& newMap)
{
  // copy!!

  map<string, SiconosVector*>::const_iterator it;
  for (it = newMap.begin(); it != newMap.end(); ++it)
  {
    parametersList[it->first] = new SimpleVector(*(it->second));
    string alloc = "parameter_for_" + it->first;
    isAllocatedIn[alloc] = true;
  }
}

void DynamicalSystem::setParameter(const SiconosVector& newValue, const string id)
{
  parametersList[id] = new SimpleVector(newValue);
  string alloc = "parameter_for_" + id;
  isAllocatedIn[alloc] = true;
}

void DynamicalSystem::setParameterPtr(SiconosVector *newPtr, const string id)
{
  parametersList[id] = newPtr;
  string alloc = "parameter_for_" + id;
  isAllocatedIn[alloc] = false;
}

void DynamicalSystem::computeF(const double time)
{
  if (isPlugin["f"])
  {
    if (computeFPtr == NULL)
      RuntimeException::selfThrow("DynamicalSystem::computeVF() f is not linked to a plugin function");
    SiconosVector* param = parametersList["f"];
    computeFPtr(n, time, &((*(x[0]))(0)) , &(*f)(0), &(*param)(0));
  }
  // else nothing!
}

void DynamicalSystem::computeJacobianXF(const double time, const bool)
{
  // second argument is useless at the time - Used in derived classes
  if (isPlugin["jacobianXF"])
  {
    if (computeJacobianXFPtr == NULL)
      RuntimeException::selfThrow("DynamicalSystem::computeJacobianXF() is not linked to a plugin function");
    SiconosVector* param = parametersList["jacobianXF"];
    computeJacobianXFPtr(n, time, &((*(x[0]))(0)), &(*jacobianXF)(0, 0), &(*param)(0));
  }
  // else nothing!
}

void DynamicalSystem::computeRhs(const double time, const bool)
{
  // second argument is useless at the time - Used in derived classes

  // compute rhs = M-1*( f + Tu + r ).

  computeF(time);
  *(x[1]) = *f;

  // Mind that at the time f = rhs (pointer equality)

  if (u != NULL)
  {
    if (isPlugin["u"]) // if u is a plug-in function
      computeU(time);
    if (T != NULL)
    {
      if (isPlugin["T"]) // if T is a plug-in function
        computeT();
      *(x[1]) += prod(*T, *u);
    }
    else
      *(x[1]) += *u;
  }

  *(x[1]) += *r; // Warning: r update is done in Interactions/Relations

  if (M != NULL)
  {
    // compute M-1 at the first call of the present function
    if (invM == NULL)
    {
      invM = new SimpleMatrix(*M);
      invM->PLUFactorizationInPlace();
      invM->PLUInverseInPlace();
    }

    *(x[1]) = prod(*invM, *(x[1]));
  }
}

void DynamicalSystem::computeJacobianXRhs(const double time, const bool)
{
  // second argument is useless at the time - Used in derived classes

  // compute jacobian of rhs according to x, = M-1(jacobianXF + jacobianX(T.u))
  // At the time, second term is set to zero.
  computeJacobianXF(time);
  *jacobianXRhs = *jacobianXF;
  // compute M-1 at the first call of the present function, if required (it may have been computed in a previous computeRhs call)
  if (M != NULL)
  {
    // compute M-1 at the first call of the present function
    if (invM == NULL)
    {
      invM = new SimpleMatrix(*M);
      invM->PLUFactorizationInPlace();
      invM->PLUInverseInPlace();
    }

    *jacobianXRhs = prod(*invM, *jacobianXRhs);
  }
}

void DynamicalSystem::computeU(const double time)
{
  if (isPlugin["u"])
  {
    if (computeUPtr == NULL)
      RuntimeException::selfThrow("DynamicalSystem::computeU() is not linked to a plugin function");
    if (u == NULL)
      RuntimeException::selfThrow("DynamicalSystem::computeU(), warning: u = NULL");

    SiconosVector* param = parametersList["u"];
    computeUPtr(uSize, n, time, &((*(x[0]))(0)), &(*u)(0), &(*param)(0));
  }
  // else nothing!
}

void DynamicalSystem::computeU(const double time, SiconosVector* xx)
{
  if (isPlugin["u"])
  {
    if (computeUPtr == NULL)
      RuntimeException::selfThrow("DynamicalSystem::computeU() is not linked to a plugin function");
    if (u == NULL)
      RuntimeException::selfThrow("DynamicalSystem::computeU(), warning: u = NULL");
    SiconosVector* param = parametersList["u"];
    computeUPtr(uSize, n, time, &(*xx)(0), &(*u)(0), &(*param)(0));
  }
  // else nothing!
}

void DynamicalSystem::computeT()
{
  if (isPlugin["T"])
  {
    if (computeTPtr == NULL)
      RuntimeException::selfThrow("DynamicalSystem::computeT() is not linked to a plugin function");
    if (T == NULL)
      RuntimeException::selfThrow("DynamicalSystem::computeT(), warning: T = NULL");
    SiconosVector* param = parametersList["T"];
    computeTPtr(uSize, n, &((*(x[0]))(0)), &(*T)(0, 0), &(*param)(0));
  }
  // else nothing!
}

// ===== XML MANAGEMENT FUNCTIONS =====

void DynamicalSystem::saveDSToXML()
{
  // --- general DS data ---
  saveDSDataToXML();

  // --- other data ---
  if (dsxml != NULL)
  {
    dsxml->setFPlugin(computeFFunctionName);
    if (computeJacobianXFFunctionName != "")
      dsxml->setJacobianXFPlugin(computeJacobianXFFunctionName);
    else
      dsxml->setJacobianXFPlugin("DefaultPlugin:jacobianXF");
  }
  else RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - The DynamicalSystemXML object doesn't exists");
}

// Save data common to each system into the xml file
void DynamicalSystem::saveDSDataToXML()
{
}

// ===== MISCELLANEOUS ====

void DynamicalSystem::display() const
{
  cout << " ===== General dynamical system display =====" << endl;
  cout << "- number : " << number << endl;
  cout << "- id : " << id << endl;
  cout << "- n (size) : " << n << endl;
  cout << "- x " << endl;
  if (x[0] != NULL) x[0]->display();
  else cout << "-> NULL" << endl;
  cout << "- x0 " << endl;
  if (x0 != NULL) x0->display();
  else cout << "-> NULL" << endl;
  cout << "- xFree " << endl;
  if (xFree != NULL) xFree->display();
  else cout << "-> NULL" << endl;
  cout << "- stepsInMemory : " << stepsInMemory << endl;
  cout << "- r " << endl;
  if (r != NULL) r->display();
  else cout << "-> NULL" << endl;
  cout << "- M: " << endl;
  if (M != NULL) M->display();
  else cout << "-> NULL" << endl;
  cout << "- u " << endl;
  if (u != NULL) u->display();
  else cout << "-> NULL" << endl;
  cout << "- T " << endl;
  if (T != NULL) T->display();
  else cout << "-> NULL" << endl;
  map<string, bool>::const_iterator it;
  cout << "Plug-in state (0: unplugged, ie constant, 1: plugged):" << endl;
  for (it = isPlugin.begin(); it != isPlugin.end(); ++it)
    cout << "Object " << it->first << " " << it->second << endl;
  cout << " ============================================" << endl;
}

double DynamicalSystem::dsConvergenceIndicator()
{
  RuntimeException::selfThrow("DynamicalSystem:dsConvergenceIndicator - not yet implemented for Dynamical system type :" + DSType);
  return 0;
}

void DynamicalSystem::resetNonSmoothPart()
{
  r->zero();
}
