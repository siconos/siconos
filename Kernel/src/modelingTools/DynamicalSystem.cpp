/* Siconos-Kernel version 2.0.0, Copyright INRIA 2005-2006.
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
  if (in) // set default, minimum required, configuration
  {
    isAllocatedIn["x0"] = true;
    isAllocatedIn["x"] = true;
    isAllocatedIn["xMemory"] = false;
    isAllocatedIn["rhs"] = true; // always allocated ??
    isAllocatedIn["jacobianXRhs"] = true;
    isAllocatedIn["xFree"] = false;
    isAllocatedIn["r"] = true;
    isAllocatedIn["rMemory"] = false;
    isAllocatedIn["f"] = false;
    isAllocatedIn["jacobianXF"] = false;
    isAllocatedIn["u"] = false;
    isAllocatedIn["T"] = false;
  }
  else // if not default, set all to false
  {
    isAllocatedIn["x0"] = false;
    isAllocatedIn["x"] = false;
    isAllocatedIn["xMemory"] = false;
    isAllocatedIn["rhs"] = false;
    isAllocatedIn["jacobianXRhs"] = false;
    isAllocatedIn["xFree"] = false;
    isAllocatedIn["r"] = false;
    isAllocatedIn["rMemory"] = false;
    isAllocatedIn["f"] = false;
    isAllocatedIn["jacobianXF"] = false;
    isAllocatedIn["u"] = false;
    isAllocatedIn["T"] = false;
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

// Default constructor (private)
DynamicalSystem::DynamicalSystem():
  DSType(NLDS), nsds(NULL), number(0), id("none"), n(0), x0(NULL), x(NULL), xMemory(NULL),
  rhs(NULL), jacobianXRhs(NULL), xFree(NULL), r(NULL), rMemory(NULL), f(NULL), jacobianXF(NULL),
  uSize(0), u(NULL), T(NULL), stepsInMemory(1), dsxml(NULL),
  computeFFunctionName("none"), computeJacobianXFFunctionName("none"), computeUFunctionName("none"),
  computeTFunctionName("none"), computeFPtr(NULL), computeJacobianXFPtr(NULL),
  computeUPtr(NULL), computeTPtr(NULL)
{
  initAllocationFlags(false);
  initPluginFlags(false);
}

// From XML file (warning: newNsds is optional, default = NULL)
DynamicalSystem::DynamicalSystem(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* newNsds):
  DSType(NLDS), nsds(newNsds), number(0), id("none"), n(0), x0(NULL), x(NULL), xMemory(NULL),
  rhs(NULL), jacobianXRhs(NULL), xFree(NULL), r(NULL), rMemory(NULL), f(NULL), jacobianXF(NULL),
  uSize(0), u(NULL), T(NULL), stepsInMemory(1), dsxml(dsXML),
  computeFFunctionName("none"), computeJacobianXFFunctionName("none"), computeUFunctionName("none"),
  computeTFunctionName("none"), computeFPtr(NULL), computeJacobianXFPtr(NULL),
  computeUPtr(NULL), computeTPtr(NULL)
{
  // --- get values in xml file ---
  if (dsXML != NULL)
  {
    number = dsxml->getNumber();
    if (dsxml->hasId()) id = dsxml->getId();
    // n is required only for first order systems (not for Lagrangian ones, replaced by ndof)
    if (dsxml->hasN()) n = dsxml->getN();
    else if (DSType == NLDS || DSType == LDS)
      RuntimeException::selfThrow("DynamicalSystem:: xml constructor, n (problem size) is a required input");

    // --- Memory allocation for vector and matrix members ---

    // x and related vectors
    if (dsxml->hasX0())
      x0 = new SimpleVector(dsxml->getX0());
    else
      RuntimeException::selfThrow("DynamicalSystem:: xml constructor, x0 is a required input");
    if (x0->size() != n)
      RuntimeException::selfThrow("DynamicalSystem:: xml constructor, x size differs from n!");

    if (dsxml->hasX()) x = new SimpleVector(dsxml->getX());
    else x = new SimpleVector(*x0);
    if (x->size() != n)
      RuntimeException::selfThrow("DynamicalSystem:: xml constructor, x size differs from n!");

    // r
    r = new SimpleVector(n);

    string plugin;
    // f and jacobianXF are required for DynamicalSystem but not for derived class.
    // Then we can not set exception if they are not given. By default, they are connected to defaultPlugin
    // and isPlugin corresponding flag has to be set in the derived class.

    // rhs and f   Warning: Only rhs is allocated, f  points to the same memory location!!!
    rhs = new SimpleVector(n); // rhs is always allocated by default.
    if (dsxml->hasF())
    {
      f = rhs ; // POINTER LINK
      if (dsxml->isFPlugin())
      {
        plugin = dsxml->getFPlugin();
        setComputeFFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));     // this set isPlugin["f"] to true.
      }
      else
      {
        if (dsxml->getFVector().size() != n)
          RuntimeException::selfThrow("DynamicalSystem:: xml constructor, f size differs from n!");
        *rhs = dsxml->getFVector(); // Usefull case?
        // warning: at the time, f is saved in rhs vector.
        isPlugin["f"] = false;
      }
    }
    else // set f plug-in to default. Usefull?
    {
      f = rhs ; // POINTER LINK
      setComputeFFunction("DefaultPlugin.so", "computeF"); // this set isPlugin["f"] to true.
    }
    //RuntimeException::selfThrow("DynamicalSystem:: xml constructor,f input is missing.");

    // jacobianXF  Warning: Only jacobianXRhs is allocated, jacobianXF points to the same memory location!!!
    if (dsxml->hasJacobianXF())
    {
      if (dsxml->isJacobianXFPlugin())
      {
        jacobianXRhs = new SimpleMatrix(n, n);
        jacobianXF = jacobianXRhs; // POINTER LINK
        plugin = dsxml->getJacobianXFPlugin();
        setComputeJacobianXFFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));     // this set isPlugin["jacobianXF"] to true
      }
      else // This means that jacobianXF is constant
      {
        jacobianXRhs = new SimpleMatrix(dsxml->getJacobianXFMatrix());
        jacobianXF = jacobianXRhs; // POINTER LINK
        if (jacobianXF->size(0) != n || jacobianXF->size(1) != n)
          RuntimeException::selfThrow("DynamicalSystem:: xml constructor, jacobianXF size differs from n!");
        isPlugin["jacobianXF"] = false;
      }
    }
    else
    {
      jacobianXRhs = new SimpleMatrix(n, n);
      jacobianXF = jacobianXRhs; // POINTER LINK
      setComputeJacobianXFFunction("DefaultPlugin.so", "computeJacobianXF"); // this set isPlugin["jacobianXF"] to true.
    }

    //RuntimeException::selfThrow("DynamicalSystem:: xml constructor,jacobianXF input is missing.");

    initAllocationFlags(); // default values

    // Memory
    isAllocatedIn["rMemory"] = false ; // rMemory
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
  }
  else
    RuntimeException::selfThrow("DynamicalSystem::DynamicalSystem - DynamicalSystemXML paramater must not be NULL");

  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// From a minimum set of data
DynamicalSystem::DynamicalSystem(const int newNumber, const unsigned int newN, const SiconosVector& newX0,
                                 const string fPlugin, const string jacobianXFPlugin):
  DSType(NLDS), nsds(NULL), number(newNumber), id("none"), n(newN), x0(NULL), x(NULL), xMemory(NULL),
  rhs(NULL), jacobianXRhs(NULL), xFree(NULL), r(NULL), rMemory(NULL), f(NULL), jacobianXF(NULL), uSize(0), u(NULL),
  T(NULL), stepsInMemory(1), dsxml(NULL),
  computeFFunctionName("none"), computeJacobianXFFunctionName("none"), computeUFunctionName("none"),
  computeTFunctionName("none"), computeFPtr(NULL), computeJacobianXFPtr(NULL),
  computeUPtr(NULL), computeTPtr(NULL)
{
  // --- Memory allocation ---
  if (n != newX0.size())
    RuntimeException::selfThrow("DynamicalSystem::constructor from data, inconsistent sizes between problem size and x0.");

  // x and related vectors
  x0 = new SimpleVector(newX0);
  x = new SimpleVector(*x0); // x is initialized with x0.

  // r
  r = new SimpleVector(n);

  // right-hand side
  rhs = new SimpleVector(n);
  f = rhs ; // No allocation for f !!! Pointer link.
  setComputeFFunction(cShared.getPluginName(fPlugin), cShared.getPluginFunctionName(fPlugin));

  // jacobianXF
  jacobianXRhs = new SimpleMatrix(n, n);
  jacobianXF = jacobianXRhs; // No allocation !!! Pointer link.
  setComputeJacobianXFFunction(cShared.getPluginName(jacobianXFPlugin), cShared.getPluginFunctionName(jacobianXFPlugin));

  initAllocationFlags();

  // u and T are optional
  isPlugin["u"] = false;
  isPlugin["T"] = false;
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// copy constructor
DynamicalSystem::DynamicalSystem(const DynamicalSystem& newDS):
  DSType(NLDS), nsds(newDS.getNonSmoothDynamicalSystemPtr()), number(-2), id("copy"), n(newDS.getN()),
  x0(NULL), x(NULL), xMemory(NULL), rhs(NULL), jacobianXRhs(NULL),
  xFree(NULL), r(NULL), rMemory(NULL), f(NULL), jacobianXF(NULL), uSize(newDS.getUSize()), u(NULL), T(NULL),
  stepsInMemory(newDS.getStepsInMemory()), dsxml(NULL),
  computeFFunctionName(newDS.getComputeFFunctionName()), computeJacobianXFFunctionName(newDS. getComputeJacobianXFFunctionName()),
  computeUFunctionName(newDS.getComputeUFunctionName()), computeTFunctionName(newDS.getComputeTFunctionName()),
  computeFPtr(NULL), computeJacobianXFPtr(NULL),
  computeUPtr(NULL), computeTPtr(NULL)
{

  cout << "!!! Warning: Dynamical System copy, do not forget to set id and number for the new system !!! " << endl;

  // Only DS or LinearDS or LinearTIDS copy are authorized. No copies from a Lagrangian.
  //
  //   if(newDS.getX0Ptr()->isBlock())
  //     {
  //       x0 = new BlockVector(newDS.getX0());
  //       x = new BlockVector(newDS.getX());
  //       rhs = new BlockVector(newDS.getXDot());
  //       xFree = new BlockVector(newDS.getXFree());
  //     }
  //   else
  //     {
  x0 = new SimpleVector(newDS.getX0());
  x = new SimpleVector(newDS.getX());
  rhs = new SimpleVector(newDS.getRhs());
  f = rhs;  // POINTER LINK
  //    }
  initAllocationFlags();

  if (newDS.getRPtr() != NULL)
    r = new SimpleVector(newDS.getR());
  else
    r = new SimpleVector(n);

  if (newDS.getJacobianXRhsPtr() != NULL)
    jacobianXRhs = new SimpleMatrix(newDS.getJacobianXRhs());
  else
    jacobianXRhs = new SimpleMatrix(n, n);

  jacobianXF = jacobianXRhs; // POINTER LINK

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
  if (isAllocatedIn["x"]) delete x;
  x = NULL;
  if (isAllocatedIn["xMemory"]) delete xMemory;
  xMemory = NULL;
  if (isAllocatedIn["rhs"]) delete rhs;
  rhs = NULL;
  if (isAllocatedIn["jacobianXRhs"]) delete jacobianXRhs;
  jacobianXRhs = NULL;
  if (isAllocatedIn["xFree"]) delete xFree;
  xFree = NULL;
  if (isAllocatedIn["r"]) delete r;
  r = NULL;
  if (isAllocatedIn["rMemory"]) delete rMemory;
  rMemory = NULL;
  if (isAllocatedIn["f"]) delete f;
  f = NULL;
  if (isAllocatedIn["jacobianXF"]) delete jacobianXF;
  jacobianXF = NULL;
  if (isAllocatedIn["u"]) delete u;
  u = NULL;
  if (isAllocatedIn["T"]) delete T;
  T = NULL;

  map<string, SimpleVector*>::iterator it;
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

  // f/rhs
  if (rhs == NULL)
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
    xFree = x;
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
  // check dimensions ...
  if (newValue.size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setX - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (x != NULL)
    *x = newValue;

  else
  {
    if (newValue.isBlock())
      x = new BlockVector(newValue);
    else
      x = new SimpleVector(newValue);
    isAllocatedIn["x"] = true;
  }
}

void DynamicalSystem::setXPtr(SiconosVector* newPtr)
{
  // check dimensions ...
  if (newPtr->size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setXPtr - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (isAllocatedIn["x"]) delete x;
  x = newPtr;
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
  // check dimensions ...
  if (newValue.size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setRhs - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (rhs != NULL)
    *rhs = newValue;

  else
  {
    if (newValue.isBlock())
      rhs = new BlockVector(newValue);
    else
      rhs = new SimpleVector(newValue);
    isAllocatedIn["rhs"] = true;
  }
}

void DynamicalSystem::setRhsPtr(SiconosVector* newPtr)
{
  // check dimensions ...
  if (newPtr->size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setRhsPtr - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (isAllocatedIn["rhs"]) delete rhs;
  rhs = newPtr;
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
  *x = *x0;

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  // compute initial values for f and jacobianXF, initialize right-hand side.
  rhs->zero();
  computeRhs(time); // this will compute, if required, f, u and T.

  jacobianXRhs->zero();
  computeJacobianXRhs(time);

}

void DynamicalSystem::update(const double time)
{

  // compute initial values for f and jacobianXF, initialize right-hand side.
  rhs->zero();
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
  xMemory->swap(x);
  rMemory->swap(r);
}

// ===== COMPUTE PLUGINS FUNCTIONS =====

void DynamicalSystem::setComputeFFunction(const string pluginPath, const string functionName)
{
  if (rhs == NULL) // Warning: f is saved in rhs.
  {
    rhs = new SimpleVector(n);
    isAllocatedIn["rhs"] = true ;
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

void DynamicalSystem::setParameters(const std::map<string, SimpleVector*>& newMap)
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

void DynamicalSystem::setParameter(const SimpleVector& newValue, const string id)
{
  parametersList[id] = new SimpleVector(newValue);
  string alloc = "parameter_for_" + id;
  isAllocatedIn[alloc] = true;
}

void DynamicalSystem::setParameterPtr(SimpleVector *newPtr, const string id)
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
    SimpleVector* param = parametersList["f"];
    computeFPtr(n, time, &(*x)(0) , &(*f)(0), &(*param)(0)); // warning: f is saved in rhs!!
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
    SimpleVector* param = parametersList["jacobianXF"];
    computeJacobianXFPtr(n, time, &(*x)(0), &(*jacobianXF)(0, 0), &(*param)(0));
  }
  // else nothing!
}

void DynamicalSystem::computeRhs(const double time, const bool)
{
  // second argument is useless at the time - Used in derived classes

  // compute rhs = f + Tu.

  computeF(time);

  // Mind that at the time f = rhs (pointer equality)

  if (u != NULL)
  {
    if (isPlugin["u"]) // if u is a plug-in function
      computeU(time);
    if (T != NULL)
    {
      if (isPlugin["T"]) // if T is a plug-in function
        computeT();
      *rhs += prod(*T, *u);
    }
    else
      *rhs += *u;
  }

  *rhs += * r; // Warning: r update is done in Interactions/Relations
}

void DynamicalSystem::computeJacobianXRhs(const double time, const bool)
{
  // second argument is useless at the time - Used in derived classes

  // compute jacobian of rhs according to x, = jacobianXF + jacobianXT.u.

  computeJacobianXF(time); // this compute f and save it into rhs vector.
  // Mind that at the time jacobianXF = jacobianXRhs (pointer equality)

  // Warning: at the time no jacobianXT is taken into account !!!


}

void DynamicalSystem::computeU(const double time)
{
  if (isPlugin["u"])
  {
    if (computeUPtr == NULL)
      RuntimeException::selfThrow("DynamicalSystem::computeU() is not linked to a plugin function");
    if (u == NULL)
      RuntimeException::selfThrow("DynamicalSystem::computeU(), warning: u = NULL");

    SimpleVector* param = parametersList["u"];
    computeUPtr(uSize, n, time, &(*x)(0), &(*u)(0), &(*param)(0));
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
    SimpleVector* param = parametersList["u"];
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
    SimpleVector* param = parametersList["T"];
    computeTPtr(uSize, n, &(*x)(0), &(*T)(0, 0), &(*param)(0));
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
    dsxml->setN(n);
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
  if (x != NULL) x->display();
  else cout << "-> NULL" << endl;
  cout << "- x0 " << endl;
  if (x0 != NULL) x0->display();
  else cout << "-> NULL" << endl;
  cout << "- xFree " << endl;
  if (xFree != NULL) xFree->display();
  else cout << "-> NULL" << endl;
  cout << "- rhs " << endl;
  if (rhs != NULL) rhs->display();
  else cout << "-> NULL" << endl;
  cout << "- stepsInMemory : " << stepsInMemory << endl;
  cout << "- r " << endl;
  if (r != NULL) r->display();
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
