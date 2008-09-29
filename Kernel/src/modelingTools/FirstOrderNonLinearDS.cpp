/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include "FirstOrderNonLinearDS.h"
#include "FirstOrderNonLinearDSXML.h"
#include "BlockVector.h"

using namespace std;

void FirstOrderNonLinearDS::initPluginFlags(bool val)
{
  isPlugin["f"] = val;
  isPlugin["jacobianXF"] = val;
}

// ===== CONSTRUCTORS =====

// Default constructor (protected)
FirstOrderNonLinearDS::FirstOrderNonLinearDS(DSTYPES type):
  DynamicalSystem(FONLDS)
{
  initPluginFlags(false);
}

// From XML file (warning: newNsds is optional, default = NULL)
#ifndef WithSmartPtr
FirstOrderNonLinearDS::FirstOrderNonLinearDS(DynamicalSystemXML * dsXML, NonSmoothSP::DynamicalSystem newNsds):
  DynamicalSystem(dsXML, newNsds), M(NULL), f(NULL), jacobianXF(NULL),
  computeFPtr(NULL), computeJacobianXFPtr(NULL), r(NULL), rMemory(NULL), invM(NULL)
#else
FirstOrderNonLinearDS::FirstOrderNonLinearDS(DynamicalSystemXMLSPtr dsXML, SP::NonSmoothDynamicalSystem newNsds):
  DynamicalSystem(FONLDS), computeFPtr(NULL), computeJacobianXFPtr(NULL)
#endif

{
  // -- FONLDS xml object --
  FirstOrderNonLinearDSXMLSPtr fonlds = boost::static_pointer_cast <FirstOrderNonLinearDSXML>(dsxml);

#ifndef WithSmartPtr
  initAllocationFlags(false);
#endif

  initPluginFlags(false);

  // === Initial conditions ===
  // Warning: n is set thanks to x0 size
  if (! fonlds->hasX0())
    RuntimeException::selfThrow("FirstOrderNonLinearDS:: xml constructor, x0 is a required input");

#ifndef WithSmartPtr
  x0 = new SimpleVector(fonlds->getX0());
  isAllocatedIn["x0"] = true;
#else
  x0.reset(new SimpleVector(fonlds->getX0()));
#endif

  n = x0->size();

  // === Current state (optional input) ===
  // x is composed of two blocks of size n, (*x)[0] = \f$ x \f$ and (*x)[1]=\f$ \dot x \f$.

#ifndef WithSmartPtr
  if (fonlds->hasX())
    x[0] = new SimpleVector(fonlds->getX());
  else // (*x)[0] initialize with x0.
    x[0] = new SimpleVector(*x0);
  isAllocatedIn["x"] = true;
#else
  if (fonlds->hasX())
    x[0].reset(new SimpleVector(fonlds->getX()));
  else // (*x)[0] initialize with x0.
    x[0].reset(new SimpleVector(*x0));
#endif

  // build and initialize right-hand side

#ifndef WithSmartPtr
  x[1] = new SimpleVector(n);
  isAllocatedIn["rhs"] = true;
#else
  x[1].reset(new SimpleVector(n));
#endif

  // r

#ifndef WithSmartPtr
  r = new SimpleVector(n);
  isAllocatedIn["r"] = true;
#else
  r.reset(new SimpleVector(n));
#endif

  string plugin;

  // M - Optional parameter supposed to be a SimpleMatrix with xml
  if (fonlds->hasM())
  {

#ifndef WithSmartPtr
    M = new SimpleMatrix(fonlds->getM());
    isAllocatedIn["M"] = true;
#else
    M.reset(new SimpleMatrix(fonlds->getM()));
#endif

  }

  // f and jacobianXF are required for DynamicalSystem but not for derived class.
  // Then we can not set exception if they are not given.
  if (fonlds->hasF())
  {
    if (fonlds->isFPlugin())
    {
      plugin = fonlds->getFPlugin();
      setComputeFFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));     // this set isPlugin["f"] to true.
    }
    else
    {
      if (fonlds->getFVector().size() != n)
        RuntimeException::selfThrow("FirstOrderNonLinearDS:: xml constructor, f size differs from n!");

#ifndef WithSmartPtr
      f = new SimpleVector(fonlds->getFVector());
      isAllocatedIn["f"] = true;
#else
      f.reset(new SimpleVector(fonlds->getFVector()));
#endif

    }
  }

  if (fonlds->hasJacobianXF())
  {
    if (fonlds->isJacobianXFPlugin())
    {
      plugin = fonlds->getJacobianXFPlugin();
      setComputeJacobianXFFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));     // this set isPlugin["jacobianXF"] to true
    }
    else // This means that jacobianXF is constant
    {
#ifndef WithSmartPtr
      jacobianXF = new SimpleMatrix(fonlds->getJacobianXFMatrix());
#else
      jacobianXF.reset(new SimpleMatrix(fonlds->getJacobianXFMatrix()));
#endif

      if (jacobianXF->size(0) != n || jacobianXF->size(1) != n)
        RuntimeException::selfThrow("FirstOrderNonLinearDS:: xml constructor, jacobianXF size differs from n!");

#ifndef WithSmartPtr
      isAllocatedIn["jacobianXF"] = true;
#endif
    }
  }

  // Memory

  if (fonlds->hasXMemory())
  {

#ifndef WithSmartPtr
    xMemory = new SiconosMemory(fonlds->getXMemoryXML());
    isAllocatedIn["xMemory"] = true;
#else
    xMemory.reset(new SiconosMemory(fonlds->getXMemoryXML()));
#endif

  }

  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// From a minimum set of data
#ifndef WithSmartPtr
FirstOrderNonLinearDS::FirstOrderNonLinearDS(int newNumber, const SiconosVector& newX0, const string& fPlugin, const string& jacobianXFPlugin):
  DynamicalSystem(FONLDS, newNumber, newX0.size()), M(NULL), f(NULL), jacobianXF(NULL),
  computeFPtr(NULL), computeJacobianXFPtr(NULL), r(NULL), rMemory(NULL), invM(NULL)
#else
FirstOrderNonLinearDS::FirstOrderNonLinearDS(int newNumber, const SiconosVector& newX0, const string& fPlugin, const string& jacobianXFPlugin):
  DynamicalSystem(FONLDS, newNumber, newX0.size()), computeFPtr(NULL), computeJacobianXFPtr(NULL)
#endif
{
  // === The dimension of the problem is given newX0.size() ===

#ifndef WithSmartPtr
  initAllocationFlags(false);
#endif
  initPluginFlags(false);

  // == Initial conditions ==

#ifndef WithSmartPtr
  x0 = new SimpleVector(newX0);
  isAllocatedIn["x0"] = true;
#else
  x0.reset(new SimpleVector(newX0));
#endif

  // == Current state ==
  // x is composed of two blocks of size n, x[0] = \f$ x \f$ and x[1]=\f$ \dot x \f$.
  // x[0] initialized with x0.

#ifndef WithSmartPtr
  x[0] = new SimpleVector(*x0);
  // x[1]
  x[1] = new SimpleVector(n);
  isAllocatedIn["x"] = true;
  isAllocatedIn["rhs"] = true;
#else
  x[0].reset(new SimpleVector(*x0));
  x[1].reset(new SimpleVector(n));
#endif

  // == r ==

#ifndef WithSmartPtr
  r = new SimpleVector(n);
  isAllocatedIn["r"] = true;
#else
  r.reset(new SimpleVector(n));
#endif

  // == f and its jacobian ==
  // Allocation and link with the plug-in
  setComputeFFunction(cShared.getPluginName(fPlugin), cShared.getPluginFunctionName(fPlugin));
  setComputeJacobianXFFunction(cShared.getPluginName(jacobianXFPlugin), cShared.getPluginFunctionName(jacobianXFPlugin));
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// --- Destructor ---
FirstOrderNonLinearDS::~FirstOrderNonLinearDS()
{

#ifndef WithSmartPtr
  if (isAllocatedIn["r"]) delete r;
  r = NULL;
  if (isAllocatedIn["rMemory"]) delete rMemory;
  rMemory = NULL;
  if (isAllocatedIn["M"]) delete M;
  M = NULL;
  if (isAllocatedIn["invM"]) delete invM;
  invM = NULL;
  if (isAllocatedIn["f"]) delete f;
  f = NULL;
  if (isAllocatedIn["jacobianXF"]) delete jacobianXF;
  jacobianXF = NULL;
#endif

}

bool FirstOrderNonLinearDS::checkDynamicalSystem()
{
  bool output = true;
  // n
  if (n == 0)
  {
    RuntimeException::selfThrow("FirstOrderNonLinearDS::checkDynamicalSystem - number of degrees of freedom is equal to 0.");
    output = false;
  }
  // x0 != NULL
  if (! x0)
  {
    RuntimeException::selfThrow("FirstOrderNonLinearDS::checkDynamicalSystem - x0 not set.");
    output = false;
  }

  return output;
}

void FirstOrderNonLinearDS::setR(const SiconosVector& newValue)
{
  // check dimensions ...
  if (newValue.size() != n)
    RuntimeException::selfThrow("FirstOrderNonLinearDS::setR - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (r)
    *r = newValue;

  else
  {

#ifndef WithSmartPtr
    r = new SimpleVector(newValue);
    isAllocatedIn["r"] = true;
#else
    r.reset(new SimpleVector(newValue));
#endif

  }
}

void FirstOrderNonLinearDS::setRPtr(SiconosVectorSPtr newPtr)
{
  // check dimensions ...
  if (newPtr->size() != n)
    RuntimeException::selfThrow("FirstOrderNonLinearDS::setRPtr - inconsistent sizes between x0 input and n - Maybe you forget to set n?");


#ifndef WithSmartPtr
  if (isAllocatedIn["r"]) delete r;
  r = newPtr;
  isAllocatedIn["r"] = false;
#else
  r = newPtr;
#endif

}

void FirstOrderNonLinearDS::setRMemory(const SiconosMemory& newValue)
{
  if (rMemory)
  {
    if (newValue.getMemorySize() != rMemory->getMemorySize())
      RuntimeException::selfThrow("FirstOrderNonLinearDS::setRMemory - inconsistent sizes between rMemory input and existing memorySize");
    else
      *rMemory = newValue;
  }
  else
  {

#ifndef WithSmartPtr
    rMemory = new SiconosMemory(newValue);
    isAllocatedIn["rMemory"] = true;
#else
    rMemory.reset(new SiconosMemory(newValue));
#endif

  }
}

void FirstOrderNonLinearDS::setRMemoryPtr(SiconosMemorySPtr newPtr)
{

#ifndef WithSmartPtr
  if (isAllocatedIn["rMemory"]) delete rMemory;
  isAllocatedIn["rMemory"] = false;
#endif
  rMemory = newPtr;


}

void FirstOrderNonLinearDS::setM(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != n || newValue.size(1) != n)
    RuntimeException::selfThrow("FirstOrderNonLinearDS::setM: inconsistent dimensions with problem size for input matrix M.");

  if (! M)
  {

#ifndef WithSmartPtr
    M = new SimpleMatrix(n, n);
    isAllocatedIn["M"] = true;
#else
    M.reset(new SimpleMatrix(n, n));
#endif

  }
  *M = newValue;
}

void FirstOrderNonLinearDS::setMPtr(SiconosMatrixSPtr newPtr)
{

#ifndef WithSmartPtr
  if (isAllocatedIn["M"]) delete M;
  isAllocatedIn["M"] = false;
#endif
  M = newPtr;

}

void FirstOrderNonLinearDS::setInvM(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != n || newValue.size(1) != n)
    RuntimeException::selfThrow("FirstOrderNonLinearDS::setInvM: inconsistent dimensions with problem size for input matrix.");

  if (! invM)
  {

#ifndef WithSmartPtr
    invM = new SimpleMatrix(n, n);
    isAllocatedIn["invM"] = true;
#else
    invM.reset(new SimpleMatrix(n, n));
#endif

  }
  *invM = newValue;
}

void FirstOrderNonLinearDS::setInvMPtr(SiconosMatrixSPtr newPtr)
{

#ifndef WithSmartPtr
  if (isAllocatedIn["invM"]) delete invM;
  isAllocatedIn["invM"] = false;
#endif

  invM = newPtr;

}

void FirstOrderNonLinearDS::setF(const SiconosVector& newValue)
{
  // check dimensions ...
  if (newValue.size() != n)
    RuntimeException::selfThrow("FirstOrderNonLinearDS::setF - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (f)
    *f = newValue;

  else
  {

#ifndef WithSmartPtr
    if (newValue.isBlock())
      f = new BlockVector(newValue);
    else
      f = new SimpleVector(newValue);
    isAllocatedIn["f"] = true;
#else
    if (newValue.isBlock())
      f.reset(new BlockVector(newValue));
    else
      f.reset(new SimpleVector(newValue));
#endif

  }
}

void FirstOrderNonLinearDS::setFPtr(SiconosVectorSPtr newPtr)
{
  // check dimensions ...
  if (newPtr->size() != n)
    RuntimeException::selfThrow("FirstOrderNonLinearDS::setFPtr - inconsistent sizes between x0 input and n - Maybe you forget to set n?");


#ifndef WithSmartPtr
  if (isAllocatedIn["f"]) delete f;
  isAllocatedIn["f"] = false;
#endif

  f = newPtr;

}

void FirstOrderNonLinearDS::setJacobianXF(const SiconosMatrix& newValue)
{
  // check dimensions ...
  if (newValue.size(0) != n || newValue.size(1) != n)
    RuntimeException::selfThrow("FirstOrderNonLinearDS::setJacobianXF - inconsistent sizes between jacobianXF input and n - Maybe you forget to set n?");

  if (jacobianXF)
    *jacobianXF = newValue;

  else
  {

#ifndef WithSmartPtr
    jacobianXF = new SimpleMatrix(newValue);
    isAllocatedIn["jacobianXF"] = true;
#else
    jacobianXF.reset(new SimpleMatrix(newValue));
#endif

  }
  isPlugin["jacobianXF"] = false;
}

void FirstOrderNonLinearDS::setJacobianXFPtr(SiconosMatrixSPtr newPtr)
{
  // check dimensions ...
  if (newPtr->size(0) != n || newPtr->size(1) != n)
    RuntimeException::selfThrow("FirstOrderNonLinearDS::setJacobianXFPtr - inconsistent sizes between jacobianXF input and n - Maybe you forget to set n?");

#ifndef WithSmartPtr
  if (isAllocatedIn["jacobianXF"]) delete jacobianXF;
  isAllocatedIn["jacobianXF"] = false;
#endif
  jacobianXF = newPtr;

  isPlugin["jacobianXF"] = false;
}

void FirstOrderNonLinearDS::initRhs(double time)
{
  // compute initial values for f and jacobianXF, initialize right-hand side.
  computeRhs(time); // this will compute, if required, f.

  if (! jacobianXRhs)  // if not allocated with a set or anything else
  {
    if (jacobianXF && ! M)  // if M is not defined, then jacobianXF = jacobianXRhs, no memory allocation for that one.
      jacobianXRhs = jacobianXF;
    else if (jacobianXF && M)
    {

#ifndef WithSmartPtr
      jacobianXRhs = new SimpleMatrix(n, n);
      isAllocatedIn["jacobianXRhs"] = true;
#else
      jacobianXRhs.reset(new SimpleMatrix(n, n));
#endif

    }
    // else no allocation, jacobian is equal to 0.
  }
  computeJacobianXRhs(time);
}

void FirstOrderNonLinearDS::initialize(const string& simulationType, double time, unsigned int sizeOfMemory)
{
  // reset x to x0 and r to zero.
  r->zero();
  *(x[0]) = *x0;

  // If z is NULL (ie has not been set), we initialize it with a null vector of size 1, since z is required in plug-in functions call.
  if (! z)
  {

#ifndef WithSmartPtr
    z = new SimpleVector(1);
    isAllocatedIn["z"] = true;
#else
    z.reset(new SimpleVector(1));
#endif

  }

  // Initialize memory vectors
  initMemory(sizeOfMemory);
  if (simulationType == "EventDriven")
  {
    // Rhs and its jacobian
    initRhs(time);
  }
}

// ===== MEMORY MANAGEMENT FUNCTIONS =====

void FirstOrderNonLinearDS::initMemory(unsigned int steps)
{
  DynamicalSystem::initMemory(steps);

  if (steps == 0)
    cout << "Warning : FirstOrderNonLinearDS::initMemory with size equal to zero" << endl;
  else
  {

#ifndef WithSmartPtr
    if (isAllocatedIn["rMemory"]) delete rMemory;
    rMemory = new SiconosMemory(steps);
    isAllocatedIn["rMemory"] = true;
#else
    rMemory.reset(new SiconosMemory(steps));
#endif

  }
}

void FirstOrderNonLinearDS::swapInMemory()
{
  xMemory->swap(x[0]);
  rMemory->swap(r);
}

// ===== COMPUTE PLUGINS FUNCTIONS =====

void FirstOrderNonLinearDS::setComputeFFunction(const string& pluginPath, const string& functionName)
{
  if (! f)
  {

#ifndef WithSmartPtr
    f = new SimpleVector(n);
    isAllocatedIn["f"] = true ;
#else
    f.reset(new SimpleVector(n));
#endif

  }

  computeFPtr = NULL;
  cShared.setFunction(&computeFPtr, pluginPath, functionName);
  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["f"] = plugin + ":" + functionName;
  isPlugin["f"] = true;
}

void FirstOrderNonLinearDS::setComputeJacobianXFFunction(const string& pluginPath, const string& functionName)
{
  if (! jacobianXF)
  {

#ifndef WithSmartPtr
    jacobianXF = new SimpleMatrix(n, n);
    isAllocatedIn["jacobianXF"] = true ;
#else
    jacobianXF.reset(new SimpleMatrix(n, n));
#endif

  }

  computeJacobianXFPtr = NULL;
  cShared.setFunction(&computeJacobianXFPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["jacobianXF"] = plugin + ":" + functionName;
  isPlugin["jacobianXF"] = true;
}

void FirstOrderNonLinearDS::computeF(double time)
{
  if (isPlugin["f"])
  {
    if (computeFPtr == NULL)
      RuntimeException::selfThrow("FirstOrderNonLinearDS::computeVF() f is not linked to a plugin function");
    computeFPtr(time, n, &((*(x[0]))(0)) , &(*f)(0), z->size(), &(*z)(0));
  }
  // else nothing!
}

void FirstOrderNonLinearDS::computeF(double time, SiconosVectorSPtr x2)
{
  if (isPlugin["f"])
  {
    if (computeFPtr == NULL)
      RuntimeException::selfThrow("FirstOrderNonLinearDS::computeVF() f is not linked to a plugin function");
    if (x2->size() != n)
      RuntimeException::selfThrow("FirstOrderNonLinearDS::computeJacobianXF(t,x) x size does not fit with the system size.");

    computeFPtr(time, n, &((*x2)(0)) , &(*f)(0), z->size(), &(*z)(0));
  }
  // else nothing!
}

void FirstOrderNonLinearDS::computeJacobianXF(double time, bool)
{
  // second argument is useless at the time - Used in derived classes
  if (isPlugin["jacobianXF"])
  {
    if (computeJacobianXFPtr == NULL)
      RuntimeException::selfThrow("FirstOrderNonLinearDS::computeJacobianXF() is not linked to a plugin function");
    computeJacobianXFPtr(time, n, &((*(x[0]))(0)), &(*jacobianXF)(0, 0), z->size(), &(*z)(0));
  }
  // else nothing!
}

void FirstOrderNonLinearDS::computeJacobianXF(double time, SP::SiconosVector x2)
{
  // second argument is useless at the time - Used in derived classes
  if (isPlugin["jacobianXF"])
  {
    if (computeJacobianXFPtr == NULL)
      RuntimeException::selfThrow("FirstOrderNonLinearDS::computeJacobianXF() is not linked to a plugin function");
    if (x2->size() != n)
      RuntimeException::selfThrow("FirstOrderNonLinearDS::computeJacobianXF(t,x) x size does not fit with the system size.");

    computeJacobianXFPtr(time, n, &((*x2)(0)), &(*jacobianXF)(0, 0), z->size(), &(*z)(0));
  }
  // else nothing!
}

void FirstOrderNonLinearDS::computeRhs(double time, bool)
{
  // second argument is useless at the time - Used in derived classes

  // compute rhs = M-1*( f + r ).

  *x[1] = *r; // Warning: p update is done in Interactions/Relations

  if (f)
  {
    computeF(time);
    *(x[1]) += *f;
  }

  if (M)
  {
    // allocate invM at the first call of the present function
    if (! invM)
    {

#ifndef WithSmartPtr
      invM = new SimpleMatrix(*M);
      isAllocatedIn["invM"] = true;
#else
      invM.reset(new SimpleMatrix(*M));
#endif
    }

    invM->PLUForwardBackwardInPlace(*x[1]);
  }
}

void FirstOrderNonLinearDS::computeJacobianXRhs(double time, bool)
{
  // second argument is useless at the time - Used in derived classes

  // compute jacobian of rhs according to x, = M-1(jacobianXF + jacobianX(T.u))
  // At the time, second term is set to zero.
  computeJacobianXF(time);
  // solve M*jacobianXRhS = jacobianXF
  if (M && jacobianXF)
  {
    *jacobianXRhs = *jacobianXF;
    // copy M into invM for LU-factorisation, at the first call of this function.
    if (! invM)
    {

#ifndef WithSmartPtr
      invM = new SimpleMatrix(*M);
      isAllocatedIn["invM"] = true;
#else
      invM.reset(new SimpleMatrix(*M));
#endif

    }

    invM->PLUForwardBackwardInPlace(*jacobianXRhs);
  }
  // else jacobianXRhs = jacobianXF, pointers equality set in initRhs

}

// ===== XML MANAGEMENT FUNCTIONS =====

void FirstOrderNonLinearDS::saveSpecificDataToXML()
{
  // -- FirstOrderNonLinearDS  xml object --
  FirstOrderNonLinearDSXMLSPtr fonlds = boost::static_pointer_cast <FirstOrderNonLinearDSXML>(dsxml);
  // --- other data ---
  if (fonlds == NULL)
    RuntimeException::selfThrow("FirstOrderNonLinearDS::saveSpecificDataToXML - The DynamicalSystemXML object doesn't exists");

  if (isPlugin["f"])
    fonlds->setFPlugin(pluginNames["f"]);
  if (isPlugin["jacobianXF"])
    fonlds->setJacobianXFPlugin(pluginNames["jacobianXF"]);
}

// ===== MISCELLANEOUS ====

void FirstOrderNonLinearDS::display() const
{
  cout << " =====> First Order Non Linear DS (number: " << number << ")." << endl;
  cout << "- n (size) : " << n << endl;
  cout << "- x " << endl;
  if (x[0]) x[0]->display();
  else cout << "-> NULL" << endl;
  cout << "- x0 " << endl;
  if (x0) x0->display();
  else cout << "-> NULL" << endl;
  cout << "- M: " << endl;
  if (M) M->display();
  else cout << "-> NULL" << endl;
  map<string, bool>::const_iterator it;
  cout << "Plug-in state (0: unplugged, ie constant, 1: plugged):" << endl;
  for (it = isPlugin.begin(); it != isPlugin.end(); ++it)
    cout << "Object " << it->first << " " << it->second << endl;
  cout << " ============================================" << endl;
}

void FirstOrderNonLinearDS::resetNonSmoothPart()
{
  r->zero();
}

double FirstOrderNonLinearDS::dsConvergenceIndicator()
{
  double dsCvgIndic;
  // Velocity is used to calculate the indicator.
  SP::SiconosVector diff(new SimpleVector(x[0]->size()));
  // Compute difference between present and previous Newton steps
  SiconosVectorSPtr valRef = workVector["NewtonCvg"];
  *diff =  *(x[0]) - *valRef;
  if (valRef->norm2() != 0)
    dsCvgIndic = diff->norm2() / (valRef->norm2());
  else
    dsCvgIndic = diff->norm2();
  return (dsCvgIndic);
}

FirstOrderNonLinearDS* FirstOrderNonLinearDS::convert(DynamicalSystem* ds)
{
  FirstOrderNonLinearDS* fonlds = dynamic_cast<FirstOrderNonLinearDS*>(ds);
  return fonlds;
}

