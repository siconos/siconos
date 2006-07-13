/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
#include "LinearDS.h"
using namespace std;

void LinearDS::initAllocationFlags(const bool in) // default in = true.
{
  // isAllocatedIn["A"]=in; useless, since A is never allocated in, only link to jacobianXF
  isAllocatedIn["Mxdot"] = in;
  isAllocatedIn["b"] = in;
}

void LinearDS::initPluginFlags(const bool val)
{
  isPlugin["f"] = val;
  isPlugin["jacobianXF"] = val;
  isPlugin["A"] = val;
  isPlugin["b"] = val;
}

// --- Constructors ---

// Default constructor
LinearDS::LinearDS():
  DynamicalSystem(), A(NULL), Mxdot(NULL), b(NULL), computeAFunctionName("none"), computeBFunctionName("none"), APtr(NULL), bPtr(NULL)
{
  DSType = LDS;
  initAllocationFlags(false);
  initPluginFlags(false);
}

// From xml file (newNsds is optional)
LinearDS::LinearDS(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* newNsds):
  DynamicalSystem(dsXML, newNsds), A(NULL), Mxdot(NULL), b(NULL), computeAFunctionName("none"), computeBFunctionName("none"), APtr(NULL), bPtr(NULL)
{
  if (dsXML != NULL)
  {
    DSType = LDS;
    // pointer to xml
    LinearDSXML * ldsxml = (static_cast <LinearDSXML*>(dsxml));

    // --- vector and matrix members memory allocation ---
    // (only those specific to LinearDS) and values loading
    string plugin;

    // Check if f is given as a plug-in in xml input file.
    if (ldsxml->hasF())
      RuntimeException::selfThrow("LinearDS - xml constructor, you give a f plug-in for a LinearDS -> set rather A (or jacobianXF) and b plug-in.");

    initAllocationFlags(false);
    initPluginFlags(false);

    // A = jacobianXF
    A = jacobianXF; // jacobianXF is allocated during DynamicalSystem constructor call

    // reject case were jacobianXF plug-in is given
    if (ldsxml->hasJacobianXF())
      RuntimeException::selfThrow("LinearDS - xml constructor, you give a plug-in for jacobianXF, set rather A.");

    // set A or A plug-in - A is a required input in xml (xml error if node not found)
    if (ldsxml->hasA())
    {
      if (ldsxml->isAPlugin())
      {
        plugin = ldsxml->getAPlugin();
        setComputeAFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      }
      else
        *A = ldsxml->getA();
    }
    else
      RuntimeException::selfThrow("LinearDS - xml constructor, no input (plug-in or matrix) find for A.");

    // Mxdot - Optional parameter supposed to be a SimpleMatrix with xml
    if (ldsxml->hasMxdot())
    {
      Mxdot = new SimpleMatrix(ldsxml->getMxdot());
      isAllocatedIn["Mxdot"] = true;
    }

    // b - Optional parameter
    if (ldsxml->hasB())
    {
      if (ldsxml->isBPlugin())
      {
        plugin = ldsxml->getBPlugin();
        setComputeBFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      }
      else
      {
        b = new SimpleVector(ldsxml->getBVector());
        isAllocatedIn["b"] = true;
      }
    }
  }
  else
    RuntimeException::selfThrow("LinearDS - xml constructor, xml file = NULL");
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// For the following constructors, only A is required. If necessary b or u can be defined thanks
// to set or setCompute, depending on they are plugins or not.

// From a minimum set of data, A and b connected to a plug-in
LinearDS::LinearDS(const int newNumber, const unsigned int newN, const SiconosVector& newX0,
                   const string APlugin, const string bPlugin):
  DynamicalSystem(newNumber, newN, newX0),
  A(NULL), Mxdot(NULL), b(NULL), computeAFunctionName("none"), computeBFunctionName("none"), APtr(NULL), bPtr(NULL)
{
  DSType = LDS;
  initAllocationFlags(false);
  initPluginFlags(false);
  A = jacobianXF; // jacobianXF is allocated during DynamicalSystem constructor call
  setComputeAFunction(cShared.getPluginName(APlugin), cShared.getPluginFunctionName(APlugin));
  setComputeBFunction(cShared.getPluginName(bPlugin), cShared.getPluginFunctionName(bPlugin));
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// From a minimum set of data, A from a given matrix
LinearDS::LinearDS(const int newNumber, const SiconosVector& newX0, const SiconosMatrix& newA):
  DynamicalSystem(newNumber, newA.size(0), newX0),
  A(NULL), Mxdot(NULL), b(NULL), computeAFunctionName("none"), computeBFunctionName("none"), APtr(NULL), bPtr(NULL)
{
  if (newA.size(0) != n || newA.size(1) != n)
    RuntimeException::selfThrow("LinearDS - constructor(number,x0,A): inconsistent dimensions with problem size for input matrix A");

  initAllocationFlags(false);
  initPluginFlags(false);
  DSType = LDS;
  A = jacobianXF; // jacobianXF is allocated during DynamicalSystem constructor call
  *A = newA;
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// From a minimum set of data, A from a given matrix
LinearDS::LinearDS(const int newNumber, const SiconosVector& newX0, const SiconosMatrix& newA, const SiconosVector& newB):
  DynamicalSystem(newNumber, newA.size(0), newX0),
  A(NULL), Mxdot(NULL), b(NULL), computeAFunctionName("none"), computeBFunctionName("none"), APtr(NULL), bPtr(NULL)
{
  if (newA.size(0) != n || newA.size(1) != n)
    RuntimeException::selfThrow("LinearDS - constructor(number,x0,A,b): inconsistent dimensions with problem size for input matrix A");
  if (newB.size() != n)
    RuntimeException::selfThrow("LinearDS - constructor(number,x0,A,b): inconsistent dimensions with problem size for input vector b.");

  initAllocationFlags(false);
  initPluginFlags(false);
  DSType = LDS;
  A = jacobianXF; // jacobianXF is allocated during DynamicalSystem constructor call
  *A = newA;
  b = new SimpleVector(newB);
  isAllocatedIn["b"] = true;
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// Copy constructor
LinearDS::LinearDS(const LinearDS & lds):
  DynamicalSystem(lds), A(NULL), Mxdot(NULL), b(NULL), computeAFunctionName("none"), computeBFunctionName("none"), APtr(NULL), bPtr(NULL)
{
  DSType = LDS;
  initAllocationFlags(false);
  initPluginFlags(false);

  A = jacobianXF; // jacobianXF is allocated during DynamicalSystem constructor call
  //*A =lds->getA(); // this may be useless, since jacobianXF copy has already been done?

  if (lds.getMxdotPtr() != NULL)
  {
    if (lds.getMxdotPtr()->isBlock())
      Mxdot = new BlockMatrix(lds.getMxdotBlock());
    else
      Mxdot = new SimpleMatrix(lds.getMxdotSimple());

    isAllocatedIn["Mxdot"] = true;
  }

  if (lds.getBPtr() != NULL)
  {
    b = new SimpleVector(lds.getB());
    isAllocatedIn["b"] = true;
  }

  // note that isPlugin has been copied during DynamicalSystem constructor call.
  string pluginPath, functionName;
  if (isPlugin["A"])
  {
    computeAFunctionName = lds.getComputeAFunctionName();
    functionName = cShared.getPluginFunctionName(computeAFunctionName);
    pluginPath  = cShared.getPluginName(computeAFunctionName);
    setComputeJacobianXFFunction(pluginPath, functionName);
  }
  if (isPlugin["b"])
  {
    computeBFunctionName = lds.getComputeBFunctionName();
    functionName = cShared.getPluginFunctionName(computeBFunctionName);
    pluginPath  = cShared.getPluginName(computeBFunctionName);
    setComputeBFunction(pluginPath, functionName);
  }
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

LinearDS::LinearDS(const DynamicalSystem & newDS):
  DynamicalSystem(newDS), A(NULL), Mxdot(NULL), b(NULL), computeAFunctionName("none"), computeBFunctionName("none"), APtr(NULL), bPtr(NULL)
{
  if (newDS.getType() != LDS || newDS.getType() != LITIDS)
    RuntimeException::selfThrow("LinearDS - copy constructor: try to copy into a LinearDS a DS of type: " + newDS.getType());

  DSType = LDS;
  initAllocationFlags(false);
  initPluginFlags(false);

  // convert newDS to linearDS by keeping const options
  const LinearDS * lds = static_cast<const LinearDS*>(&newDS);

  A = jacobianXF; // jacobianXF is allocated during DynamicalSystem constructor call
  //*A =lds->getA(); // this may be useless, since jacobianXF copy has already been done?

  if (lds->getMxdotPtr() != NULL)
  {
    if (lds->getMxdotPtr()->isBlock())
      Mxdot = new BlockMatrix(lds->getMxdotBlock());
    else
      Mxdot = new SimpleMatrix(lds->getMxdotSimple());

    isAllocatedIn["Mxdot"] = true;
  }

  if (lds->getBPtr() != NULL)
  {
    b = new SimpleVector(lds->getB());
    isAllocatedIn["b"] = true;
  }

  string pluginPath, functionName;
  // note that isPlugin has been copied during DynamicalSystem constructor call.
  if (isPlugin["A"])
  {
    computeAFunctionName = lds->getComputeAFunctionName();
    functionName = cShared.getPluginFunctionName(computeAFunctionName);
    pluginPath  = cShared.getPluginName(computeAFunctionName);
    setComputeJacobianXFFunction(pluginPath, functionName);
  }
  if (isPlugin["b"])
  {
    computeBFunctionName = lds->getComputeBFunctionName();
    functionName = cShared.getPluginFunctionName(computeBFunctionName);
    pluginPath  = cShared.getPluginName(computeBFunctionName);
    setComputeBFunction(pluginPath, functionName);
  }
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

LinearDS::~LinearDS()
{
  A = NULL ;
  if (isAllocatedIn["Mxdot"]) delete Mxdot;
  Mxdot = NULL;
  if (isAllocatedIn["b"]) delete b;
  b = NULL;
}

bool LinearDS::checkDynamicalSystem()
{
  bool output = true;
  // n
  if (n == 0)
  {
    RuntimeException::selfThrow("LinearDS::checkDynamicalSystem - number of degrees of freedom is equal to 0.");
    output = false;
  }
  // x0 != NULL
  if (x0 == NULL)
  {
    RuntimeException::selfThrow("LinearDS::checkDynamicalSystem - x0 not set.");
    output = false;
  }

  // A
  if (A == NULL)
  {
    RuntimeException::selfThrow("LinearDS::checkDynamicalSystem - A not set.");
    output = false;
  }

  return output;
}

void LinearDS::initialize(const double time, const unsigned int sizeOfMemory)
{
  // reset x to x0, xFree and r to zero.
  *x = *x0;
  xFree->zero();
  r->zero();

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  computeRhs(time); // If necessary, this will also compute A, b, u and T

}


void LinearDS::update(const double time)
{
  computeRhs(time); // If necessary, this will also compute A, b, u and T
}

void LinearDS::setA(const SiconosMatrix& newValue)
{
  setJacobianXF(newValue);
  isPlugin["A"] = false;
}

void LinearDS::setAPtr(SiconosMatrix *newPtr)
{
  setJacobianXFPtr(newPtr);
  isPlugin["A"] = false;
}

void LinearDS::setJacobianXF(const SiconosMatrix& newValue)
{
  DynamicalSystem::setJacobianXF(newValue);
  A = jacobianXF;
  isPlugin["jacobianXF"] = false;
  isPlugin["A"] = false;
}

void LinearDS::setJacobianXFPtr(SiconosMatrix *newPtr)
{
  DynamicalSystem::setJacobianXFPtr(newPtr);
  A = jacobianXF;
  isPlugin["jacobianXF"] = false;
  isPlugin["A"] = false;
}

void LinearDS::setMxdot(const SimpleMatrix& newValue)
{
  if (newValue.size(0) != n || newValue.size(1) != n)
    RuntimeException::selfThrow("LinearDS - setMxdot: inconsistent dimensions with problem size for input matrix Mxdot");

  if (Mxdot == NULL)
  {
    Mxdot = new SimpleMatrix(n, n);
    isAllocatedIn["Mxdot"] = true;
  }
  *Mxdot = newValue;
}

void LinearDS::setMxdotPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn["Mxdot"]) delete Mxdot;
  Mxdot = newPtr;
  isAllocatedIn["Mxdot"] = false;
}

void LinearDS::setB(const SimpleVector& newValue)
{
  if (newValue.size() != n)
    RuntimeException::selfThrow("LinearDS - setB: inconsistent dimensions with problem size for input vector b");

  if (b == NULL)
  {
    b = new SimpleVector(n);
    isAllocatedIn["b"] = true;
  }
  *b = newValue;
  isPlugin["b"] = false;
}

void LinearDS::setBPtr(SimpleVector *newPtr)
{
  if (isAllocatedIn["b"]) delete b;
  b = newPtr;
  isAllocatedIn["b"] = false;
  isPlugin["b"] = false;
}

void LinearDS::setComputeAFunction(const string pluginPath, const string functionName)
{
  if (jacobianXF == NULL)
  {
    jacobianXF = new SimpleMatrix(n, n);
    isAllocatedIn["jacobianXF"] = true ;
    A = jacobianXF;
  }

  initParameter("A");

  APtr = NULL;
  cShared.setFunction(&APtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  computeAFunctionName = plugin + ":" + functionName;
  isPlugin["A"] = true;
}

void LinearDS::setComputeBFunction(const string pluginPath, const string functionName)
{
  if (b == NULL)
  {
    b = new SimpleVector(n);
    isAllocatedIn["b"] = true;
  }
  cShared.setFunction(&bPtr, pluginPath, functionName);

  initParameter("b");

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  computeBFunctionName = plugin + ":" + functionName;
  isPlugin["b"] = true;
}

void LinearDS::computeF(const double time)// is it really necessary to reimplement this?
{
  if (isPlugin["A"])
    computeA(time);
}

void LinearDS::computeJacobianXF(const double time, const bool) // is it really necessary to reimplement this?
{
  if (isPlugin["A"])
    computeA(time);
}

void LinearDS::computeA(const double time)
{
  if (isPlugin["A"])
  {
    if (APtr == NULL)
      RuntimeException::selfThrow("computeA() is not linked to a plugin function");
    SimpleVector* param = parametersList["A"];
    APtr(n, &time, &(*A)(0, 0), &(*param)(0));
  }
  // else nothing
}

void LinearDS::computeB(const double time)
{
  if (isPlugin["b"])
  {
    if (bPtr == NULL)
      RuntimeException::selfThrow("computeB() is not linked to a plugin function");
    SimpleVector* param = parametersList["b"];
    bPtr(n, &time, &(*b)(0), &(*param)(0));
  }
  // else nothing
}

void LinearDS::computeRhs(const double time, const bool)
{
  // second argument is useless at the time - Used in derived classes
  // compute A=jacobianXF
  if (isPlugin["A"])
    computeA(time);

  // compute right-hand side
  *rhs = *A * *x;

  // compute and add b if required
  if (b != NULL)
  {
    if (isPlugin["b"])
      computeB(time);
    *rhs += *b;
  }

  // compute and add Tu if required
  if (u != NULL)
  {
    if (isPlugin["u"]) // if u is a plug-in function
      computeU(time);
    if (T != NULL)
    {
      if (isPlugin["T"]) // if T is a plug-in function
        computeT();
      *rhs += *T ** u;
    }
    else
      *rhs += * u;
  }

  *rhs += * r; // Warning: r update is done in Interactions/Relations

}

void LinearDS::computeJacobianXRhs(const double time, const bool)
{
  if (isPlugin["A"])
    computeA(time);
}

void LinearDS::display() const
{
  DynamicalSystem::display();
  cout << "=== Linear system display ===" << endl;
  cout << "- A " << endl;
  if (A != NULL) A->display();
  else cout << "-> NULL" << endl;
  cout << "- Mxdot " << endl;
  if (Mxdot != NULL) Mxdot->display();
  else cout << "-> NULL" << endl;
  cout << "- b " << endl;
  if (b != NULL) b->display();
  else cout << "-> NULL" << endl;
  cout << "=============================" << endl;
}

void LinearDS::saveDSToXML()
{
  //--- Common data ---
  saveDSDataToXML();
  // --- other data ---
  if (dsxml != NULL)
  {
    dsxml->setN(n);
    static_cast<LinearDSXML*>(dsxml)->setA(*A);

    if (Mxdot != NULL) static_cast<LinearDSXML*>(dsxml)->setMxdot(*Mxdot);
    // b
    if (b != NULL)
    {
      if (!(static_cast <LinearDSXML*>(dsxml))->isBPlugin())
        static_cast<LinearDSXML*>(dsxml)->setB(*b);
    }
  }
  else RuntimeException::selfThrow("LinearDS::saveDSToXML - The DynamicalSystemXML object doesn't exists");
}

LinearDS* LinearDS::convert(DynamicalSystem* ds)
{
  LinearDS* lsds = dynamic_cast<LinearDS*>(ds);
  return lsds;
}
