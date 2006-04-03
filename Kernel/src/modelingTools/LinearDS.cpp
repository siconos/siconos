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
#include "LinearDS.h"
using namespace std;

// --- Constructors ---

// Default constructor
LinearDS::LinearDS():
  DynamicalSystem(NULL), A(NULL), b(NULL), bFunctionName("none"), computeBPtr(NULL)
{
  DSType = LDS;
  isAllocatedIn["b"] = false;
  isPlugin["b"] = false;
}

// From xml file (newNsds is optional)
LinearDS::LinearDS(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* newNsds):
  DynamicalSystem(dsXML, newNsds), A(NULL), b(NULL), bFunctionName("none"), computeBPtr(NULL)
{
  if (dsXML != NULL)
  {
    DSType = LDS;
    // pointer to xml
    LinearDSXML * ldsxml = (static_cast <LinearDSXML*>(dsxml));

    // --- vector and matrix members memory allocation ---
    // (only those specific to LinearDS) and values loading
    string plugin;

    // Check if vectorField is not given as a plug-in in xml input file.
    if (ldsxml->hasVectorField())
      RuntimeException::selfThrow("LinearDS - xml constructor, you give a vectorField plug-in for a LinearDS -> set rather A (or jacobianX) and b plug-in.");

    isPlugin["vectorField"] = false;

    // A = jacobianX
    A = jacobianX; // jacobianX is allocated during DynamicalSystem constructor call

    // reject case were jacobianX plug-in is given
    if (ldsxml->hasJacobianX())
      RuntimeException::selfThrow("LinearDS - xml constructor, you give a plug-in for jacobianX, set rather A.");

    // set A or A plug-in (ie jacobianX) - A is a required input in xml (xml error if node not found)
    if (ldsxml->hasA())
    {
      if (ldsxml->isAPlugin())
      {
        plugin = ldsxml->getAPlugin();
        setComputeJacobianXFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      }
      else
      {
        *A = ldsxml->getA();
        isPlugin["jacobianX"] = false;
      }
    }
    else
      RuntimeException::selfThrow("LinearDS - xml constructor, no input (plug-in or matrix) find for A.");

    // b - Optional parameter
    isAllocatedIn["b"] = false;
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
        isPlugin["b"] = false;
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
LinearDS::LinearDS(const int& newNumber, const unsigned int& newN, const SiconosVector& newX0,
                   const string& APlugin, const string& bPlugin):
  DynamicalSystem(newNumber, newN, newX0),
  A(NULL), b(NULL), bFunctionName("none"), computeBPtr(NULL)
{
  DSType = LDS;
  isAllocatedIn["b"] = false;
  isPlugin["vectorField"] = false;
  A = jacobianX; // jacobianX is allocated during DynamicalSystem constructor call
  setComputeJacobianXFunction(cShared.getPluginName(APlugin), cShared.getPluginFunctionName(APlugin));
  setComputeBFunction(cShared.getPluginName(bPlugin), cShared.getPluginFunctionName(bPlugin));
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// From a minimum set of data, A from a given matrix
LinearDS::LinearDS(const int& newNumber, const SiconosVector& newX0, const SiconosMatrix& newA):
  DynamicalSystem(newNumber, newA.size(0), newX0),
  A(NULL), b(NULL), bFunctionName("none"), computeBPtr(NULL)
{
  if (newA.size(0) != n || newA.size(1) != n)
    RuntimeException::selfThrow("LinearDS - constructor(number,x0,A): inconsistent dimensions with problem size for input matrix A");

  isAllocatedIn["b"] = false;
  isPlugin["b"] = false;
  isPlugin["vectorField"] = false;
  isPlugin["jacobianX"] = false; // A
  DSType = LDS;
  A = jacobianX; // jacobianX is allocated during DynamicalSystem constructor call
  *A = newA;
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// From a minimum set of data, A from a given matrix
LinearDS::LinearDS(const int& newNumber, const SiconosVector& newX0, const SiconosMatrix& newA, const SiconosVector& newB):
  DynamicalSystem(newNumber, newA.size(0), newX0),
  A(NULL), b(NULL), bFunctionName("none"), computeBPtr(NULL)
{
  if (newA.size(0) != n || newA.size(1) != n)
    RuntimeException::selfThrow("LinearDS - constructor(number,x0,A,b): inconsistent dimensions with problem size for input matrix A");
  if (newB.size() != n)
    RuntimeException::selfThrow("LinearDS - constructor(number,x0,A,b): inconsistent dimensions with problem size for input vector b.");

  isPlugin["vectorField"] = false;
  isPlugin["jacobianX"] = false; // A
  DSType = LDS;
  A = jacobianX; // jacobianX is allocated during DynamicalSystem constructor call
  *A = newA;
  b = new SimpleVector(newB);
  isAllocatedIn["b"] = true;
  isPlugin["b"] = false;
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// Copy constructor
LinearDS::LinearDS(const LinearDS & lds):
  DynamicalSystem(lds), A(NULL), b(NULL), bFunctionName("none"), computeBPtr(NULL)
{
  DSType = LDS;

  A = jacobianX; // jacobianX is allocated during DynamicalSystem constructor call
  //*A =lds->getA(); // this may be useless, since jacobianX copy has already been done?

  if (lds.getBPtr() != NULL)
  {
    b = new SimpleVector(lds.getB());
    isAllocatedIn["b"] = true;
  }
  else isAllocatedIn["b"] = false;

  // note that isPlugin has been copied during DynamicalSystem constructor call.
  isPlugin["vectorField"] = false;
  string pluginPath, functionName;
  if (isPlugin["jacobianX"])
  {
    string AFunctionName = lds.getComputeJacobianXFunctionName();
    functionName = cShared.getPluginFunctionName(AFunctionName);
    pluginPath  = cShared.getPluginName(AFunctionName);
    setComputeJacobianXFunction(pluginPath, functionName);
  }
  if (isPlugin["b"])
  {
    bFunctionName = lds.getBFunctionName();
    functionName = cShared.getPluginFunctionName(bFunctionName);
    pluginPath  = cShared.getPluginName(bFunctionName);
    setComputeBFunction(pluginPath, functionName);
  }
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

LinearDS::LinearDS(const DynamicalSystem & newDS):
  DynamicalSystem(newDS), A(NULL), b(NULL), bFunctionName("none"), computeBPtr(NULL)
{
  if (newDS.getType() != LDS || newDS.getType() != LITIDS)
    RuntimeException::selfThrow("LinearDS - copy constructor: try to copy into a LinearDS a DS of type: " + newDS.getType());

  DSType = LDS;

  // convert newDS to linearDS by keeping const options
  const LinearDS * lds = static_cast<const LinearDS*>(&newDS);

  A = jacobianX; // jacobianX is allocated during DynamicalSystem constructor call
  //*A =lds->getA(); // this may be useless, since jacobianX copy has already been done?

  if (lds->getBPtr() != NULL)
  {
    b = new SimpleVector(lds->getB());
    isAllocatedIn["b"] = true;
  }
  else isAllocatedIn["b"] = false;

  string pluginPath, functionName;
  // note that isPlugin has been copied during DynamicalSystem constructor call.
  isPlugin["vectorField"] = false;
  if (isPlugin["jacobianX"])
  {
    string AFunctionName = lds->getComputeJacobianXFunctionName();
    functionName = cShared.getPluginFunctionName(AFunctionName);
    pluginPath  = cShared.getPluginName(AFunctionName);
    setComputeJacobianXFunction(pluginPath, functionName);
  }
  if (isPlugin["b"])
  {
    bFunctionName = lds->getBFunctionName();
    functionName = cShared.getPluginFunctionName(bFunctionName);
    pluginPath  = cShared.getPluginName(bFunctionName);
    setComputeBFunction(pluginPath, functionName);
  }
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

LinearDS::~LinearDS()
{
  A = NULL ;
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

  // if vectorField not constant (plugged) => jacobianX required
  if (isPlugin["vectorField"])
  {
    RuntimeException::selfThrow("LinearDS::checkDynamicalSystem - vectorField is plugged and should not be.");
    output = false;
  }
  return output;
}

void LinearDS::initialize(const double& time, const unsigned int& sizeOfMemory)
{
  // reset x to x0, xFree and r to zero.
  *x = *x0;
  xFree->zero();
  r->zero();

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  computeVectorField(time); // If necessary, this will also compute A, b, u and T

}

void LinearDS::setA(const SiconosMatrix& newValue)
{
  setJacobianX(newValue);
  isPlugin["jacobianX"] = false;
}

void LinearDS::setAPtr(SiconosMatrix *newPtr)
{
  setJacobianXPtr(newPtr);
  isPlugin["jacobianX"] = false;
}

void LinearDS::setJacobianX(const SiconosMatrix& newValue)
{
  DynamicalSystem::setJacobianX(newValue);
  A = jacobianX;
  isPlugin["jacobianX"] = false;
}

void LinearDS::setJacobianXPtr(SiconosMatrix *newPtr)
{
  DynamicalSystem::setJacobianXPtr(newPtr);
  A = jacobianX;
  isPlugin["jacobianX"] = false;
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

void LinearDS::setVectorFieldFunction(const string& pluginPath, const string& functionName)
{
  cout << " /!\\ LinearDS setVectorFieldFunction: useless function call. Set A (or jacobianX) and b plug-in /!\\ ." << endl;
}

void LinearDS::setComputeJacobianXFunction(const string& pluginPath, const string& functionName)
{
  DynamicalSystem::setComputeJacobianXFunction(pluginPath, functionName);
  A = jacobianX;
  isPlugin["jacobianX"] = true;
}

void LinearDS::setComputeAFunction(const string& pluginPath, const string& functionName)
{
  setComputeJacobianXFunction(pluginPath, functionName);
}

void LinearDS::setComputeBFunction(const string& pluginPath, const string& functionName)
{
  if (b == NULL)
  {
    b = new SimpleVector(n);
    isAllocatedIn["b"] = true;
  }
  cShared.setFunction(&computeBPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  bFunctionName = plugin + ":" + functionName;
  isPlugin["b"] = true;
}

void LinearDS::computeVectorField(const double& time)
{
  // compute A=jacobianX
  if (isPlugin["jacobianX"])
    computeJacobianX(time);

  // compute xDot = vectorField
  *xDot = *A * *x;

  // compute and add b if required
  if (b != NULL)
  {
    if (isPlugin["b"])
      computeB(time);
    *xDot += *b;
  }

  // compute and add Tu if required
  if (u != NULL && T != NULL)
  {
    if (isPlugin["u"]) // if u is a plug-in function
      computeU(time);
    if (isPlugin["T"]) // if T is a plug-in function
      computeT();
    *xDot += *T ** u;
  }
  else if (u != NULL && T == NULL)
  {
    if (isPlugin["u"]) // if u is a plug-in function
      computeU(time);
    *xDot += * u;
  }
}

void LinearDS::computeA(const double& time)
{
  DynamicalSystem::computeJacobianX(time);
}

void LinearDS::computeB(const double& time)
{
  if (isPlugin["b"])
  {
    if (computeBPtr == NULL)
      RuntimeException::selfThrow("computeB() is not linked to a plugin function");
    SimpleVector* param = parametersList0[0]; // Since vectorField plug-in is not used in LinearDS, we use parameter(0) which corresponds to vectorField, for b.
    computeBPtr(n, &time, &(*b)(0), &(*param)(0));
  }
  // else nothing
}

void LinearDS::display() const
{
  DynamicalSystem::display();
  cout << "=== Linear system display ===" << endl;
  cout << "- A " << endl;
  if (A != NULL) A->display();
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

    // b
    if (!(static_cast <LinearDSXML*>(dsxml))->isBPlugin())
    {
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
