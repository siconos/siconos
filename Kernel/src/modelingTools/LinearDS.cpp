/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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
  DynamicalSystem(NULL), A(NULL), b(NULL), bFunctionName("none"), computeBPtr(NULL), isBAllocatedIn(false)
{
  DSType = LDS;
}

// From xml file (newNsds is optional)
LinearDS::LinearDS(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* newNsds):
  DynamicalSystem(dsXML, newNsds), A(NULL), b(NULL), bFunctionName("none"), computeBPtr(NULL), isBAllocatedIn(false)
{
  if (dsXML != NULL)
  {
    DSType = LDS;

    // pointer to xml
    LinearDSXML * ldsxml = (static_cast <LinearDSXML*>(dsxml));

    // --- vector and matrix members memory allocation ---
    // (only those specific to LinearDS) and values loading
    string plugin;
    isLDSPlugin.resize(2, false);

    // Check if vectorField is not given as a plug-in in xml input file.
    if (ldsxml->hasVectorFieldPlugin())
      RuntimeException::selfThrow("LinearDS - xml constructor, you give a vectorField plug-in for a LinearDS -> set rather A (or jacobianX) and b plug-in.");

    // A = jacobianX
    A = jacobianX; // jacobianX is allocated during DynamicalSystem constructor call

    // reject case were jacobianX plug-in is given
    if (ldsxml->hasComputeJacobianXPlugin())
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
        *A = ldsxml->getA();
    }
    else
      RuntimeException::selfThrow("LinearDS - xml constructor, no input (plug-in or matrix) find for A.");

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
        isBAllocatedIn = true;
      }
    }
  }
  else
    RuntimeException::selfThrow("LinearDS - xml constructor, xml file = NULL");
}

// For the following constructors, only A is required. If necessary b or u can be defined thanks
// to set or setCompute, depending on they are plugins or not.

// From a minimum set of data, A from a plugin
LinearDS::LinearDS(const int& newNumber, const unsigned int& newN, const SiconosVector& newX0,
                   const string& APlugin, const string& bPlugin):
  DynamicalSystem(newNumber, newN, newX0),
  A(NULL), b(NULL), bFunctionName("none"), computeBPtr(NULL), isBAllocatedIn(false)
{
  DSType = LDS;
  isLDSPlugin.resize(2, false);

  A = jacobianX; // jacobianX is allocated during DynamicalSystem constructor call
  setComputeJacobianXFunction(cShared.getPluginName(APlugin), cShared.getPluginFunctionName(APlugin));
  setComputeBFunction(cShared.getPluginName(bPlugin), cShared.getPluginFunctionName(bPlugin));
}

// From a minimum set of data, A from a given matrix
LinearDS::LinearDS(const int& newNumber, const SiconosVector& newX0, const SiconosMatrix& newA):
  DynamicalSystem(newNumber, newA.size(0), newX0),
  A(NULL), b(NULL), bFunctionName("none"), computeBPtr(NULL), isBAllocatedIn(false)
{
  if (newA.size(0) != n || newA.size(1) != n)
    RuntimeException::selfThrow("LinearDS - constructor(3): inconsistent dimensions with problem size for input matrix A");

  DSType = LDS;
  A = jacobianX; // jacobianX is allocated during DynamicalSystem constructor call
  *A = newA;
  isLDSPlugin.resize(2, false);
}

// Copy constructor
LinearDS::LinearDS(const LinearDS & lds):
  DynamicalSystem(lds), A(NULL), b(NULL), bFunctionName("none"), computeBPtr(NULL), isBAllocatedIn(false)
{
  DSType = LDS;

  A = jacobianX; // jacobianX is allocated during DynamicalSystem constructor call
  //*A =lds->getA(); // this may be useless, since jacobianX copy has already been done?

  if (lds.getBPtr() != NULL)
  {
    b = new SimpleVector(lds.getB());
    isBAllocatedIn = true;
  }

  isLDSPlugin = lds.getIsLDSPlugin();
  string pluginPath, functionName;
  if (isLDSPlugin[0])
  {
    string AFunctionName = lds.getComputeJacobianXFunctionName();
    functionName = cShared.getPluginFunctionName(AFunctionName);
    pluginPath  = cShared.getPluginName(AFunctionName);
    setComputeJacobianXFunction(pluginPath, functionName);
  }
  if (isLDSPlugin[1])
  {
    bFunctionName = lds.getBFunctionName();
    functionName = cShared.getPluginFunctionName(bFunctionName);
    pluginPath  = cShared.getPluginName(bFunctionName);
    setComputeBFunction(pluginPath, functionName);
  }
}

LinearDS::LinearDS(const DynamicalSystem & newDS):
  DynamicalSystem(newDS), A(NULL), b(NULL), bFunctionName("none"), computeBPtr(NULL), isBAllocatedIn(false)
{
  if (newDS.getType() != LDS)
    RuntimeException::selfThrow("LinearDS - copy constructor: try to copy into a LinearDS a DS of type: " + newDS.getType());

  DSType = LDS;

  // convert newDS to linearDS by keeping const options
  const LinearDS * lds = static_cast<const LinearDS*>(&newDS);

  A = jacobianX; // jacobianX is allocated during DynamicalSystem constructor call
  //*A =lds->getA(); // this may be useless, since jacobianX copy has already been done?

  if (lds->getBPtr() != NULL)
  {
    b = new SimpleVector(lds->getB());
    isBAllocatedIn = true;
  }

  isLDSPlugin = lds->getIsLDSPlugin();
  string pluginPath, functionName;
  if (isLDSPlugin[0])
  {
    string AFunctionName = lds->getComputeJacobianXFunctionName();
    functionName = cShared.getPluginFunctionName(AFunctionName);
    pluginPath  = cShared.getPluginName(AFunctionName);
    setComputeJacobianXFunction(pluginPath, functionName);
  }
  if (isLDSPlugin[1])
  {
    bFunctionName = lds->getBFunctionName();
    functionName = cShared.getPluginFunctionName(bFunctionName);
    pluginPath  = cShared.getPluginName(bFunctionName);
    setComputeBFunction(pluginPath, functionName);
  }
}

LinearDS::~LinearDS()
{
  A = NULL ;
  if (isBAllocatedIn)
  {
    delete b;
    b = NULL;
  }
}

void LinearDS::initialize(const double& time, const unsigned int& sizeOfMemory)
{
  // reset x to x0, xFree and r to zero.
  *x = *x0;
  xFree->zero();
  r->zero();

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  SimpleVector* param;
  // compute A if it is a plug-in
  if (isLDSPlugin[0])
  {
    param = parametersList0[1];
    computeJacobianXPtr(n, &time, &(*x)(0), &(*jacobianX)(0, 0), &(*param)(0));
  }

  // compute xDot = vectorField
  *xDot = *A * *x;

  if (b != NULL) // if b is given
  {
    if (isLDSPlugin[1] && computeBPtr == NULL) // if a plug-in is given for b
    {
      parametersList0[0]; // Since vectorField plug-in is not used in LinearDS, we use parameter(0) which corresponds to vectorField, for b.
      computeBPtr(n, &time, &(*b)(0), &(*param)(0));
    }
    *xDot += *b;
  }

  // If they are plugged, initialize u and T, and then complete xDot
  if (u != NULL && T != NULL)
  {
    if (isPlugin[0] && computeUPtr != NULL) // if u is a plug-in function
    {
      param = parametersList0[2];
      computeUPtr(uSize, n, &time, &(*x)(0), &(*u)(0), &(*param)(0));
    }
    if (isPlugin[1] && computeTPtr != NULL) // if T is a plug-in function
    {
      param = parametersList0[3];
      computeTPtr(uSize, n, &(*x)(0), &(*T)(0, 0), &(*param)(0));
    }

    *xDot += *T ** u;
  }
  else if (u != NULL && T == NULL)
  {
    if (isPlugin[0] && computeUPtr != NULL) // if u is a plug-in function
    {
      param = parametersList0[2];
      computeUPtr(uSize, n, &time, &(*x)(0), &(*u)(0), &(*param)(0));
    }
    *xDot += * u;
  }
}

void LinearDS::setA(const SiconosMatrix& newValue)
{
  setJacobianX(newValue);
  isLDSPlugin[0] = false;
}

void LinearDS::setAPtr(SiconosMatrix *newPtr)
{
  setJacobianXPtr(newPtr);
  isLDSPlugin[0] = false;
}

void LinearDS::setJacobianX(const SiconosMatrix& newValue)
{
  DynamicalSystem::setJacobianX(newValue);
  A = jacobianX;
  isLDSPlugin[0] = false;
}

void LinearDS::setJacobianXPtr(SiconosMatrix *newPtr)
{
  DynamicalSystem::setJacobianXPtr(newPtr);
  A = jacobianX;
  isLDSPlugin[0] = false;
}

void LinearDS::setB(const SimpleVector& newValue)
{
  if (newValue.size() != n)
    RuntimeException::selfThrow("LinearDS - setB: inconsistent dimensions with problem size for input vector b");

  if (b == NULL)
  {
    b = new SimpleVector(n);
    isBAllocatedIn = true;
  }
  *b = newValue;
  isLDSPlugin[1] = false;
}

void LinearDS::setBPtr(SimpleVector *newPtr)
{
  if (isBAllocatedIn) delete b;
  b = newPtr;
  isBAllocatedIn = false;
  isLDSPlugin[1] = false;
}

void LinearDS::setVectorFieldFunction(const string& pluginPath, const string& functionName)
{
  cout << " /!\\ LinearDS setVectorFieldFunction: useless function call. Set A (or jacobianX) and b plug-in /!\\ ." << endl;
}

void LinearDS::setComputeJacobianXFunction(const string& pluginPath, const string& functionName)
{
  DynamicalSystem::setComputeJacobianXFunction(pluginPath, functionName);
  A = jacobianX;
  isLDSPlugin[0] = true;
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
    isBAllocatedIn = true ;
  }
  cShared.setFunction(&computeBPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  bFunctionName = plugin + ":" + functionName;
  isLDSPlugin[1] = true;
}

void LinearDS::computeVectorField(const double& time)
{
  // compute A=jacobianX
  computeJacobianX(time);

  // compute xDot = vectorField
  *xDot = *A * *x;

  // compute and add b if required
  if (b != NULL)
  {
    computeB(time);
    *xDot += *b;
  }

  // compute and add Tu if required
  if (u != NULL && T != NULL)
  {
    if (isPlugin[0]) // if u is a plug-in function
      computeU(time);
    if (isPlugin[1]) // if T is a plug-in function
      computeT();
    *xDot += *T ** u;
  }
  else if (u != NULL && T == NULL)
  {
    if (isPlugin[0]) // if u is a plug-in function
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
  if (computeBPtr == NULL)
    RuntimeException::selfThrow("computeB() is not linked to a plugin function");

  unsigned int size = b->size();
  SimpleVector* param = parametersList0[0]; // Since vectorField plug-in is not used in LinearDS, we use parameter(0) which corresponds to vectorField, for b.
  computeBPtr(size, &time, &(*b)(0), &(*param)(0));
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
      static_cast<LinearDSXML*>(dsxml)->setBVector(*b);
    }
  }
  else RuntimeException::selfThrow("LinearDS::saveDSToXML - The DynamicalSystemXML object doesn't exists");
}

LinearDS* LinearDS::convert(DynamicalSystem* ds)
{
  LinearDS* lsds = dynamic_cast<LinearDS*>(ds);
  return lsds;
}
