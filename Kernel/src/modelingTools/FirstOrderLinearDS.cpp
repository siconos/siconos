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
#include "FirstOrderLinearDS.h"
#include "FirstOrderLinearDSXML.h"

using namespace std;

void FirstOrderLinearDS::initAllocationFlags(bool in) // default in = true.
{
  isAllocatedIn["A"] = in;
  isAllocatedIn["b"] = in;
}

void FirstOrderLinearDS::initPluginFlags(bool val)
{
  isPlugin["A"] = val;
  isPlugin["b"] = val;
}

// --- Constructors ---

// Default constructor
FirstOrderLinearDS::FirstOrderLinearDS():
  FirstOrderNonLinearDS(FOLDS), A(NULL), b(NULL), APtr(NULL), bPtr(NULL)
{
  initAllocationFlags(false);
  initPluginFlags(false);
}

// From xml file (newNsds is optional)
FirstOrderLinearDS::FirstOrderLinearDS(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* newNsds):
  FirstOrderNonLinearDS(dsXML, newNsds), A(NULL), b(NULL), APtr(NULL), bPtr(NULL)
{
  // pointer to xml
  FirstOrderLinearDSXML * foldsxml = (static_cast <FirstOrderLinearDSXML*>(dsxml));

  // Check if f is given as a plug-in in xml input file.
  if (foldsxml->hasF() || foldsxml->hasJacobianXF())
    RuntimeException::selfThrow("FirstOrderLinearDS - xml constructor, you give a f or its jacobian as a plug-in for a FirstOrderLinearDS -> set rather A and b plug-in.");

  initAllocationFlags(false);
  initPluginFlags(false);
  string plugin;
  // A
  if (foldsxml->hasA())
  {
    if (foldsxml->isAPlugin())
    {
      plugin = foldsxml->getAPlugin();
      setComputeAFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {
      A = new SimpleMatrix(foldsxml->getA());
      isAllocatedIn["A"] = true;
    }
  }

  // b
  if (foldsxml->hasB())
  {
    if (foldsxml->isBPlugin())
    {
      plugin = foldsxml->getBPlugin();
      setComputeBFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {
      b = new SimpleVector(foldsxml->getBVector());
      isAllocatedIn["b"] = true;
    }
  }
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// From a minimum set of data, A and b connected to a plug-in
FirstOrderLinearDS::FirstOrderLinearDS(int newNumber, const SiconosVector& newX0,
                                       const string& APlugin, const string& bPlugin):
  FirstOrderNonLinearDS(newNumber, newX0),
  A(NULL), b(NULL), APtr(NULL), bPtr(NULL)
{
  DSType = FOLDS;
  initAllocationFlags(false);
  initPluginFlags(false);
  setComputeAFunction(cShared.getPluginName(APlugin), cShared.getPluginFunctionName(APlugin));
  setComputeBFunction(cShared.getPluginName(bPlugin), cShared.getPluginFunctionName(bPlugin));
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// From a minimum set of data, A from a given matrix
FirstOrderLinearDS::FirstOrderLinearDS(int newNumber, const SiconosVector& newX0, const SiconosMatrix& newA):
  FirstOrderNonLinearDS(newNumber, newX0),
  A(NULL), b(NULL), APtr(NULL), bPtr(NULL)
{
  DSType = FOLDS;
  if (newA.size(0) != n || newA.size(1) != n)
    RuntimeException::selfThrow("FirstOrderLinearDS - constructor(number,x0,A): inconsistent dimensions with problem size for input matrix A");

  initAllocationFlags(false);
  initPluginFlags(false);
  A = new SimpleMatrix(newA);
  isAllocatedIn["A"] = true;

  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// From a minimum set of data, A from a given matrix
FirstOrderLinearDS::FirstOrderLinearDS(const int newNumber, const SiconosVector& newX0, const SiconosMatrix& newA, const SiconosVector& newB):
  FirstOrderNonLinearDS(newNumber, newX0),
  A(NULL), b(NULL), APtr(NULL), bPtr(NULL)
{
  DSType = FOLDS;
  if (newA.size(0) != n || newA.size(1) != n)
    RuntimeException::selfThrow("FirstOrderLinearDS - constructor(number,x0,A,b): inconsistent dimensions with problem size for input matrix A");
  if (newB.size() != n)
    RuntimeException::selfThrow("FirstOrderLinearDS - constructor(number,x0,A,b): inconsistent dimensions with problem size for input vector b.");

  initAllocationFlags(false);
  initPluginFlags(false);
  A = new SimpleMatrix(newA);
  isAllocatedIn["A"] = true;
  b = new SimpleVector(newB);
  isAllocatedIn["b"] = true;
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

FirstOrderLinearDS::~FirstOrderLinearDS()
{
  if (isAllocatedIn["A"]) delete A;
  A = NULL;
  if (isAllocatedIn["b"]) delete b;
  b = NULL;
}

bool FirstOrderLinearDS::checkDynamicalSystem() // useless ...
{
  bool output = true;
  // n
  if (n == 0)
  {
    RuntimeException::selfThrow("FirstOrderLinearDS::checkDynamicalSystem - number of degrees of freedom is equal to 0.");
    output = false;
  }
  // x0 != NULL
  if (x0 == NULL)
  {
    RuntimeException::selfThrow("FirstOrderLinearDS::checkDynamicalSystem - x0 not set.");
    output = false;
  }

  return output;
}

void FirstOrderLinearDS::initRhs(double time)
{
  computeRhs(time); // If necessary, this will also compute A and b.
  if (jacobianXRhs == NULL) // if not allocated with a set or anything else
  {
    if (A != NULL && M == NULL) // if M is not defined, then A = jacobianXRhs, no memory allocation for that one.
      jacobianXRhs = A;
    else if (A != NULL && M != NULL)
    {
      jacobianXRhs = new SimpleMatrix(n, n);
      isAllocatedIn["jacobianXRhs"] = true;
    }
    // else no allocation, jacobian is equal to 0.
  }
  computeJacobianXRhs(time);
}

void FirstOrderLinearDS::setA(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != n || newValue.size(1) != n)
    RuntimeException::selfThrow("FirstOrderLinearDS - setA: inconsistent dimensions with problem size for input matrix A");

  if (A == NULL)
  {
    A = new SimpleMatrix(newValue);
    isAllocatedIn["A"] = true;
  }
  else
    *A = newValue;
  isPlugin["A"] = false;
}

void FirstOrderLinearDS::setAPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn["A"]) delete A;
  A = newPtr;
  isAllocatedIn["A"] = false;
  isPlugin["A"] = false;
}

void FirstOrderLinearDS::setB(const SimpleVector& newValue)
{
  if (newValue.size() != n)
    RuntimeException::selfThrow("FirstOrderLinearDS - setB: inconsistent dimensions with problem size for input vector b");

  if (b == NULL)
  {
    b = new SimpleVector(newValue);
    isAllocatedIn["b"] = true;
  }
  else
    *b = newValue;
  isPlugin["b"] = false;
}

void FirstOrderLinearDS::setBPtr(SimpleVector *newPtr)
{
  if (isAllocatedIn["b"]) delete b;
  b = newPtr;
  isAllocatedIn["b"] = false;
  isPlugin["b"] = false;
}

void FirstOrderLinearDS::setComputeAFunction(const string& pluginPath, const string& functionName)
{
  if (A == NULL)
  {
    A = new SimpleMatrix(n, n);
    isAllocatedIn["A"] = true ;
  }

  APtr = NULL;
  cShared.setFunction(&APtr, pluginPath, functionName);

  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["A"] = plugin + ":" + functionName;
  isPlugin["A"] = true;
}

void FirstOrderLinearDS::setComputeBFunction(const string& pluginPath, const string& functionName)
{
  if (b == NULL)
  {
    b = new SimpleVector(n);
    isAllocatedIn["b"] = true;
  }
  cShared.setFunction(&bPtr, pluginPath, functionName);

  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["b"] = plugin + ":" + functionName;
  isPlugin["b"] = true;
}
void   FirstOrderLinearDS::setComputeBFunction(bPtrFunction fct)
{
  if (b == NULL)
  {
    b = new SimpleVector(n);
    isAllocatedIn["b"] = true;
  }
  bPtr = fct;
  isPlugin["b"] = true;
}


void FirstOrderLinearDS::computeA(const double time)
{
  if (isPlugin["A"])
  {
    if (APtr == NULL)
      RuntimeException::selfThrow("computeA() is not linked to a plugin function");
    APtr(time, n, &(*A)(0, 0), z->size(), &(*z)(0));
  }
  // else nothing
}

void FirstOrderLinearDS::computeB(const double time)
{
  if (isPlugin["b"])
  {
    if (bPtr == NULL)
      RuntimeException::selfThrow("computeB() is not linked to a plugin function");
    bPtr(time, n, &(*b)(0), z->size(), &(*z)(0));
  }
  // else nothing
}

void FirstOrderLinearDS::computeRhs(const double time, const bool)
{
  // second argument is useless at the time - Used in derived classes
  // compute A=jacobianXF

  *x[1] = * r;

  if (A != NULL)
  {
    computeA(time);
    prod(*A, *x[0], *x[1], false);
  }

  // compute and add b if required
  if (b != NULL)
  {
    computeB(time);
    *x[1] += *b;
  }

  if (M != NULL)
  {
    // allocate invM at the first call of the present function
    if (invM == NULL)
    {
      invM = new SimpleMatrix(*M);
      isAllocatedIn["invM"] = true;
    }

    invM->PLUForwardBackwardInPlace(*x[1]);
  }
}

void FirstOrderLinearDS::computeJacobianXRhs(const double time, const bool)
{
  if (isPlugin["A"])
    computeA(time);

  if (M != NULL && A != NULL)
  {
    *jacobianXRhs = *A;
    // copy M into invM for LU-factorisation, at the first call of this function.
    if (invM == NULL)
    {
      invM = new SimpleMatrix(*M);
      isAllocatedIn["invM"] = true;
    }
    // solve MjacobianXRhs = A
    invM->PLUForwardBackwardInPlace(*jacobianXRhs);
  }
  // else jacobianXRhs = A, pointers equality.

}

void FirstOrderLinearDS::display() const
{
  cout << "=== Linear system display, " << number << ", id: " << id << ")." << endl;
  NamesConstIterator it;
  cout << "The following operators are linked to plug-in: " << endl;
  for (it = pluginNames.begin(); it != pluginNames.end(); ++it)
    cout << (*it).first << " plugged to:" << (*it).second << endl;
  cout << "=============================" << endl;
}

void FirstOrderLinearDS::saveSpecificDataToXML()
{
  if (dsxml == NULL)
    RuntimeException::selfThrow("FirstOrderLinearDS::saveDSToXML - The DynamicalSystemXML object doesn't exists");
  static_cast<FirstOrderLinearDSXML*>(dsxml)->setA(*A);

  // b
  if (b != NULL)
  {
    if (!(static_cast <FirstOrderLinearDSXML*>(dsxml))->isBPlugin())
      static_cast<FirstOrderLinearDSXML*>(dsxml)->setB(*b);
  }

  else RuntimeException::selfThrow("FirstOrderLinearDS::saveDSToXML - The DynamicalSystemXML object doesn't exists");
}

FirstOrderLinearDS* FirstOrderLinearDS::convert(DynamicalSystem* ds)
{
  FirstOrderLinearDS* lsds = dynamic_cast<FirstOrderLinearDS*>(ds);
  return lsds;
}
