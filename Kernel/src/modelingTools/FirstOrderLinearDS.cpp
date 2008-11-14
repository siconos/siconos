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

// --- Constructors ---

// From xml file
FirstOrderLinearDS::FirstOrderLinearDS(SP::DynamicalSystemXML dsXML): FirstOrderNonLinearDS(dsXML)
{
  // pointer to xml
  SP::FirstOrderLinearDSXML foldsxml = (boost::static_pointer_cast <FirstOrderLinearDSXML>(dsxml));

  // Check if f is given as a plug-in in xml input file.
  if (foldsxml->hasF() || foldsxml->hasJacobianXF())
    RuntimeException::selfThrow("FirstOrderLinearDS - xml constructor, you give a f or its jacobian as a plug-in for a FirstOrderLinearDS -> set rather A and b plug-in.");

  string plugin;
  // A
  if (foldsxml->hasA())
  {
    if (foldsxml->isAPlugin())
    {
      plugin = foldsxml->getAPlugin();
      setComputeAFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      A.reset(new Plugged_Matrix_FTime(foldsxml->getA()));
  }

  // b
  if (foldsxml->hasB())
  {
    if (foldsxml->isBPlugin())
    {
      plugin = foldsxml->getBPlugin();
      setComputeBFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      b.reset(new Plugged_Vector_FTime(foldsxml->getBVector()));
  }

  checkDynamicalSystem();
}

// From a minimum set of data, A and b connected to a plug-in
FirstOrderLinearDS::FirstOrderLinearDS(const SiconosVector& newX0, const string& APlugin, const string& bPlugin):
  FirstOrderNonLinearDS(newX0)
{
  DSType = DS::FOLDS;
  setComputeAFunction(SSL::getPluginName(APlugin), SSL::getPluginFunctionName(APlugin));
  setComputeBFunction(SSL::getPluginName(bPlugin), SSL::getPluginFunctionName(bPlugin));

  checkDynamicalSystem();
}

// From a minimum set of data, A from a given matrix
FirstOrderLinearDS::FirstOrderLinearDS(const SiconosVector& newX0, const SiconosMatrix& newA):
  FirstOrderNonLinearDS(newX0)
{
  DSType = DS::FOLDS;
  assert(((newA.size(0) == n) && (newA.size(1) == n)) &&
         "FirstOrderLinearDS - constructor(number,x0,A): inconsistent dimensions with problem size for input matrix A");

  A.reset(new Plugged_Matrix_FTime(newA));
  checkDynamicalSystem();
}

// From a minimum set of data, A from a given matrix
FirstOrderLinearDS::FirstOrderLinearDS(const SiconosVector& newX0, const SiconosMatrix& newA, const SiconosVector& newB):
  FirstOrderNonLinearDS(newX0)
{
  DSType = DS::FOLDS;
  assert(((newA.size(0) == n) && (newA.size(1) == n)) &&
         "FirstOrderLinearDS - constructor(x0,A,b): inconsistent dimensions with problem size for input matrix A");
  assert(newB.size() == n &&
         "FirstOrderLinearDS - constructor(x0,A,b): inconsistent dimensions with problem size for input vector b ");

  A.reset(new Plugged_Matrix_FTime(newA));
  b.reset(new Plugged_Vector_FTime(newB));

  checkDynamicalSystem();
}

bool FirstOrderLinearDS::checkDynamicalSystem() // useless ...?
{
  bool output = DynamicalSystem::checkDynamicalSystem();
  if (!output) cout << "FirstOrderLinearDS Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
  return output;
}

void FirstOrderLinearDS::initRhs(double time)
{
  computeRhs(time); // If necessary, this will also compute A and b.
  if (! jacobianXRhs)  // if not allocated with a set or anything else
  {
    if (A && ! M)  // if M is not defined, then A = jacobianXRhs, no memory allocation for that one.
      jacobianXRhs = A;
    else if (A && M)
      jacobianXRhs.reset(new SimpleMatrix(n, n));
    // else no allocation, jacobian is equal to 0.
  }
  computeJacobianXRhs(time);
}

void FirstOrderLinearDS::updatePlugins(double time)
{
  computeA(time);
  computeB(time);
}

void FirstOrderLinearDS::setA(const Plugged_Matrix_FTime& newValue)
{
  assert(newValue.size(0) == n && "FirstOrderLinearDS - setA: inconsistent dimensions with problem size for input matrix A.");
  assert(newValue.size(1) == n && "FirstOrderLinearDS - setA: inconsistent dimensions with problem size for input matrix A.");

  if (! A)
    A.reset(new Plugged_Matrix_FTime(newValue));
  else
    *A = newValue;
}

void FirstOrderLinearDS::setB(const Plugged_Vector_FTime& newValue)
{
  assert(newValue.size() == n && "FirstOrderLinearDS - setB: inconsistent dimensions with problem size for input vector b");

  if (! b)
    b.reset(new Plugged_Vector_FTime(newValue));
  else
    *b = newValue;
}

void FirstOrderLinearDS::setComputeAFunction(const string& pluginPath, const string& functionName)
{
  if (!A)
    A.reset(new Plugged_Matrix_FTime(n, n));
  A->setComputeFunction(pluginPath, functionName);
}

void FirstOrderLinearDS::setComputeAFunction(MatrixFunctionOfTime fct)
{
  if (!A)
    A.reset(new Plugged_Matrix_FTime(n, n));
  A->setComputeFunction(fct);
}
void FirstOrderLinearDS::setComputeBFunction(const string& pluginPath, const string& functionName)
{
  if (! b)
    b.reset(new Plugged_Vector_FTime(n));
  b->setComputeFunction(pluginPath, functionName);
}

void FirstOrderLinearDS::setComputeBFunction(VectorFunctionOfTime fct)
{
  if (! b)
    b.reset(new Plugged_Vector_FTime(n));
  b->setComputeFunction(fct);
}


void FirstOrderLinearDS::computeA(const double time)
{
  if (A->isPlugged())
  {
    if (!A->fPtr)
      RuntimeException::selfThrow("computeA() is not linked to a plugin function");
    (A->fPtr)(time, n, n, &(*A)(0, 0), z->size(), &(*z)(0));
  }
  // else nothing
}

void FirstOrderLinearDS::computeB(const double time)
{
  if (b->isPlugged())
  {
    if (!b->fPtr)
      RuntimeException::selfThrow("computeB() is not linked to a plugin function");
    (b->fPtr)(time, n, &(*b)(0), z->size(), &(*z)(0));
  }
  // else nothing
}

void FirstOrderLinearDS::computeRhs(const double time, const bool)
{
  // second argument is useless at the time - Used in derived classes
  // compute A=jacobianXF

  *x[1] = * r;

  if (A)
  {
    computeA(time);
    prod(*A, *x[0], *x[1], false);
  }

  // compute and add b if required
  if (b)
  {
    computeB(time);
    *x[1] += *b;
  }

  if (M)
  {
    // allocate invM at the first call of the present function
    if (! invM)
      invM.reset(new SimpleMatrix(*M));

    invM->PLUForwardBackwardInPlace(*x[1]);
  }
}

void FirstOrderLinearDS::computeJacobianXRhs(const double time, const bool)
{
  computeA(time);

  if (M && A)
  {
    *jacobianXRhs = *A;
    // copy M into invM for LU-factorisation, at the first call of this function.
    if (! invM)
      invM.reset(new SimpleMatrix(*M));
    // solve MjacobianXRhs = A
    invM->PLUForwardBackwardInPlace(*jacobianXRhs);
  }
  // else jacobianXRhs = A, pointers equality.

}

void FirstOrderLinearDS::display() const
{
  cout << "=== Linear system display, " << number << endl;
  cout << "=============================" << endl;
}

void FirstOrderLinearDS::saveSpecificDataToXML()
{
  if (!dsxml)
    RuntimeException::selfThrow("FirstOrderLinearDS::saveDSToXML - The DynamicalSystemXML object doesn't exists");
  boost::static_pointer_cast<FirstOrderLinearDSXML>(dsxml)->setA(*A);

  // b
  if (b)
  {
    if (!(boost::static_pointer_cast <FirstOrderLinearDSXML>(dsxml))->isBPlugin())
      boost::static_pointer_cast<FirstOrderLinearDSXML>(dsxml)->setB(*b);
  }

  else RuntimeException::selfThrow("FirstOrderLinearDS::saveDSToXML - The DynamicalSystemXML object doesn't exists");
}

FirstOrderLinearDS* FirstOrderLinearDS::convert(DynamicalSystem* ds)
{
  FirstOrderLinearDS* lsds = dynamic_cast<FirstOrderLinearDS*>(ds);
  return lsds;
}
