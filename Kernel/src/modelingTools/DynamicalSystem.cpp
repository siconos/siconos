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
#include "DynamicalSystem.h"
#include "DynamicalSystemXML.h"
#include "BlockVector.h"

using namespace std;
using namespace DS;

unsigned int DynamicalSystem::count = 0;

// ===== CONSTRUCTORS =====

// Default constructor (protected)
DynamicalSystem::DynamicalSystem(DS::TYPES type):
  DSType(type), number(count++), n(0), stepsInMemory(1)
{
  mNormRef = 1;
  x.resize(2);
  workV.resize(sizeWorkV);
}

// From XML file
DynamicalSystem::DynamicalSystem(SP::DynamicalSystemXML dsXML):
  DSType(dsXML->getType()), number(dsXML->getNumber()), n(0), stepsInMemory(1), dsxml(dsXML)
{
  mNormRef = 1;
  assert(dsXML && "DynamicalSystem::DynamicalSystem - DynamicalSystemXML paramater must not be NULL");

  // Update count: must be at least equal to number for future DS creation
  count = number;

  // Only the following data are set in this general constructor:
  //  - DSTye
  //  - number
  //  - z
  //  - stepsInMemory
  // All others are dependent of the derived class type.

  // === Initial conditions ===
  // Warning: n is set thanks to vector of initial conditions size. That will be set in derived classes constructor.

  x.resize(2);

  // z - Optional parameter.
  if (dsxml->hasZ())
    z.reset(new SimpleVector(dsxml->getZ()));

  if (dsxml->hasStepsInMemory()) stepsInMemory = dsxml->getStepsInMemory();
  workV.resize(sizeWorkV);
}

// From a minimum set of data
DynamicalSystem::DynamicalSystem(DS::TYPES type, unsigned int newN):
  DSType(type), number(count++), n(newN), stepsInMemory(1)
{
  mNormRef = 1;
  x.resize(2);
  workV.resize(sizeWorkV);
  mResiduFree.reset(new SimpleVector(getDim()));
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
  if (!x0)
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

  if (x0)
    *x0 = newValue;

  else
  {
    if (newValue.isBlock())
      x0.reset(new BlockVector(newValue));
    else
      x0.reset(new SimpleVector(newValue));
  }
  mNormRef = x0->norm2() + 1;
}

void DynamicalSystem::setX0Ptr(SP::SiconosVector newPtr)
{
  // check dimensions ...
  assert(newPtr->size() == n && "DynamicalSystem::setX0Ptr - inconsistent sizes between x0 input and n - Maybe you forget to set n?");
  x0 = newPtr;
  mNormRef = x0->norm2() + 1;
}

void DynamicalSystem::setX(const SiconosVector& newValue)
{
  // Warning: this only sets the value of x[0]
  // We suppose that both x and (*x)[0] are properly allocated.

  // check dimensions ...
  if (newValue.size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setX - inconsistent sizes between x input and n - Maybe you forget to set n?");

  if (! x[0])
    x[0].reset(new SimpleVector(newValue));
  else
    *(x[0]) = newValue;
}

void DynamicalSystem::setXPtr(SP::SiconosVector newPtr)
{
  // Warning: this only sets the pointer x[0]

  // check dimensions ...
  if (newPtr->size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setXPtr - inconsistent sizes between x input and n - Maybe you forget to set n?");

  x[0] = newPtr;
}

void DynamicalSystem::setRhs(const SiconosVector& newValue)
{
  // Warning: this only sets the value of x[1]

  // check dimensions ...
  if (newValue.size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setRhs - inconsistent sizes between x input and n - Maybe you forget to set n?");

  if (! x[1])
    x[1].reset(new SimpleVector(newValue));
  else
    *(x[1]) = newValue;
}

void DynamicalSystem::setRhsPtr(SP::SiconosVector newPtr)
{
  // Warning: this only sets the pointer (*x)[1]

  // check dimensions ...
  if (newPtr->size() != n)
    RuntimeException::selfThrow("DynamicalSystem::setRhsPtr - inconsistent sizes between x input and n - Maybe you forget to set n?");

  x[1] = newPtr;
}

void DynamicalSystem::setJacobianXRhs(const SiconosMatrix& newValue)
{
  // check dimensions ...
  if (newValue.size(0) != n || newValue.size(1) != n)
    RuntimeException::selfThrow("DynamicalSystem::setJacobianXRhs - inconsistent sizes between jacobianXRhs input and n - Maybe you forget to set n?");

  if (jacobianXRhs)
    *jacobianXRhs = newValue;

  else
    jacobianXRhs.reset(new SimpleMatrix(newValue));
}

void DynamicalSystem::setJacobianXRhsPtr(SP::SiconosMatrix newPtr)
{
  // check dimensions ...
  if (newPtr->size(0) != n || newPtr->size(1) != n)
    RuntimeException::selfThrow("DynamicalSystem::setJacobianXRhsPtr - inconsistent sizes between jacobianXRhs input and n - Maybe you forget to set n?");

  jacobianXRhs = newPtr;
}

void DynamicalSystem::setZ(const SiconosVector& newValue)
{
  if (z)
  {
    if (newValue.size() != z->size())
      RuntimeException::selfThrow("DynamicalSystem::setZ - inconsistent sizes between input and existing z - To change z size use setZPtr.");
    *z = newValue;
  }
  else
  {

    if (newValue.isBlock())
      z.reset(new BlockVector(newValue));
    else
      z.reset(new SimpleVector(newValue));
  }
}

void DynamicalSystem::setZPtr(SP::SiconosVector newPtr)
{
  z = newPtr;
}

void DynamicalSystem::setG(const PVFint& newValue)
{
  assert(newValue.size() == n && "DynamicalSystem - setG: inconsistent dimensions with problem size for input vector g");

  if (!g)
    g.reset(new PVFint(newValue));
  else
    *g = newValue;
}

void DynamicalSystem::setJacobianG(unsigned int i, const PMFint& newValue)
{
  assert(newValue.size(0) == n && "DynamicalSystem - setJacobianG: inconsistent dimensions with problem size for matrix jacobianG.");
  assert(newValue.size(1) == n && "DynamicalSystem - setJacobianG: inconsistent dimensions with problem size for matrix jacobianG.");

  if (!jacobianG [i])
    jacobianG[i].reset(new PMFint(newValue));
  else
    *jacobianG[i] = newValue;
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
    xMemory.reset(new SiconosMemory(steps));

  }
}

void DynamicalSystem::setComputeGFunction(const string& pluginPath, const string& functionName)
{
  if (! g)
    g.reset(new PVFint(n));
  g->setComputeFunction(pluginPath, functionName);
}

void DynamicalSystem::setComputeGFunction(FPtr6 fct)
{
  if (! g)
    g.reset(new PVFint(n));
  g->setComputeFunction(fct);
}

void DynamicalSystem::setComputeJacobianGFunction(unsigned int i, const string& pluginPath, const string& functionName)
{
  if (! jacobianG[i])
    jacobianG[i].reset(new PMFint(n, n));
  jacobianG[i]->setComputeFunction(pluginPath, functionName);
}

void DynamicalSystem::setComputeJacobianGFunction(unsigned int i, FPtr6 fct)
{
  if (! jacobianG[i])
    jacobianG[i].reset(new PMFint(n, n));
  jacobianG[i]->setComputeFunction(fct);
}

void DynamicalSystem::computeG(double time)
{
  if (g->isPlugged())
  {
    if (!(g->fPtr))
      RuntimeException::selfThrow("DynamicalSystem::computeG() is not linked to a plugin function");

    (g->fPtr)(time, n, &(*x[0])(0), &(*x[1])(0), &(*g)(0), z->size(), &(*z)(0));
  }
  else RuntimeException::selfThrow("DynamicalSystem::computeG - Not yet implemented for DS of type " + DSType);
}

void DynamicalSystem::computeJacobianG(unsigned int i, double time)
{
  if (jacobianG[i]->isPlugged())
  {
    if (!(jacobianG[i]->fPtr))
      RuntimeException::selfThrow("computeJacobianG(i,time) is not linked to a plugin function. i=" + i);

    (jacobianG[i]->fPtr)(time, n, &(*x[0])(0), &(*x[1])(0), &(*jacobianG[i])(0, 0), z->size(), &(*z)(0));
  }
  else RuntimeException::selfThrow("DynamicalSystem::computeJacobianG - Not yet implemented for DS of type " + DSType);
}

// ===== XML MANAGEMENT FUNCTIONS =====

void DynamicalSystem::saveDSToXML()
{
  RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - Not yet tested for DS: do not use it.");

  if (!dsxml)
    RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - The DynamicalSystemXML object doesn't exists");

  // --- Specific (depending on derived class) DS data ---
  saveSpecificDataToXML();

  // --- other data ---
}

// ===== MISCELLANEOUS ====

double DynamicalSystem::dsConvergenceIndicator()
{
  RuntimeException::selfThrow("DynamicalSystem:dsConvergenceIndicator - not yet implemented for Dynamical system type :" + DSType);
  return 1.0;
}

