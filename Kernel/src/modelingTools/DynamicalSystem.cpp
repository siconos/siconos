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
#include "DynamicalSystem.hpp"
#include "DynamicalSystemXML.hpp"
#include "BlockVector.hpp"
#include "Plugin.hpp"

using namespace std;
using namespace DS;

unsigned int DynamicalSystem::count = 0;

// ===== CONSTRUCTORS =====

// Default constructor (protected)
DynamicalSystem::DynamicalSystem(DS::TYPES type):
  _DSType(type), _number(count++), _n(0), _stepsInMemory(1)
{
  zeroPlugin();
  _normRef = 1;
  _x.resize(2);
  _workV.resize(sizeWorkV);
}

// From XML file
DynamicalSystem::DynamicalSystem(SP::DynamicalSystemXML dsXML):
  _DSType(dsXML->getType()), _number(dsXML->number()), _n(0), _stepsInMemory(1), _dsxml(dsXML)
{
  _normRef = 1;
  assert(dsXML && "DynamicalSystem::DynamicalSystem - DynamicalSystemXML paramater must not be NULL");

  // Update count: must be at least equal to number for future DS creation
  count = number();
  zeroPlugin();
  // Only the following data are set in this general constructor:
  //  - DSTye
  //  - _number
  //  - z
  //  - stepsInMemory
  // All others are dependent of the derived class type.

  // === Initial conditions ===
  // Warning: n is set thanks to vector of initial conditions size. That will be set in derived classes constructor.

  _x.resize(2);

  // z - Optional parameter.
  if (_dsxml->hasZ())
    _z.reset(new SimpleVector(_dsxml->getZ()));

  if (_dsxml->hasStepsInMemory()) _stepsInMemory = _dsxml->getStepsInMemory();
  _workV.resize(sizeWorkV);
}

// From a minimum set of data
DynamicalSystem::DynamicalSystem(DS::TYPES type, unsigned int newN):
  _DSType(type), _number(count++), _n(newN), _stepsInMemory(1)
{
  _normRef = 1;
  _x.resize(2);
  _workV.resize(sizeWorkV);
  _residuFree.reset(new SimpleVector(getDim()));
  _r.reset(new SimpleVector(getDim()));
}


bool DynamicalSystem::checkDynamicalSystem()
{
  bool output = true;
  // n
  if (_n == 0)
  {
    RuntimeException::selfThrow("DynamicalSystem::checkDynamicalSystem - number of degrees of freedom is equal to 0.");
    output = false;
  }
  if (!_x0)
  {
    RuntimeException::selfThrow("DynamicalSystem::checkDynamicalSystem - x0 not set.");
    output = false;
  }
  return output;
}
void DynamicalSystem::zeroPlugin()
{
  _pluginJacgx.reset(new PluggedObject());
  _pluginJacxDotG.reset(new PluggedObject());
  _pluging.reset(new PluggedObject());
}
// Setters

void DynamicalSystem::setX0(const SiconosVector& newValue)
{
  // check dimensions ...
  if (newValue.size() != _n)
    RuntimeException::selfThrow("DynamicalSystem::setX0 - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (_x0)
    *_x0 = newValue;

  else
  {
    if (newValue.isBlock())
      _x0.reset(new BlockVector(newValue));
    else
      _x0.reset(new SimpleVector(newValue));
  }
  _normRef = _x0->norm2() + 1;
}

void DynamicalSystem::setX0Ptr(SP::SiconosVector newPtr)
{
  // check dimensions ...
  assert(newPtr->size() == _n && "DynamicalSystem::setX0Ptr - inconsistent sizes between x0 input and n - Maybe you forget to set n?");
  _x0 = newPtr;
  _normRef = _x0->norm2() + 1;
}

void DynamicalSystem::setX(const SiconosVector& newValue)
{
  // Warning: this only sets the value of x[0]
  // We suppose that both x and (*x)[0] are properly allocated.

  // check dimensions ...
  if (newValue.size() != _n)
    RuntimeException::selfThrow("DynamicalSystem::setX - inconsistent sizes between x input and n - Maybe you forget to set n?");

  if (! _x[0])
    _x[0].reset(new SimpleVector(newValue));
  else
    *(_x[0]) = newValue;
}

void DynamicalSystem::setXPtr(SP::SiconosVector newPtr)
{
  // Warning: this only sets the pointer x[0]

  // check dimensions ...
  if (newPtr->size() != _n)
    RuntimeException::selfThrow("DynamicalSystem::setXPtr - inconsistent sizes between x input and n - Maybe you forget to set n?");

  _x[0] = newPtr;
}

void DynamicalSystem::setRhs(const SiconosVector& newValue)
{
  // Warning: this only sets the value of x[1]

  // check dimensions ...
  if (newValue.size() != _n)
    RuntimeException::selfThrow("DynamicalSystem::setRhs - inconsistent sizes between x input and n - Maybe you forget to set n?");

  if (! _x[1])
    _x[1].reset(new SimpleVector(newValue));
  else
    *(_x[1]) = newValue;
}

void DynamicalSystem::setRhsPtr(SP::SiconosVector newPtr)
{
  // Warning: this only sets the pointer (*x)[1]

  // check dimensions ...
  if (newPtr->size() != _n)
    RuntimeException::selfThrow("DynamicalSystem::setRhsPtr - inconsistent sizes between x input and n - Maybe you forget to set n?");

  _x[1] = newPtr;
}
void DynamicalSystem::setR(const SiconosVector& newValue)
{
  // check dimensions ...
  if (newValue.size() != _n)
    RuntimeException::selfThrow("DynamicalSystem::setR - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  if (_r)
    *_r = newValue;

  else
    _r.reset(new SimpleVector(newValue));
}

void DynamicalSystem::setRPtr(SP::SiconosVector newPtr)
{
  // check dimensions ...
  if (newPtr->size() != _n)
    RuntimeException::selfThrow("DynamicalSystem::setRPtr - inconsistent sizes between x0 input and n - Maybe you forget to set n?");

  _r = newPtr;

}

void DynamicalSystem::setJacobianRhsx(const SiconosMatrix& newValue)
{
  // check dimensions ...
  if (newValue.size(0) != _n || newValue.size(1) != _n)
    RuntimeException::selfThrow("DynamicalSystem::setJacobianRhsx - inconsistent sizes between jacobianRhsx input and n - Maybe you forget to set n?");

  if (_jacxRhs)
    *_jacxRhs = newValue;

  else
    _jacxRhs.reset(new SimpleMatrix(newValue));
}

void DynamicalSystem::setJacobianRhsxPtr(SP::SiconosMatrix newPtr)
{
  // check dimensions ...
  if (newPtr->size(0) != _n || newPtr->size(1) != _n)
    RuntimeException::selfThrow("DynamicalSystem::setJacobianRhsxPtr - inconsistent sizes between _jacxRhs input and n - Maybe you forget to set n?");

  _jacxRhs = newPtr;
}

void DynamicalSystem::setZ(const SiconosVector& newValue)
{
  if (_z)
  {
    if (newValue.size() != _z->size())
      RuntimeException::selfThrow("DynamicalSystem::setZ - inconsistent sizes between input and existing z - To change z size use setZPtr.");
    *_z = newValue;
  }
  else
  {

    if (newValue.isBlock())
      _z.reset(new BlockVector(newValue));
    else
      _z.reset(new SimpleVector(newValue));
  }
}

void DynamicalSystem::setZPtr(SP::SiconosVector newPtr)
{
  _z = newPtr;
}
/*
void DynamicalSystem::setG(const PVFInt& newValue)
{
  assert(newValue.size()==n&&"DynamicalSystem - setG: inconsistent dimensions with problem size for input vector g");

  if( !g  )
    g.reset(new PVFInt(newValue));
  else
    *g = newValue;
}
*/
/*
void DynamicalSystem::setJacobianG(unsigned int i, const PMFInt& newValue)
{
  assert(newValue.size(0)==n&&"DynamicalSystem - setJacobianG: inconsistent dimensions with problem size for matrix jacobianG.");
  assert(newValue.size(1)==n&&"DynamicalSystem - setJacobianG: inconsistent dimensions with problem size for matrix jacobianG.");

  if( !jacobianG [i] )
    jacobianG[i].reset(new PMFInt(newValue));
  else
    *jacobianG[i] = newValue;
}
*/

void DynamicalSystem::update(double time)
{
  computeRhs(time);
  computeJacobianRhsx(time);
}

// ===== MEMORY MANAGEMENT FUNCTIONS =====

void DynamicalSystem::initMemory(unsigned int steps)
{
  if (steps == 0)
    cout << "Warning : FirstOrderNonLinearDS::initMemory with size equal to zero" << endl;
  else
  {
    _stepsInMemory = steps;
    _xMemory.reset(new SiconosMemory(steps));

  }
}

void DynamicalSystem::setComputegFunction(const string& pluginPath, const string& functionName)
{
  _pluging->setComputeFunction(pluginPath, functionName);
}

void DynamicalSystem::setComputegFunction(FPtr6 fct)
{
  _pluging->setComputeFunction((void *)fct);
}

void DynamicalSystem::setComputeJacobianXGFunction(const string& pluginPath, const string& functionName)
{
  _pluginJacgx->setComputeFunction(pluginPath, functionName);
}
void DynamicalSystem::setComputeJacobianDotXGFunction(const string& pluginPath, const string& functionName)
{
  _pluginJacxDotG->setComputeFunction(pluginPath, functionName);
}
// void DynamicalSystem::setComputeJacobianZGFunction( const string& pluginPath, const string& functionName){
//   Plugin::setFunction(&pluginJacobianZGPtr, pluginPath,functionName);
// }

void DynamicalSystem::computeg(double time)
{
  if (_pluging->fPtr)
    ((FPtr6)(_pluging->fPtr))(time, _n, &(*_x[0])(0), &(*_x[1])(0), &(*_g)(0), _z->size(), &(*_z)(0));
}

void DynamicalSystem::computeJacobianXG(double time)
{
  if (_pluginJacgx->fPtr)
    ((FPtr6) _pluginJacgx->fPtr)(time, _n, &(*_x[0])(0), &(*_x[1])(0), &(*_jacgx)(0, 0), _z->size(), &(*_z)(0));
}
void DynamicalSystem::computeJacobianDotXG(double time)
{
  if (_pluginJacxDotG->fPtr)
    ((FPtr6)(_pluginJacxDotG->fPtr))(time, _n, &(*_x[0])(0), &(*_x[1])(0), &(*_jacxDotG)(0, 0), _z->size(), &(*_z)(0));
}
// void DynamicalSystem::computeJacobianZG(double time){
//   if (pluginJacobianXGPtr)
//     pluginJacobianZGPtr(time, n, &(*x[0])(0), &(*x[1])(0), &(*jacobianG[i])(0,0), z->size(), &(*z)(0));
// }

// ===== XML MANAGEMENT FUNCTIONS =====

void DynamicalSystem::saveDSToXML()
{
  RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - Not yet tested for DS: do not use it.");

  if (!_dsxml)
    RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - The DynamicalSystemXML object doesn't exists");

  // --- Specific (depending on derived class) DS data ---
  saveSpecificDataToXML();

  // --- other data ---
}

// ===== MISCELLANEOUS ====

double DynamicalSystem::dsConvergenceIndicator()
{
  RuntimeException::selfThrow
  ("DynamicalSystem:dsConvergenceIndicator - not yet implemented for Dynamical system type :"
   + _DSType);
  return 1.0;
}

