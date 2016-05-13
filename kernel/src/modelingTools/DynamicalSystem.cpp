/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include "DynamicalSystem.hpp"

#include <iostream>


unsigned int DynamicalSystem::count = 0;

// ===== CONSTRUCTORS =====

// Default constructor (protected)
DynamicalSystem::DynamicalSystem():
  _number(count++), _n(0), _stepsInMemory(1)
{
  zeroPlugin();
  _normRef = 1;
  _x.resize(2);
  _workspace.resize(sizeWorkV);
  _z.reset(new SiconosVector(1));
}

// From a minimum set of data
DynamicalSystem::DynamicalSystem(unsigned int newN):
  _number(count++), _n(newN), _stepsInMemory(1)
{
  zeroPlugin();
  _normRef = 1;
  _x.resize(2);
  _workspace.resize(sizeWorkV);
  _workspace[freeresidu].reset(new SiconosVector(getDim()));
  _r.reset(new SiconosVector(getDim()));
  _z.reset(new SiconosVector(1));
}

// Copy constructor
DynamicalSystem::DynamicalSystem(const DynamicalSystem & ds):
  _number(count++)
{
  // The following data should always be initialize
  _n = ds.getN();
  _normRef = ds.normRef();
  _x0.reset(new SiconosVector(*(ds.x0())));
  _workspace.resize(sizeWorkV);
  _workspace[freeresidu].reset(new SiconosVector(*(ds.workspace(freeresidu))));
  _r.reset(new SiconosVector(*(ds.r())));
  _x.resize(2);
  _x[0].reset(new SiconosVector(*(ds.x())));
  _x[1].reset(new SiconosVector(*(ds.rhs())));

  // These  were not always initialised
  if (ds.jacobianRhsx())
    _jacxRhs.reset(new SimpleMatrix(*(ds.jacobianRhsx())));
  //  if (ds.jacobianXG())
  //    _jacgx.reset(new SimpleMatrix(*(ds.jacobianXG())));
  //  if (ds.jacobianXDotG())
  //    _jacxDotG.reset(new SimpleMatrix(*(ds.jacobianXDotG())));
  if (ds.z())
    _z.reset(new SiconosVector(*(ds.z())));

  if (_pluging)
    _pluging.reset(new PluggedObject(*(ds.getPluginG())));
  if (_pluginJacgx)
    _pluginJacgx.reset(new PluggedObject(*(ds.getPluginJacGX())));
  if (_pluginJacxDotG)
    _pluginJacxDotG.reset(new PluggedObject(*(ds.getPluginJacXDotG())));

  if (ds.xMemory())
    _xMemory.reset(new SiconosMemory(*(ds.xMemory())));
  _stepsInMemory = ds.getStepsInMemory();

  _workspace.resize(sizeWorkV);

  if (ds.workspace(local_buffer))
    _workspace[local_buffer].reset(new SiconosVector(*(ds.workspace(local_buffer))));
  if (ds.workspace(free))
    _workspace[free].reset(new SiconosVector(*(ds.workspace(free))));
  //  _workspace[sizeWorkV].reset(new SiconosVector(*(ds.workspace(sizeWorkV))));
  // XXX See how to implement the copy of _workMatrix

//  _workFree.reset(new SiconosVector(*(ds.workspace(DynamicalSystem::free))));
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
    _x0.reset(new SiconosVector(newValue));
  }
  _normRef = _x0->norm2() + 1;
}

void DynamicalSystem::setX0Ptr(SP::SiconosVector newPtr)
{
  // check dimensions ...
  if (newPtr->size() != _n)
    RuntimeException::selfThrow("DynamicalSystem::setX0Ptr - inconsistent sizes between x0 input and n - Maybe you forget to set n?");
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
    _x[0].reset(new SiconosVector(newValue));
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
    _x[1].reset(new SiconosVector(newValue));
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
    _r.reset(new SiconosVector(newValue));
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

void DynamicalSystem::setz(const SiconosVector& newValue)
{
  if (_z)
  {
    if (newValue.size() != _z->size())
      RuntimeException::selfThrow("DynamicalSystem::setz - inconsistent sizes between input and existing z - To change z size use setzPtr.");
    *_z = newValue;
  }
  else
  {
    _z.reset(new SiconosVector(newValue));
  }
}

void DynamicalSystem::setzPtr(SP::SiconosVector newPtr)
{
  _z = newPtr;
}

void DynamicalSystem::update(double time)
{
  computeRhs(time);
  computeJacobianRhsx(time);
}

// ===== MEMORY MANAGEMENT FUNCTIONS =====

void DynamicalSystem::initMemory(unsigned int steps)
{

  if (steps == 0)
    std::cout << "Warning : initMemory with size equal to zero" << std::endl;
  else
  {
    _stepsInMemory = steps;
    _xMemory.reset(new SiconosMemory(steps, _n));
  }

}

void DynamicalSystem::setComputegFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluging->setComputeFunction(pluginPath, functionName);
}

void DynamicalSystem::setComputegFunction(FPtr6 fct)
{
  _pluging->setComputeFunction((void *)fct);
}

void DynamicalSystem::setComputeJacobianXGFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginJacgx->setComputeFunction(pluginPath, functionName);
}
void DynamicalSystem::setComputeJacobianDotXGFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginJacxDotG->setComputeFunction(pluginPath, functionName);
}
// void DynamicalSystem::setComputeJacobianZGFunction( const std::string& pluginPath, const std::string& functionName){
//   Plugin::setFunction(&pluginJacobianZGPtr, pluginPath,functionName);
// }

//void DynamicalSystem::computeg(double time)
//{
//  if (_pluging->fPtr)
//    ((FPtr6)(_pluging->fPtr))(time, _n, &(*_x[0])(0), &(*_x[1])(0), &(*_g)(0), _z->size(), &(*_z)(0));
//}

//void DynamicalSystem::computeJacobianXG(double time){
//  if (_pluginJacgx->fPtr)
//    ((FPtr6) _pluginJacgx->fPtr)(time, _n, &(*_x[0])(0), &(*_x[1])(0), &(*_jacgx)(0,0), _z->size(), &(*_z)(0));
//}
//void DynamicalSystem::computeJacobianDotXG(double time){
//  if (_pluginJacxDotG->fPtr)
//    ((FPtr6) (_pluginJacxDotG->fPtr))(time, _n, &(*_x[0])(0), &(*_x[1])(0), &(*_jacxDotG)(0,0), _z->size(), &(*_z)(0));
//}
// void DynamicalSystem::computeJacobianZG(double time){
//   if (pluginJacobianXGPtr)
//     pluginJacobianZGPtr(time, n, &(*x[0])(0), &(*x[1])(0), &(*jacobianG[i])(0,0), z->size(), &(*z)(0));
// }

// ===== MISCELLANEOUS ====

double DynamicalSystem::dsConvergenceIndicator()
{
  RuntimeException::selfThrow
  ("DynamicalSystem:dsConvergenceIndicator - not yet implemented for this Dynamical system type");
  return 1.0;
}

