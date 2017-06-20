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

// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include "debug.h"

#include <iostream>


unsigned int DynamicalSystem::__count = 0;

void DynamicalSystem::_init()
{
  DEBUG_PRINT("internal _init from DynamicalSystem\n");

  // No memory allocation, only resize for containers.
  // Everything should be done in derived class init for required operators
  // and variables and in 'set'-like methods for optional
  // components.
  
  _stepsInMemory = 1;
  _x.resize(2);
  _z.reset(new SiconosVector(1));
}

// ===== CONSTRUCTORS =====

// Default constructor (protected)
DynamicalSystem::DynamicalSystem():
  _number(__count++), _n(0)
{
  _init();
}

// From a minimum set of data
DynamicalSystem::DynamicalSystem(unsigned int dimension):
  _number(__count++), _n(dimension)
{
  _init();
}

// Copy constructor
DynamicalSystem::DynamicalSystem(const DynamicalSystem & ds):
  _number(__count++), _n(ds.n()), _stepsInMemory(ds.stepsInMemory())
{
  // The following data should always be initialize
  if(ds.x0())
    _x0.reset(new SiconosVector(*(ds.x0())));
  if(ds.r())
    _r.reset(new SiconosVector(*(ds.r())));
  _x.resize(2);
  if(ds.x())
     _x[0].reset(new SiconosVector(*(ds.x())));
  if(ds.rhs())
    _x[1].reset(new SiconosVector(*(ds.rhs())));
  if (ds.jacobianRhsx())
    _jacxRhs.reset(new SimpleMatrix(*(ds.jacobianRhsx())));

  _z.reset(new SiconosVector(*(ds.z())));

  if (ds.xMemory())
    _xMemory.reset(new SiconosMemory(*(ds.xMemory())));
  _stepsInMemory = ds.stepsInMemory();
}

void DynamicalSystem::resetToInitialState()
{
  if(_x0)
    {
      *(_x[0]) = *_x0;
    }
  else
    RuntimeException::selfThrow("DynamicalSystem::resetToInitialState() - initial state _x0 is null");
}


// Setters

void DynamicalSystem::setX0(const SiconosVector& newValue)
{
  // check dimensions ...
  if (newValue.size() != _n)
    RuntimeException::selfThrow("DynamicalSystem::setX0 - inconsistent sizes between x0 input and system dimension.");
  if (_x0)
    *_x0 = newValue;

  else
  {
    _x0.reset(new SiconosVector(newValue));
  }
}

void DynamicalSystem::setX0Ptr(SP::SiconosVector newPtr)
{
  // check dimensions ...
  if (newPtr->size() != _n)
    RuntimeException::selfThrow("DynamicalSystem::setX0Ptr - inconsistent sizes between x0 input and system dimension.");
  _x0 = newPtr;
}

void DynamicalSystem::setX(const SiconosVector& newValue)
{
  // Warning: this only sets the value of x[0]
  // We suppose that both x and (*x)[0] are properly allocated.

  // check dimensions ...
  if (newValue.size() != _n)
    RuntimeException::selfThrow("DynamicalSystem::setX - inconsistent sizes between x input and system dimension.");

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
    RuntimeException::selfThrow("DynamicalSystem::setXPtr - inconsistent sizes between x input and system dimension.");

  _x[0] = newPtr;
}

void DynamicalSystem::setRhs(const SiconosVector& newValue)
{
  // Warning: this only sets the value of x[1]

  // check dimensions ...
  if (newValue.size() != _n)
    RuntimeException::selfThrow("DynamicalSystem::setRhs - inconsistent sizes between rhs input and system dimension.");

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
    RuntimeException::selfThrow("DynamicalSystem::setRhsPtr - inconsistent sizes between rhs input and system dimension.");

  _x[1] = newPtr;
}
void DynamicalSystem::setR(const SiconosVector& newValue)
{
  // check dimensions ...
  if (newValue.size() != _n)
    RuntimeException::selfThrow("DynamicalSystem::setR - inconsistent sizes between input and system dimension.");

  if (_r)
    *_r = newValue;

  else
    _r.reset(new SiconosVector(newValue));
}

void DynamicalSystem::setRPtr(SP::SiconosVector newPtr)
{
  // check dimensions ...
  if (newPtr->size() != _n)
    RuntimeException::selfThrow("DynamicalSystem::setRPtr - inconsistent sizes between input and system dimension.");

  _r = newPtr;

}

void DynamicalSystem::setJacobianRhsx(const SiconosMatrix& newValue)
{
  // check dimensions ...
  if (newValue.size(0) != _n || newValue.size(1) != _n)
    RuntimeException::selfThrow("DynamicalSystem::setJacobianRhsx - inconsistent sizes between and system dimension.");

  if (_jacxRhs)
    *_jacxRhs = newValue;

  else
    _jacxRhs.reset(new SimpleMatrix(newValue));
}

void DynamicalSystem::setJacobianRhsxPtr(SP::SiconosMatrix newPtr)
{
  // check dimensions ...
  if (newPtr->size(0) != _n || newPtr->size(1) != _n)
    RuntimeException::selfThrow("DynamicalSystem::setJacobianRhsxPtr - inconsistent sizes between and system dimension.");

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

