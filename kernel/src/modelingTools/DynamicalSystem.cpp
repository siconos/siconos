/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include "SiconosVector.hpp"
#include "SiconosException.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosMemory.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include "siconos_debug.h"

#include <iostream>

  
size_t siconos::modeling::DynamicalSystem::__count = 0;

// void siconos::modeling::siconos::modeling::DynamicalSystem::_init()
// {
//   DEBUG_PRINT("internal _init from DynamicalSystem\n");

//   // No memory allocation, only resize for containers.
//   // Everything should be done in derived class init for required operators
//   // and variables and in 'set'-like methods for optional
//   // components.

//   _z = std::make_shared<siconos::algebra::SiconosVector>(1);
// }

// ===== CONSTRUCTORS =====

// From a minimum set of data
siconos::modeling::DynamicalSystem::DynamicalSystem(unsigned int dimension) :  _n(dimension)
{
  _z = std::make_shared<siconos::algebra::SiconosVector>(1);
}

siconos::modeling::DynamicalSystem::DynamicalSystem(const DynamicalSystem &ds)
    :  _n(ds.n()), _stepsInMemory(ds.stepsInMemory())
{
  // The following data should always be initialize
  if (ds.x0())
    _x0 = std::make_shared<siconos::algebra::SiconosVector>(*(ds.x0()));
  if (ds.r())
    _r = std::make_shared<siconos::algebra::SiconosVector>(*(ds.r()));
  _x.resize(2);
  if (ds.x())
    _x[0] = std::make_shared<siconos::algebra::SiconosVector>(*(ds.x()));
  if (ds.rhs())
    _x[1] = std::make_shared<siconos::algebra::SiconosVector>(*(ds.rhs()));
  if (ds.jacobianRhsx())
    _jacxRhs = std::make_shared<siconos::algebra::SimpleMatrix>(*(ds.jacobianRhsx()));

  _z = std::make_shared<siconos::algebra::SiconosVector>(*(ds.z()));

  _xMemory = ds.xMemory();
  _stepsInMemory = ds.stepsInMemory();
}

void siconos::modeling::DynamicalSystem::resetToInitialState()
{
  if (_x0) {
    *(_x[0]) = *_x0;
  }
  else
    THROW_EXCEPTION("siconos::modeling::DynamicalSystem::resetToInitialState() - initial state _x0 is null");
}

// Setters

void siconos::modeling::DynamicalSystem::setX0(const siconos::algebra::SiconosVector &newValue)
{
  // check dimensions ...
  if (newValue.size() != _n)
    THROW_EXCEPTION(
        "siconos::modeling::DynamicalSystem::setX0 - inconsistent sizes between x0 input and system dimension.");
  if (_x0)
    *_x0 = newValue;

  else {
    _x0 = std::make_shared<siconos::algebra::SiconosVector>(newValue);
  }
}

void siconos::modeling::DynamicalSystem::setX0Ptr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr)
{
  // check dimensions ...
  if (newPtr->size() != _n)
    THROW_EXCEPTION("siconos::modeling::DynamicalSystem::setX0Ptr - inconsistent sizes between x0 input and "
                    "system dimension.");
  _x0 = newPtr;
}

void siconos::modeling::DynamicalSystem::setX(const siconos::algebra::SiconosVector &newValue)
{
  // Warning: this only sets the value of x[0]
  // We suppose that both x and (*x)[0] are properly allocated.

  // check dimensions ...
  if (newValue.size() != _n)
    THROW_EXCEPTION(
        "siconos::modeling::DynamicalSystem::setX - inconsistent sizes between x input and system dimension.");

  if (!_x[0])
    _x[0] = std::make_shared<siconos::algebra::SiconosVector>(newValue);
  else
    *(_x[0]) = newValue;
}

void siconos::modeling::DynamicalSystem::setXPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr)
{
  // Warning: this only sets the pointer x[0]

  // check dimensions ...
  if (newPtr->size() != _n)
    THROW_EXCEPTION(
        "siconos::modeling::DynamicalSystem::setXPtr - inconsistent sizes between x input and system dimension.");

  _x[0] = newPtr;
}

void siconos::modeling::DynamicalSystem::setRhs(const siconos::algebra::SiconosVector &newValue)
{
  // Warning: this only sets the value of x[1]

  // check dimensions ...
  if (newValue.size() != _n)
    THROW_EXCEPTION("siconos::modeling::DynamicalSystem::setRhs - inconsistent sizes between rhs input and "
                    "system dimension.");

  if (!_x[1])
    _x[1] = std::make_shared<siconos::algebra::SiconosVector>(newValue);
  else
    *(_x[1]) = newValue;
}

void siconos::modeling::DynamicalSystem::setRhsPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr)
{
  // Warning: this only sets the pointer (*x)[1]

  // check dimensions ...
  if (newPtr->size() != _n)
    THROW_EXCEPTION("siconos::modeling::DynamicalSystem::setRhsPtr - inconsistent sizes between rhs input and "
                    "system dimension.");

  _x[1] = newPtr;
}
void siconos::modeling::DynamicalSystem::setR(const siconos::algebra::SiconosVector &newValue)
{
  // check dimensions ...
  if (newValue.size() != _n)
    THROW_EXCEPTION(
        "siconos::modeling::DynamicalSystem::setR - inconsistent sizes between input and system dimension.");

  if (_r)
    *_r = newValue;

  else
    _r = std::make_shared<siconos::algebra::SiconosVector>(newValue);
}

void siconos::modeling::DynamicalSystem::setRPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr)
{
  // check dimensions ...
  if (newPtr->size() != _n)
    THROW_EXCEPTION(
        "siconos::modeling::DynamicalSystem::setRPtr - inconsistent sizes between input and system dimension.");

  _r = newPtr;
}

void siconos::modeling::DynamicalSystem::setJacobianRhsx(const siconos::algebra::SiconosMatrix &newValue)
{
  // check dimensions ...
  if (newValue.size(0) != _n || newValue.size(1) != _n)
    THROW_EXCEPTION(
        "siconos::modeling::DynamicalSystem::setJacobianRhsx - inconsistent sizes between and system dimension.");

  if (_jacxRhs)
    *_jacxRhs = newValue;

  else
    _jacxRhs = std::make_shared<siconos::algebra::SimpleMatrix>(newValue);
}

void siconos::modeling::DynamicalSystem::setJacobianRhsxPtr(std::shared_ptr<siconos::algebra::SiconosMatrix> newPtr)
{
  // check dimensions ...
  if (newPtr->size(0) != _n || newPtr->size(1) != _n)
    THROW_EXCEPTION("siconos::modeling::DynamicalSystem::setJacobianRhsxPtr - inconsistent sizes between and "
                    "system dimension.");

  _jacxRhs = newPtr;
}

void siconos::modeling::DynamicalSystem::setz(const siconos::algebra::SiconosVector &newValue)
{
  if (_z) {
    if (newValue.size() != _z->size())
      THROW_EXCEPTION("siconos::modeling::DynamicalSystem::setz - inconsistent sizes between input and existing "
                      "z - To change z size use setzPtr.");
    *_z = newValue;
  }
  else {
    _z = std::make_shared<siconos::algebra::SiconosVector>(newValue);
  }
}

void siconos::modeling::DynamicalSystem::setzPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) { _z = newPtr; }

void siconos::modeling::DynamicalSystem::update(double time)
{
  computeRhs(time);
  computeJacobianRhsx(time);
}

// ===== MEMORY MANAGEMENT FUNCTIONS =====

void siconos::modeling::DynamicalSystem::initMemory(unsigned int steps)
{
  DEBUG_BEGIN("void siconos::modeling::DynamicalSystem::initMemory(unsigned int steps)\n");
  if (steps > _xMemory.size()) {
    if (steps == 0)
      std::cout << "Warning : initMemory with size equal to zero" << std::endl;
    else {
      _stepsInMemory = steps;
      _xMemory.setMemorySize(steps, _n);
    }
  }
  DEBUG_EXPR(_xMemory.display(););

  DEBUG_END("void siconos::modeling::DynamicalSystem::initMemory(unsigned int steps)\n");
}

