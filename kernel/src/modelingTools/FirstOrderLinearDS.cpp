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
#include "FirstOrderLinearDS.hpp"
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include "debug.h"
#include <iostream>

typedef void (*computeAfct)(double, unsigned int, unsigned int, double*, unsigned int, double*);

// --- Constructors ---
// From a minimum set of data, A and b connected to a plug-in
FirstOrderLinearDS::FirstOrderLinearDS(SP::SiconosVector newX0, const std::string& APlugin, const std::string& bPlugin):
  FirstOrderNonLinearDS(newX0)
{
  _zeroPlugin();
  setComputeAFunction(SSLH::getPluginName(APlugin), SSLH::getPluginFunctionName(APlugin));
  setComputebFunction(SSLH::getPluginName(bPlugin), SSLH::getPluginFunctionName(bPlugin));
  // dot x = A(t)x + b(t)
}

// From a minimum set of data, A from a given matrix
FirstOrderLinearDS::FirstOrderLinearDS(SP::SiconosVector newX0, SP::SiconosMatrix newA):
  FirstOrderNonLinearDS(newX0)
{
  _zeroPlugin();
  if ((newA->size(0) != _n) || (newA->size(1) != _n))
    RuntimeException::selfThrow("FirstOrderLinearDS - constructor(number,x0,A): inconsistent dimensions with problem size for input matrix A");
  _A = newA;
}

FirstOrderLinearDS::FirstOrderLinearDS(SP::SiconosVector newX0):
  FirstOrderNonLinearDS(newX0)
{
  _zeroPlugin();
}

// From a minimum set of data, A from a given matrix
FirstOrderLinearDS::FirstOrderLinearDS(SP::SiconosVector newX0, SP::SiconosMatrix newA, SP::SiconosVector newB):
  FirstOrderNonLinearDS(newX0)
{
  _zeroPlugin();
    
  if ((newA->size(0) != _n) || (newA->size(1) != _n))
    RuntimeException::selfThrow("FirstOrderLinearDS - constructor(x0,A,b): inconsistent dimensions with problem size for input matrix A");

  if (newB->size() != _n)
    RuntimeException::selfThrow("FirstOrderLinearDS - constructor(x0,A,b): inconsistent dimensions with problem size for input vector b ");

  _A = newA;
  _b = newB;
}

// Copy constructor
FirstOrderLinearDS::FirstOrderLinearDS(const FirstOrderLinearDS & FOLDS): FirstOrderNonLinearDS(FOLDS)
{
  _zeroPlugin();
  if(FOLDS.A())
    _A.reset(new SimpleMatrix(*(FOLDS.A())));
  if(FOLDS.b())
    _b.reset(new SiconosVector(*(FOLDS.b())));
  
  if (Type::value(FOLDS) == Type::FirstOrderLinearDS)
  {
    _pluginA.reset(new PluggedObject(*(FOLDS.getPluginA())));
    _pluginb.reset(new PluggedObject(*(FOLDS.getPluginb())));
  }
}

void FirstOrderLinearDS::initRhs(double time)
{

  DEBUG_PRINT("init Rhs in FirstOrderLinearDS");
  computeRhs(time); // If necessary, this will also compute A and b.
  if (! _jacxRhs)  // if not allocated with a set or anything else
  {
    if (_A && ! _M)  // if M is not defined, then A = jacobianRhsx, no memory allocation for that one.
      _jacxRhs = _A;
    else if (_A && _M)
      _jacxRhs.reset(new SimpleMatrix(_n, _n));
    // else no allocation, jacobian is equal to 0.
  }
  computeJacobianRhsx(time);
}

void FirstOrderLinearDS::updatePlugins(double time)
{
  if(_M)
    computeM(time);
  if(_A)
    computeA(time);
  if (_b)
    computeb(time);
}

void FirstOrderLinearDS::setComputeAFunction(const std::string& pluginPath, const std::string& functionName)
{
  if(!_A)
    _A.reset(new SimpleMatrix(_n, _n));
  _pluginA->setComputeFunction(pluginPath, functionName);
}

void FirstOrderLinearDS::setComputeAFunction(LDSPtrFunction fct)
{
  if(!_A)
    _A.reset(new SimpleMatrix(_n, _n));
  _pluginA->setComputeFunction((void*)fct);
}

void FirstOrderLinearDS::setComputebFunction(const std::string& pluginPath, const std::string& functionName)
{
  if (!_b)
    _b.reset(new SiconosVector(_n));
  _pluginb->setComputeFunction(pluginPath, functionName);
}

void FirstOrderLinearDS::setComputebFunction(LDSPtrFunction fct)
{
  if(!_b)
    _b.reset(new SiconosVector(_n));
  _pluginb->setComputeFunction((void*)fct);
}

void FirstOrderLinearDS::computeA(double time)
{
  if (_A && _pluginA->fPtr)
  {
    ((computeAfct)_pluginA->fPtr)(time, _n, _n, &(*_A)(0, 0), _z->size(), &(*_z)(0));
  }
}

void FirstOrderLinearDS::computeb(double time)
{
  if (_b && _pluginb->fPtr)
    ((LDSPtrFunction)_pluginb->fPtr)(time, _n, &(*_b)(0), _z->size(), &(*_z)(0));
}
/*This function is called only by LsodarOSI and eventDriven*/
void FirstOrderLinearDS::computeRhs(double time, bool isDSup)
{
  // second argument is useless at the time - Used in derived classes
  // compute A=jacobianfx

  *_x[1] = * _r;

  if (_A)
  {
    computeA(time);
    prod(*_A, *_x[0], *_x[1], false);
  }

  // compute and add b if required
  if (_b)
  {
    computeb(time);
    *_x[1] += *_b;
  }

  if (_M)
  {
    computeM(time);
    // allocate invM at the first call of the present function
    if (! _invM)
      _invM.reset(new SimpleMatrix(*_M));

    _invM->PLUForwardBackwardInPlace(*_x[1]);
  }
}

void FirstOrderLinearDS::computeJacobianRhsx(double time, bool isDSup)
{
  if(_A)
    {
      computeA(time);
      if (_M)
	{
	  computeM(time);
	  *_jacxRhs = *_A;
	  if (! _invM)
	    _invM.reset(new SimpleMatrix(*_M));
	  else if(_pluginM->fPtr) // if M is plugged, invM must be updated
	    *_invM = *_M;
	  // solve MjacobianRhsx = A
	  _invM->PLUForwardBackwardInPlace(*_jacxRhs);
	}
    }
  // else 0
}

void FirstOrderLinearDS::display() const
{
  std::cout << "=== Linear system display, " << _number << std::endl;
  std::cout << "=============================" << std::endl;
}

void FirstOrderLinearDS::setA(const SiconosMatrix& newA)
{
  if (_A)
    *_A = newA;
  else
    _A.reset(new SimpleMatrix(newA));
}

void FirstOrderLinearDS::_zeroPlugin()
{
  _pluginA.reset(new PluggedObject());
  _pluginb.reset(new PluggedObject());
}
