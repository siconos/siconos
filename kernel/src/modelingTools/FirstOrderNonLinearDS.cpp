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

#include "FirstOrderNonLinearDS.hpp"
#include "PluginTypes.hpp"

// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include "debug.h"
#include <iostream>

// ===== CONSTRUCTORS =====

void FirstOrderNonLinearDS::_init(SP::SiconosVector initial_state)
{
  DEBUG_PRINT("internal _init from FirstOrderNonLinearDS\n");

  // Memory allocation only for required parts of the DS:
  // state (initial and current). All other operators are optional and
  // allocated with 'set'-like methods.
  assert(_n > 0 && "dynamical system dimension should be greater than 0.");
  // Set initial conditions
  _x0 = initial_state;

  // == Current state ==
  // x is composed of two blocks of size n, x[0] = \f$ x \f$ and x[1]=\f$ \dot x \f$.
  // x[0] initialized with x0.
  // _x.resize(2); done in base class constructor.
  _x[0].reset(new SiconosVector(*_x0));
  _x[1].reset(new SiconosVector(_n));
  _r.reset(new SiconosVector(_n)); // FP: move this to initializeNonSmoothInput?
  _zeroPlugin();
}

// From a minimum set of data
FirstOrderNonLinearDS::FirstOrderNonLinearDS(SP::SiconosVector initial_state):
  DynamicalSystem(initial_state->size())
{
  _init(initial_state);
  // dot x = r
}

// From a minimum set of data
FirstOrderNonLinearDS::FirstOrderNonLinearDS(SP::SiconosVector initial_state, const std::string& fPlugin, const std::string& jacobianfxPlugin):
  DynamicalSystem(initial_state->size())
{
  _init(initial_state);
  // == f and its jacobian ==
  // Allocation and link with the plug-in
  setComputeFFunction(SSLH::getPluginName(fPlugin), SSLH::getPluginFunctionName(fPlugin));
  setComputeJacobianfxFunction(SSLH::getPluginName(jacobianfxPlugin), SSLH::getPluginFunctionName(jacobianfxPlugin));
  // dot x  = f(x, z , t) + r
}

// Copy constructor
FirstOrderNonLinearDS::FirstOrderNonLinearDS(const FirstOrderNonLinearDS & FONLDS): DynamicalSystem(FONLDS)
{
  _zeroPlugin();
  
  if (FONLDS.M())
    _M.reset(new SimpleMatrix(*(FONLDS.M())));
  if (FONLDS.f())
    _f.reset(new SiconosVector(*(FONLDS.f())));
  if (FONLDS.jacobianfx())
    _jacobianfx.reset(new SimpleMatrix(*(FONLDS.jacobianfx())));
  if (FONLDS.b())
    _b.reset(new SiconosVector(*(FONLDS.b())));
  if (FONLDS.getPluginF())
    _pluginf.reset(new PluggedObject(*(FONLDS.getPluginF())));
  if (FONLDS.getPluginJacxf())
    _pluginJacxf.reset(new PluggedObject(*(FONLDS.getPluginJacxf())));
  if (FONLDS.getPluginM())
    _pluginM.reset(new PluggedObject(*(FONLDS.getPluginM())));
  if (FONLDS.invM())
    _invM.reset(new SimpleMatrix(*(FONLDS.invM())));

  // Memory stuff to me moved to graph/osi
  if(FONLDS.fold())
    _fold.reset(new SiconosVector(*(FONLDS.fold())));
  if(FONLDS.rMemory())
    _rMemory.reset(new SiconosMemory(*(FONLDS.rMemory())));
}


void FirstOrderNonLinearDS::_zeroPlugin()
{
  _pluginf.reset(new PluggedObject());
  _pluginJacxf.reset(new PluggedObject());
  _pluginM.reset(new PluggedObject());
}

void FirstOrderNonLinearDS::initRhs(double time)
{
  computeRhs(time);


  // !! jacxRhs must always be allocated (we must check this?)!!
  if (! _jacxRhs)  // if not allocated with a set or anything else
  {
    if (_jacobianfx && ! _M)  // if M is not defined, then jacobianfx = jacobianRhsx, no memory allocation for that one.
      _jacxRhs = _jacobianfx;
    else//  if (_jacobianfx && _M) or if(!jacobianRhsx)
      _jacxRhs.reset(new SimpleMatrix(_n, _n));

    // else no allocation, jacobian is equal to 0.
  }
  computeJacobianRhsx(time);
}

void FirstOrderNonLinearDS::updatePlugins(double time)
{
  if (_M)
    computeM(time);
  if(_f)
    {
      computef(time, _x[0]);
      computeJacobianfx(time, _x[0]);
    }
}

void FirstOrderNonLinearDS::initializeNonSmoothInput(unsigned int level)
{
  /**\warning V.A. _r should be initialized here and not in  the constructor
   * The level should also be used if we need more that one _r
   */
  if (!_r)
    _r.reset(new SiconosVector(_n));
}



// ===== MEMORY MANAGEMENT FUNCTIONS =====

void FirstOrderNonLinearDS::initMemory(unsigned int steps)
{
  DynamicalSystem::initMemory(steps);

  if(_f && !_fold)
    _fold.reset(new SiconosVector(_n));

  if (steps == 0)
    std::cout << "Warning : FirstOrderNonLinearDS::initMemory with size equal to zero" <<std::endl;
  else
    _rMemory.reset(new SiconosMemory(steps, _n));
}

void FirstOrderNonLinearDS::swapInMemory()
{
  _xMemory->swap(*_x[0]);
  if(_rMemory && _r)
    _rMemory->swap(*_r);
  if(_f && _fold)
    *_fold = *_f;
}

// ===== COMPUTE PLUGINS FUNCTIONS =====

void FirstOrderNonLinearDS::setComputeMFunction(const std::string& pluginPath, const std::string& functionName)
{
  if(!_M)
    _M.reset(new SimpleMatrix(_n, _n));

  _pluginM->setComputeFunction(pluginPath, functionName);
}

void FirstOrderNonLinearDS::setComputeMFunction(FPtr1 fct)
{
  if(!_M)
    _M.reset(new SimpleMatrix(_n, _n));

  _pluginM->setComputeFunction((void *)fct);
}

void FirstOrderNonLinearDS::setComputeFFunction(const std::string& pluginPath, const std::string& functionName)
{
  if(!_f)
    _f.reset(new SiconosVector(_n));
  
  _pluginf->setComputeFunction(pluginPath, functionName);
}

void FirstOrderNonLinearDS::setComputeFFunction(FPtr1 fct)
{
  if(!_f)
    _f.reset(new SiconosVector(_n));
  _pluginf->setComputeFunction((void *)fct);
}

void FirstOrderNonLinearDS::setComputeJacobianfxFunction(const std::string& pluginPath, const std::string& functionName)
{
  if(!_jacobianfx)
    _jacobianfx.reset(new SimpleMatrix(_n, _n));
  _pluginJacxf->setComputeFunction(pluginPath, functionName);
}

void FirstOrderNonLinearDS::setComputeJacobianfxFunction(FPtr1 fct)
{
  if(!_jacobianfx)
    _jacobianfx.reset(new SimpleMatrix(_n, _n));
  _pluginJacxf->setComputeFunction((void *)fct);
}

void FirstOrderNonLinearDS::computeM(double time)
{
  if (_pluginM->fPtr && _M) 
  {
    ((FNLDSPtrfct)_pluginM->fPtr)(time, _n, &((*(_x[0]))(0)), &(*_M)(0, 0), _z->size(), &(*_z)(0));
  }
}

void FirstOrderNonLinearDS::computef(double time, SP::SiconosVector state)
{
  if (_f && _pluginf->fPtr)
    ((FNLDSPtrfct)_pluginf->fPtr)(time, _n, &(*state)(0) , &(*_f)(0), _z->size(), &(*_z)(0));
}

void FirstOrderNonLinearDS::computeJacobianfx(double time, SP::SiconosVector state)
{
  if (_jacobianfx && _pluginJacxf->fPtr)
    ((FNLDSPtrfct)_pluginJacxf->fPtr)(time, _n, state->getArray(), &(*_jacobianfx)(0, 0), _z->size(), _z->getArray());
}

void FirstOrderNonLinearDS::computeRhs(double time, bool isDSUp)
{
  // second argument is useless at the time - Used in derived classes

  // compute rhs = M-1*( f + r ).

  *_x[1] = *_r; // Warning: p update is done in Interactions/Relations

  if (_f)
  {
    computef(time, _x[0]);
    *(_x[1]) += *_f;
  }

  if (_M)
  {
    computeM(time);
    // allocate invM at the first call of the present function
    if (! _invM)
      _invM.reset(new SimpleMatrix(*_M));
    else if(_pluginM->fPtr) // if M is plugged, invM must be updated
      *_invM = *_M;
    
    _invM->PLUForwardBackwardInPlace(*_x[1]);
  }
}

void FirstOrderNonLinearDS::computeJacobianRhsx(double time, bool isDSUp)
{
  // second argument is useless at the time - Used in derived classes

  // compute jacobian of rhs according to x, = M-1(jacobianfx + jacobianX(T.u))
  // At the time, second term is set to zero.
  //assert(!_pluginJacxf->fPtr && "FirstOrderNonLinearDS::computeJacobianRhsx: there is no plugin to compute the jacobian of f");

  computeJacobianfx(time, _x[0]);
  // solve M*jacobianXRhS = jacobianfx
  if (_M && _jacobianfx)
  {
    *_jacxRhs = *_jacobianfx;
    // copy _M into _invM for LU-factorisation, at the first call of this function.

    computeM(time);
    
    if (! _invM)
      _invM.reset(new SimpleMatrix(*_M));
    else if(_pluginM->fPtr) // if M is plugged, invM must be updated
      *_invM = *_M;
    _invM->PLUForwardBackwardInPlace(*_jacxRhs);
  }
  // else jacobianRhsx = jacobianfx, pointers equality set in initRhs

}

// ===== MISCELLANEOUS ====

void FirstOrderNonLinearDS::display() const
{
  std::cout << " =====> First Order Non Linear DS (number: " << _number << ")." <<std::endl;
  std::cout << "- n (size) : " << _n <<std::endl;
  std::cout << "- x " <<std::endl;
  if (_x[0]) _x[0]->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- x0 " <<std::endl;
  if (_x0) _x0->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- M: " <<std::endl;
  if (_M) _M->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << " ============================================" <<std::endl;
}

void FirstOrderNonLinearDS::resetAllNonSmoothParts()
{
  _r->zero();
}

void FirstOrderNonLinearDS::resetNonSmoothPart(unsigned int level)
{
  // V.A. 28/05/2012:  for the moment various level are not used for First Order systems
  //assert(0);
  _r->zero();
}

void FirstOrderNonLinearDS::setb(const SiconosVector& b)
{
  if (_b)
    *_b = b;
  else
    _b.reset(new SiconosVector(b));
}
