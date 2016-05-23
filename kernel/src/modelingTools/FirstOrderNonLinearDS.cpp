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

// From a minimum set of data
FirstOrderNonLinearDS::FirstOrderNonLinearDS(SP::SiconosVector newX0):
  DynamicalSystem(newX0->size())
{
  DEBUG_PRINT("FirstOrderNonLinearDS::FirstOrderNonLinearDS(SP::SiconosVector newX0)\n");
  DEBUG_EXPR( newX0->display(););


  zeroPlugin();
  // == Initial conditions ==
  _x0 = newX0;

  // == Current state ==
  // x is composed of two blocks of size n, x[0] = \f$ x \f$ and x[1]=\f$ \dot x \f$.
  // x[0] initialized with x0.

  _x[0].reset(new SiconosVector(*_x0));
  _x[1].reset(new SiconosVector(_n));

  //mG
  _workspace[free].reset(new SiconosVector(dimension()));
  _fold.reset(new SiconosVector(dimension()));
  _f.reset(new SiconosVector(_n));
  _b.reset(new SiconosVector(dimension()));
  _jacobianfx.reset(new SimpleMatrix(_n, _n));
  // == r ==

  _r.reset(new SiconosVector(_n));

  checkDynamicalSystem();
}

// From a minimum set of data
FirstOrderNonLinearDS::FirstOrderNonLinearDS(SP::SiconosVector newX0, const std::string& fPlugin, const std::string& jacobianfxPlugin):
  DynamicalSystem(newX0->size())
{
  zeroPlugin();
  // == Initial conditions ==
  _x0 = newX0;

  // == Current state ==
  // x is composed of two blocks of size n, x[0] = \f$ x \f$ and x[1]=\f$ \dot x \f$.
  // x[0] initialized with x0.

  _x[0].reset(new SiconosVector(*_x0));
  _x[1].reset(new SiconosVector(_n));
  _f.reset(new SiconosVector(_n));
  _b.reset(new SiconosVector(dimension()));
  _jacobianfx.reset(new SimpleMatrix(_n, _n));
  _workspace[free].reset(new SiconosVector(dimension()));
  _r.reset(new SiconosVector(dimension()));
  _fold.reset(new SiconosVector(dimension()));

  // == r ==

  _r.reset(new SiconosVector(_n));

  // == f and its jacobian ==
  // Allocation and link with the plug-in
  _pluginf->setComputeFunction(fPlugin);
  _pluginJacxf->setComputeFunction(jacobianfxPlugin);
  //  Plugin::setFunction(&computeJacobianfxPtr, SSLH::getPluginName( jacobianfxPlugin ),SSLH::getPluginFunctionName( jacobianfxPlugin ));
  //  _pluginNameComputeFPtr = fPlugin;
  //  pluginNameComputeJacobianfxPtr = jacobianfxPlugin;

  checkDynamicalSystem();
}

// Copy constructor
FirstOrderNonLinearDS::FirstOrderNonLinearDS(const FirstOrderNonLinearDS & FONLDS): DynamicalSystem(FONLDS)
{
  // Always initialized
  _fold.reset(new SiconosVector(*(FONLDS.fold())));
  _rMemory.reset(new SiconosMemory(*(FONLDS.rMemory())));

  // Not always initialized
  if (FONLDS.getPluginF())
    _pluginf.reset(new PluggedObject(*(FONLDS.getPluginF())));
  if (FONLDS.getPluginJacxf())
    _pluginJacxf.reset(new PluggedObject(*(FONLDS.getPluginJacxf())));
  if (FONLDS.jacobianfx())
    _jacobianfx.reset(new SimpleMatrix(*(FONLDS.jacobianfx())));

  if (FONLDS.M())
    _M.reset(new SimpleMatrix(*(FONLDS.M())));
  if (FONLDS.f())
    _f.reset(new SiconosVector(*(FONLDS.f())));

  if (FONLDS.b())
    _b.reset(new SiconosVector(*(FONLDS.b())));

  if (FONLDS.getPluginM())
    _pluginM.reset(new PluggedObject(*(FONLDS.getPluginM())));

  // data - not always initialized
  if (FONLDS.invM())
    _invM.reset(new SimpleMatrix(*(FONLDS.invM())));
}


void FirstOrderNonLinearDS::zeroPlugin()
{
  // DynamicalSystem::zeroPlugin();
  _pluginf.reset(new PluggedObject());
  _pluginJacxf.reset(new PluggedObject());
  _pluginM.reset(new PluggedObject());
}

bool FirstOrderNonLinearDS::checkDynamicalSystem()
{
  DynamicalSystem::checkDynamicalSystem();
  bool output = DynamicalSystem::checkDynamicalSystem();
  if (!output) std::cout << "FirstOrderNonLinearDS Warning: your dynamical system seems to be uncomplete (check = false)" <<std::endl;
  return output;
}


/*void FirstOrderNonLinearDS::setM(const PMJF& newValue)
{
  assert(newValue.size(0)==n&&"FirstOrderNonLinearDS - setM: inconsistent dimensions with problem size for input matrix M.");
  assert(newValue.size(1)==n&&"FirstOrderNonLinearDS - setM: inconsistent dimensions with problem size for input matrix M.");

  if( ! M )
    M.reset(new PMJF(newValue));
  else
    *M = newValue;
    }*/

void FirstOrderNonLinearDS::setInvM(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != _n || newValue.size(1) != _n)
    RuntimeException::selfThrow("FirstOrderNonLinearDS::setInvM: inconsistent dimensions with problem size for input matrix.");

  if (! _invM)
    _invM.reset(new SimpleMatrix(_n, _n));
  *_invM = newValue;
}

void FirstOrderNonLinearDS::setInvMPtr(SP::SiconosMatrix newPtr)
{
  _invM = newPtr;
}

void FirstOrderNonLinearDS::initRhs(double time)
{
  // compute initial values for f and jacobianfx, initialize right-hand side.
  computeRhs(time); // this will compute, if required, f and M.

  if (! _jacxRhs)  // if not allocated with a set or anything else
  {
    if (_jacobianfx && ! _M)  // if M is not defined, then jacobianfx = jacobianRhsx, no memory allocation for that one.
      _jacxRhs = _jacobianfx;
    else if (_jacobianfx && _M)
      _jacxRhs.reset(new SimpleMatrix(_n, _n));

    // else no allocation, jacobian is equal to 0.
  }
}

void FirstOrderNonLinearDS::updatePlugins(double time)
{
  if (_M)
    computeM(time);

  computef(time);
  computeJacobianfx(time);
}

void FirstOrderNonLinearDS::initialize(double time, unsigned int sizeOfMemory)
{
  DEBUG_PRINT("FirstOrderNonLinearDS::initialize(double time, unsigned int sizeOfMemory)" );
  DEBUG_EXPR( _x0->display(););
  DEBUG_EXPR( _x[0]->display(););


  // reset x to x0.
  *(_x[0]) = *_x0;

  // If z has not been set, we initialize it with a null vector of size 1, since z is required in plug-in functions call.
  if (! _z)
    _z.reset(new SiconosVector(1));

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  updatePlugins(time);
  if (_f)
    *_fold = *_f;

  //   if (simulationType == "EventDriven"){
  //     // Rhs and its jacobian ==> the right is to put in initOSNA of EventDriven
  //     initRhs(time);
  //   }
}

void FirstOrderNonLinearDS::initializeNonSmoothInput(unsigned int level)
{

  /**\warning V.A. _r should be initialized here and not in  the constructor
   * The level should also be used if we need more thatn one _r
   */

  // reset  r to zero.
  _r->zero();
}



// ===== MEMORY MANAGEMENT FUNCTIONS =====

void FirstOrderNonLinearDS::initMemory(unsigned int steps)
{
  DynamicalSystem::initMemory(steps);

  if (steps == 0)
    std::cout << "Warning : FirstOrderNonLinearDS::initMemory with size equal to zero" <<std::endl;
  else
    _rMemory.reset(new SiconosMemory(steps, _n));
}

void FirstOrderNonLinearDS::swapInMemory()
{
  _xMemory->swap(*_x[0]);
  _rMemory->swap(*_r);
  *_fold = *_f;
}

// ===== COMPUTE PLUGINS FUNCTIONS =====

void FirstOrderNonLinearDS::setComputeMFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginM->setComputeFunction(pluginPath, functionName);
}

void FirstOrderNonLinearDS::setComputeMFunction(FPtr1 fct)
{
  _pluginM->setComputeFunction((void *)fct);
}

void FirstOrderNonLinearDS::setComputeFFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginf->setComputeFunction(pluginPath, functionName);
}

void FirstOrderNonLinearDS::setComputeFFunction(FPtr1 fct)
{
  _pluginf->setComputeFunction((void *)fct);
}

void FirstOrderNonLinearDS::setComputeJacobianfxFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginJacxf->setComputeFunction(pluginPath, functionName);
}

void FirstOrderNonLinearDS::setComputeJacobianfxFunction(FPtr1 fct)
{
  _pluginJacxf->setComputeFunction((void *)fct);
}

void FirstOrderNonLinearDS::computeM(double time)
{
  // second argument is useless at the time - Used in derived classes
  if (_pluginM->fPtr)
  {
    ((FNLDSPtrfct)_pluginM->fPtr)(time, _n, &((*(_x[0]))(0)), &(*_M)(0, 0), _z->size(), &(*_z)(0));
  }
}

void FirstOrderNonLinearDS::computef(double time)
{
  if (_pluginf->fPtr)
    ((FNLDSPtrfct)_pluginf->fPtr)(time, _n, _x[0]->getArray(), _f->getArray(), _z->size(), &(*_z)(0));
}

void FirstOrderNonLinearDS::computef(double time, SiconosVector& x2)
{
  if (_pluginf->fPtr)
    ((FNLDSPtrfct)_pluginf->fPtr)(time, _n, &((x2)(0)) , &(*_f)(0), _z->size(), &(*_z)(0));
  // else nothing!
}

void FirstOrderNonLinearDS::computeJacobianfx(double time, bool isDSUp)
{
  // second argument is useless at the time - Used in derived classes
  if (_pluginJacxf->fPtr)
    ((FNLDSPtrfct)_pluginJacxf->fPtr)(time, _n, &((*(_x[0]))(0)), &(*_jacobianfx)(0, 0), _z->size(), &(*_z)(0));
}

void FirstOrderNonLinearDS::computeJacobianfx(double time, const SiconosVector& x2)
{
  // second argument is useless at the time - Used in derived classes
  if (_pluginJacxf->fPtr)
    ((FNLDSPtrfct)_pluginJacxf->fPtr)(time, _n, x2.getArray(), &(*_jacobianfx)(0, 0), _z->size(), _z->getArray());
}

void FirstOrderNonLinearDS::computeRhs(double time, bool isDSUp)
{
  // second argument is useless at the time - Used in derived classes

  // compute rhs = M-1*( f + r ).

  *_x[1] = *_r; // Warning: p update is done in Interactions/Relations

  if (_f)
  {
    computef(time);
    *(_x[1]) += *_f;
  }

  if (_M)
  {
    // allocate invM at the first call of the present function
    if (! _invM)
      _invM.reset(new SimpleMatrix(*_M));
    _invM->PLUForwardBackwardInPlace(*_x[1]);
  }
}

void FirstOrderNonLinearDS::computeJacobianRhsx(double time, bool isDSUp)
{
  // second argument is useless at the time - Used in derived classes

  // compute jacobian of rhs according to x, = M-1(jacobianfx + jacobianX(T.u))
  // At the time, second term is set to zero.
  assert(!_pluginJacxf->fPtr && "FirstOrderNonLinearDS::computeJacobianRhsx: there is no plugin to compute the jacobian of f");

  computeJacobianfx(time);
  // solve M*jacobianXRhS = jacobianfx
  if (_M && _jacobianfx)
  {
    *_jacxRhs = *_jacobianfx;
    // copy _M into _invM for LU-factorisation, at the first call of this function.
    if (! _invM)
      _invM.reset(new SimpleMatrix(*_M));

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

void FirstOrderNonLinearDS::resetAllNonSmoothPart()
{
  _r->zero();
}

void FirstOrderNonLinearDS::resetNonSmoothPart(unsigned int level)
{
  // V.A. 28/05/2012:  for the moment various level are not used for First Order systems
  //assert(0);
  _r->zero();
}
void FirstOrderNonLinearDS::initWorkSpace(VectorOfVectors& workVector, VectorOfMatrices& workMatrices)
{
  workVector.resize(FirstOrderDS::sizeWorkV);
  workVector[FirstOrderDS::residu].reset(new SiconosVector(_n));
  workVector[FirstOrderDS::residuFree].reset(new SiconosVector(_n));
  workVector[FirstOrderDS::xfree].reset(new SiconosVector(_n));
  workVector[FirstOrderDS::xPartialNS].reset(new SiconosVector(_n));
  workVector[FirstOrderDS::deltaxForRelation].reset(new SiconosVector(_n));
  workVector[FirstOrderDS::xBuffer].reset(new SiconosVector(_n));
}

