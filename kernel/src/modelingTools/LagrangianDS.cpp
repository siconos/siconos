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
#include "LagrangianDS.hpp"
#include "BlockVector.hpp"
#include "BlockMatrix.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"
#include <iostream>

void LagrangianDS::_init(SP::SiconosVector position, SP::SiconosVector velocity)
{
  assert(_ndof > 0 && "lagrangian dynamical system dimension should be greater than 0.");

  // Set initial conditions
  _q0 = position;
  _velocity0 = velocity;

  // -- Memory allocation for vector and matrix members --
  _q.resize(3);
  _q[0].reset(new SiconosVector(*_q0));
  _q[1].reset(new SiconosVector(*_velocity0));

  /** \todo lazy Memory allocation */
  _p.resize(3);
  _p[1].reset(new SiconosVector(_ndof));

  _zeroPlugin();
}


// Build from initial state only
LagrangianDS::LagrangianDS(SP::SiconosVector q0, SP::SiconosVector v0):
  DynamicalSystem(2 * q0->size()), _ndof(v0->size()),
  _hasConstantMass(true), _hasConstantFExt(true)
{
  // Initial conditions
  _init(q0, v0);
}

// From initial state and constant mass matrix, \f$ M\ddot q = p \f$
LagrangianDS::LagrangianDS(SP::SiconosVector q0, SP::SiconosVector v0, SP::SiconosMatrix newMass):
  DynamicalSystem(2 * q0->size()), _ndof(v0->size()),
  _hasConstantMass(true), _hasConstantFExt(true)

{
  _init(q0, v0);
  // Mass matrix
  _mass = newMass;
}

// From a set of data - Mass loaded from a plugin
// This constructor leads to the minimum Lagrangian System form: \f$ M(q)\ddot q = p \f$
LagrangianDS::LagrangianDS(SP::SiconosVector q0, SP::SiconosVector v0, const std::string& massName):
  DynamicalSystem(), _ndof(q0->size()),
  _hasConstantMass(false), _hasConstantFExt(true)
{
  _init(q0, v0);
  // Mass
  _mass.reset(new SimpleMatrix(_ndof, _ndof));
  setComputeMassFunction(SSLH::getPluginName(massName), SSLH::getPluginFunctionName(massName));
}

void LagrangianDS::_zeroPlugin()
{
  _pluginMass.reset(new PluggedObject());
  _pluginFInt.reset(new PluggedObject());
  _pluginFExt.reset(new PluggedObject());
  _pluginFGyr.reset(new PluggedObject());
  _pluginJacqFInt.reset(new PluggedObject());
  _pluginJacqDotFInt.reset(new PluggedObject());
  _pluginJacqFGyr.reset(new PluggedObject());
  _pluginJacqDotFGyr.reset(new PluggedObject());
}

void LagrangianDS::initializeNonSmoothInput(unsigned int level)
{
  if(!_p[level])
    _p[level].reset(new SiconosVector(_ndof));
}

void LagrangianDS::resetToInitialState()
{
  if(_q0)
  {
    *(_q[0]) = *_q0;
  }
  else
    RuntimeException::selfThrow("LagrangianDS::resetToInitialState - initial position _q0 is null");
  if(_velocity0)
  {
    *(_q[1]) = *_velocity0;
  }
  else
    RuntimeException::selfThrow("LagrangianDS::resetToInitialState - initial velocity _velocity0 is null");
}

void LagrangianDS::init_generalized_coordinates(unsigned int level)
{
  assert(level>1);
  if(!_q[level])
    _q[level].reset(new SiconosVector(_ndof));
}


void LagrangianDS::init_inverse_mass()
{
  if(_mass && !_inverseMass)
  {
    computeMass();
    _inverseMass.reset(new SimpleMatrix(*_mass));
  }
}

void LagrangianDS::update_inverse_mass()
{
  if(_mass && _inverseMass && !_hasConstantMass)
  {
    computeMass();
    *_inverseMass = *_mass;
  }
}

void LagrangianDS::init_forces()
{
  // Allocate memory for forces and its jacobians.
  // Needed only for integrators with first-order formulation.
  if(_fInt || _fExt || _fGyr)
  {
    if(!_forces)
      _forces.reset(new SiconosVector(_ndof));
  }

  if(_fInt || _fGyr)
  {
    if(!_jacobianqForces)
      _jacobianqForces.reset(new SimpleMatrix(_ndof, _ndof));
    if(!_jacobianqDotForces)
      _jacobianqDotForces.reset(new SimpleMatrix(_ndof, _ndof));
  }
}

void LagrangianDS::initRhs(double time)
{
  // dim
  _n = 2 * _ndof;

  // All links between DS and LagrangianDS class members are pointer links, which means
  // that no useless memory is allocated when connection is established.
  // One exception: zero and identity matrices, used to filled in M and jacobianfx.

  // Initial conditions and state

  // WARNING : this function is supposed to be called
  // by the OneStepIntegrator, and maybe several times for the same DS
  // if the system is involved in more than one interaction. So, we must check
  // if p2 and q2 already exist to be sure that DSlink won't be lost.

  _x0.reset(new SiconosVector(*_q0, *_velocity0));

  _x[0].reset(new SiconosVector(*_q[0], *_q[1]));

  if(!_q[2])
    _q[2].reset(new SiconosVector(_ndof));

  _x[1].reset(new SiconosVector(*_q[1], *_q[2]));

  // Everything concerning rhs and its jacobian is handled in initRhs and computeXXX related functions.
  _rhsMatrices.resize(numberOfRhsMatrices);

  if(!_p[2])
    _p[2].reset(new SiconosVector(_ndof));

  init_forces();
  init_inverse_mass();

  computeRhs(time);

  bool flag1 = false, flag2 = false;
  if(_jacobianqForces)
  {
    // Solve MjacobianX(1,0) = jacobianFL[0]
    computeJacobianqForces(time);

    _rhsMatrices[jacobianXBloc10].reset(new SimpleMatrix(*_jacobianqForces));
    _inverseMass->PLUForwardBackwardInPlace(*_rhsMatrices[jacobianXBloc10]);
    flag1 = true;
  }

  if(_jacobianqDotForces)
  {
    // Solve MjacobianX(1,1) = jacobianFL[1]
    computeJacobianqDotForces(time);
    _rhsMatrices[jacobianXBloc11].reset(new SimpleMatrix(*_jacobianqDotForces));
    _inverseMass->PLUForwardBackwardInPlace(*_rhsMatrices[jacobianXBloc11]);
    flag2 = true;
  }

  if(!_rhsMatrices[zeroMatrix])
    _rhsMatrices[zeroMatrix].reset(new SimpleMatrix(_ndof, _ndof, Siconos::ZERO));
  if(!_rhsMatrices[idMatrix])
    _rhsMatrices[idMatrix].reset(new SimpleMatrix(_ndof, _ndof, Siconos::IDENTITY));

  if(flag1 && flag2)
    _jacxRhs.reset(new BlockMatrix(_rhsMatrices[zeroMatrix], _rhsMatrices[idMatrix],
                                   _rhsMatrices[jacobianXBloc10], _rhsMatrices[jacobianXBloc11]));
  else if(flag1)  // flag2 = false
    _jacxRhs.reset(new BlockMatrix(_rhsMatrices[zeroMatrix], _rhsMatrices[idMatrix],
                                   _rhsMatrices[jacobianXBloc10], _rhsMatrices[zeroMatrix]));
  else if(flag2)  // flag1 = false
    _jacxRhs.reset(new BlockMatrix(_rhsMatrices[zeroMatrix], _rhsMatrices[idMatrix],
                                   _rhsMatrices[zeroMatrix], _rhsMatrices[jacobianXBloc11]));
  else
    _jacxRhs.reset(new BlockMatrix(_rhsMatrices[zeroMatrix], _rhsMatrices[idMatrix],
                                   _rhsMatrices[zeroMatrix], _rhsMatrices[zeroMatrix]));
}

// --- GETTERS/SETTERS ---

void LagrangianDS::setQ(const SiconosVector& newValue)
{
  if(newValue.size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ: inconsistent input vector size ");

  if(! _q[0])
    _q[0].reset(new SiconosVector(newValue));
  else
    *_q[0] = newValue;
}

void LagrangianDS::setQPtr(SP::SiconosVector newPtr)
{
  if(newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setQPtr: inconsistent input vector size ");
  _q[0] = newPtr;

}

void LagrangianDS::setQ0(const SiconosVector& newValue)
{
  if(newValue.size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ0: inconsistent input vector size ");

  if(! _q0)
    _q0.reset(new SiconosVector(newValue));
  else
    *_q0 = newValue;
}

void LagrangianDS::setQ0Ptr(SP::SiconosVector newPtr)
{
  if(newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ0Ptr: inconsistent input vector size ");
  _q0 = newPtr;
}

void LagrangianDS::setVelocity0(const SiconosVector& newValue)
{
  if(newValue.size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity0: inconsistent input vector size ");

  if(! _velocity0)
    _velocity0.reset(new SiconosVector(newValue));
  else
    *_velocity0 = newValue;
}

void LagrangianDS::setVelocity(const SiconosVector& newValue)
{
  if(newValue.size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity: inconsistent input vector size ");

  if(! _q[1])
    _q[1].reset(new SiconosVector(newValue));
  else
    *_q[1] = newValue;
}

void LagrangianDS::setVelocityPtr(SP::SiconosVector newPtr)
{
  if(newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocityPtr: inconsistent input vector size ");
  _q[1] = newPtr;
}


void LagrangianDS::setVelocity0Ptr(SP::SiconosVector newPtr)
{
  if(newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity0Ptr: inconsistent input vector size ");
  _velocity0 = newPtr;
}

void LagrangianDS::setP(const SiconosVector& newValue, unsigned int level)
{
  if(newValue.size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setP: inconsistent input vector size ");

  if(! _p[level])
    _p[level].reset(new SiconosVector(newValue));
  else
    *(_p[level]) = newValue;
}

void LagrangianDS::setPPtr(SP::SiconosVector newPtr, unsigned int level)
{

  if(newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setPPtr: inconsistent input vector size ");
  _p[level] = newPtr;
}

void LagrangianDS::computeMass()
{
  DEBUG_PRINT("LagrangianDS::computeMass()\n");
  DEBUG_EXPR(_q[0]->display());
  if(!_hasConstantMass)
  {
    if(_mass && _pluginMass->fPtr)
    {
      ((FPtr7)_pluginMass->fPtr)(_ndof, &(*_q[0])(0), &(*_mass)(0, 0), _z->size(), &(*_z)(0));
      _mass->resetLU();
    }
  }
  DEBUG_EXPR(_mass->display());
}

void LagrangianDS::computeMass(SP::SiconosVector position)
{
  if(_mass && !_hasConstantMass && _pluginMass->fPtr)
  {
    ((FPtr7)_pluginMass->fPtr)(_ndof, &(*position)(0), &(*_mass)(0, 0), _z->size(), &(*_z)(0));
    _mass->resetLU();
  }
}

void LagrangianDS::computeFInt(double time)
{
  if(_fInt && _pluginFInt->fPtr)
    ((FPtr6)_pluginFInt->fPtr)(time, _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_fInt)(0), _z->size(), &(*_z)(0));
}
void LagrangianDS::computeFInt(double time, SP::SiconosVector position, SP::SiconosVector velocity)
{
  if(_fInt && _pluginFInt->fPtr)
    ((FPtr6)_pluginFInt->fPtr)(time, _ndof, &(*position)(0), &(*velocity)(0), &(*_fInt)(0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeFExt(double time)
{
  if(!_hasConstantFExt)
  {
    if(_fExt && _pluginFExt->fPtr)
      ((VectorFunctionOfTime)_pluginFExt->fPtr)(time, _ndof, &(*_fExt)(0), _z->size(), &(*_z)(0));
  }

}
void LagrangianDS::computeFGyr()
{
  if(_fGyr && _pluginFGyr->fPtr)
    ((FPtr5)_pluginFGyr->fPtr)(_ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_fGyr)(0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeFGyr(SP::SiconosVector position, SP::SiconosVector velocity)
{
  if(_fGyr && _pluginFGyr->fPtr)
    ((FPtr5)_pluginFGyr->fPtr)(_ndof, &(*position)(0), &(*velocity)(0), &(*_fGyr)(0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeJacobianFIntq(double time)
{
  if(_jacobianFIntq&& _pluginJacqFInt->fPtr)
    ((FPtr6)_pluginJacqFInt->fPtr)(time, _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianFIntq)(0, 0), _z->size(), &(*_z)(0));
}
void LagrangianDS::computeJacobianFIntqDot(double time)
{
  if(_jacobianFIntqDot && _pluginJacqDotFInt->fPtr)
    ((FPtr6)_pluginJacqDotFInt->fPtr)(time, _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianFIntqDot)(0, 0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeJacobianFIntq(double time, SP::SiconosVector position, SP::SiconosVector velocity)
{
  if(_jacobianFIntq && _pluginJacqFInt->fPtr)
    ((FPtr6)_pluginJacqFInt->fPtr)(time, _ndof, &(*position)(0), &(*velocity)(0), &(*_jacobianFIntq)(0, 0), _z->size(), &(*_z)(0));
}
void LagrangianDS::computeJacobianFIntqDot(double time, SP::SiconosVector position, SP::SiconosVector velocity)
{
  if(_jacobianFIntqDot && _pluginJacqDotFInt->fPtr)
    ((FPtr6)_pluginJacqDotFInt->fPtr)(time, _ndof, &(*position)(0), &(*velocity)(0), &(*_jacobianFIntqDot)(0, 0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeJacobianFGyrq()
{
  if(_pluginJacqFGyr->fPtr)
    ((FPtr5)_pluginJacqFGyr->fPtr)(_ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianFGyrq)(0, 0), _z->size(), &(*_z)(0));
}
void LagrangianDS::computeJacobianFGyrqDot()
{
  if(_jacobianFGyrqDot && _pluginJacqDotFGyr->fPtr)
    ((FPtr5)_pluginJacqDotFGyr->fPtr)(_ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianFGyrqDot)(0, 0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeJacobianFGyrq(SP::SiconosVector position, SP::SiconosVector velocity)
{
  if(_jacobianFGyrq && _pluginJacqFGyr->fPtr)
    ((FPtr5)_pluginJacqFGyr->fPtr)(_ndof, &(*position)(0), &(*velocity)(0), &(*_jacobianFGyrq)(0, 0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeJacobianFGyrqDot(SP::SiconosVector position, SP::SiconosVector velocity)
{
  if(_jacobianFGyrqDot && _pluginJacqDotFGyr->fPtr)
    ((FPtr5)_pluginJacqDotFGyr->fPtr)(_ndof, &(*position)(0), &(*velocity)(0), &(*_jacobianFGyrqDot)(0, 0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeRhs(double time, bool isDSup)
{
  // if isDSup == true, this means that there is no need to re-compute mass ...
  *_q[2] = *(_p[2]); // Warning: r/p update is done in Interactions/Relations

  // if(_forces)
  //   {
  computeForces(time, _q[0], _q[1]);
  *_q[2] += *_forces;
  //#  }

  // Computes q[2] = inv(mass)*(fL+p) by solving Mq[2]=fL+p.
  // -- Case 1: if mass is constant, then a copy of imass is LU-factorized during initialization and saved into _inverseMasse
  // -- Case 2: mass is not constant, we copy it into _inverseMass
  // Then we proceed with PLUForwardBackward.
  // mass and inv(mass) computation
  if(_mass && !isDSup && !_hasConstantMass)  // if it is necessary to re-compute mass, FInt ..., ie if they have not been compiled during the present time step
  {
    computeMass();
    *_inverseMass = *_mass;
  }

  //  if(mass->isPlugged()) : mass may be not plugged in LagrangianDS children
  if(_inverseMass)
    _inverseMass->PLUForwardBackwardInPlace(*_q[2]);

  _x[1]->setBlock(0, *_q[1]);
  _x[1]->setBlock(_ndof, *_q[2]);
}

void LagrangianDS::computeJacobianRhsx(double time, bool isDSup)
{
  // if isDSup == true, this means that there is no need to re-compute mass ...

  if(!isDSup && !_hasConstantMass)
    computeMass();

  //  if(mass->isPlugged()) : mass may b not plugged in LagrangianDS children

  if(_jacobianqForces || _jacobianqDotForces)
  {
    if(!_hasConstantMass) // else inverseMass is already uptodate
      *_inverseMass = *_mass;
  }

  if(_jacobianqForces)
  {
    SP::SiconosMatrix bloc10 = _jacxRhs->block(1, 0);
    computeJacobianqForces(time);
    *bloc10 = *_jacobianqForces;
    _inverseMass->PLUForwardBackwardInPlace(*bloc10);
  }

  if(_jacobianqDotForces)
  {
    SP::SiconosMatrix bloc11 = _jacxRhs->block(1, 1);
    computeJacobianqDotForces(time);
    *bloc11 = *_jacobianqDotForces;
    _inverseMass->PLUForwardBackwardInPlace(*bloc11);
  }
}

void LagrangianDS::computeForces(double time, SP::SiconosVector position, SP::SiconosVector velocity)
{
  if(!_forces)
  {
    _forces.reset(new SiconosVector(_ndof));
  }
  else
    _forces->zero();

  // 1 - Computes the required function
  computeFInt(time, position, velocity);
  computeFExt(time);
  computeFGyr(position, velocity);

  // seems ok.
  // if (_forces.use_count() == 1)
  // {
  //   //if not that means that fL is already (pointer-)connected with
  //   // either fInt, FGyr OR fExt.
  //_forces->zero();

  if(_fInt)
    *_forces -= *_fInt;

  if(_fExt)
    *_forces += *_fExt;

  if(_fGyr)
    *_forces -= *_fGyr;
  // }

}

void LagrangianDS::computeJacobianqForces(double time)
{
  if(_jacobianqForces)
  {
    computeJacobianFIntq(time);
    computeJacobianFGyrq();

    // not true!
    // if( jacobianFL[i].use_count() == 1 )
    {
      //if not that means that jacobianFL[i] is already (pointer-)connected with
      // either jacobianFInt or jacobianFGyr
      _jacobianqForces->zero();
      if(_jacobianFIntq)
        *_jacobianqForces -= *_jacobianFIntq;
      if(_jacobianFGyrq)
        *_jacobianqForces -= *_jacobianFGyrq;
    }
  }
  //else nothing.
}
void LagrangianDS::computeJacobianqDotForces(double time)
{
  if(_jacobianqDotForces)
  {
    computeJacobianFIntqDot(time);
    computeJacobianFGyrqDot();

    // not true!
    // if( jacobianFL[i].use_count() == 1 )
    {
      //if not that means that jacobianFL[i] is already (pointer-)connected with
      // either jacobianFInt or jacobianFGyr
      _jacobianqDotForces->zero();
      if(_jacobianFIntqDot)
        *_jacobianqDotForces -= *_jacobianFIntqDot;
      if(_jacobianFGyrqDot)
        *_jacobianqDotForces -= *_jacobianFGyrqDot;
    }
  }
  //else nothing.
}
// void LagrangianDS::computeJacobianZFL( double time){
//    RuntimeException::selfThrow("LagrangianDS::computeJacobianZFL - not implemented");
// }

void LagrangianDS::display() const
{
  std::cout << "=====> Lagrangian System display (number: " << _number << ")." <<std::endl;
  std::cout << "- _ndof : " << _ndof <<std::endl;
  std::cout << "- q " <<std::endl;
  if(_q[0]) _q[0]->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- q0 " <<std::endl;
  if(_q0) _q0->display();
  std::cout << "- velocity " <<std::endl;
  if(_q[1]) _q[1]->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- acceleration " <<std::endl;
  if(_q[2]) _q[2]->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- v0 " <<std::endl;
  if(_velocity0) _velocity0->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- p[0] " <<std::endl;
  if(_p[0]) _p[0]->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- p[1] " <<std::endl;
  if(_p[1]) _p[1]->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- p[2] " <<std::endl;
  if(_p[2]) _p[2]->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "===================================== " <<std::endl;
}

// --- Functions for memory handling ---
void LagrangianDS::initMemory(unsigned int steps)
{
  DEBUG_PRINTF("LagrangianDS::initMemory(unsigned int steps) with steps = %i\n", steps);
  if(steps == 0)
    std::cout << "Warning : LagragianDS::initMemory with size equal to zero" <<std::endl;
  else
  {
    _qMemory.reset(new SiconosMemory(steps, _ndof));
    _velocityMemory.reset(new SiconosMemory(steps, _ndof));
    _forcesMemory.reset(new SiconosMemory(steps, _ndof));
    _pMemory.resize(3);

    //TODO : initMemory in graph + f(OSI/level)
    for(unsigned int level=0; level<3; ++level)
    {
      if(!_pMemory[level])
        _pMemory[level].reset(new SiconosMemory(steps, _ndof));
    }

    //swapInMemory();
  }
}

void LagrangianDS::swapInMemory()
{
  _qMemory->swap(*_q[0]);
  _velocityMemory->swap(*_q[1]);
  if(_forces && _forcesMemory)
    _forcesMemory->swap(*_forces);
  if(_p[0] && _pMemory[0])
  {
    _pMemory[0]->swap(*_p[0]);
  }
  if(_p[1] && _pMemory[1])
  {
    _pMemory[1]->swap(*_p[1]);
  }
  if(_p[2] && _pMemory[2])
  {
    _pMemory[2]->swap(*_p[2]);
  }
  if(_x[0] && _xMemory)
    _xMemory->swap(*_x[0]);
}

void LagrangianDS::resetAllNonSmoothParts()
{
  if(_p[0])
    _p[0]->zero();
  if(_p[1])
    _p[1]->zero();
  if(_p[2])
    _p[2]->zero();
}

void LagrangianDS::resetNonSmoothPart(unsigned int level)
{
  if(level < LEVELMAX)
    if(_p[level])
      _p[level]->zero();
}

void LagrangianDS::computePostImpactVelocity()
{
  // When this function is call, q[1] is supposed to be pre-impact velocity.
  // We solve M(v+ - v-) = p - The result is saved in(place of) p[1].
  DEBUG_BEGIN("LagrangianDS::computePostImpactV()\n");
  SiconosVector tmp(*_p[1]);
  if(_inverseMass)
    _inverseMass->PLUForwardBackwardInPlace(tmp);
  *_q[1] += tmp;  // v+ = v- + p
  DEBUG_BEGIN("LagrangianDS::computePostImpactV() END \n");
}

void LagrangianDS::setComputeFGyrFunction(const std::string& pluginPath, const std::string&  functionName)
{
  _pluginFGyr->setComputeFunction(pluginPath, functionName);
  if(!_fGyr)
    _fGyr.reset(new SiconosVector(_ndof));
  init_forces();
}

void LagrangianDS::setComputeFGyrFunction(FPtr5 fct)
{
  _pluginFGyr->setComputeFunction((void *)fct);
  if(!_fGyr)
    _fGyr.reset(new SiconosVector(_ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFIntqFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  _pluginJacqFInt->setComputeFunction(pluginPath, functionName);
  if(!_jacobianFIntq)
    _jacobianFIntq.reset(new SimpleMatrix(_ndof, _ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFIntqDotFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  _pluginJacqDotFInt->setComputeFunction(pluginPath, functionName);
  if(!_jacobianFIntqDot)
    _jacobianFIntqDot.reset(new SimpleMatrix(_ndof, _ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFIntqFunction(FPtr6 fct)
{
  _pluginJacqFInt->setComputeFunction((void *)fct);
  if(!_jacobianFIntq)
    _jacobianFIntq.reset(new SimpleMatrix(_ndof, _ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFIntqDotFunction(FPtr6 fct)
{
  _pluginJacqDotFInt->setComputeFunction((void *)fct);
  if(!_jacobianFIntqDot)
    _jacobianFIntqDot.reset(new SimpleMatrix(_ndof, _ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFGyrqFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  _pluginJacqFGyr->setComputeFunction(pluginPath, functionName);
  if(!_jacobianFGyrq)
    _jacobianFGyrq.reset(new SimpleMatrix(_ndof, _ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFGyrqDotFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  _pluginJacqDotFGyr->setComputeFunction(pluginPath, functionName);
  if(!_jacobianFGyrqDot)
    _jacobianFGyrqDot.reset(new SimpleMatrix(_ndof, _ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFGyrqFunction(FPtr5 fct)
{
  _pluginJacqFGyr->setComputeFunction((void *)fct);
  if(!_jacobianFGyrq)
    _jacobianFGyrq.reset(new SimpleMatrix(_ndof, _ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFGyrqDotFunction(FPtr5 fct)
{
  _pluginJacqDotFGyr->setComputeFunction((void *)fct);
  if(!_jacobianFGyrqDot)
    _jacobianFGyrqDot.reset(new SimpleMatrix(_ndof, _ndof));
}

double LagrangianDS::computeKineticEnergy()
{
  DEBUG_BEGIN("NewtonEulerDS::computeKineticEnergy()\n");
  SP::SiconosVector velo = velocity();
  assert(velo);
  DEBUG_EXPR(velo->display());

  SP::SiconosVector tmp(new SiconosVector(*velo));
  if(_mass)
    prod(*_mass, *velo, *tmp, true);

  double K =0.5*inner_prod(*tmp,*velo);

  DEBUG_PRINTF("Kinetic Energy = %e\n", K);
  DEBUG_END("LagrangianDS::computeKineticEnergy()\n");
  return K;
}

void LagrangianDS::setBoundaryConditions(SP::BoundaryCondition newbd)
{
  if(!_boundaryConditions)
  {
    std::cout << "Warning : LagrangianDS::setBoundaryConditions. old boundary conditions were pre-existing" <<std::endl;
  }
  _boundaryConditions = newbd;
  _reactionToBoundaryConditions.reset(new SiconosVector(_boundaryConditions->velocityIndices()->size()));
};
