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
#include "LagrangianLinearDiagonalDS.hpp"
#include "BlockMatrix.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

#include <iostream>

// --- Constructor for the complete system
LagrangianLinearDiagonalDS::LagrangianLinearDiagonalDS(SP::SiconosVector q0, SP::SiconosVector velocity0,
                                                       SP::SiconosVector stiffness, SP::SiconosVector damping,
                                                       SP::SiconosVector mass):
  LagrangianDS(q0, velocity0)
{
  _mass.reset(new SimpleMatrix(dimension(), dimension(), Siconos::BANDED, 0, 0));
  for(unsigned int i = 0;i<dimension(); ++i)
    (*_mass)(i, i) = (*mass)(i);
  _stiffness = stiffness;
  _damping = damping;
}

// --- Constructor for the complete system with identity mass matrix
LagrangianLinearDiagonalDS::LagrangianLinearDiagonalDS(SP::SiconosVector q0, SP::SiconosVector velocity0,
                                                       SP::SiconosVector stiffness, SP::SiconosVector damping):
  LagrangianDS(q0, velocity0)
{
  _stiffness = stiffness;
  _damping = damping;
}

// --- Constructor for the undamped system with identity mass matrix
LagrangianLinearDiagonalDS::LagrangianLinearDiagonalDS(SP::SiconosVector q0, SP::SiconosVector velocity0,
                                                       SP::SiconosVector stiffness):
  LagrangianDS(q0, velocity0)
{
  _stiffness = stiffness;
}

void LagrangianLinearDiagonalDS::initRhs(double time)
{
  RuntimeException::selfThrow("LagrangianLinearDiagonalDS::initRhs - not yet implemented for LagrangianLinearDiagonalDS.");

}

void LagrangianLinearDiagonalDS::computeForces(double time, SP::SiconosVector q2, SP::SiconosVector v2)
{
  DEBUG_PRINT("LagrangianLinearTIDS::computeForces(double time, SP::SiconosVector q2, SP::SiconosVector v2) \n");

  if (!_forces)
  {
    _forces.reset(new SiconosVector(_ndof));
  }
  else
    _forces->zero();

  if (_fExt)
  {
    computeFExt(time);
    *_forces += *_fExt;
  }

  if (_stiffness)
    for(unsigned int i=0;i<_ndof;++i)
      (*_forces)(i) -= (*_stiffness)(i) * (*q2)(i);
  if (_damping)
    for(unsigned int i=0;i<_ndof;++i)
      (*_forces)(i) -= (*_damping)(i) * (*v2)(i);

  // if (_stiffness)
  // {
  //   DenseVect & stiff = *_stiffness->dense();
  //   noalias(*force->dense()) -= element_prod(stiff, *q2->dense());
  // }
  // if (_damping)
  // {
  //   DenseVect & damp = *_damping->dense();
  //   noalias(*force->dense()) -= element_prod(damp, *v2->dense());
  // }

  // Test this:
  // DenseVect& omega = *_stiffness.dense();
  // DenseVect& sigma = *_damping.dense();
  // DenseVect& res = *_forces.dense();
  // DenseVect& qvect = *q2.dense();
  // DenseVect& vvect = *v2.dense();
  // if (_stiffness)
  //   std::transform(omega.begin(), omega.end(),
  //                  qvect.begin(), res.begin(), std::multiplies<double>() );
  // if (_damping)
  //   std::transform(sigma.begin(), sigma.end(),
  //                  vvect.begin(), res.begin(), std::multiplies<double>() );
  // *forces *= -1.;
}




void LagrangianLinearDiagonalDS::display() const
{
  LagrangianDS::display();
  std::cout << "===== Lagrangian Linear Diagonal System display ===== " <<std::endl;
  std::cout << "- Mass Matrix M : " <<std::endl;
  if (_mass) _mass->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- Stiffness Matrix K : " <<std::endl;
  if (_stiffness) _stiffness->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- Viscosity Matrix C : " <<std::endl;
  if (_damping) _damping->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "=========================================================== " <<std::endl;
}

