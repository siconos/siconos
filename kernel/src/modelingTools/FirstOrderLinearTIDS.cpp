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
#include "FirstOrderLinearTIDS.hpp"

#include <iostream>

void FirstOrderLinearTIDS::initRhs(double time)
{
  if (_M && !_invM)
    _invM.reset(new SimpleMatrix(*_M));

  computeRhs(time);

  if (! _jacxRhs)  // if not allocated with a set or anything else
  {
    if (_A && ! _M)  // if M is not defined, then A = _jacxRhs, no memory allocation for that one.
      _jacxRhs = _A;
    else if (_A && _M)
    {
      _jacxRhs.reset(new SimpleMatrix(*_A)); // Copy A into _jacxRhs
      // Solve M_jacxRhs = A
      _invM->PLUForwardBackwardInPlace(*_jacxRhs);
    }
    // else no allocation, jacobian is equal to 0.
  }
}

void FirstOrderLinearTIDS::computeRhs(double time, const bool isDSup)
{

  *_x[1] = * _r; // Warning: r update is done in Interactions/Relations

  if (_A)
    prod(*_A, *_x[0], *_x[1], false);

  // compute and add b if required
  if (_b)
    *_x[1] += *_b;

  if (_M)
    {
      // allocate invM at the first call of the present function
      if (! _invM)
	_invM.reset(new SimpleMatrix(*_M));
      _invM->PLUForwardBackwardInPlace(*_x[1]);
    }
}

void FirstOrderLinearTIDS::computeJacobianRhsx(double time, const bool isDSup)
{
  // Nothing to be done: _jacxRhs is constant and computed during initialize. But this function is required to avoid call to base class function.
}

void FirstOrderLinearTIDS::display() const
{
  std::cout << "===> Linear Time-invariant First Order System display, " << _number << ")." <<std::endl;
  std::cout << "- A " <<std::endl;
  if (_A) _A->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- b " <<std::endl;
  if (_b) _b->display();
  else std::cout << "-> NULL" <<std::endl;

  std::cout << "- M: " <<std::endl;
  if (_M) _M->display();
  else std::cout << "-> NULL" <<std::endl;

  std::cout << "============================================" <<std::endl;
}
