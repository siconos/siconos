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
#include "Equality.hpp"
#include "Simulation.hpp"
#include "OSNSMatrix.hpp"
#include <NumericsMatrix.h>

using namespace RELATION;

int Equality::compute(double time)
{
  int info = 0;
  // --- Prepare data for EQUALITY computing ---
  bool cont = preCompute(time);
  if (!cont)
    return info;


  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)

  if (_sizeOutput != 0)
  {
    double* q_ = q()->getArray();
    double* z_ =  _z->getArray();
    for (size_t i = 0; i < _sizeOutput; ++i) z_[i] = -q_[i];
    info = NM_gesv(&*_M->numericsMatrix(), z_, true);
    // --- Recovering of the desired variables from EQUALITY output ---
    postCompute();

  }

  return info;
}

void Equality::initialize(SP::Simulation sim)
{
  // General initialize for LinearOSNS
  LinearOSNS::initialize(sim);
  //SP::InteractionsGraph indexSet = simulation()->indexSet(levelMin());
  //_M.reset(new OSNSMatrix(indexSet,_numericsMatrixStorageType));
}

void Equality::updateM()
{
  assert(0);
  // Get index set from Simulation
  InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());

  if (!_M)
  {
    // Creates and fills M using Interactionof indexSet
    _M.reset(new OSNSMatrix(indexSet, _numericsMatrixStorageType));
  }
  else
  {
    _M->setStorageType(_numericsMatrixStorageType);
    _M->fillW(indexSet);

  }
  _sizeOutput = _M->size();
}


void Equality::display() const
{
  std::cout << "======= EQUALITY of size " << _sizeOutput << " with: " <<std::endl;
  LinearOSNS::display();
}
