/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include "Equality.hpp"
#include "Simulation.hpp"
using namespace std;
using namespace RELATION;

int Equality::compute(double time)
{
  // --- Prepare data for EQUALITY computing ---
  preCompute(time);

  int info = 0;
  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)

  if (_sizeOutput != 0)
  {
    // The EQUALITY in Numerics format
    // Call EQUALITY Driver
    _numerics_problem.q = q()->getArray();
    _numerics_problem.size = _sizeOutput;
    //      displayLS(&_numerics_problem);
    info = LinearSystem_driver(&_numerics_problem, _z->getArray() , _w->getArray() , 0);

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
  //_M.reset(new OSNSMatrix(indexSet,_MStorageType));
  _numerics_problem.M = &*_M->getNumericsMatrix();
}

void Equality::updateM()
{
  assert(0);
  // Get index set from Simulation
  SP::InteractionsGraph indexSet = simulation()->indexSet(levelMin());

  if (!_M)
  {
    // Creates and fills M using Interactionof indexSet
    _M.reset(new OSNSMatrix(indexSet, _MStorageType));
    _numerics_problem.M = &*_M->getNumericsMatrix();
  }
  else
  {
    _M->setStorageType(_MStorageType);
    _M->fill(indexSet);

  }
  _sizeOutput = _M->size();
}


void Equality::display() const
{
  cout << "======= EQUALITY of size " << _sizeOutput << " with: " << endl;
  LinearOSNS::display();
}

Equality* Equality::convert(OneStepNSProblem* osnsp)
{
  Equality* equality = dynamic_cast<Equality*>(osnsp);
  return equality;
}


