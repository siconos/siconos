/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#include "Relay.h"


using namespace std;
using namespace RELATION;


void Relay::initialize(SP::Simulation sim)
{
  LinearOSNS::initialize(sim);

}

int Relay::compute(double time)
{
  // --- Prepare data for Relay computing ---
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
    // The Relay in Numerics format
    Relay_Problem numerics_problem;
    numerics_problem.M = &*_M->getNumericsMatrix();
    numerics_problem.q = _q->getArray();
    //numerics_problem.lb = _lb->getArray();
    //numerics_problem.ub = _ub->getArray();
    numerics_problem.size = _sizeOutput;

    int nbSolvers = 1;
    // Call Relay Driver
    info = relay_driver(&numerics_problem, _z->getArray() , _w->getArray() ,
                        &*_solver->numericsSolverOptions(), &*_numerics_options);

    // --- Recovering of the desired variables from Relay output ---
    postCompute();

  }

  return info;
}

void Relay::display() const
{
  cout << "======= Relay of size " << _sizeOutput << " with: " << endl;
  LinearOSNS::display();
}

Relay* Relay::convert(OneStepNSProblem* osnsp)
{
  Relay* lcp = dynamic_cast<Relay*>(osnsp);
  return lcp;
}


