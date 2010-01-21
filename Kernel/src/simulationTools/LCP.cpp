/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
#include "LCP.hpp"

using namespace std;
using namespace RELATION;

int LCP::compute(double time)
{
  // --- Prepare data for LCP computing ---
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
    // The LCP in Numerics format
    LinearComplementarity_Problem numerics_problem;
    numerics_problem.M = &*_M->getNumericsMatrix();
    numerics_problem.q = _q->getArray();
    numerics_problem.size = _sizeOutput;
    int nbSolvers = 1;
    const char * name = &*_solver->numericsSolverOptions()->solverName;
    if ((strcmp(name , "ENUM") == 0))
    {
      lcp_enum_init(&numerics_problem, &*_solver->numericsSolverOptions(), 1);


    }

    // Call LCP Driver
    info = lcp_driver(&numerics_problem, _z->getArray() , _w->getArray() ,
                      &*_solver->numericsSolverOptions(), nbSolvers, &*_numerics_options);

    if ((strcmp(name , "ENUM") == 0))
    {
      lcp_enum_reset(&numerics_problem, &*_solver->numericsSolverOptions(), 1);


    }

    // --- Recovering of the desired variables from LCP output ---
    postCompute();

  }

  return info;
}

void LCP::display() const
{
  cout << "======= LCP of size " << _sizeOutput << " with: " << endl;
  LinearOSNS::display();
}

LCP* LCP::convert(OneStepNSProblem* osnsp)
{
  LCP* lcp = dynamic_cast<LCP*>(osnsp);
  return lcp;
}

void LCP::initialize(SP::Simulation sim)
{
  // General initialize for LinearOSNS
  LinearOSNS::initialize(sim);

  // Initialization of the NonSmoothSolver
  _solver->initialize(this);

}
