/* Siconos-Kernel, Copyright INRIA 2005-2010.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY ory FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include "GenericMechanical.hpp"
#include "FrictionContactXML.hpp"
#include "Topology.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"

using namespace std;
using namespace RELATION;


GenericMechanical::GenericMechanical():
  LinearOSNS()
{
  _numerics_problem.reset(new  GenericMechanicalProblem);
}


void GenericMechanical::initialize(SP::Simulation sim)
{
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  LinearOSNS::initialize(sim);




}

int GenericMechanical::compute(double time)
{
  int info = 0;
  // --- Prepare data for GenericMechanical computing ---
  preCompute(time);



  /*
    La matrice _M est construite.
    Ici, il faut construire les sous-problemes, c'est a dire completer les champs des NumericsProblem (_mu, _e, _en, les dimentions...).
    Il faut aussi remplir la sous matrice M du sous-probleme.
    Pour cela, on peut boucler sur les interactions et completer le membres _numerics_problem.problems[i] and _numerics_problem.problemsType[i].
   */

  //......


  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)
  if (_sizeOutput != 0)
  {
    // The GenericMechanical Problem in Numerics format
    GenericMechanicalProblem numerics_problem;
    numerics_problem.M = &*_M->getNumericsMatrix();
    numerics_problem.q = &*_q->getArray();

    // Call Numerics Driver for GenericMechanical
    info = genericMechanical_driver(&numerics_problem,
                                    &*_z->getArray() ,
                                    &*_w->getArray() ,
                                    &*_numerics_solver_options);
    postCompute();

  }

  return info;
}

void GenericMechanical::display() const
{
  cout << "===== " << "Generic mechanical Problem " << endl;
  LinearOSNS::display();
}

GenericMechanical* GenericMechanical::convert(OneStepNSProblem* osnsp)
{
  GenericMechanical* fc2d = dynamic_cast<GenericMechanical*>(osnsp);
  return fc2d;
}

GenericMechanical::~GenericMechanical()
{
  deleteSolverOptions(&*_numerics_solver_options);
}


