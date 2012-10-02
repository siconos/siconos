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
#include "LCP.hpp"

#include "debug.h"
#define DEBUG_MESSAGES 1

using namespace std;
using namespace RELATION;

LCP::LCP(SP::OneStepNSProblemXML onestepnspbxml) :
  LinearOSNS(onestepnspbxml, "LCP")
{

  if (onestepnspbxml->hasNumericsSolverName())
    _numerics_solver_id = nameToId((char *)onestepnspbxml->getNumericsSolverName().c_str());
  else
    _numerics_solver_id = SICONOS_LCP_LEMKE;

  _numerics_problem.reset(new LinearComplementarityProblem);

  linearComplementarity_setDefaultSolverOptions(NULL, &*_numerics_solver_options, _numerics_solver_id);
}

LCP::LCP(const int newNewNumericsSolverId , const std::string& newId):
  LinearOSNS(newNewNumericsSolverId, "LCP", newId)
{
  _numerics_problem.reset(new  LinearComplementarityProblem);


  linearComplementarity_setDefaultSolverOptions(NULL, &*_numerics_solver_options, _numerics_solver_id);

}


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
    DEBUG_PRINTF("LCP : sizeOutput=%d\n", _sizeOutput);

    // The LCP in Numerics format
    _numerics_problem->M = &*_M->getNumericsMatrix();
    _numerics_problem->q = _q->getArray();
    _numerics_problem->size = _sizeOutput;

    //const char * name = &*_numerics_solver_options->solverName;
    if (_numerics_solver_options->solverId == SICONOS_LCP_ENUM)
    {
      lcp_enum_init(&*_numerics_problem, &*_numerics_solver_options, 1);


    }
    info = linearComplementarity_driver(&*_numerics_problem, _z->getArray() , _w->getArray() ,
                                        &*_numerics_solver_options, &*_numerics_options);

    if (_numerics_solver_options->solverId == SICONOS_LCP_ENUM)
    {
      lcp_enum_reset(&*_numerics_problem, &*_numerics_solver_options, 1);


    }



    // --- Recovering of the desired variables from LCP output ---
    postCompute();


    DEBUG_EXPR(display());

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
}

LCP::~LCP()
{
  deleteSolverOptions(&*_numerics_solver_options);
}
