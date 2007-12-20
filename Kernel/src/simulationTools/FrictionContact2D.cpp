/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
#include "FrictionContact2D.h"
#include "Interaction.h"
#include "Simulation.h"
#include "Model.h"
#include "NonSmoothDynamicalSystem.h"

#include "NewtonImpactFrictionNSL.h"

using namespace std;

// xml constructor
FrictionContact2D::FrictionContact2D(OneStepNSProblemXML* osNsPbXml, Simulation* newSimu):
  FrictionContact("FrictionContact2D", osNsPbXml, newSimu)
{}

// Constructor from a set of data
FrictionContact2D::FrictionContact2D(Simulation* newSimu, const string newId, const string solverName,
                                     const unsigned int MaxIter, const double  Tolerance, const unsigned int  Verbose, const string  NormType,
                                     const double  SearchDirection): FrictionContact("FrictionContact2D", newSimu, newId)
{
  // set solver:
  solver = new Solver(nspbType, solverName, MaxIter, Tolerance, Verbose, NormType, SearchDirection);
  isSolverAllocatedIn = true;
}

// Constructor from a set of data
FrictionContact2D::FrictionContact2D(Solver*  newSolver, Simulation* newSimu, const string newId):
  FrictionContact("FrictionContact2D", newSolver, newSimu, newId)
{}

// destructor
FrictionContact2D::~FrictionContact2D()
{}

int FrictionContact2D::compute(double time)
{
  // --- Prepare data for FrictionContact2D computing ---
  preCompute(time);

  int info = 0;
  // --- Call Numerics solver ---
  if (sizeOutput != 0)
  {

    clock_t startSolve;
    // get solving method.
    method solvingMethod = *(solver->getSolvingMethodPtr());
    startSolve = clock();
    if (solver->useBlocks()) // Use solver block
      //  info = pfc_2D_solver_block((int)sizeOutput/2, Mspbl , q->getArray() , &solvingMethod , z->getArray() , w->getArray() ,mu->getArray());
      RuntimeException::selfThrow("FrictionContact2D::compute - Solver block not yet available.");
    else // Use classical solver
      info = pfc_2D_solver((int)sizeOutput / 2, M->getArray(), q->getArray(), &solvingMethod  , z->getArray(), w->getArray(), mu->getArray());
    CPUtime += (clock() - startSolve);
    nbIter++;

    // Remark : the output result, info, from solver is treated when this function is call from Simulation
    // --- Recover the desired variables from FrictionContact2D output ---
    postCompute();
  }
  return info;
}

FrictionContact2D* FrictionContact2D::convert(OneStepNSProblem* osnsp)
{
  FrictionContact2D* fc2d = dynamic_cast<FrictionContact2D*>(osnsp);
  return fc2d;
}


