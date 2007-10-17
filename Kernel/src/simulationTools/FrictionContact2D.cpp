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

void FrictionContact2D::compute(const double time)
{
  // --- Prepare data for FrictionContact2D computing ---
  preCompute(time);

  // --- Call Numerics solver ---
  if (sizeOutput != 0)
  {

    int info;
    int SizeOutput = (int)sizeOutput;
    // get solving method and friction coefficient value.
    method solvingMethod = *(solver->getSolvingMethodPtr());
    Interaction * currentInteraction = simulation->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0);
    // call Numerics method for 2D or 3D problem:

    solvingMethod.pfc_2D.mu = static_cast<NewtonImpactFrictionNSL*>(currentInteraction->getNonSmoothLawPtr())->getMu();
    info = pfc_2D_solver(M->getArray(), q->getArray(), &SizeOutput, &solvingMethod  , z->getArray(), w->getArray());

    check_solver(info);
    // --- Recover the desired variables from FrictionContact2D output ---
    postCompute();
  }
}

FrictionContact2D* FrictionContact2D::convert(OneStepNSProblem* osnsp)
{
  FrictionContact2D* fc2d = dynamic_cast<FrictionContact2D*>(osnsp);
  return fc2d;
}


