/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
#include "FrictionContact3D.h"
#include "Interaction.h"
#include "Simulation.h"
#include "Model.h"
#include "NonSmoothDynamicalSystem.h"

#include "NewtonImpactFrictionNSL.h"

using namespace std;

// Default constructor
FrictionContact3D::FrictionContact3D(): FrictionContact("FrictionContact3D")
{}

// xml constructor (Simulation is optional)
FrictionContact3D::FrictionContact3D(OneStepNSProblemXML* osNsPbXml, Simulation* newSimu):
  FrictionContact("FrictionContact3D", osNsPbXml, newSimu)
{}

// From data (the only required argument is the simulation)
FrictionContact3D::FrictionContact3D(Simulation * newSimu, const string newId, const string newSolver, const unsigned int MaxIter,
                                     const double  Tolerance, const unsigned int Verbose,  const string  NormType,
                                     const double  SearchDirection): FrictionContact("FrictionContact3D", newSimu, newId)
{
  // set solver:
  solver = new Solver(nspbType, newSolver, MaxIter, Tolerance, Verbose, NormType, SearchDirection);
  isSolverAllocatedIn = true;
}

// Constructor from a set of data
FrictionContact3D::FrictionContact3D(Solver*  newSolver, Simulation* newSimu, const string newId):
  FrictionContact("FrictionContact3D", newSolver, newSimu, newId)
{}

// destructor
FrictionContact3D::~FrictionContact3D()
{}

void FrictionContact3D::compute(const double time)
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

    SizeOutput = SizeOutput / 3; // in pfc_3D, SizeOutput is the number of contact points.
    solvingMethod.pfc_3D.mu = static_cast<NewtonImpactFrictionNSL*>(currentInteraction->getNonSmoothLawPtr())->getMu();
    info = pfc_3D_solver(M->getArray(), q->getArray(), &SizeOutput, &solvingMethod  , z->getArray(), w->getArray());

    check_solver(info);
    // --- Recover the desired variables from FrictionContact2D output ---
    postCompute(w, z);
  }
}

FrictionContact3D* FrictionContact3D::convert(OneStepNSProblem* osnsp)
{
  FrictionContact3D* fc3d = dynamic_cast<FrictionContact3D*>(osnsp);
  return fc3d;
}


