/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#include "FrictionContact.h"
#include "FrictionContactXML.hpp"
#include "Topology.h"
#include "Simulation.h"
#include "Model.h"
#include "NonSmoothDynamicalSystem.h"
#include "NewtonImpactFrictionNSL.h"

using namespace std;
using namespace RELATION;

// xml constructor
FrictionContact::FrictionContact(SP::OneStepNSProblemXML osNsPbXml):
  LinearOSNS(osNsPbXml, "FrictionContact"), contactProblemDim(3)
{
  SP::FrictionContactXML xmlFC = boost::static_pointer_cast<FrictionContactXML>(osNsPbXml);

  // Read dimension of the problem (required parameter)
  if (!xmlFC->hasProblemDim())
    RuntimeException::selfThrow("FrictionContact: xml constructor failed, attribute for dimension of the problem (2D or 3D) is missing.");

  contactProblemDim = xmlFC->getProblemDim();

}

void FrictionContact::initialize(SP::Simulation sim)
{
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  LinearOSNS::initialize(sim);

  // Connect to the right function according to dim. of the problem
  if (contactProblemDim == 2)
    frictionContact_driver = &pfc_2D_driver;
  else // if(contactProblemDim == 3)
    frictionContact_driver = &frictionContact3D_driver;

  // get topology
  SP::Topology topology = simulation->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();

  // Note that unitaryBlocks is up to date since updateUnitaryBlocks has been called during OneStepNSProblem::initialize()

  // Fill vector of friction coefficients
  int sizeMu = simulation->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->getIndexSet0Ptr()->size();
  mu.reset(new MuStorage());
  mu->reserve(sizeMu);

  // If the topology is TimeInvariant ie if M structure does not change during simulation:
  if (topology->isTimeInvariant() &&   !OSNSInteractions->isEmpty())
  {
    // Get index set from Simulation
    SP::UnitaryRelationsSet indexSet = simulation->getIndexSetPtr(levelMin);
    for (ConstUnitaryRelationsIterator itUR = indexSet->begin(); itUR != indexSet->end(); ++itUR)
    {
      mu->push_back(boost::static_pointer_cast<NewtonImpactFrictionNSL> ((*itUR)->getInteractionPtr()->getNonSmoothLawPtr())->getMu());
    }
  }
}

int FrictionContact::compute(double time)
{
  int info = 0;
  // --- Prepare data for FrictionContact computing ---
  preCompute(time);

  // Update mu
  mu->clear();

  SP::UnitaryRelationsSet indexSet = simulation->getIndexSetPtr(levelMin);
  for (ConstUnitaryRelationsIterator itUR = indexSet->begin(); itUR != indexSet->end(); ++itUR)
  {
    mu->push_back(boost::static_pointer_cast<NewtonImpactFrictionNSL>((*itUR)->getInteractionPtr()->getNonSmoothLawPtr())->getMu());
  }
  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)
  if (sizeOutput != 0)
  {
    // The FrictionContact Problem in Numerics format
    FrictionContact_Problem numerics_problem;
    numerics_problem.M = &*M->getNumericsMatrix();
    numerics_problem.q = &*q->getArray();
    numerics_problem.numberOfContacts = sizeOutput / contactProblemDim;
    numerics_problem.isComplete = 1;
    numerics_problem.mu = &((*mu)[0]);
    // Call Numerics Driver for FrictionContact
    info = (*frictionContact_driver)(&numerics_problem, &*_z->getArray() , &*_w->getArray() , (solver->getNumericsSolverOptionsPtr()).get(), &*numerics_options);
    postCompute();

  }

  return info;
}

void FrictionContact::display() const
{
  cout << "===== " << contactProblemDim << "D Friction Contact Problem " << endl;
  cout << "of size " << sizeOutput << "(ie " << sizeOutput / contactProblemDim << " contacts)." << endl;
  LinearOSNS::display();
}

FrictionContact* FrictionContact::convert(OneStepNSProblem* osnsp)
{
  FrictionContact* fc2d = dynamic_cast<FrictionContact*>(osnsp);
  return fc2d;
}


