/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include "FrictionContact.hpp"
#include "Topology.hpp"
#include "Simulation.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "OSNSMatrix.hpp"
#include "NonSmoothDrivers.h" // from numerics, for fcX_driver
#include <fc2d_Solvers.h>
#include <fc3d_Solvers.h>

using namespace RELATION;


FrictionContact::FrictionContact(int dimPb, int numericsSolverId):
  FrictionContact(dimPb, SP::SolverOptions(solver_options_create(numericsSolverId),
                  solver_options_delete))
{}

FrictionContact::FrictionContact(int dimPb, SP::SolverOptions options):
  LinearOSNS(options), _contactProblemDim(dimPb)
{
  if(dimPb == 2 && options->solverId == SICONOS_FRICTION_3D_NSGS)
  {
    _numerics_solver_options.reset(solver_options_create(SICONOS_FRICTION_2D_NSGS),
                                   solver_options_delete);
  }

  if(dimPb == 2)
  {
    _frictionContact_driver = &fc2d_driver;
  }
  else if(dimPb == 3)
  {
    _frictionContact_driver = &fc3d_driver;
  }
  else
    THROW_EXCEPTION("Wrong dimension value (must be 2 or 3) for FrictionContact constructor.");

  _mu.reset(new MuStorage());
}

void FrictionContact::initialize(SP::Simulation sim)
{
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  LinearOSNS::initialize(sim);

  // Connect to the right function according to dim. of the problem

  // get topology
  SP::Topology topology =
    simulation()->nonSmoothDynamicalSystem()->topology();

  // Note that interactionBlocks is up to date since updateInteractionBlocks
  // has been called during OneStepNSProblem::initialize()

  // Fill vector of friction coefficients
  int sizeMu = simulation()->nonSmoothDynamicalSystem()
               ->topology()->indexSet(0)->size();
  _mu->reserve(sizeMu);

  // If the topology is TimeInvariant ie if M structure does not
  // change during simulation:

  if(topology->indexSet0()->size()>0)
  {
    // Get index set from Simulation
    SP::InteractionsGraph indexSet =
      simulation()->indexSet(indexSetLevel());
    InteractionsGraph::VIterator ui, uiend;
    for(std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
    {
      _mu->push_back(std::static_pointer_cast<NewtonImpactFrictionNSL>
                     (indexSet->bundle(*ui)->nonSmoothLaw())->mu());
    }
  }
}

void FrictionContact::updateMu()
{
  _mu->clear();
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  InteractionsGraph::VIterator ui, uiend;
  for(std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    _mu->push_back(std::static_pointer_cast<NewtonImpactFrictionNSL>
                   (indexSet->bundle(*ui)->nonSmoothLaw())->mu());
  }
}

SP::FrictionContactProblem FrictionContact::frictionContactProblem()
{
  SP::FrictionContactProblem numerics_problem(new FrictionContactProblem());
  numerics_problem->dimension = _contactProblemDim;
  numerics_problem->numberOfContacts = _sizeOutput / _contactProblemDim;
  numerics_problem->M = &*_M->numericsMatrix();
  numerics_problem->q = &*_q->getArray();
  numerics_problem->mu = _mu->data();
  return numerics_problem;
}

FrictionContactProblem *FrictionContact::frictionContactProblemPtr()
{
  FrictionContactProblem *numerics_problem = &_numerics_problem;
  numerics_problem->dimension = _contactProblemDim;
  numerics_problem->numberOfContacts = _sizeOutput / _contactProblemDim;
  numerics_problem->M = &*_M->numericsMatrix();
  numerics_problem->q = &*_q->getArray();
  numerics_problem->mu = _mu->data();
  return numerics_problem;
}

int FrictionContact::solve(SP::FrictionContactProblem problem)
{
  if(!problem)
  {
    problem = frictionContactProblem();
  }

  return (*_frictionContact_driver)(&*problem,
                                    &*_z->getArray(),
                                    &*_w->getArray(),
                                    &*_numerics_solver_options);
}




bool FrictionContact::checkCompatibleNSLaw(NonSmoothLaw& nslaw)
{

  float type_number= (float) (Type::value(nslaw) + 0.1 * nslaw.size());
  _nslawtype.insert(type_number);

  if (Type::value(nslaw) != Type::NewtonImpactFrictionNSL)
  {
    THROW_EXCEPTION("\nFrictionContact::checkCompatibleNSLaw -  \n\
                      The chosen nonsmooth law is not compatible with FrictionalContact one step nonsmooth problem. \n\
                      Compatible NonSmoothLaw are: NewtonImpactFrictionNSL (2D or 3D) \n");
    return false;
  }
  if (_nslawtype.size() > 1)
  {
    THROW_EXCEPTION("\nFrictionContact::checkCompatibleNSLaw -  \n\
                     Compatible NonSmoothLaw are: NewtonImpactFrictionNSL (2D or 3D), but you cannot mix them \n");
    return false;
  }

  return true;
}

int FrictionContact::compute(double time)
{
  int info = 0;
  // --- Prepare data for FrictionContact computing ---
  bool cont = preCompute(time);
  if(!cont)
  {
    return info;
  }
  // nothing to do
  if(indexSetLevel() == LEVELMAX)
  {
    return info;
  }

  updateMu();

  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)
  if(_sizeOutput != 0)
  {
    // Call Numerics Driver for FrictionContact
    info = solve();
    postCompute();
  }

  return info;
}

void FrictionContact::display() const
{
  std::cout << "===== " << _contactProblemDim << "D Friction Contact Problem " <<std::endl;
  std::cout << "of size " << _sizeOutput << "(ie " << _sizeOutput / _contactProblemDim << " contacts)." <<std::endl;
  LinearOSNS::display();
}
