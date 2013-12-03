/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
#include "FrictionContact.hpp"
#include "FrictionContactXML.hpp"
#include "Topology.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "OSNSMatrix.hpp"

using namespace RELATION;


FrictionContact::FrictionContact(int dimPb, const int newNumericsSolverId,
                                 const std::string& newId):
  LinearOSNS(newNumericsSolverId, "FrictionContact", newId), _contactProblemDim(dimPb)
{
  if (dimPb == 2 && newNumericsSolverId == SICONOS_FRICTION_3D_NSGS)
    _numerics_solver_id = SICONOS_FRICTION_2D_NSGS;
  _numerics_problem.reset(new FrictionContactProblem);


  if (dimPb == 2)
    frictionContact2D_setDefaultSolverOptions(&*_numerics_solver_options, _numerics_solver_id);
  else if (dimPb == 3)
    frictionContact3D_setDefaultSolverOptions(&*_numerics_solver_options, _numerics_solver_id);
  else
    RuntimeException::selfThrow("cannot set defaults solver options for other problem dimension than 2 or 3");
}

// xml constructor
FrictionContact::FrictionContact(SP::OneStepNSProblemXML osNsPbXml):
  LinearOSNS(osNsPbXml, "FrictionContact"), _contactProblemDim(3)
{
  SP::FrictionContactXML xmlFC = std11::static_pointer_cast<FrictionContactXML>(osNsPbXml);


  if (osNsPbXml->hasNumericsSolverName())
    _numerics_solver_id = nameToId((char *)osNsPbXml->getNumericsSolverName().c_str());
  else
    _numerics_solver_id = SICONOS_FRICTION_3D_NSGS;

  _numerics_problem.reset(new  FrictionContactProblem);

  // Read dimension of the problem (required parameter)
  if (!xmlFC->hasProblemDim())
    RuntimeException::selfThrow("FrictionContact: xml constructor failed, attribute for dimension of the problem (2D or 3D) is missing.");

  _contactProblemDim = xmlFC->getProblemDim();
  if (_contactProblemDim == 2 && _numerics_solver_id == SICONOS_FRICTION_3D_NSGS) _numerics_solver_id = SICONOS_FRICTION_2D_NSGS;

  // initialize the _numerics_solver_options


  if (_contactProblemDim == 2)
    frictionContact2D_setDefaultSolverOptions(&*_numerics_solver_options, _numerics_solver_id);
  else // if(_contactProblemDim == 3)
    frictionContact3D_setDefaultSolverOptions(&*_numerics_solver_options, _numerics_solver_id);

}

void FrictionContact::initialize(SP::Simulation sim)
{
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  LinearOSNS::initialize(sim);

  // Connect to the right function according to dim. of the problem
  if (_contactProblemDim == 2)
    _frictionContact_driver = &frictionContact2D_driver;
  else // if(_contactProblemDim == 3)
    _frictionContact_driver = &frictionContact3D_driver;

  // get topology
  SP::Topology topology =
    simulation()->model()->nonSmoothDynamicalSystem()->topology();

  // Note that interactionBlocks is up to date since updateInteractionBlocks
  // has been called during OneStepNSProblem::initialize()

  // Fill vector of friction coefficients
  int sizeMu = simulation()->model()->nonSmoothDynamicalSystem()
               ->topology()->indexSet(0)->size();
  _mu.reset(new MuStorage());
  _mu->reserve(sizeMu);

  // If the topology is TimeInvariant ie if M structure does not
  // change during simulation:

  //if( !topology->hasChanged() &&   !interactions()->isEmpty())
  if (!topology->interactions()->isEmpty())
  {
    // Get index set from Simulation
    SP::InteractionsGraph indexSet =
      simulation()->indexSet(indexSetLevel());
    InteractionsGraph::VIterator ui, uiend;
    for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
    {
      _mu->push_back(std11::static_pointer_cast<NewtonImpactFrictionNSL>
                     (indexSet->bundle(*ui)->nonSmoothLaw())->mu());
    }
  }
}

void FrictionContact::updateMu()
{
  _mu->clear();
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    _mu->push_back(std11::static_pointer_cast<NewtonImpactFrictionNSL>
                   (indexSet->bundle(*ui)->nonSmoothLaw())->mu());
  }
}

SP::FrictionContactProblem FrictionContact::frictionContactProblem()
{
  SP::FrictionContactProblem numerics_problem(new FrictionContactProblem());
  numerics_problem->dimension = _contactProblemDim;
  numerics_problem->numberOfContacts = _sizeOutput / _contactProblemDim;
  numerics_problem->M = &*_M->getNumericsMatrix();
  numerics_problem->q = &*_q->getArray();
  numerics_problem->mu = &(_mu->at(0));
  return numerics_problem;
}

int FrictionContact::solve(SP::FrictionContactProblem problem)
{
  if (!problem)
  {
    problem = frictionContactProblem();
  }
  
  return (*_frictionContact_driver)(&*problem,
                                    &*_z->getArray(),
                                    &*_w->getArray(),
                                    &*_numerics_solver_options,
                                    &*_numerics_options);
}


int FrictionContact::compute(double time)
{
  int info = 0;
  // --- Prepare data for FrictionContact computing ---
  bool cont = preCompute(time);
  if (!cont)
  {
    return info;
  }
  // nothing to do
  if (indexSetLevel() == LEVELMAX)
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
  if (_sizeOutput != 0)
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

FrictionContact::~FrictionContact()
{
  deleteSolverOptions(&*_numerics_solver_options);
}


