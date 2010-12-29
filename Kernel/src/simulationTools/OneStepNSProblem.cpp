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
#include "OneStepNSProblem.hpp"
#include "OneStepNSProblemXML.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "UnitaryRelation.hpp"
#include "Interaction.hpp"
#include "Topology.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "Moreau.hpp"
#include "LagrangianDS.hpp"
#include "NewtonEulerDS.hpp"

using namespace std;
OneStepNSProblem::OneStepNSProblem():
  _levelMin(0), _levelMax(0), _maxSize(0), _CPUtime(0), _nbIter(0)
{
  _numerics_solver_options.reset(new SolverOptions);
  _numerics_solver_options->iWork = NULL;
  _numerics_solver_options->dWork = NULL;
}
// --- CONSTRUCTORS/DESTRUCTOR ---
// xml constructor
OneStepNSProblem::OneStepNSProblem(const string& pbType,
                                   SP::OneStepNSProblemXML osnspbxml):
/*_nspbType(pbType),*/ _id(DEFAULT_OSNS_NAME), _sizeOutput(0),
  _onestepnspbxml(osnspbxml), _levelMin(0), _levelMax(0), _maxSize(0), _CPUtime(0), _nbIter(0)
{
  if (!_onestepnspbxml)
    RuntimeException::selfThrow("OneStepNSProblem::xml constructor, xml file == NULL");

  // === get dimension of the problem ===
  if (_onestepnspbxml->hasDim())
    _sizeOutput = _onestepnspbxml->getDimNSProblem();

  // === get Id ===

  if (_onestepnspbxml->hasId())
    _id = _onestepnspbxml->getId();

  if (_onestepnspbxml->hasNumericsSolverName())
    _id = _onestepnspbxml->getNumericsSolverName();


  // Numerics general options
  _numerics_options.reset(new NumericsOptions());
  _numerics_options->verboseMode = 0; // turn verbose mode to off by default

  _numerics_solver_options.reset(new SolverOptions);
  _numerics_solver_options->iWork = NULL;
  _numerics_solver_options->dWork = NULL;

  printf("OneStepNSProblem::OneStepNSProblem 1: Depressed inertface, first parameter ignored\n");
}
OneStepNSProblem::OneStepNSProblem(SP::OneStepNSProblemXML osnspbxml):
  _id(DEFAULT_OSNS_NAME), _sizeOutput(0),
  _onestepnspbxml(osnspbxml), _levelMin(0), _levelMax(0), _maxSize(0), _CPUtime(0), _nbIter(0)
{
  if (!_onestepnspbxml)
    RuntimeException::selfThrow("OneStepNSProblem::xml constructor, xml file == NULL");

  // === get dimension of the problem ===
  if (_onestepnspbxml->hasDim())
    _sizeOutput = _onestepnspbxml->getDimNSProblem();

  // === get Id ===

  if (_onestepnspbxml->hasId())
    _id = _onestepnspbxml->getId();

  if (_onestepnspbxml->hasNumericsSolverName())
    _id = _onestepnspbxml->getNumericsSolverName();


  // Numerics general options
  _numerics_options.reset(new NumericsOptions());
  _numerics_options->verboseMode = 0; // turn verbose mode to off by default

  _numerics_solver_options.reset(new SolverOptions);
  _numerics_solver_options->iWork = NULL;
  _numerics_solver_options->dWork = NULL;

}
// Constructor with given simulation and a pointer on Solver (Warning, solver is an optional argument)
OneStepNSProblem::OneStepNSProblem(const string& pbType, const string& newId, const int newNumericsSolverId):
  _numerics_solver_id(newNumericsSolverId),/*_nspbType(pbType),*/ _id(newId), _sizeOutput(0), _levelMin(0), _levelMax(0), _maxSize(0), _CPUtime(0), _nbIter(0)
{

  // Numerics general options
  _numerics_options.reset(new NumericsOptions());
  _numerics_options->verboseMode = 0; // turn verbose mode to off by default

  _numerics_solver_options.reset(new SolverOptions);
  _numerics_solver_options->iWork = NULL;
  _numerics_solver_options->dWork = NULL;
  _numerics_solver_options->solverId = newNumericsSolverId;
  printf("OneStepNSProblem::OneStepNSProblem 2: Depressed inertface, first parameter ignored, removed it.\n");
}

OneStepNSProblem::OneStepNSProblem(const string& newId, const int newNumericsSolverId):
  _numerics_solver_id(newNumericsSolverId), _id(newId), _sizeOutput(0), _levelMin(0), _levelMax(0), _maxSize(0), _CPUtime(0), _nbIter(0)
{

  // Numerics general options
  _numerics_options.reset(new NumericsOptions());
  _numerics_options->verboseMode = 0; // turn verbose mode to off by default

  _numerics_solver_options.reset(new SolverOptions);
  _numerics_solver_options->iWork = NULL;
  _numerics_solver_options->dWork = NULL;

}
OneStepNSProblem::OneStepNSProblem(const int newNumericsSolverId):
  _numerics_solver_id(newNumericsSolverId), _sizeOutput(0), _levelMin(0), _levelMax(0), _maxSize(0), _CPUtime(0), _nbIter(0)
{

  // Numerics general options
  _numerics_options.reset(new NumericsOptions());
  _numerics_options->verboseMode = 0; // turn verbose mode to off by default

  _numerics_solver_options.reset(new SolverOptions);
  _numerics_solver_options->iWork = NULL;
  _numerics_solver_options->dWork = NULL;

}

SP::SiconosMatrix OneStepNSProblem::dSBlock(SP::DynamicalSystem DS1) const
{

  ConstMatIterator itDS = _DSBlocks.find(DS1);
  return itDS->second;

}


void OneStepNSProblem::setDSBlocks(const MapOfDSMatrices& newMap)
{
  RuntimeException::selfThrow("OneStepNSProblem::setDSBlocks - Not implemented: forbidden operation.");
}

void OneStepNSProblem::updateUnitaryBlocks()
{
  // The present functions checks various conditions and possibly
  // compute unitaryBlocks matrices.
  //
  // Let URi and URj be two Unitary Relations.
  //
  // Things to be checked are:
  //  1 - is the topology time invariant?
  //  2 - does unitaryBlocks[URi][URj] already exists (ie has been
  //  computed in a previous time step)?
  //  3 - do we need to compute this unitaryBlock? A unitaryBlock is
  //  to be computed if URi and URj are in IndexSet1 AND if URi and
  //  URj have common DynamicalSystems.
  //
  // The possible cases are:
  //
  //  - If 1 and 2 are true then it does nothing. 3 is not checked.
  //  - If 1 == true, 2 == false, 3 == false, it does nothing.
  //  - If 1 == true, 2 == false, 3 == true, it computes the
  //    unitaryBlock.
  //  - If 1==false, 2 is not checked, and the unitaryBlock is
  //    computed if 3==true.
  //

  // Get index set from Simulation
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(_levelMin);


  bool isTimeInvariant = simulation()->model()->
                         nonSmoothDynamicalSystem()->topology()->isTimeInvariant();
  bool isLinear = simulation()->model()->nonSmoothDynamicalSystem()->isLinear();

  // we put diagonal informations on vertices
  // self loops with bgl are a *nightmare* at the moment
  // (patch 65198 on standard boost install)

  if (simulation()->indexSet(_levelMin)->properties().symmetric)
  {
    UnitaryRelationsGraph::VIterator vi, viend;
    for (boost::tie(vi, viend) = indexSet->vertices();
         vi != viend; ++vi)
    {

      if (!isTimeInvariant || !isLinear)
      {
        computeDiagonalUnitaryBlock(*vi);
      }
      else
      {
        if (! indexSet->properties(*vi).block)
        {
          computeDiagonalUnitaryBlock(*vi);
        }
      }
    }

    UnitaryRelationsGraph::EIterator ei, eiend;
    for (boost::tie(ei, eiend) = indexSet->edges();
         ei != eiend; ++ei)
    {

      if (!isTimeInvariant || !isLinear)
      {
        computeUnitaryBlock(*ei);

        unsigned int isrc = indexSet->index(indexSet->source(*ei));
        unsigned int itar = indexSet->index(indexSet->target(*ei));

        // allocation for transposed block
        // should be avoided

        if (itar > isrc) // upper block has been computed
        {
          if (!indexSet->properties(*ei).lower_block)
          {
            indexSet->properties(*ei).lower_block.
            reset(new SimpleMatrix(indexSet->properties(*ei).upper_block->size(1),
                                   indexSet->properties(*ei).upper_block->size(0)));
          }
          indexSet->properties(*ei).lower_block->trans(*indexSet->properties(*ei).upper_block);
        }
        else
        {
          assert(itar < isrc);    // lower block has been computed
          if (!indexSet->properties(*ei).upper_block)
          {
            indexSet->properties(*ei).upper_block.
            reset(new SimpleMatrix(indexSet->properties(*ei).lower_block->size(1),
                                   indexSet->properties(*ei).lower_block->size(0)));
          }
          indexSet->properties(*ei).upper_block->trans(*indexSet->properties(*ei).lower_block);
        }
      }
      else // isTimeInvariant or isLinear
      {
        // symmetric case : lower or upper can be tested
        if (! indexSet->properties(*ei).lower_block)
        {
          computeUnitaryBlock(*ei);

          // if lower block has been computed
          if (! indexSet->properties(*ei).upper_block)
          {
            assert(indexSet->properties(*ei).lower_block);
            indexSet->properties(*ei).upper_block.
            reset(new SimpleMatrix(indexSet->properties(*ei).lower_block->size(1),
                                   indexSet->properties(*ei).lower_block->size(0)));
            indexSet->properties(*ei).upper_block->trans(*indexSet->properties(*ei).lower_block);
          }
          else // upper block has been computed
          {
            assert(indexSet->properties(*ei).upper_block);
            indexSet->properties(*ei).lower_block.
            reset(new SimpleMatrix(indexSet->properties(*ei).upper_block->size(1),
                                   indexSet->properties(*ei).upper_block->size(0)));
            indexSet->properties(*ei).lower_block->trans(*indexSet->properties(*ei).upper_block);
          }
        }
      }
    }
  }
  else // not symmetric => follow out_edges for each vertices
  {
    UnitaryRelationsGraph::VIterator vi, viend;
    for (boost::tie(vi, viend) = indexSet->vertices();
         vi != viend; ++vi)
    {

      if (!isTimeInvariant || !isLinear)
      {
        computeDiagonalUnitaryBlock(*vi);
      }
      else
      {
        if (! indexSet->properties(*vi).block)
        {
          computeDiagonalUnitaryBlock(*vi);
        }
      }

      /* on a undirected graph, out_edges gives all incident edges */
      UnitaryRelationsGraph::OEIterator oei, oeiend;
      for (boost::tie(oei, oeiend) = indexSet->out_edges(*vi);
           oei != oeiend; ++oei)
      {

        if (!isTimeInvariant || !isLinear)
        {
          computeUnitaryBlock(*oei);
        }
        else
        {
          /* upper or lower ? */
          unsigned int isrc = indexSet->index(indexSet->source(*oei));
          unsigned int itar = indexSet->index(indexSet->target(*oei));

          assert(isrc != itar);

          if (((itar > isrc) && !indexSet->properties(*oei).upper_block) ||
              ((isrc < itar) && !indexSet->properties(*oei).lower_block))
          {
            computeUnitaryBlock(*oei);
          }
        }
      }
    }
  }
}

void OneStepNSProblem::computeAllUnitaryBlocks()
{
  assert(0);
}

void OneStepNSProblem::updateDSBlocks()
{
  // The present functions checks various conditions and possibly compute DSBlocks matrices.
  //
  // Let URi and URj be two Unitary Relations.
  //
  // Things to be checked are:
  //  1 - is the topology time invariant?
  //  2 - does DSBlocks[DSi] already exists (ie has been computed in a previous time step)?
  //  3 - do we need to compute this DSBlock? A DSBlock is to be computed if DSi is concerned by a UR in IndexSet1
  //
  // The possible cases are:
  //
  //  - If 1 and 2 are true then it does nothing. 3 is not checked.
  //  - If 1 == true, 2 == false, 3 == false, it does nothing.
  //  - If 1 == true, 2 == false, 3 == true, it computes the unitaryBlock.
  //  - If 1==false, 2 is not checked, and the unitaryBlock is computed if 3==true.
  //

  // \warning We decided to include all dynamical systems test 3 is not satisfied


  //   bool isTimeInvariant= simulation()->model()->nonSmoothDynamicalSystem()
  //     ->topology()->isTimeInvariant();
  //   SP::DynamicalSystemsSet allDS = simulation()->model()->nonSmoothDynamicalSystem()->dynamicalSystems();

  //   DSIterator itDS;
  //   for(itDS = allDS->begin(); itDS!=allDS->end();++itDS)
  //     {
  //       if(!isTimeInvariant)
  //  computeDSBlock(*itDS);
  //       else // if(isTimeInvariant)
  //  {
  //    if( (DSBlocks.find(*itDS)) != DSBlocks.end())  // if unitaryBlocks[UR1] exists
  //      {
  //        ; // do nothing
  //      }
  //    else computeDSBlock(*itDS);
  //  }
  //     }

}

void OneStepNSProblem::computeAllDSBlocks()
{
  assert(0);
  //  SP::DynamicalSystemsSet allDS;
  //   DSIterator itDS;
  //   allDS = simulation()->model()->nonSmoothDynamicalSystem()->dynamicalSystems();

  //   for(itDS = allDS->begin(); itDS!=allDS->end();++itDS)
  //     computeDSBlock(*itDS);
}

void OneStepNSProblem::computeDSBlock(SP::DynamicalSystem)
{
  RuntimeException::selfThrow
  ("OneStepNSProblem::computeDSBlock - not yet implemented for problem type ="
  );
}


void OneStepNSProblem::initialize(SP::Simulation sim)
{
  // Link with the simulation that owns this osnsp

  assert(sim && "OneStepNSProblem::initialize(sim), sim is null.");

  _simulation = sim;

  bool isTimeInvariant = simulation()->model()->nonSmoothDynamicalSystem()->topology()->isTimeInvariant();
  bool isLinear = simulation()->model()->nonSmoothDynamicalSystem()->isLinear();

  // === Link to the Interactions of the Non Smooth Dynamical System
  // (through the Simulation) === Warning: this means that all
  // Interactions of the NSProblem are included in the OSNS !!
  _OSNSInteractions = simulation()->model()->nonSmoothDynamicalSystem()->interactions();

  // === Adds this in the simulation set of OneStepNSProblem === First
  // checks the id if required.  An id is required if there is more
  // than one OneStepNSProblem in the simulation

  //    if( !(simulation()->oneStepNSProblems())->empty() && _id ==
  //      DEFAULT_OSNS_NAME)
  //      RuntimeException::selfThrow("OneStepNSProblem::constructor(...). Since
  //      the simulation has several one step non smooth problem, an
  //      id is required for each of them.");

  // Checks that the set of Interactions is not empty - Empty set is
  // not forbidden, then we just display a warning message.
  if (!_OSNSInteractions->isEmpty())
    if (isTimeInvariant || !isLinear) // if time variant it is done in
      // precompute
      updateUnitaryBlocks();

  // The maximum size of the problem (for example, the dim. of M in
  // LCP or Friction problems).  Set to the number of possible scalar
  // constraints declared in the topology.
  if (_maxSize == 0) // if maxSize not set explicitely by user before
    // initialize
    _maxSize = simulation()->model()->
               nonSmoothDynamicalSystem()->topology()->numberOfConstraints();


}

void OneStepNSProblem::saveInMemory()
{
  assert(_OSNSInteractions);
  std::for_each(_OSNSInteractions->begin(), _OSNSInteractions->end(),
                boost::bind(&Interaction::swapInMemory, _1));
}

void OneStepNSProblem::saveNSProblemToXML()
{
  // OUT OF DATE - TO BE REVIEWED

  RuntimeException::selfThrow("OneStepNSProblem::saveNSProblemToXML - Not yet implemented");
}

void OneStepNSProblem::getOSIMaps(SP::UnitaryRelation UR, MapOfDSMatrices& centralUnitaryBlocks)
{
  // === OSI = MOREAU : gets W matrices ===
  // === OSI = LSODAR : gets M matrices of each DS concerned by the UnitaryRelation ===

  SP::OneStepIntegrator Osi;
  OSI::TYPES osiType; // type of the current one step integrator
  Type::Siconos dsType; // type of the current Dynamical System
  DSIterator itDS = UR->interaction()->dynamicalSystemsBegin();
  while (itDS != (UR->interaction()->dynamicalSystemsEnd()))
  {
    Osi = simulation()->integratorOfDS(*itDS); // get OneStepIntegrator of current dynamical system
    osiType = Osi->getType();
    if (osiType == OSI::MOREAU)
    {
      dsType = Type::value(**itDS);
      if (dsType != Type::NewtonEulerDS)
        centralUnitaryBlocks[*itDS] = (boost::static_pointer_cast<Moreau> (Osi))->W(*itDS); // get its W matrix ( pointer link!)
      else
        centralUnitaryBlocks[*itDS] = (boost::static_pointer_cast<NewtonEulerDS> (*itDS))->luM(); // get its W matrix ( pointer link!)
    }
    else if (osiType == OSI::LSODAR) // Warning: LagrangianDS only at the time !!!
    {
      dsType = Type::value(**itDS);
      if (dsType != Type::LagrangianDS && dsType != Type::LagrangianLinearTIDS)
        RuntimeException::selfThrow("OneStepNSProblem::getOSIMaps not yet implemented for Lsodar Integrator with dynamical system of type " + dsType);

      // get lu-factorized mass
      centralUnitaryBlocks[*itDS] =
        (boost::static_pointer_cast<LagrangianDS>(*itDS))->massLU();

    }
    else
      RuntimeException::selfThrow("OneStepNSProblem::getOSIMaps not yet implemented for Integrator of type " + osiType);
    ++itDS;
  }
}

void OneStepNSProblem::printStat()
{
  cout << " CPU time for solving : " << _CPUtime / (double)CLOCKS_PER_SEC << endl;
  cout << " Number of iterations done: " << _nbIter << endl;
}

void OneStepNSProblem::clear()
{

  _DSBlocks.clear();
  if (_OSNSInteractions)
    _OSNSInteractions->clear();
}

OneStepNSProblem::~OneStepNSProblem()
{
  clear();
};


