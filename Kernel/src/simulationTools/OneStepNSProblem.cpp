/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
//#include "Interaction.hpp"
#include "Interaction.hpp"
#include "Topology.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "Moreau.hpp"
#include "LagrangianDS.hpp"
#include "NewtonEulerDS.hpp"
#include "ZeroOrderHold.hpp"
//#define OSNS_DEBUG
using namespace std;
OneStepNSProblem::OneStepNSProblem():
  _levelMin(0), _levelMax(0), _maxSize(0), _CPUtime(0), _nbIter(0), _hasBeenUpdated(false)
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
  _onestepnspbxml(osnspbxml), _levelMin(0), _levelMax(0), _maxSize(0), _CPUtime(0), _nbIter(0), _hasBeenUpdated(false)
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
  _onestepnspbxml(osnspbxml), _levelMin(0), _levelMax(0), _maxSize(0), _CPUtime(0), _nbIter(0), _hasBeenUpdated(false)
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
  _numerics_solver_id(newNumericsSolverId),/*_nspbType(pbType),*/ _id(newId), _sizeOutput(0), _levelMin(0), _levelMax(0), _maxSize(0), _CPUtime(0), _nbIter(0), _hasBeenUpdated(false)
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
  _numerics_solver_id(newNumericsSolverId), _id(newId), _sizeOutput(0), _levelMin(0), _levelMax(0), _maxSize(0), _CPUtime(0), _nbIter(0), _hasBeenUpdated(false)
{

  // Numerics general options
  _numerics_options.reset(new NumericsOptions());
  _numerics_options->verboseMode = 0; // turn verbose mode to off by default

  _numerics_solver_options.reset(new SolverOptions);
  _numerics_solver_options->iWork = NULL;
  _numerics_solver_options->dWork = NULL;

}
OneStepNSProblem::OneStepNSProblem(const int newNumericsSolverId):
  _numerics_solver_id(newNumericsSolverId), _sizeOutput(0), _levelMin(0), _levelMax(0), _maxSize(0), _CPUtime(0), _nbIter(0), _hasBeenUpdated(false)
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

  ConstMatIterator itDS = _DSBlocks.find(DS1->number());
  return itDS->second;

}


void OneStepNSProblem::setDSBlocks(const MapOfDSMatrices& newMap)
{
  RuntimeException::selfThrow("OneStepNSProblem::setDSBlocks - Not implemented: forbidden operation.");
}

void OneStepNSProblem::updateInteractionBlocks()
{
  // The present functions checks various conditions and possibly
  // compute interactionBlocks matrices.
  //
  // Let interi and interj be two Interactions.
  //
  // Things to be checked are:
  //  1 - is the topology time invariant?
  //  2 - does interactionBlocks[interi][interj] already exists (ie has been
  //  computed in a previous time step)?
  //  3 - do we need to compute this interactionBlock? A interactionBlock is
  //  to be computed if interi and interj are in IndexSet1 AND if interi and
  //  interj have common DynamicalSystems.
  //
  // The possible cases are:
  //
  //  - If 1 and 2 are true then it does nothing. 3 is not checked.
  //  - If 1 == true, 2 == false, 3 == false, it does nothing.
  //  - If 1 == true, 2 == false, 3 == true, it computes the
  //    interactionBlock.
  //  - If 1==false, 2 is not checked, and the interactionBlock is
  //    computed if 3==true.
  //

  // Get index set from Simulation
  SP::InteractionsGraph indexSet = simulation()->indexSet(_levelMin);


  bool isLinear = simulation()->model()->nonSmoothDynamicalSystem()->isLinear();

  // we put diagonal informations on vertices
  // self loops with bgl are a *nightmare* at the moment
  // (patch 65198 on standard boost install)

  if (indexSet->properties().symmetric)
  {
    InteractionsGraph::VIterator vi, viend;
    for (boost::tie(vi, viend) = indexSet->vertices();
         vi != viend; ++vi)
    {
      SP::Interaction inter = indexSet->bundle(*vi);
      unsigned int nslawSize = inter->getNonSmoothLawSize();
      if (! indexSet->properties(*vi).block)
      {
        indexSet->properties(*vi).block.reset(new SimpleMatrix(nslawSize, nslawSize));
      }

      if (!isLinear || !_hasBeenUpdated)
      {
        computeDiagonalInteractionBlock(*vi);
      }
    }

    /* interactionBlock must be zeroed at init */
    std::vector<bool> initialized;
    initialized.resize(indexSet->edges_number());
    std::fill(initialized.begin(), initialized.end(), false);

    InteractionsGraph::EIterator ei, eiend;
    for (boost::tie(ei, eiend) = indexSet->edges();
         ei != eiend; ++ei)
    {
      SP::Interaction inter1 = indexSet->bundle(indexSet->source(*ei));
      SP::Interaction inter2 = indexSet->bundle(indexSet->target(*ei));

      /* on adjoint graph there is at most 2 edges between source and target */
      InteractionsGraph::EDescriptor ed1, ed2;
      boost::tie(ed1, ed2) = indexSet->edges(indexSet->source(*ei), indexSet->target(*ei));

      assert(*ei == ed1 || *ei == ed2);

      /* the first edge as the lower index */
      assert(indexSet->index(ed1) <= indexSet->index(ed2));

      // Memory allocation if needed
      unsigned int nslawSize1 = inter1->getNonSmoothLawSize();
      unsigned int nslawSize2 = inter2->getNonSmoothLawSize();
      unsigned int isrc = indexSet->index(indexSet->source(*ei));
      unsigned int itar = indexSet->index(indexSet->target(*ei));

      SP::SiconosMatrix currentInteractionBlock;

      if (itar > isrc) // upper block
      {
        if (! indexSet->properties(ed1).upper_block)
        {
          indexSet->properties(ed1).upper_block.reset(new SimpleMatrix(nslawSize1, nslawSize2));
          if (ed2 != ed1)
            indexSet->properties(ed2).upper_block = indexSet->properties(ed1).upper_block;
        }
        currentInteractionBlock = indexSet->properties(ed1).upper_block;
      }
      else  // lower block
      {
        if (! indexSet->properties(ed1).lower_block)
        {
          indexSet->properties(ed1).lower_block.reset(new SimpleMatrix(nslawSize1, nslawSize2));
          if (ed2 != ed1)
            indexSet->properties(ed2).lower_block = indexSet->properties(ed1).lower_block;
        }
        currentInteractionBlock = indexSet->properties(ed1).lower_block;
      }

      if (!initialized[indexSet->index(ed1)])
      {
        initialized[indexSet->index(ed1)] = true;
        currentInteractionBlock->zero();
      }

      if (!isLinear || !_hasBeenUpdated)
      {
        computeInteractionBlock(*ei);


        // allocation for transposed block
        // should be avoided

        if (itar > isrc) // upper block has been computed
        {
          if (!indexSet->properties(ed1).lower_block)
          {
            indexSet->properties(ed1).lower_block.
            reset(new SimpleMatrix(indexSet->properties(ed1).upper_block->size(1),
                                   indexSet->properties(ed1).upper_block->size(0)));
          }
          indexSet->properties(ed1).lower_block->trans(*indexSet->properties(ed1).upper_block);
          indexSet->properties(ed2).lower_block = indexSet->properties(ed1).lower_block;
        }
        else
        {
          assert(itar < isrc);    // lower block has been computed
          if (!indexSet->properties(ed1).upper_block)
          {
            indexSet->properties(ed1).upper_block.
            reset(new SimpleMatrix(indexSet->properties(ed1).lower_block->size(1),
                                   indexSet->properties(ed1).lower_block->size(0)));
          }
          indexSet->properties(ed1).upper_block->trans(*indexSet->properties(ed1).lower_block);
          indexSet->properties(ed2).upper_block = indexSet->properties(ed1).upper_block;
        }
      }
    }
  }
  else // not symmetric => follow out_edges for each vertices
  {

    InteractionsGraph::VIterator vi, viend;
    for (boost::tie(vi, viend) = indexSet->vertices();
         vi != viend; ++vi)
    {
      SP::Interaction inter = indexSet->bundle(*vi);
      unsigned int nslawSize = inter->getNonSmoothLawSize();
      if (! indexSet->properties(*vi).block)
      {
        indexSet->properties(*vi).block.reset(new SimpleMatrix(nslawSize, nslawSize));
      }

      if (!isLinear || !_hasBeenUpdated)
      {
        computeDiagonalInteractionBlock(*vi);
      }

      /* interactionBlock must be zeroed at init */
      std::vector<bool> initialized;
      initialized.resize(indexSet->edges_number());
      std::fill(initialized.begin(), initialized.end(), false);

      /* on a undirected graph, out_edges gives all incident edges */
      InteractionsGraph::OEIterator oei, oeiend;
      for (boost::tie(oei, oeiend) = indexSet->out_edges(*vi);
           oei != oeiend; ++oei)
      {

        /* on adjoint graph there is at most 2 edges between source and target */
        InteractionsGraph::EDescriptor ed1, ed2;
        boost::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));

        assert(*oei == ed1 || *oei == ed2);

        /* the first edge as the lower index */
        assert(indexSet->index(ed1) == indexSet->index(ed2));

        SP::Interaction inter1 = indexSet->bundle(indexSet->source(*oei));
        SP::Interaction inter2 = indexSet->bundle(indexSet->target(*oei));

        // Memory allocation if needed
        unsigned int nslawSize1 = inter1->getNonSmoothLawSize();
        unsigned int nslawSize2 = inter2->getNonSmoothLawSize();
        unsigned int isrc = indexSet->index(indexSet->source(*oei));
        unsigned int itar = indexSet->index(indexSet->target(*oei));

        SP::SiconosMatrix currentInteractionBlock;

        if (itar > isrc) // upper block
        {
          if (! indexSet->properties(ed1).upper_block)
          {
            indexSet->properties(ed1).upper_block.reset(new SimpleMatrix(nslawSize1, nslawSize2));
            if (ed2 != ed1)
              indexSet->properties(ed2).upper_block = indexSet->properties(ed1).upper_block;
          }
          currentInteractionBlock = indexSet->properties(ed1).upper_block;

        }
        else  // lower block
        {
          if (! indexSet->properties(ed1).lower_block)
          {
            indexSet->properties(ed1).lower_block.reset(new SimpleMatrix(nslawSize1, nslawSize2));
            if (ed2 != ed1)
              indexSet->properties(ed2).lower_block = indexSet->properties(ed1).lower_block;
          }
          currentInteractionBlock = indexSet->properties(ed1).lower_block;
        }


        if (!initialized[indexSet->index(ed1)])
        {
          initialized[indexSet->index(ed1)] = true;
          currentInteractionBlock->zero();
        }


        if (!isLinear || !_hasBeenUpdated)
        {
          if (isrc != itar)
            computeInteractionBlock(*oei);
        }

      }
    }
  }

#ifdef OSNS_DEBUG
  displayBlocks(indexSet);
#endif


}

void OneStepNSProblem::displayBlocks(SP::InteractionsGraph indexSet)
{

  std::cout <<  "OneStepNSProblem::displayBlocks(SP::InteractionsGraph indexSet) " << std::endl;
  InteractionsGraph::VIterator vi, viend;
  for (boost::tie(vi, viend) = indexSet->vertices();
       vi != viend; ++vi)
  {
    SP::Interaction inter = indexSet->bundle(*vi);
    if (indexSet->properties(*vi).block)
    {
      indexSet->properties(*vi).block->display();
    }

    InteractionsGraph::OEIterator oei, oeiend;
    for (boost::tie(oei, oeiend) = indexSet->out_edges(*vi);
         oei != oeiend; ++oei)
    {
      InteractionsGraph::EDescriptor ed1, ed2;
      boost::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));

      if (indexSet->properties(ed1).upper_block)
      {
        indexSet->properties(ed1).upper_block->display();
      }
      if (indexSet->properties(ed1).lower_block)
      {
        indexSet->properties(ed1).lower_block->display();
      }
      if (indexSet->properties(ed2).upper_block)
      {
        indexSet->properties(ed2).upper_block->display();
      }
      if (indexSet->properties(ed2).lower_block)
      {
        indexSet->properties(ed2).lower_block->display();
      }
    }

  }
}

void OneStepNSProblem::computeAllInteractionBlocks()
{
  assert(0);
}

void OneStepNSProblem::updateDSBlocks()
{
  // The present functions checks various conditions and possibly compute DSBlocks matrices.
  //
  // Let interi and interj be two Interactions.
  //
  // Things to be checked are:
  //  1 - is the topology time invariant?
  //  2 - does DSBlocks[DSi] already exists (ie has been computed in a previous time step)?
  //  3 - do we need to compute this DSBlock? A DSBlock is to be computed if DSi is concerned by a Interactionin IndexSet1
  //
  // The possible cases are:
  //
  //  - If 1 and 2 are true then it does nothing. 3 is not checked.
  //  - If 1 == true, 2 == false, 3 == false, it does nothing.
  //  - If 1 == true, 2 == false, 3 == true, it computes the interactionBlock.
  //  - If 1==false, 2 is not checked, and the interactionBlock is computed if 3==true.
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
  //    if( (DSBlocks.find(*itDS)) != DSBlocks.end())  // if interactionBlocks[inter1] exists
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
void OneStepNSProblem::saveTimeStepInMemory()
{
  assert(_OSNSInteractions);
  std::for_each(_OSNSInteractions->begin(), _OSNSInteractions->end(),
                boost::bind(&Interaction::swapTimeStepInMemory, _1));
}

void OneStepNSProblem::saveNSProblemToXML()
{
  // OUT OF DATE - TO BE REVIEWED

  RuntimeException::selfThrow("OneStepNSProblem::saveNSProblemToXML - Not yet implemented");
}

void OneStepNSProblem::getOSIMaps(SP::Interaction inter, MapOfDSMatrices& centralInteractionBlocks)
{
  // === OSI = MOREAU : gets W matrices ===
  // === OSI = LSODAR : gets M matrices of each DS concerned by the Interaction ===

  SP::OneStepIntegrator Osi;
  OSI::TYPES osiType; // type of the current one step integrator
  Type::Siconos dsType; // type of the current Dynamical System
  DSIterator itDS = inter->dynamicalSystemsBegin();
  while (itDS != (inter->dynamicalSystemsEnd()))
  {
    Osi = simulation()->integratorOfDS(*itDS); // get OneStepIntegrator of current dynamical system
    osiType = Osi->getType();
    unsigned int itN = (*itDS)->number();
    if (osiType == OSI::MOREAU || osiType == OSI::SCHATZMANPAOLI)
    {
      dsType = Type::value(**itDS);
      if (dsType != Type::NewtonEulerDS)
        centralInteractionBlocks[itN] = (boost::static_pointer_cast<Moreau> (Osi))->W(*itDS); // get its W matrix ( pointer link!)
      else
        centralInteractionBlocks[itN] = (boost::static_pointer_cast<NewtonEulerDS> (*itDS))->luW(); // get its W matrix ( pointer link!)
    }
    else if (osiType == OSI::LSODAR) // Warning: LagrangianDS only at the time !!!
    {
      dsType = Type::value(**itDS);
      if (dsType != Type::LagrangianDS && dsType != Type::LagrangianLinearTIDS)
        RuntimeException::selfThrow("OneStepNSProblem::getOSIMaps not yet implemented for Lsodar Integrator with dynamical system of type " + dsType);

      // get lu-factorized mass
      centralInteractionBlocks[itN] =
        (boost::static_pointer_cast<LagrangianDS>(*itDS))->massLU();

    }
    else if (osiType == OSI::D1MINUSLINEAR)
    {
      dsType = Type::value(**itDS);
      if (dsType != Type::LagrangianDS && dsType != Type::LagrangianLinearTIDS)
        RuntimeException::selfThrow("OneStepNSProblem::getOSIMaps not yet implemented for D1MinusLinear integrator with dynamical system of type " + dsType);

      centralInteractionBlocks[itN].reset(new SimpleMatrix(*((boost::static_pointer_cast<LagrangianDS>(*itDS))->mass())));
    }
    else if (osiType == OSI::ZOH)
    {
      centralInteractionBlocks[itN] = (boost::static_pointer_cast<ZeroOrderHold>(Osi))->Phi(**itDS);
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


