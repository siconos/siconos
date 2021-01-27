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
#include "MLCP.hpp"
#include "MixedComplementarityConditionNSL.hpp"
#include "EqualityConditionNSL.hpp"
#include "Simulation.hpp"
#include "OSNSMatrix.hpp"
#include "NonSmoothDynamicalSystem.hpp"

// --- Numerics headers ---
#include "NonSmoothDrivers.h"
#include "MLCP_Solvers.h"
#include "SiconosCompat.h"

using namespace RELATION;
// #define DEBUG_NCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

// Constructor from a set of data, use delegated constructor
MLCP::MLCP(int numericsSolverId):
  MLCP(SP::SolverOptions(solver_options_create(numericsSolverId),
                         solver_options_delete)) {}

// Constructor from a set of data
MLCP::MLCP(SP::SolverOptions options):
  LinearOSNS(options)
{
  _numerics_problem.reset(mixedLinearComplementarity_new());

  _numerics_problem->blocksRows = (int*)malloc(MLCP_NB_BLOCKS_MAX * sizeof(int));
  _numerics_problem->blocksIsComp = (int*)malloc(MLCP_NB_BLOCKS_MAX * sizeof(int));
  _numerics_problem->blocksRows[0] = 0;

  // The storage with only one matrix M is chosen
  _numerics_problem->isStorageType1 = 1;
  _numerics_problem->isStorageType2 = 0;
}

void  MLCP::reset()
{
  if(_numerics_problem->blocksRows)
    free(_numerics_problem->blocksRows);
  _numerics_problem->blocksRows = nullptr;
  if(_numerics_problem->blocksIsComp)
    free(_numerics_problem->blocksIsComp);
  _numerics_problem->blocksIsComp = nullptr;
  mlcp_driver_reset(&*_numerics_problem, &*_numerics_solver_options);
  _numerics_solver_options.reset();
}

void MLCP::computeOptions(SP::Interaction inter1, SP::Interaction inter2)
{
  DEBUG_BEGIN("MLCP::computeOptions(SP::Interaction inter1, SP::Interaction inter2)\n");
  // Get dimension of the NonSmoothLaw (ie dim of the interactionBlock)
  unsigned int nslawSize1 = inter1->nonSmoothLaw()->size();
  //  unsigned int nslawSize2 = inter2->nonSmoothLaw()->size();

  unsigned int equalitySize1 =  0;
  //unsigned int equalitySize2 =  0;
  if(Type::value(*(inter1->nonSmoothLaw()))
      == Type::MixedComplementarityConditionNSL)
    equalitySize1 = std::static_pointer_cast<MixedComplementarityConditionNSL>(inter1->nonSmoothLaw())->equalitySize();
  else if(Type::value(*(inter1->nonSmoothLaw()))
          == Type::EqualityConditionNSL)
    equalitySize1 = nslawSize1;

  if(inter1 == inter2)
  {
    //inter1->getExtraInteractionBlock(currentInteractionBlock);
    _m += nslawSize1 - equalitySize1;
    _n += equalitySize1;
    if(_curBlock > MLCP_NB_BLOCKS_MAX - 2)
      printf("MLCP.cpp : number of block to small, memory crach below!!!\n");
    /*add an equality block.*/
    if(equalitySize1 > 0)
    {
      DEBUG_PRINT("add an equality block.\n");
      _numerics_problem->blocksRows[_curBlock + 1] = _numerics_problem->blocksRows[_curBlock] + equalitySize1;
      _numerics_problem->blocksIsComp[_curBlock] = 0;
      _curBlock++;
    }
    /*add a complementarity block.*/
    if(nslawSize1 - equalitySize1 > 0)
    {
      DEBUG_PRINT("add a complementarity block.\n");
      _numerics_problem->blocksRows[_curBlock + 1] = _numerics_problem->blocksRows[_curBlock] + nslawSize1 - equalitySize1;
      _numerics_problem->blocksIsComp[_curBlock] = 1;
      _curBlock++;
    }
  }
  DEBUG_END("MLCP::computeOptions(SP::Interaction inter1, SP::Interaction inter2)\n");
}

bool MLCP::checkCompatibleNSLaw(NonSmoothLaw& nslaw)
{
  float type_number= (float) (Type::value(nslaw));
  _nslawtype.insert(type_number);

  if (not (Type::value(nslaw) == Type::MixedComplementarityConditionNSL ||
           Type::value(nslaw) == Type::ComplementarityConditionNSL ||
           Type::value(nslaw) == Type::NewtonImpactNSL ||
           Type::value(nslaw) == Type::EqualityConditionNSL)
    )
  {
    THROW_EXCEPTION("\nMLCP::checkCompatibleNSLaw -  \n\
                      The chosen nonsmooth law is not compatible with LCP one step nonsmooth problem. \n \
                      Compatible NonSmoothLaw are: ComplementarityConditionNSL, EqualityConditionNSL or NewtonImpactNSL\n");
    return false;
  }

  return true;
}



void MLCP::computeInteractionBlock(const InteractionsGraph::EDescriptor& ed)
{
  DEBUG_BEGIN("MLCP::computeInteractionBlock(const InteractionsGraph::EDescriptor& ed)\n")
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  SP::Interaction inter1 = indexSet->bundle(indexSet->source(ed));
  SP::Interaction inter2 = indexSet->bundle(indexSet->target(ed));

  assert(inter1 != inter2);
  bool isLinear = simulation()->nonSmoothDynamicalSystem()->isLinear();

  if(!_hasBeenUpdated || !isLinear)
    LinearOSNS::computeInteractionBlock(ed);

  DEBUG_END("MLCP::computeInteractionBlock(const InteractionsGraph::EDescriptor& ed)\n")
}

void MLCP::computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor& vd)
{

  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  SP::DynamicalSystem DS1 = indexSet->properties(vd).source;
  SP::DynamicalSystem DS2 = indexSet->properties(vd).target;
  SP::Interaction inter = indexSet->bundle(vd);

  // commonDS here...
  if(!_hasBeenUpdated)
    computeOptions(inter, inter);
  LinearOSNS::computeDiagonalInteractionBlock(vd);
}

bool MLCP::preCompute(double time)
{
  bool res = LinearOSNS::preCompute(time);
  _numerics_problem->n = _n;
  _numerics_problem->m = _m;
  return res;
}

int MLCP::compute(double time)
{
  DEBUG_BEGIN("MLCP::compute(double time)\n");
  int info = 0;
  // --- Prepare data for MLCP computing ---
  bool cont = preCompute(time);
  if(!cont)
    return info;
  // cf GenericMechanical for the explanation of this line commented
  // _hasBeenUpdated=true;
  DEBUG_PRINTF("MLCP::compute m n :%d,%d\n", _n, _m);

  /*If user has not allocted the working memory, do it. */
  mlcp_driver_init(&*_numerics_problem, &*_numerics_solver_options);
  // After the first call the mlcp_direct_init must not reset the previous guess
  // But as the problem may change the MLCP update flag is raised
  _numerics_solver_options->iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED] = 1;

  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)

  if(_sizeOutput != 0)
  {
    _numerics_problem->q = _q->getArray();

    // Call MLCP Driver
    DEBUG_PRINT("MLCP display");
    //printf("n %d m %d",n,m);
    //displayNM(_numerics_problem->M);
    //      exit(1);
    //mlcpDefaultSolver *pSolver = new mlcpDefaultSolver(m,n);
    DEBUG_EXPR(display(););

    try
    {
      info = mlcp_driver(&*_numerics_problem, _z->getArray(), _w->getArray(),
                         &*_numerics_solver_options);
    }
    catch(...)
    {
      std::cout << "exception caught" <<std::endl;
      info = 1;
    }

    // --- Recovering of the desired variables from MLCP output ---
    if(!info)
      postCompute();
    else
      printf("[kernel] MLCP::compute -- MLCP solver failed\n");

  }
  else
  {
    DEBUG_PRINT("MLCP::compute : sizeoutput is null\n");
  }
  DEBUG_END("MLCP::compute(double time)\n");
  return info;
}

void MLCP::display() const
{
  std::cout << "======= MLCP of size " << _sizeOutput << " with: " <<std::endl;
  std::cout << " m (number of inequality constraints)" << _m <<std::endl;
  std::cout << " n (number of equality constraints)  " << _n <<std::endl;
  LinearOSNS::display();
}

void MLCP::initialize(SP::Simulation sim)
{
  // General initialize for LinearOSNS
  LinearOSNS::initialize(sim);

  _numerics_problem->M = &*_M->numericsMatrix();
}
void  MLCP::updateInteractionBlocks()
{
  if(!_hasBeenUpdated)
  {
    _curBlock = 0;
    _m = 0;
    _n = 0;
  }
  LinearOSNS::updateInteractionBlocks();
}
