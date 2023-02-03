/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "Interaction.hpp"
#include "NumericsSolversNamespace.h"  // solver_options stuff
#include "OSNSMatrix.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include "MixedComplementarityConditionNSL.hpp"
#include "TypeName.hpp"
#include "SiconosException.hpp"
// #define DEBUG_NCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

// Constructor from a set of data, use delegated constructor
siconos::nonsmooth_formulations::MLCP::MLCP(int numericsSolverId)
    : MLCP(std::shared_ptr<SolverOptions>(
          solver_options_create(numericsSolverId),
          solver_options_delete))
{
}

// Constructor from a set of data
siconos::nonsmooth_formulations::MLCP::MLCP(std::shared_ptr<SolverOptions> options)
    : LinearOSNS(options)
{
  _numerics_problem = std::make_shared<MixedLinearComplementarityProblem>();

  _numerics_problem->blocksRows = (int*)malloc(MLCP_NB_BLOCKS_MAX * sizeof(int));
  _numerics_problem->blocksIsComp = (int*)malloc(MLCP_NB_BLOCKS_MAX * sizeof(int));
  _numerics_problem->blocksRows[0] = 0;

  // The storage with only one matrix M is chosen
  _numerics_problem->isStorageType1 = 1;
  _numerics_problem->isStorageType2 = 0;
}

siconos::nonsmooth_formulations::MLCP::~MLCP() noexcept
{
  if (_numerics_problem->blocksRows) free(_numerics_problem->blocksRows);
  _numerics_problem->blocksRows = nullptr;
  if (_numerics_problem->blocksIsComp) free(_numerics_problem->blocksIsComp);
  _numerics_problem->blocksIsComp = nullptr;
  mlcp_driver_reset(&*_numerics_problem, &*_numerics_solver_options);
  _numerics_solver_options.reset();
}

void siconos::nonsmooth_formulations::MLCP::computeOptions(
    std::shared_ptr<siconos::modeling::Interaction> inter1,
    std::shared_ptr<siconos::modeling::Interaction> inter2)
{
  DEBUG_BEGIN(
      "siconos::nonsmooth_formulations::MLCP::computeOptions(std::shared_ptr<siconos::modeling::"
      "Interaction> inter1, "
      "std::shared_ptr<siconos::modeling::Interaction> inter2)\n");
  // Get dimension of the siconos::modeling::NonSmoothLaw (ie dim of the interactionBlock)
  auto nslawSize1 = inter1->nonSmoothLaw()->size();
  //  unsigned int nslawSize2 = inter2->nonSmoothLaw()->size();

  auto equalitySize1 = 0;
  // unsigned int equalitySize2 =  0;
  if (siconos::types::type_value(*(inter1->nonSmoothLaw())) ==
      siconos::modeling::Type::MixedComplementarityConditionNSL)
    equalitySize1 =
        std::static_pointer_cast<siconos::modeling::MixedComplementarityConditionNSL>(inter1->nonSmoothLaw())
            ->equalitySize();
  else if (siconos::types::type_value(*(inter1->nonSmoothLaw())) ==
           siconos::modeling::Type::EqualityConditionNSL)
    equalitySize1 = nslawSize1;

  if (inter1 == inter2) {
    // inter1->getExtraInteractionBlock(currentInteractionBlock);
    _m += nslawSize1 - equalitySize1;
    _n += equalitySize1;
    if (_curBlock > MLCP_NB_BLOCKS_MAX - 2)
      printf("MLCP.cpp : number of block to small, memory crach below!!!\n");
    /*add an equality block.*/
    if (equalitySize1 > 0) {
      DEBUG_PRINT("add an equality block.\n");
      _numerics_problem->blocksRows[_curBlock + 1] =
          _numerics_problem->blocksRows[_curBlock] + equalitySize1;
      _numerics_problem->blocksIsComp[_curBlock] = 0;
      _curBlock++;
    }
    /*add a complementarity block.*/
    if (nslawSize1 - equalitySize1 > 0) {
      DEBUG_PRINT("add a complementarity block.\n");
      _numerics_problem->blocksRows[_curBlock + 1] =
          _numerics_problem->blocksRows[_curBlock] + nslawSize1 - equalitySize1;
      _numerics_problem->blocksIsComp[_curBlock] = 1;
      _curBlock++;
    }
  }
  DEBUG_END(
      "siconos::nonsmooth_formulations::MLCP::computeOptions(std::shared_ptr<siconos::modeling::"
      "Interaction> inter1, "
      "std::shared_ptr<siconos::modeling::Interaction> inter2)\n");
}

bool siconos::nonsmooth_formulations::MLCP::checkCompatibleNSLaw(siconos::modeling::NonSmoothLaw& nslaw)
{
  float type_number = static_cast<float>(siconos::types::type_value(nslaw));
  _nslawtype.insert(type_number);

  if (not(siconos::types::type_value(nslaw) ==
              siconos::modeling::Type::MixedComplementarityConditionNSL ||
          siconos::types::type_value(nslaw) ==
              siconos::modeling::Type::ComplementarityConditionNSL ||
          siconos::types::type_value(nslaw) == siconos::modeling::Type::NewtonImpactNSL ||
          siconos::types::type_value(nslaw) ==
              siconos::modeling::Type::EqualityConditionNSL)) {
    THROW_EXCEPTION(
        "\nsiconos::nonsmooth_formulations::MLCP::checkCompatibleNSLaw -  \n\
                      The chosen nonsmooth law is not compatible with LCP one step nonsmooth problem. \n \
                      Compatible siconos::modeling::NonSmoothLaw are: ComplementarityConditionNSL, EqualityConditionNSL or NewtonImpactNSL\n");
    return false;
  }

  return true;
}

void siconos::nonsmooth_formulations::MLCP::computeInteractionBlock(
    const siconos::graphs::InteractionsGraph::EDescriptor& ed)
{
  DEBUG_BEGIN(
      "siconos::nonsmooth_formulations::MLCP::computeInteractionBlock(const "
      "siconos::graphs::InteractionsGraph::EDescriptor& "
      "ed)\n")
  auto indexSet = simulation()->indexSet(indexSetLevel());
  auto inter1 = indexSet->bundle(indexSet->source(ed));
  auto inter2 = indexSet->bundle(indexSet->target(ed));

  assert(inter1 != inter2);
  bool isLinear = simulation()->nonSmoothDynamicalSystem()->isLinear();

  if (!_hasBeenUpdated || !isLinear) LinearOSNS::computeInteractionBlock(ed);

  DEBUG_END(
      "siconos::nonsmooth_formulations::MLCP::computeInteractionBlock(const "
      "siconos::graphs::InteractionsGraph::EDescriptor& "
      "ed)\n")
}

void siconos::nonsmooth_formulations::MLCP::computeDiagonalInteractionBlock(
    const siconos::graphs::InteractionsGraph::VDescriptor& vd)
{
  auto indexSet = simulation()->indexSet(indexSetLevel());
  auto DS1 = indexSet->properties(vd).source;
  auto DS2 = indexSet->properties(vd).target;
  auto inter = indexSet->bundle(vd);

  // commonDS here...
  if (!_hasBeenUpdated) computeOptions(inter, inter);
  LinearOSNS::computeDiagonalInteractionBlock(vd);
}

int siconos::nonsmooth_formulations::MLCP::solve()
{
  // Note FP : wrap call to numerics solver inside this function
  // for python API (e.g. to allow profiling without C struct handling)

  _numerics_problem->q = _q->getArray();
  _numerics_problem->M = &*_M->numericsMatrix();
  _numerics_problem->n = _n;
  _numerics_problem->m = _m;

  // After the first call the mlcp_direct_init must not reset the previous guess
  // But as the problem may change the MLCP update flag is raised
  _numerics_solver_options->iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED] = 1;

  /*If user has not allocted the working memory, do it. */
  mlcp_driver_init(&*_numerics_problem, &*_numerics_solver_options);

  DEBUG_PRINT("MLCP display");
  // printf("n %d m %d",n,m);
  // displayNM(_numerics_problem->M);
  //       exit(1);
  // mlcpDefaultSolver *pSolver = new mlcpDefaultSolver(m,n);
  DEBUG_EXPR(display(););

  // Call MLCP Driver
  int info = mlcp_driver(&*_numerics_problem, _z->getArray(),
                                            _w->getArray(), &*_numerics_solver_options);

  return info;
}

int siconos::nonsmooth_formulations::MLCP::compute(double time)
{
  DEBUG_BEGIN("siconos::nonsmooth_formulations::MLCP::compute(double time)\n");
  int info = 0;
  // --- Prepare data for MLCP computing ---
  bool not_empty = preCompute(time);
  if (!not_empty) return info;
  // cf GenericMechanical for the explanation of this line commented
  // _hasBeenUpdated=true;
  DEBUG_PRINTF("siconos::nonsmooth_formulations::MLCP::compute m n :%d,%d\n", _n, _m);

  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)

  if (_sizeOutput != 0) {
    info = solve();
    // --- Recovering of the desired variables from MLCP output ---
    if (!info)
      postCompute();
    else
      printf("[kernel] siconos::nonsmooth_formulations::MLCP::compute -- MLCP solver failed\n");
  }
  else {
    DEBUG_PRINT("siconos::nonsmooth_formulations::MLCP::compute : sizeoutput is null\n");
  }
  DEBUG_END("siconos::nonsmooth_formulations::MLCP::compute(double time)\n");
  return info;
}

void siconos::nonsmooth_formulations::MLCP::display() const
{
  std::cout << "======= MLCP of size " << _sizeOutput << " with: \n";
  std::cout << " m (number of inequality constraints)" << _m << std::endl;
  std::cout << " n (number of equality constraints)  " << _n << std::endl;
  LinearOSNS::display();
}

void siconos::nonsmooth_formulations::MLCP::updateInteractionBlocks()
{
  if (!_hasBeenUpdated) {
    _curBlock = 0;
    _m = 0;
    _n = 0;
  }
  LinearOSNS::updateInteractionBlocks();
}
