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
#include "GenericMechanical.hpp"
#include "Topology.hpp"
#include "Simulation.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "RelayNSL.hpp"
#include "OSNSMatrix.hpp"
#include "GenericMechanicalProblem.h" // from numerics, for GM problem struct
#include "FrictionContactProblem.h" // from numerics, for GM problem struct
#include "RelayProblem.h" // from numerics, for GM problem struct
#include "GenericMechanical_Solvers.h"


#include "SolverOptions.h"
#include "numerics_verbose.h"
using namespace RELATION;
// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #define DEBUG_WHERE_MESSAGES
#include "GenericMechanical_cst.h"          // for SICONOS_GENERIC_MECHANICA...
#include <debug.h>


GenericMechanical::GenericMechanical(int FC3D_Solver_Id):
  GenericMechanical(
    SP::SolverOptions(solver_options_create(SICONOS_GENERIC_MECHANICAL_NSGS),
                      solver_options_delete))
{
  solver_options_update_internal(_numerics_solver_options.get(), 1, FC3D_Solver_Id);
}


GenericMechanical::GenericMechanical(SP::SolverOptions options):
  LinearOSNS(options)
{
  DEBUG_BEGIN("GenericMechanical::GenericMechanical(SP::SolverOptions options)\n");
  //assert(options->solverId == SICONOS_GENERIC_MECHANICAL_NSGS); this will be checked in the driver
  _numericsMatrixStorageType = NM_SPARSE_BLOCK;
  _pnumerics_GMP = genericMechanicalProblem_new();
  DEBUG_END("GenericMechanical::GenericMechanical(SP::SolverOptions options)\n");
}


void GenericMechanical::initialize(SP::Simulation sim)
{
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  LinearOSNS::initialize(sim);
}
bool GenericMechanical::checkCompatibleNSLaw(NonSmoothLaw& nslaw)
{
  // do nothoing since it is check in GenericMechanical::computeDiagonalInteractionBlock
  return true;
}

void GenericMechanical::computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor& vd)
{
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  //bool isTimeInvariant = simulation()->nonSmoothDynamicalSystem()->topology()->isTimeInvariant();

  /*Build the corresponding numerics problems*/

  SP::DynamicalSystem DS1 = indexSet->properties(vd).source;
  SP::DynamicalSystem DS2 = indexSet->properties(vd).target;
  SP::Interaction inter = indexSet->bundle(vd);

  DEBUG_PRINT("GenericMechanical::computeInteractionBlock: add problem of type ");

  if(!_hasBeenUpdated)
  {
    int size = inter->nonSmoothLaw()->size();
    if(Type::value(*(inter->nonSmoothLaw()))
        == Type::EqualityConditionNSL)
    {
      gmp_add(_pnumerics_GMP, SICONOS_NUMERICS_PROBLEM_EQUALITY, size);
      DEBUG_PRINT("Type::EqualityConditionNSL\n");
      //pAux->size= inter->nonSmoothLaw()->size();
    }
    else if(Type::value(*(inter->nonSmoothLaw()))
            == Type::NewtonImpactNSL)
    {
      gmp_add(_pnumerics_GMP, SICONOS_NUMERICS_PROBLEM_LCP, size);
      DEBUG_PRINT(" Type::NewtonImpactNSL\n");
    }
    else if(Type::value(*(inter->nonSmoothLaw()))
            == Type::RelayNSL)
    {
      RelayProblem * pAux =
        (RelayProblem *)gmp_add(_pnumerics_GMP, SICONOS_NUMERICS_PROBLEM_RELAY, size);
      SP::RelayNSL nsLaw =
        std::static_pointer_cast<RelayNSL> (inter->nonSmoothLaw());
      for(int i=0; i<size; i++)
      {
        pAux->lb[i] = nsLaw->lb();
        pAux->ub[i] = nsLaw->ub();
      }
      DEBUG_PRINT(" Type::RelayNSL\n");
    }
    else if(Type::value(*(inter->nonSmoothLaw()))
            == Type::NewtonImpactFrictionNSL)
    {
      if (size ==3)
      {
        FrictionContactProblem * pAux =
          (FrictionContactProblem *)gmp_add(_pnumerics_GMP, SICONOS_NUMERICS_PROBLEM_FC3D, size);
        SP::NewtonImpactFrictionNSL nsLaw =
          std::static_pointer_cast<NewtonImpactFrictionNSL> (inter->nonSmoothLaw());
        pAux->dimension = 3;
        pAux->numberOfContacts = 1;
        *(pAux->mu) = nsLaw->mu();

        DEBUG_PRINT(" Type::NewtonImpactFrictionNSL\n");
      }
      else if (size ==2)
      {
        FrictionContactProblem * pAux =
          (FrictionContactProblem *)gmp_add(_pnumerics_GMP, SICONOS_NUMERICS_PROBLEM_FC2D, size);
        SP::NewtonImpactFrictionNSL nsLaw =
          std::static_pointer_cast<NewtonImpactFrictionNSL> (inter->nonSmoothLaw());
        pAux->dimension = 2;
        pAux->numberOfContacts = 1;
        *(pAux->mu) = nsLaw->mu();
        
        DEBUG_PRINT(" Type::NewtonImpactFrictionNSL\n");
      }
        
    }
    else
    {
      THROW_EXCEPTION("GenericMechanical::computeDiagonalInteractionBlock- not yet implemented for that NSLAW type");
    }
  }

  LinearOSNS::computeDiagonalInteractionBlock(vd);
}

void GenericMechanical::computeInteractionBlock(const InteractionsGraph::EDescriptor& ed)
{
  LinearOSNS::computeInteractionBlock(ed);
}

int GenericMechanical::compute(double time)
{
  DEBUG_BEGIN("GenericMechanical::compute(double time)\n");
  int info = 0;
  // --- Prepare data for GenericMechanical computing ---
  bool cont = preCompute(time);

  if(!cont)
    return info;
  // MB: if _hasBeenUpdated is set true then :
  // LinearOSNS.cpp:602
  // position unitialized, pos get a wrong value then :
  // computeqBlock(inter, pos) -> SEGFAULT
  // so I comment this:
  // _hasBeenUpdated = true;
  /*
    La matrice _M est construite.  Ici, il faut construire les
    sous-problemes, c'est a dire completer les champs des
    NumericsProblem (_mu, _e, _en, les dimentions...).  Il faut aussi
    remplir la sous matrice M du sous-probleme.  Pour cela, on peut
    boucler sur les interactions et completer le membres
    _numerics_problem.problems[i] and
    _numerics_problem.problemsType[i].
   */

  //......


  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)

  if(_sizeOutput != 0)
  {
    // The GenericMechanical Problem in Numerics format

    _pnumerics_GMP->M = &*_M->numericsMatrix();
    _pnumerics_GMP->q = &*_q->getArray();
    DEBUG_EXPR(display(););
    // Call Numerics Driver for GenericMechanical
    //    display();
    info = gmp_driver(_pnumerics_GMP,
                      &*_z->getArray(),
                      &*_w->getArray(),
                      &*_numerics_solver_options);
    //printf("GenericMechanical::compute : R:\n");
    //_z->display();
    postCompute();

  }
  else
  {

    DEBUG_PRINT("GenericMechanical::compute : sizeoutput is null\n");

  }
  DEBUG_END("GenericMechanical::compute(double time)\n");
  return info;
}

void GenericMechanical::display() const
{
  std::cout << "===== " << "Generic mechanical Problem " <<std::endl;
  LinearOSNS::display();
}


void  GenericMechanical::updateInteractionBlocks()
{
  if(!_hasBeenUpdated)
  {
    //    printf("GenericMechanical::updateInteractionBlocks : must be updated\n");
    genericMechanicalProblem_free(_pnumerics_GMP, NUMERICS_GMP_FREE_GMP);
    _pnumerics_GMP = genericMechanicalProblem_new();
  }
  LinearOSNS::updateInteractionBlocks();
}

GenericMechanical::~GenericMechanical()
{
  genericMechanicalProblem_free(_pnumerics_GMP, NUMERICS_GMP_FREE_GMP);
  _pnumerics_GMP = nullptr;
}
