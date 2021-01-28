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
#include "AVI.hpp"
#include "NumericsMatrix.h"
#include <assert.h>
#include "Simulation.hpp"
#include "NormalConeNSL.hpp"
#include "OSNSMatrix.hpp"
#include "SiconosSets.h" // from numerics, for polyhedron
//#include <AVI_Solvers.h>
#include <NonSmoothDrivers.h>
#include "SolverOptions.h" // for solver_options_create ...

#include <limits>
#include <cstdlib>
#include <iostream>
#include <cassert>

/*****************************************************
 * START visitor for nslaw
*/
#if 0
struct AVI::_BoundsNSLEffect : public SiconosVisitor
{

  using SiconosVisitor::visit;

  AVI* _parent;
  SP::Interaction _inter;
  unsigned int _pos;


  _BoundsNSLEffect(AVI *p, SP::Interaction inter, unsigned int pos) :
    _parent(p), _inter(inter), _pos(pos) {};

  void visit(const NormalConeNSL& nslaw)
  {
    if(_pos > 0)
    {
      S
    }
    // take the
    SiconosVector& K = nslaw.K();
    SimpleMatrix& H = nslaw.H();
    _numerics_problem->size = nslaw.size();
    _numerics_problem->d = nullptr;
    _numerics_problem->poly->id = SICONOS_SET_POLYHEDRON;
    _numerics_problem->poly->size_ineq = K.size();
    _numerics_problem->poly->size_eq = 0;
    _numerics_problem->poly->H = H.getArray();
    _numerics_problem->poly->K = K.getArray();
    _numerics_problem->poly->Heq = nullptr;
    _numerics_problem->poly->Keq = nullptr;
  }

  void visit(const RelayNSL& nslaw)
  {
    Siconos

  }

  void visit(const ComplementarityConditionNSL& nslaw)
  {
  }

};
#endif
/*****************************************************
 * END visitor for nslaw
*/


AVI::AVI(int numericsSolverId):
  AVI(SP::SolverOptions(solver_options_create(numericsSolverId),
                        solver_options_delete))
{}

AVI::AVI(SP::SolverOptions options):
  LinearOSNS(options), _numerics_problem(new AffineVariationalInequalities)
{
  _numerics_problem->poly.split = new polyhedron;
}


void AVI::initialize(SP::Simulation sim)
{
  LinearOSNS::initialize(sim);

  // right now we support only one (1) NonsmoothLaw associated with this AVI
  // It is not clear whether having multiple NonsmoothLaw would be beneficial given the exponential complexity of most solvers
  // TODO We should support RelayNSL with generic rectangles -- xhub
  InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());
  InteractionsGraph::VIterator ui, uiend;
  unsigned nbInter = 0;
  for(std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
  {
    NormalConeNSL& nc = static_cast<NormalConeNSL&>(*indexSet.bundle(*ui)->nonSmoothLaw());
    assert(Type::value(nc) == Type::NormalConeNSL &&
           "AVI::initialize :: found a NonSmoothLaw that is not of the NormalConeNSL type! This is currently not supported");
    SiconosVector& K = nc.K();
    SimpleMatrix& H = nc.H();
    _numerics_problem->size = nc.size();
    _numerics_problem->d = nullptr;
    _numerics_problem->poly.split->id = SICONOS_SET_POLYHEDRON;
    _numerics_problem->poly.split->size_ineq = K.size();
    _numerics_problem->poly.split->size_eq = 0;
    _numerics_problem->poly.split->H = NM_create_from_data(NM_DENSE, K.size(), nc.size(), H.getArray());
    _numerics_problem->poly.split->K = K.getArray();
    _numerics_problem->poly.split->Heq = nullptr;
    _numerics_problem->poly.split->Keq= nullptr;

    // we do not support more than one interaction
    if(!(nbInter++ == 0))
      THROW_EXCEPTION("AVI::initialize :: more than one Interactions for this OneStepNSProblem is not support ATM!");
  }

}
bool AVI::checkCompatibleNSLaw(NonSmoothLaw& nslaw)
{
  float type_number= (float) (Type::value(nslaw));
  _nslawtype.insert(type_number);

  if (not (Type::value(nslaw) == Type::NormalConeNSL))
  {
    THROW_EXCEPTION("\nAVI::checkCompatibleNSLaw -  \n\
                      The chosen nonsmooth law is not compatible with AVVI one step nonsmooth problem. \n \
                      Compatible NonSmoothLaw are: NormalConeNSL\n");
    return false;
  }
  return true;
}


int AVI::compute(double time)
{
  int info = 0;
  // --- Prepare data for AVI computing ---
  bool cont = preCompute(time);
  if(!cont)
    return info;

  if(_numerics_problem->size != _sizeOutput)
  {
    THROW_EXCEPTION("AVI::compute - size mismatch between AVI size and and the current size");
  }

  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)

  if(_sizeOutput != 0)
  {
    // The AVI in Numerics format
    _numerics_problem->M = _M->numericsMatrix().get();
    _numerics_problem->q = _q->getArray();

    info = avi_driver(_numerics_problem.get(), _z->getArray(), _w->getArray(),
                      _numerics_solver_options.get());

    if(info != 0)
    {
      std::cout << "Warning : Problem in AVI resolution" <<std::endl;
    }

    // --- Recovering of the desired variables from AVI output ---
    postCompute();
  }

  return info;
}

void AVI::display() const
{
  std::cout << "======= AVI of size " << _sizeOutput << " with: " <<std::endl;
  LinearOSNS::display();
}

AVI::~AVI()
{
  delete _numerics_problem->poly.split;
}

