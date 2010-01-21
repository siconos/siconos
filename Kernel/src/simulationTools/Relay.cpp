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
#include "Relay.hpp"
#include <iostream>
#include <assert.h>
#include "Tools.hpp"


#include "Simulation.hpp"
#include "RelayNSL.hpp"


using namespace std;
using namespace RELATION;

/* nslaw dispatch on bounds */

struct Relay::_BoundsNSLEffect : public SiconosVisitor
{
  Relay *parent;
  unsigned int pos;
  SP::UnitaryRelation UR;

  _BoundsNSLEffect(Relay *p, SP::UnitaryRelation UR, unsigned int pos) :
    parent(p), UR(UR), pos(pos) {};

  void visit(RelayNSL& nslaw)
  {

    // cout << "Relay::_BoundsNSLEffect visit"  << endl;
    for (int i = 0; i <  UR->getNonSmoothLawSize(); i++)
    {
      (*(parent->lb()))(pos + i) =
        nslaw.lb();
      (*(parent->ub()))(pos + i) =
        nslaw.ub();
    }
  }

  void visit(ComplementarityConditionNSL& nslaw)
  {
    for (int i = 0; i <  UR->getNonSmoothLawSize(); i++)
    {
      (*(parent->lb()))(pos + i) = 0.0;
      (*(parent->ub()))(pos + i) = 1e+24;
    }


    // \warning Introduce the infinte double symbol rather 1e+24

  }

};

void Relay::initialize(SP::Simulation sim)
{
  LinearOSNS::initialize(sim);
  //cout << "Relay::initialize" << endl;


  // initialize memory for _lb and _ub
  if (! _lb)
    _lb.reset(new SimpleVector(maxSize()));
  else
  {
    if (_lb->size() != maxSize())
      _lb->resize(maxSize());
  }
  if (! _ub)
    _ub.reset(new SimpleVector(maxSize()));
  else
  {
    if (_ub->size() != maxSize())
      _ub->resize(maxSize());
  }

  // fill _lb and _ub wiht the value of the NonSmooth Law

  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(levelMin());

  //cout << " _sizeOutput =" <<_sizeOutput << endl;
  if (_lb->size() != _sizeOutput)
  {
    _lb->resize(_sizeOutput, false);
    _lb->zero();
  }
  if (_ub->size() != _sizeOutput)
  {
    _ub->resize(_sizeOutput, false);
    _ub->zero();
  }

  unsigned int pos = 0;
  UnitaryRelationsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::UnitaryRelation ur = indexSet->bundle(*ui);

    // Compute q, this depends on the type of non smooth problem, on
    // the relation type and on the non smooth law
    pos = _M->getPositionOfUnitaryBlock(ur);
    SP::SiconosVisitor NSLEffect(new _BoundsNSLEffect(this, ur, pos));
    ur->interaction()->nonSmoothLaw()->accept(*NSLEffect);
  }
  // Initialization of the NonSmoothSolver
  _solver->initialize(this) ;
}


int Relay::compute(double time)
{
  // --- Prepare data for Relay computing ---
  preCompute(time);

  int info = 0;
  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)

  if (_sizeOutput != 0)
  {
    // The Relay in Numerics format
    Relay_Problem numerics_problem;
    numerics_problem.M = &*_M->getNumericsMatrix();
    numerics_problem.q = _q->getArray();
    numerics_problem.lb = _lb->getArray();
    numerics_problem.ub = _ub->getArray();
    numerics_problem.size = _sizeOutput;

    int nbSolvers = 1;
    // Call Relay Driver

    //      Relay_display(&numerics_problem);

    info = relay_driver(&numerics_problem, _z->getArray() , _w->getArray() ,
                        &*_solver->numericsSolverOptions(), nbSolvers, &*_numerics_options);

    if (info != 0)
    {
      cout << "Warning : Problem in Relay resolution" << endl;
    }

    // --- Recovering of the desired variables from Relay output ---
    postCompute();

  }

  return info;
}

void Relay::display() const
{
  cout << "======= Relay of size " << _sizeOutput << " with: " << endl;
  LinearOSNS::display();
}

Relay* Relay::convert(OneStepNSProblem* osnsp)
{
  Relay* lcp = dynamic_cast<Relay*>(osnsp);
  return lcp;
}


