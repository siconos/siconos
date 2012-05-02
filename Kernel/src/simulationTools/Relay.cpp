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
#include "Relay.hpp"
#include <iostream>
#include <assert.h>
#include "Tools.hpp"


#include "Simulation.hpp"
#include "RelayNSL.hpp"


using namespace std;
using namespace RELATION;


Relay::Relay(const int newNumericsSolverId , const std::string& newId):
  LinearOSNS(newNumericsSolverId, "Relay", newId)
{
  _numerics_problem.reset(new  RelayProblem);

  relay_setDefaultSolverOptions(NULL, &*_numerics_solver_options, newNumericsSolverId);


};
/* nslaw dispatch on bounds */

struct Relay::_BoundsNSLEffect : public SiconosVisitor
{
  Relay* _parent;
  SP::Interaction _inter;
  unsigned int _pos;


  _BoundsNSLEffect(Relay *p, SP::Interaction inter, unsigned int pos) :
    _parent(p), _inter(inter), _pos(pos) {};

  void visit(const RelayNSL& nslaw)
  {

    // cout << "Relay::_BoundsNSLEffect visit"  << endl;
    for (unsigned int i = 0; i <  _inter->getNonSmoothLawSize(); i++)
    {
      (*(_parent->lb()))(_pos + i) =
        nslaw.lb();
      (*(_parent->ub()))(_pos + i) =
        nslaw.ub();
    }
  }

  void visit(const ComplementarityConditionNSL& nslaw)
  {
    for (unsigned int i = 0; i <  _inter->getNonSmoothLawSize(); i++)
    {
      (*(_parent->lb()))(_pos + i) = 0.0;
      (*(_parent->ub()))(_pos + i) = 1e+24;
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
}


int Relay::compute(double time)
{
  // --- Prepare data for Relay computing ---
  preCompute(time);

  // fill _lb and _ub wiht the value of the NonSmooth Law

  SP::InteractionsGraph indexSet = simulation()->indexSet(levelMin());

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
  InteractionsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet->bundle(*ui);

    // Compute q, this depends on the type of non smooth problem, on
    // the relation type and on the non smooth law
    pos = _M->getPositionOfInteractionBlock(inter);
    SP::SiconosVisitor NSLEffect(new _BoundsNSLEffect(this, inter, pos));
    inter->nonSmoothLaw()->accept(*NSLEffect);
  }











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
    RelayProblem numerics_problem;
    numerics_problem.M = &*_M->getNumericsMatrix();
    numerics_problem.q = _q->getArray();
    numerics_problem.lb = _lb->getArray();
    numerics_problem.ub = _ub->getArray();
    numerics_problem.size = _sizeOutput;

    //int nbSolvers = 1;
    // Call Relay Driver

    //      Relay_display(&numerics_problem);

    info = relay_driver(&numerics_problem, _z->getArray() , _w->getArray() ,
                        &*_numerics_solver_options, &*_numerics_options);

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


Relay::~Relay()
{
  deleteSolverOptions(&*_numerics_solver_options);
}

