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
#include "Topology.hpp"
#include "NonSmoothLaw.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Interaction.hpp"
//#include "Interaction.hpp"
#include "EqualityConditionNSL.hpp"

#include <boost/bind.hpp>
#include <algorithm>
#include <limits>

//#define DEBUG_MESSAGES 1
#include <debug.h>

#define MAX_RELATIVE_DEGREE 999

// --- CONSTRUCTORS/DESTRUCTOR ---

// default
Topology::Topology(): _isTopologyUpToDate(false), _hasChanged(true),
  _numberOfConstraints(0), _symmetric(false)
{
  _IG.resize(1);
  _DSG.resize(1);

  _IG[0].reset(new InteractionsGraph());
  _DSG[0].reset(new DynamicalSystemsGraph());

  _allInteractions.reset(new InteractionsSet());
}

Topology::Topology(SP::InteractionsSet newInteractions) :
  _isTopologyUpToDate(false), _hasChanged(true),
  _numberOfConstraints(0), _symmetric(false)
{

  _IG.resize(1);
  _DSG.resize(1);

  _IG[0].reset(new InteractionsGraph());
  _DSG[0].reset(new DynamicalSystemsGraph());
  setProperties();

  _allInteractions.reset(new InteractionsSet());

  for (InteractionsIterator it = newInteractions->begin();
       it != newInteractions->end(); ++it)
  {
    insertInteraction(*it);
  }

  _isTopologyUpToDate = false;
}


// a constructor with a DS set : when some DS may not be in interactions
Topology::Topology(SP::DynamicalSystemsSet newDSset, SP::InteractionsSet newInteractions) :
  _isTopologyUpToDate(false), _hasChanged(true),
  _numberOfConstraints(0), _symmetric(false)
{

  _IG.resize(1);
  _DSG.resize(1);

  _IG[0].reset(new InteractionsGraph());
  _DSG[0].reset(new DynamicalSystemsGraph());
  setProperties();

  _allInteractions.reset(new InteractionsSet());

  for (InteractionsIterator it = newInteractions->begin();
       it != newInteractions->end(); ++it)
  {
    insertInteraction(*it);
  }

  for (DSIterator ids = newDSset->begin(); ids != newDSset->end() ; ++ids)
  {
    _DSG[0]->add_vertex(*ids);
  }
  _isTopologyUpToDate = false;
}

// destructor
Topology::~Topology()
{
  clear();
}

InteractionsGraph::VDescriptor Topology::addInteractionInIndexSet(SP::Interaction inter)
{
  // Private function
  //
  // Creates Interactions corresponding to inter and add them into
  // _IG

  // First, we get the number of relations in the interaction.  This
  // corresponds to inter->getNumberOfRelations but since Interaction
  // has not been initialized yet, this value is not set and we need
  // to get interaction size and nsLaw size.
  unsigned int nsLawSize = inter->nonSmoothLaw()->size();
  unsigned int m = inter->getSizeOfY() / nsLawSize;

  if (m > 1)
    RuntimeException::selfThrow("Topology::addInteractionInIndexSet - m > 1. Obsolete !");

  // vector of the Interaction
  InteractionsGraph::EDescriptor inter_current_edge;

  InteractionsGraph::EDescriptor ds_current_edge;

  SP::DynamicalSystemsSet systems = inter->dynamicalSystems();

  _numberOfConstraints += m * nsLawSize;

  // _DSG is the hyper forest : (vertices : dynamical systems, edges :
  // Interactions)
  //
  // _IG is the hyper graph : (vertices : Interactions, edges :
  // dynamical systems)

  // _IG = L(_DSG),  L is the line graph transformation


  // for all couples of ds in the interaction
  DEBUG_PRINTF("addInteractionInIndexSet systems->size() : %d\n", systems->size());


  InteractionsGraph::VDescriptor ig_new_ve;

  // only one or two ds! (otherwise we need a hyper-graph &
  // SiconosGraph does not provide it yet)
  for (DSIterator i1ds = systems->begin(); i1ds != systems->end(); ++i1ds)
  {
    for (DSIterator i2ds = i1ds; i2ds != systems->end(); ++i2ds)
    {

      // build the Interactions
      //      std::vector<SP::Interaction> current_interSet;
      //      for (unsigned int i=0, pos=0; i<m; ++i, pos += nsLawSize)
      //      {
      //        SP::Interaction inter(new Interaction(inter,pos,i));
      //        current_interSet.push_back(inter);
      //      };

      // one DS in the interaction is a special case : a self branch
      if ((i1ds == i2ds) && inter->dynamicalSystems()->size() == 1)
      {
        DynamicalSystemsGraph::VDescriptor dsgv;
        dsgv = _DSG[0]->add_vertex(*i1ds);

        // this may be a multi edges graph
        //        for (std::vector<SP::Interaction>::iterator itI = current_interSet.begin();
        //             itI != current_interSet.end(); ++itI)
        {
          assert(!_DSG[0]->is_edge(dsgv, dsgv, inter));
          assert(!_IG[0]->is_vertex(inter));

          DynamicalSystemsGraph::EDescriptor new_ed;
          boost::tie(new_ed, ig_new_ve) = _DSG[0]->add_edge(dsgv, dsgv, inter, *_IG[0]);

          // add self branches in vertex properties
          // note : boost graph SEGFAULT on self branch removal
          // see https://svn.boost.org/trac/boost/ticket/4622
          _IG[0]->properties(ig_new_ve).source = *i1ds;
          _IG[0]->properties(ig_new_ve).target = *i1ds;

          assert(_IG[0]->bundle(ig_new_ve) == inter);
          assert(_IG[0]->is_vertex(inter));
          assert(_DSG[0]->is_edge(dsgv, dsgv, inter));
        }
      }
      else
        // multiples DS in the interaction : no self branch
        if (i1ds != i2ds)
        {
          DynamicalSystemsGraph::VDescriptor dsgv1, dsgv2;

          dsgv1 = _DSG[0]->add_vertex(*i1ds);
          dsgv2 = _DSG[0]->add_vertex(*i2ds);

          // this may be a multi edges graph
          //          for (std::vector<SP::Interaction>::iterator itI = current_interSet.begin();
          //               itI != current_interSet.end(); ++itI)
          {
            assert(!_DSG[0]->is_edge(dsgv1, dsgv2, inter));
            assert(!_IG[0]->is_vertex(inter));

            DynamicalSystemsGraph::EDescriptor new_ed;
            boost::tie(new_ed, ig_new_ve) = _DSG[0]->add_edge(dsgv1, dsgv2, inter, *_IG[0]);

            // add self branches in vertex properties
            // note : boost graph SEGFAULT on self branch removal
            // see https://svn.boost.org/trac/boost/ticket/4622
            _IG[0]->properties(ig_new_ve).source = *i1ds;
            _IG[0]->properties(ig_new_ve).target = *i2ds;

            assert(_IG[0]->bundle(ig_new_ve) == inter);
            assert(_IG[0]->is_vertex(inter));
            assert(_DSG[0]->is_edge(dsgv1, dsgv2, inter));
          }
        }
    }
  }

  // note: only one or two ds => only one vertex in IG
  return ig_new_ve;

};

/* an edge is removed from _DSG graph if the corresponding vertex is
   removed from the adjoint graph (_IG)
*/
struct VertexIsRemoved
{
  VertexIsRemoved(SP::Interaction I,
                  SP::DynamicalSystemsGraph sg, SP::InteractionsGraph asg) :
    _I(I), __DSG(sg), __IG(asg) {};
  bool operator()(DynamicalSystemsGraph::EDescriptor ed)
  {

    if (__IG->is_vertex(__DSG->bundle(ed)))
    {
      InteractionsGraph::VDescriptor ivd = __IG->descriptor(__DSG->bundle(ed));

      if (__IG->bundle(ivd) == _I)
      {
        __IG->remove_vertex(__DSG->bundle(ed));

        assert(__IG->size() == __DSG->edges_number() - 1);

        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return true;
    }
  }
  SP::Interaction _I;
  SP::DynamicalSystemsGraph __DSG;
  SP::InteractionsGraph __IG;
};


/* remove an interaction : remove edges (Interaction) from _DSG if
   corresponding vertices are removed from _IG */
void Topology::removeInteractionFromIndexSet(SP::Interaction inter)
{

  for (DSIterator ids = inter->dynamicalSystems()->begin();
       ids != inter->dynamicalSystems()->end();
       ++ids)
  {
    _DSG[0]->remove_out_edge_if
    (_DSG[0]->descriptor(*ids),
     VertexIsRemoved(inter, _DSG[0], _IG[0]));
  };
};


void Topology::insertDynamicalSystem(SP::DynamicalSystem ds)
{
  _DSG[0]->add_vertex(ds);
};

InteractionsGraph::VDescriptor Topology::insertInteraction(SP::Interaction inter)
{
  assert(_allInteractions);
  assert(_DSG[0]->edges_number() == _IG[0]->size());

  _allInteractions->insert(inter);
  InteractionsGraph::VDescriptor ig_new_ve = addInteractionInIndexSet(inter);

  assert(_DSG[0]->edges_number() == _IG[0]->size());

  return ig_new_ve;

};

void Topology::removeInteraction(SP::Interaction inter)
{
  assert(_allInteractions);
  assert(_DSG[0]->edges_number() == _IG[0]->size());

  _allInteractions->erase(inter);
  removeInteractionFromIndexSet(inter);

  assert(_DSG[0]->edges_number() == _IG[0]->size());
};

void Topology::removeDynamicalSystem(SP::DynamicalSystem ds)
{
  RuntimeException::selfThrow("remove dynamical system not implemented");
};


void Topology::link(SP::Interaction inter, SP::DynamicalSystem ds)
{
  // interactions should not know linked dynamical systems in the
  // future
  InteractionsIterator it = _allInteractions->find(inter);
  if (it != _allInteractions->end())
  {
    _allInteractions->erase(*it);
    removeInteractionFromIndexSet(inter);
  }

  inter->insert(ds);
  insertInteraction(inter);

};




bool Topology::hasInteraction(SP::Interaction inter) const
{
  return _allInteractions->isIn(inter);
}

void Topology::initialize()
{
  _isTopologyUpToDate = true;
}

void Topology::setProperties()
{
  _IG[0]->properties().reset(new InteractionsGraphProperties(_IG[0]));

  for (unsigned int i = 1; i < _IG.size(); ++i)
  {
    _IG[i]->properties().reset(new InteractionsGraphProperties(_IG[i]));
    // .. global properties may be defined here with
    // InteractionsSubGraphProperties(), see SiconosProperties.hpp
    // VertexSubProperties or EdgeSubProperties and the macros
    // INSTALL_GRAPH_PROPERTIES

    _IG[i]->properties()->symmetric = _symmetric;
  }

  _DSG[0]->properties().reset(new DynamicalSystemsGraphProperties(_DSG[0]));
  _DSG[0]->properties()->symmetric = _symmetric;
}

void Topology::clear()
{
  _allInteractions->clear();

  _IG.clear();
  _DSG.clear();

  _isTopologyUpToDate = false;
}

