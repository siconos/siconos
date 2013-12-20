/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
#include "EqualityConditionNSL.hpp"

#if (__cplusplus >= 201103L) && !defined(USE_BOOST_FOR_CXX11)
#include <functional>
#else
#include <boost/bind.hpp>
#endif


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

  _IG[0]->update_vertices_indices();
  _IG[0]->update_edges_indices();
}

// destructor
Topology::~Topology()
{
  clear();
}

std::pair<DynamicalSystemsGraph::EDescriptor, InteractionsGraph::VDescriptor> 
Topology::addInteractionInIndexSet(SP::Interaction inter, SP::DynamicalSystem ds1, SP::DynamicalSystem ds2)
{
  // !! Private function !!
  //
  // Add inter and ds into IG/DSG

  // Compute number of constraints
  unsigned int nsLawSize = inter->nonSmoothLaw()->size();
  unsigned int m = inter->getSizeOfY() / nsLawSize;
  if (m > 1)
    RuntimeException::selfThrow("Topology::addInteractionInIndexSet - m > 1. Obsolete !");

  _numberOfConstraints += nsLawSize;

  // _DSG is the hyper forest : (vertices : dynamical systems, edges :
  // Interactions)
  //
  // _IG is the hyper graph : (vertices : Interactions, edges :
  // dynamical systems)
  assert(_DSG[0]->edges_number() == _IG[0]->size());

  // _IG = L(_DSG),  L is the line graph transformation
  // vector of the Interaction
  DynamicalSystemsGraph::VDescriptor dsgv1, dsgv2;
  dsgv1 = _DSG[0]->add_vertex(ds1);
  if(ds2)
    dsgv2 = _DSG[0]->add_vertex(ds2);
  else 
    dsgv2 = dsgv1;
  
  // this may be a multi edges graph
  assert(!_DSG[0]->is_edge(dsgv1, dsgv2, inter));
  assert(!_IG[0]->is_vertex(inter));
  InteractionsGraph::VDescriptor ig_new_ve;
  DynamicalSystemsGraph::EDescriptor new_ed;
  std11::tie(new_ed, ig_new_ve) = _DSG[0]->add_edge(dsgv1, dsgv2, inter, *_IG[0]);
  
  // add self branches in vertex properties
  // note : boost graph SEGFAULT on self branch removal
  // see https://svn.boost.org/trac/boost/ticket/4622
  _IG[0]->properties(ig_new_ve).source = ds1;
  _IG[0]->properties(ig_new_ve).source_pos = 0;
  if(!ds2)
  {
    _IG[0]->properties(ig_new_ve).target = ds1;
    _IG[0]->properties(ig_new_ve).target_pos = 0;
  }
  else
  {
    _IG[0]->properties(ig_new_ve).target = ds2;
    _IG[0]->properties(ig_new_ve).target_pos = ds1->getDim();
  }

  assert(_IG[0]->bundle(ig_new_ve) == inter);
  assert(_IG[0]->is_vertex(inter));
  assert(_DSG[0]->is_edge(dsgv1, dsgv2, inter));
  assert(_DSG[0]->edges_number() == _IG[0]->size());
  
  return std::pair<DynamicalSystemsGraph::EDescriptor, InteractionsGraph::VDescriptor>(new_ed, ig_new_ve);
}

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
  
  SP::DynamicalSystem ds1 = _IG[0]->properties(_IG[0]->descriptor(inter)).source;
  SP::DynamicalSystem ds2 = _IG[0]->properties(_IG[0]->descriptor(inter)).target;
  _DSG[0]->remove_out_edge_if(_DSG[0]->descriptor(ds1), VertexIsRemoved(inter, _DSG[0], _IG[0]));
  if (ds1 != ds2)
    _DSG[0]->remove_out_edge_if(_DSG[0]->descriptor(ds2), VertexIsRemoved(inter, _DSG[0], _IG[0]));
}


void Topology::insertDynamicalSystem(SP::DynamicalSystem ds)
{
  _DSG[0]->add_vertex(ds);
}

void Topology::setControlProperty(SP::Interaction inter,
  const bool isControlInteraction)
{
  InteractionsGraph::VDescriptor ivd = _IG[0]->descriptor(inter);
  DynamicalSystemsGraph::VDescriptor dvd = _DSG[0]->descriptor(_IG[0]->properties(ivd).target);
  SP::Interaction interG;
  DynamicalSystemsGraph::OEIterator oei, oeiend;
  for (std11::tie(oei, oeiend) = _DSG[0]->out_edges(dvd); oei != oeiend; ++oei)
  {
    interG = _DSG[0]->bundle(*oei);
    if (inter == interG)
    {
      _DSG[0]->properties(*oei).forControl = isControlInteraction;
      break;
    }
  }
  _IG[0]->properties(ivd).forControl = isControlInteraction;
}

void Topology::removeInteraction(SP::Interaction inter)
{
  assert(_DSG[0]->edges_number() == _IG[0]->size());
  removeInteractionFromIndexSet(inter);
  assert(_DSG[0]->edges_number() == _IG[0]->size());
}

void Topology::removeDynamicalSystem(SP::DynamicalSystem ds)
{
  RuntimeException::selfThrow("remove dynamical system not implemented");
}

std::pair<DynamicalSystemsGraph::EDescriptor, InteractionsGraph::VDescriptor> 
Topology::link(SP::Interaction inter, SP::DynamicalSystem ds, SP::DynamicalSystem ds2)
{
  if (indexSet0()->is_vertex(inter))
  {
    removeInteractionFromIndexSet(inter);
  }

  DEBUG_PRINT("Topology::link(SP::Interaction inter, SP::DynamicalSystem ds)");
  unsigned int sumOfDSSizes = 0, sumOfZSizes = 0;
  
  sumOfDSSizes += ds->getDim();
  if(ds->z())
    sumOfZSizes += ds->z()->size();

  if(ds2)
  {
    sumOfDSSizes += ds->getDim();
    if(ds->z())
      sumOfZSizes += ds->z()->size();
    inter->setHas2Bodies(true);
  }

  inter->setDSSizes(sumOfDSSizes, sumOfZSizes);

  return addInteractionInIndexSet(inter, ds, ds2);
}

bool Topology::hasInteraction(SP::Interaction inter) const
{
  return indexSet0()->is_vertex(inter);
}

void Topology::initialize()
{
  _isTopologyUpToDate = true;
}

void Topology::setProperties()
{
  _IG[0]->update_vertices_indices();
  _IG[0]->update_edges_indices();

  for (unsigned int i = 1; i < _IG.size(); ++i)
  {
    // .. global properties may be defined here with
    // InteractionsSubGraphProperties(), see SiconosProperties.hpp
    // VertexSubProperties or EdgeSubProperties and the macros
    // INSTALL_GRAPH_PROPERTIES

    _IG[i]->properties().symmetric = _symmetric;

    _IG[i]->update_vertices_indices();
    _IG[i]->update_edges_indices();
  }

  _DSG[0]->properties().symmetric = _symmetric;

  _DSG[0]->update_vertices_indices();
  _DSG[0]->update_edges_indices();

}

void Topology::clear()
{
  _IG.clear();
  _DSG.clear();

  _isTopologyUpToDate = false;
}

SP::DynamicalSystem Topology::getDynamicalSystem(unsigned int requiredNumber)
{
  DynamicalSystemsGraph::VIterator vi, vdend;
  SP::DynamicalSystem ds;
  unsigned int currentNumber;
  for (std11::tie(vi, vdend) = _DSG[0]->vertices(); vi != vdend; ++vi)
  {
    ds = _DSG[0]->bundle(*vi);
    currentNumber = ds->number();
    if (currentNumber == requiredNumber)
      return ds;
  }

  RuntimeException::selfThrow("Topology::getDynamicalSystem(n) ds not found.");

  return ds;
}


