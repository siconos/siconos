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
#include "Topology.hpp"

#include "DynamicalSystem.hpp"
#include "Interaction.hpp"
#include "NonSmoothLaw.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES 1
#include "siconos_debug.h"

// default
siconos::simulation::Topology::Topology()
{
  _IG.resize(1);
  _DSG.resize(1);

  _IG[0] = std::make_shared<siconos::graphs::InteractionsGraph>();
  _DSG[0] = std::make_shared<siconos::graphs::DynamicalSystemsGraph>();

  _IG[0]->update_vertices_indices();
  _IG[0]->update_edges_indices();
}

// destructor
siconos::simulation::Topology::~Topology() noexcept { clear(); }

std::pair<siconos::graphs::DynamicalSystemsGraph::EDescriptor,
          siconos::graphs::InteractionsGraph::VDescriptor>
siconos::simulation::Topology::__addInteractionInIndexSet0(
    std::shared_ptr<siconos::modeling::Interaction> inter,
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds1,
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds2)
{
  // !! Private function !!
  //
  // This function must
  // - insert interaction and ds into IG/DSG
  // - update graph properties related to modeling (DSlink ...)
  //
  // Expected result : all interaction methods can be safely called after a call
  // to this function (whatever the simulation is, if it exists or not).

  // Update total number of constraints
  _numberOfConstraints += inter->nonSmoothLaw()->size();

  auto ds2_ = ds2;
  // _DSG is the hyper forest : (vertices : dynamical systems, edges :
  // Interactions)
  //
  // _IG is the hyper graph : (vertices : Interactions, edges :
  // dynamical systems)
  assert(_DSG[0]->edges_number() == _IG[0]->size());

  // _IG = L(_DSG),  L is the line graph transformation
  // vector of the Interaction
  siconos::graphs::DynamicalSystemsGraph::VDescriptor dsgv1, dsgv2;
  dsgv1 = _DSG[0]->add_vertex(ds1);

  if (ds2) {
    dsgv2 = _DSG[0]->add_vertex(ds2);
  }
  else {
    dsgv2 = dsgv1;
    ds2_ = ds1;
  }

  // this may be a multi edges graph
  assert(!_DSG[0]->is_edge(dsgv1, dsgv2, inter));
  assert(!_IG[0]->is_vertex(inter));
  siconos::graphs::InteractionsGraph::VDescriptor ig_new_ve;
  siconos::graphs::DynamicalSystemsGraph::EDescriptor new_ed;
  std::tie(new_ed, ig_new_ve) = _DSG[0]->add_edge(dsgv1, dsgv2, inter, *_IG[0]);

  inter->initializeLinkToDsVariables(*ds1, *ds2_);

  // add self branches in vertex properties
  // note : boost graph SEGFAULT on self branch removal
  // see https://svn.boost.org/trac/boost/ticket/4622
  _IG[0]->properties(ig_new_ve).source = ds1;
  _IG[0]->properties(ig_new_ve).source_pos = 0;

  if (!ds2)  // self loop in ds graph, source=target in interaction graph
  {
    _IG[0]->properties(ig_new_ve).target = ds1;
    _IG[0]->properties(ig_new_ve).target_pos = 0;
  }
  else {
    _IG[0]->properties(ig_new_ve).target = ds2;
    _IG[0]->properties(ig_new_ve).target_pos = ds1->dimension();
  }

  assert(_IG[0]->bundle(ig_new_ve) == inter);
  assert(_IG[0]->is_vertex(inter));
  assert(_DSG[0]->is_edge(dsgv1, dsgv2, inter));
  assert(_DSG[0]->edges_number() == _IG[0]->size());

  return std::pair<siconos::graphs::DynamicalSystemsGraph::EDescriptor,
                   siconos::graphs::InteractionsGraph::VDescriptor>(new_ed, ig_new_ve);
}

/* an edge is removed from _DSG graph if the corresponding vertex is
   removed from the adjoint graph (_IG)
*/
struct VertexIsRemoved {
  VertexIsRemoved(std::shared_ptr<siconos::modeling::Interaction> I,
                  std::shared_ptr<siconos::graphs::DynamicalSystemsGraph> sg,
                  std::shared_ptr<siconos::graphs::InteractionsGraph> asg)
      : _I(I), __DSG(sg), __IG(asg){};
  bool operator()(siconos::graphs::DynamicalSystemsGraph::EDescriptor ed)
  {
    if (__IG->is_vertex(__DSG->bundle(ed))) {
      auto ivd = __IG->descriptor(__DSG->bundle(ed));

      if (__IG->bundle(ivd) == _I) {
        __IG->remove_vertex(__DSG->bundle(ed));

        assert(__IG->size() == __DSG->edges_number() - 1);

        return true;
      }
      else {
        return false;
      }
    }
    else {
      return true;
    }
  }
  std::shared_ptr<siconos::modeling::Interaction> _I{nullptr};
  std::shared_ptr<siconos::graphs::DynamicalSystemsGraph> __DSG{nullptr};
  std::shared_ptr<siconos::graphs::InteractionsGraph> __IG{nullptr};
};

/* an edge is removed from _DSG graph if the corresponding vertex is
   removed from the adjoint graph (_IG)
*/
struct VertexIsRemovedDS {
  VertexIsRemovedDS(std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
                    std::shared_ptr<siconos::graphs::DynamicalSystemsGraph> sg,
                    std::shared_ptr<siconos::graphs::InteractionsGraph> asg)
      : _ds(ds), __DSG(sg), __IG(asg){};
  bool operator()(siconos::graphs::DynamicalSystemsGraph::EDescriptor ed)
  {
    if (__IG->is_vertex(__DSG->bundle(ed))) {
      auto ivd = __IG->descriptor(__DSG->bundle(ed));

      if (__IG->properties(ivd).source == _ds || __IG->properties(ivd).target == _ds) {
        __IG->remove_vertex(__DSG->bundle(ed));

        assert(__IG->size() == __DSG->edges_number() - 1);

        return true;
      }
      else {
        return false;
      }
    }
    else {
      return true;
    }
  }
  std::shared_ptr<siconos::modeling::DynamicalSystem> _ds{nullptr};
  std::shared_ptr<siconos::graphs::DynamicalSystemsGraph> __DSG{nullptr};
  std::shared_ptr<siconos::graphs::InteractionsGraph> __IG{nullptr};
};

/* remove an interaction : remove edges (Interaction) from _DSG if
   corresponding vertices are removed from _IG */
void siconos::simulation::Topology::__removeInteractionFromIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter)
{
  std::shared_ptr<siconos::modeling::DynamicalSystem> ds1 =
      _IG[0]->properties(_IG[0]->descriptor(inter)).source;
  std::shared_ptr<siconos::modeling::DynamicalSystem> ds2 =
      _IG[0]->properties(_IG[0]->descriptor(inter)).target;
  _DSG[0]->remove_out_edge_if(_DSG[0]->descriptor(ds1),
                              VertexIsRemoved(inter, _DSG[0], _IG[0]));
  if (ds1 != ds2)
    _DSG[0]->remove_out_edge_if(_DSG[0]->descriptor(ds2),
                                VertexIsRemoved(inter, _DSG[0], _IG[0]));
}

void siconos::simulation::Topology::insertDynamicalSystem(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  _DSG[0]->add_vertex(ds);
  setHasChanged(true);
}

/* remove a dynamical system : remove edges (DynamicalSystem) from _IG if
   corresponding vertices are removed from _DSG */
void siconos::simulation::Topology::__removeDynamicalSystemFromIndexSet(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  _DSG[0]->remove_edge_if(_DSG[0]->descriptor(ds), VertexIsRemovedDS(ds, _DSG[0], _IG[0]));

  // note: remove_vertex also calls clear_vertex and removes all in/out edges
  _DSG[0]->remove_vertex(ds);
}

void siconos::simulation::Topology::setName(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds, const std::string& name)
{
  auto dsgv = _DSG[0]->descriptor(ds);
  _DSG[0]->name.insert(dsgv, name);
}

std::string siconos::simulation::Topology::name(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  auto dsgv = _DSG[0]->descriptor(ds);
  if (dsgv)
    return _DSG[0]->name.at(dsgv);
  else
    return "";
}

void siconos::simulation::Topology::setName(
    std::shared_ptr<siconos::modeling::Interaction> inter, const std::string& name)
{
  auto igv = _IG[0]->descriptor(inter);
  _IG[0]->name.insert(igv, name);
}

std::string siconos::simulation::Topology::name(
    std::shared_ptr<siconos::modeling::Interaction> inter)
{
  auto igv = _IG[0]->descriptor(inter);
  if (igv)
    return _IG[0]->name.at(igv);
  else
    return "";
}

void siconos::simulation::Topology::setOSI(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
    std::shared_ptr<siconos::integrators::OneStepIntegrator> OSI)
{
  _DSG[0]->properties(_DSG[0]->descriptor(ds)).osi = OSI;
}

void siconos::simulation::Topology::setControlProperty(
    std::shared_ptr<siconos::modeling::Interaction> inter, const bool isControlInteraction)
{
  auto ivd = _IG[0]->descriptor(inter);
  auto dvd = _DSG[0]->descriptor(_IG[0]->properties(ivd).target);
  std::shared_ptr<siconos::modeling::Interaction> interG;
  siconos::graphs::DynamicalSystemsGraph::OEIterator oei, oeiend;
  for (std::tie(oei, oeiend) = _DSG[0]->out_edges(dvd); oei != oeiend; ++oei) {
    interG = _DSG[0]->bundle(*oei);
    if (inter == interG) {
      _DSG[0]->properties(*oei).forControl = isControlInteraction;
      break;
    }
  }
  _IG[0]->properties(ivd).forControl = isControlInteraction;
}

std::shared_ptr<siconos::graphs::InteractionsGraph> siconos::simulation::Topology::indexSet(
    unsigned int num) const
{
  if (num >= _IG.size()) {
    THROW_EXCEPTION("siconos::simulation::Topology::indexSet: indexSet does not exist");
  }
  assert(num < _IG.size());
  return _IG[num];
};

void siconos::simulation::Topology::removeInteraction(
    std::shared_ptr<siconos::modeling::Interaction> inter)
{
  DEBUG_PRINTF("removeInteraction : %p\n", &*inter);

  assert(_DSG[0]->edges_number() == _IG[0]->size());
  __removeInteractionFromIndexSet(inter);
  assert(_DSG[0]->edges_number() == _IG[0]->size());
  setHasChanged(true);
}

void siconos::simulation::Topology::removeDynamicalSystem(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  DEBUG_PRINTF("removeDynamicalSystem : %p\n", &*ds);

  assert(_DSG[0]->edges_number() == _IG[0]->size() && "1");
  __removeDynamicalSystemFromIndexSet(ds);
  assert(_DSG[0]->edges_number() == _IG[0]->size() && "2");
  setHasChanged(true);
}

std::pair<siconos::graphs::DynamicalSystemsGraph::EDescriptor,
          siconos::graphs::InteractionsGraph::VDescriptor>
siconos::simulation::Topology::link(std::shared_ptr<siconos::modeling::Interaction> inter,
                                    std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
                                    std::shared_ptr<siconos::modeling::DynamicalSystem> ds2)
{
  DEBUG_PRINTF("siconos::simulation::Topology::link : inter %p, ds1 %p, ds2 %p\n", &*inter,
               &*ds, &*ds2);

  // If the interaction is already in the graph remove it
  if (indexSet0()->is_vertex(inter)) {
    __removeInteractionFromIndexSet(inter);
  }

  // Compute interaction dimension (sum of involved dynamical systems sizes)
  unsigned int sumOfDSSizes = 0;
  sumOfDSSizes += ds->dimension();
  if (ds2) {
    sumOfDSSizes += ds2->dimension();
    inter->setHas2Bodies(true);
  }
  inter->setDSSizes(sumOfDSSizes);

  return __addInteractionInIndexSet0(inter, ds, ds2);
}

bool siconos::simulation::Topology::hasInteraction(
    std::shared_ptr<siconos::modeling::Interaction> inter) const
{
  return indexSet0()->is_vertex(inter);
}

bool siconos::simulation::Topology::hasDynamicalSystem(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) const
{
  return _DSG[0]->is_vertex(ds);
}

void siconos::simulation::Topology::setProperties()
{
  _IG[0]->update_vertices_indices();
  _IG[0]->update_edges_indices();

  for (unsigned int i = 1; i < _IG.size(); ++i) {
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

void siconos::simulation::Topology::clear()
{
  _IG.clear();
  _DSG.clear();
}

std::shared_ptr<siconos::modeling::DynamicalSystem>
siconos::simulation::Topology::getDynamicalSystem(unsigned int requiredNumber) const
{
  siconos::graphs::DynamicalSystemsGraph::VIterator vi, vdend;
  std::shared_ptr<siconos::modeling::DynamicalSystem> ds;
  unsigned int currentNumber;
  for (std::tie(vi, vdend) = _DSG[0]->vertices(); vi != vdend; ++vi) {
    ds = _DSG[0]->bundle(*vi);
    currentNumber = ds->number();
    if (currentNumber == requiredNumber) return ds;
  }

  THROW_EXCEPTION("siconos::simulation::Topology::getDynamicalSystem(n) ds not found.");

  return ds;
}

void siconos::simulation::Topology::displayDynamicalSystems() const
{
  siconos::graphs::DynamicalSystemsGraph::VIterator vi, vdend;
  std::shared_ptr<siconos::modeling::DynamicalSystem> ds;
  unsigned int currentNumber;
  for (std::tie(vi, vdend) = _DSG[0]->vertices(); vi != vdend; ++vi) {
    ds = _DSG[0]->bundle(*vi);
    currentNumber = ds->number();
    std::cout << "Dynamical system number : " << currentNumber << std::endl;
    ds->display();
  }
}

std::shared_ptr<siconos::modeling::DynamicalSystem>
siconos::simulation::Topology::getDynamicalSystem(std::string name) const
{
  siconos::graphs::DynamicalSystemsGraph::VIterator vi, vdend;
  for (std::tie(vi, vdend) = _DSG[0]->vertices(); vi != vdend; ++vi) {
    if (name == _DSG[0]->name.at(*vi)) return _DSG[0]->bundle(*vi);
  }

  THROW_EXCEPTION("siconos::simulation::Topology::getDynamicalSystem() ds not found.");

  return std::shared_ptr<siconos::modeling::DynamicalSystem>();
}

unsigned int siconos::simulation::Topology::numberOfInvolvedDS(unsigned int inumber)
{
  if (inumber >= _IG.size()) {
    THROW_EXCEPTION(
        "siconos::simulation::Topology::numberOfInvolvedDS :: index number must be smaller "
        "than the number of indexSets");
  }

  /* on an adjoint graph a dynamical system may be on several edges */
  std::map<std::shared_ptr<siconos::modeling::DynamicalSystem>, bool> flag;

  unsigned int return_value = 0;

  auto indexSet = _IG[inumber];

  siconos::graphs::InteractionsGraph::VIterator vi, viend;
  for (std::tie(vi, viend) = indexSet->vertices(); vi != viend; ++vi) {
    if (indexSet->properties(*vi).source) {
      if (flag.find(indexSet->properties(*vi).source) == flag.end()) {
        flag[indexSet->properties(*vi).source] = true;
        return_value++;
      }
    }
    if (indexSet->properties(*vi).target) {
      if (flag.find(indexSet->properties(*vi).target) == flag.end()) {
        flag[indexSet->properties(*vi).target] = true;
        return_value++;
      }
    }
  }

  siconos::graphs::InteractionsGraph::EIterator ei, eiend;

  for (std::tie(ei, eiend) = indexSet->edges(); ei != eiend; ++ei) {
    if (flag.find(indexSet->bundle(*ei)) == flag.end()) {
      flag[indexSet->bundle(*ei)] = true;
      return_value++;
    }
  }

  return return_value;
}

std::shared_ptr<siconos::modeling::Interaction> siconos::simulation::Topology::getInteraction(
    unsigned int requiredNumber) const
{
  siconos::graphs::InteractionsGraph::VIterator vi, vdend;
  std::shared_ptr<siconos::modeling::Interaction> inter;
  unsigned int currentNumber;
  for (std::tie(vi, vdend) = _IG[0]->vertices(); vi != vdend; ++vi) {
    inter = _IG[0]->bundle(*vi);
    currentNumber = inter->number();
    if (currentNumber == requiredNumber) return inter;
  }

  return inter;
}

std::shared_ptr<siconos::modeling::Interaction> siconos::simulation::Topology::getInteraction(
    std::string name) const
{
  siconos::graphs::DynamicalSystemsGraph::VIterator vi, vdend;
  for (std::tie(vi, vdend) = _IG[0]->vertices(); vi != vdend; ++vi) {
    if (name == _IG[0]->name.at(*vi)) return _IG[0]->bundle(*vi);
  }

  return std::shared_ptr<siconos::modeling::Interaction>();
}

std::vector<std::shared_ptr<siconos::modeling::Interaction>>
siconos::simulation::Topology::interactionsForDS(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) const
{
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  std::shared_ptr<siconos::modeling::Interaction> inter;
  std::vector<std::shared_ptr<siconos::modeling::Interaction>> result;
  if (!ds) return result;
  for (std::tie(ui, uiend) = _IG[0]->vertices(); ui != uiend; ++ui) {
    inter = _IG[0]->bundle(*ui);
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds1 =
        _IG[0]->properties(_IG[0]->descriptor(inter)).source;
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds2 =
        _IG[0]->properties(_IG[0]->descriptor(inter)).target;
    if (ds == ds1 || ds == ds2) result.push_back(inter);
  }
  return result;
}

std::vector<std::shared_ptr<siconos::modeling::Interaction>>
siconos::simulation::Topology::interactionsForPairOfDS(
    std::shared_ptr<siconos::modeling::DynamicalSystem> dsA,
    std::shared_ptr<siconos::modeling::DynamicalSystem> dsB) const
{
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  std::shared_ptr<siconos::modeling::Interaction> inter;
  std::vector<std::shared_ptr<siconos::modeling::Interaction>> result;
  if (!dsA && !dsB) return result;
  for (std::tie(ui, uiend) = _IG[0]->vertices(); ui != uiend; ++ui) {
    inter = _IG[0]->bundle(*ui);
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds1 =
        _IG[0]->properties(_IG[0]->descriptor(inter)).source;
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds2 =
        _IG[0]->properties(_IG[0]->descriptor(inter)).target;
    int found = 0;
    if (dsA == ds1)
      found = 1;
    else if (dsA == ds2)
      found = 2;
    if (found == 2 && dsB != ds1)
      found = 0;
    else if (found == 1 && dsB == ds2)
      found = 0;
    if (found) result.push_back(inter);
  }
  return result;
}

std::vector<std::shared_ptr<siconos::modeling::DynamicalSystem>>
siconos::simulation::Topology::dynamicalSystemsForInteraction(
    std::shared_ptr<siconos::modeling::Interaction> inter) const
{
  std::shared_ptr<siconos::modeling::DynamicalSystem> ds1 =
      _IG[0]->properties(_IG[0]->descriptor(inter)).source;
  std::shared_ptr<siconos::modeling::DynamicalSystem> ds2 =
      _IG[0]->properties(_IG[0]->descriptor(inter)).target;
  std::vector<std::shared_ptr<siconos::modeling::DynamicalSystem>> result;
  if (ds1) result.push_back(ds1);
  if (ds2 && ds1 != ds2) result.push_back(ds2);
  return result;
}
