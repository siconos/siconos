/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "NonSmoothLaw.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Interaction.hpp"
#include "EqualityConditionNSL.hpp"


#include "MoreauJeanOSI.hpp"
#include "MoreauJeanGOSI.hpp"
#include "EulerMoreauOSI.hpp"
#include "SchatzmanPaoliOSI.hpp"

#define VISITOR_CLASSES() \
  REGISTER(MoreauJeanOSI) \
  REGISTER(MoreauJeanGOSI) \
  REGISTER(EulerMoreauOSI) \
  REGISTER(SchatzmanPaoliOSI)

#include <VisitorMaker.hpp>
using namespace Experimental;

// to be removed, once the mess with allOSI and OSIDynamicalSystems has been cleaned -- xhub
#include "OneStepIntegrator.hpp"

#if (__cplusplus >= 201103L) && !defined(USE_BOOST_FOR_CXX11)
#include <functional>
#else
#include <boost/bind.hpp>
#endif


#include <algorithm>
#include <limits>

//#define DEBUG_STDOUT
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
Topology::addInteractionInIndexSet0(SP::Interaction inter, SP::DynamicalSystem ds1, SP::DynamicalSystem ds2)
{
  // !! Private function !!
  //
  // Add inter and ds into IG/DSG

  // Compute number of constraints
  unsigned int nsLawSize = inter->nonSmoothLaw()->size();
  unsigned int m = inter->getSizeOfY() / nsLawSize;
  if (m > 1)
    RuntimeException::selfThrow("Topology::addInteractionInIndexSet0 - m > 1. Obsolete !");

  _numberOfConstraints += nsLawSize;

  SP::DynamicalSystem ds2_ = ds2;
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

  SP::VectorOfVectors workVds1 = _DSG[0]->properties(dsgv1).workVectors;

  SP::VectorOfVectors workVds2;
  if (!workVds1)
  {
    // V.A. 210/06/2017 Could we defer  this initialization ?
    workVds1.reset(new VectorOfVectors());
    _DSG[0]->properties(dsgv1).workMatrices.reset(new VectorOfMatrices());
    ds1->initializeWorkSpace(*workVds1, *_DSG[0]->properties(dsgv1).workMatrices);
  }
  if(ds2)
  {
    dsgv2 = _DSG[0]->add_vertex(ds2);
    workVds2 = _DSG[0]->properties(dsgv2).workVectors;
    if (!workVds2)
    {
      // V.A. 210/06/2017 Could we defer  this initialization ?
      workVds2.reset(new VectorOfVectors());
      _DSG[0]->properties(dsgv2).workMatrices.reset(new VectorOfMatrices());
      ds2->initializeWorkSpace(*workVds2, *_DSG[0]->properties(dsgv2).workMatrices);
    }
  }
  else
  {
    dsgv2 = dsgv1;
    ds2_ = ds1;
    workVds2 = workVds1;
  }
  // this may be a multi edges graph
  assert(!_DSG[0]->is_edge(dsgv1, dsgv2, inter));
  assert(!_IG[0]->is_vertex(inter));
  InteractionsGraph::VDescriptor ig_new_ve;
  DynamicalSystemsGraph::EDescriptor new_ed;
  std11::tie(new_ed, ig_new_ve) = _DSG[0]->add_edge(dsgv1, dsgv2, inter, *_IG[0]);
  InteractionProperties& interProp = _IG[0]->properties(ig_new_ve);

  // V.A. 210/06/2017 Could we defer  this initialization ?
  interProp.DSlink.reset(new VectorOfBlockVectors);
  interProp.workVectors.reset(new VectorOfVectors);
  interProp.workMatrices.reset(new VectorOfSMatrices);


  unsigned int nslawSize = inter->nonSmoothLaw()->size();
  interProp.block.reset(new SimpleMatrix(nslawSize, nslawSize));
  inter->setDSLinkAndWorkspace(interProp, *ds1, *workVds1, *ds2_, *workVds2);


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
    _IG[0]->properties(ig_new_ve).target_pos = ds1->dimension();
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
  DynamicalSystemsGraph::VDescriptor dsgv = _DSG[0]->add_vertex(ds);
  _DSG[0]->properties(dsgv).workVectors.reset(new VectorOfVectors());
  _DSG[0]->properties(dsgv).workMatrices.reset(new VectorOfMatrices());
  ds->initializeWorkSpace(*_DSG[0]->properties(dsgv).workVectors, *_DSG[0]->properties(dsgv).workMatrices);
}

void Topology::setName(SP::DynamicalSystem ds, const std::string& name)
{
  DynamicalSystemsGraph::VDescriptor dsgv = _DSG[0]->descriptor(ds);
  _DSG[0]->name.insert(dsgv, name);
}

void Topology::setName(SP::Interaction inter, const std::string& name)
{
  InteractionsGraph::VDescriptor igv = _IG[0]->descriptor(inter);
  _IG[0]->name.insert(igv, name);
}

/* initializeIterationMatrixW is not a member of OneStepIntegrator (should it be ?),
   so we visit some integrators which provide this initialization.
*/

/* first, a generic visitor is defined. */
struct CallInitDS : public SiconosVisitor
{
  double time;
  SP::DynamicalSystem ds;
  SP::Model m ; 
  template<typename T>
  void operator()(const T& osi)
  {
    const_cast<T*>(&osi)->initializeDynamicalSystem(*(this->m), this->time, this->ds);
  }
};

/* the visit is made on classes which provide the function initializeIterationMatrixW */
// typedef Visitor < Classes < MoreauJeanOSI >,
//                   CallInitDS >::Make InitDynamicalSystem;

typedef Visitor < Classes < MoreauJeanOSI,
                            MoreauJeanGOSI,
                            EulerMoreauOSI,
                            SchatzmanPaoliOSI >,
                  CallInitDS >::Make InitDynamicalSystem;

void Topology::initDS(SP::Model m, double time, SP::DynamicalSystem ds, SP::OneStepIntegrator OSI)
{
  InitDynamicalSystem initDynamicalSystem;
  initDynamicalSystem.time = 0;
  initDynamicalSystem.ds = ds;
  initDynamicalSystem.m = m;
  OSI->accept(initDynamicalSystem);
}

void Topology::setOSI(SP::DynamicalSystem ds, SP::OneStepIntegrator OSI)
{
  _DSG[0]->properties(_DSG[0]->descriptor(ds)).osi = OSI;
  ds->initMemory(OSI->getSizeMem());
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
  DEBUG_PRINTF("removeInteraction : %p\n", &*inter);

  assert(_DSG[0]->edges_number() == _IG[0]->size());
  removeInteractionFromIndexSet(inter);
  assert(_DSG[0]->edges_number() == _IG[0]->size());
}


std::pair<DynamicalSystemsGraph::EDescriptor, InteractionsGraph::VDescriptor>
Topology::link(SP::Interaction inter, SP::DynamicalSystem ds, SP::DynamicalSystem ds2)
{
  DEBUG_PRINTF("Topology::link : inter %p, ds1 %p, ds2 %p\n", &*inter, &*ds,
               &*ds2);
  if (indexSet0()->is_vertex(inter))
  {
    removeInteractionFromIndexSet(inter);
  }

  unsigned int sumOfDSSizes = 0, sumOfZSizes = 0;

  sumOfDSSizes += ds->dimension();
  if(ds->z())
    sumOfZSizes += ds->z()->size();

  if(ds2)
  {
    sumOfDSSizes += ds2->dimension();
    if(ds->z())
      sumOfZSizes += ds2->z()->size();
    inter->setHas2Bodies(true);
  }
  DEBUG_PRINTF("sumOfDSSizes = %i\t, sumOfZSizes = %i\n ", sumOfDSSizes, sumOfZSizes);

  inter->setDSSizes(sumOfDSSizes, sumOfZSizes);

  return addInteractionInIndexSet0(inter, ds, ds2);
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


SP::DynamicalSystem Topology::getDynamicalSystem(std::string name)
{
  DynamicalSystemsGraph::VIterator vi, vdend;
  for (std11::tie(vi, vdend) = _DSG[0]->vertices(); vi != vdend; ++vi)
  {
    if (name == _DSG[0]->name.at(*vi))
      return _DSG[0]->bundle(*vi);
  }

  RuntimeException::selfThrow("Topology::getDynamicalSystem() ds not found.");

  return SP::DynamicalSystem();
}


unsigned int Topology::numberOfInvolvedDS(unsigned int inumber)
{
  if (inumber >= _IG.size())
  {
    RuntimeException::selfThrow("Topology::numberOfInvolvedDS :: index number must be smaller than the number of indexSets");
  }

  /* on an adjoint graph a dynamical system may be on several edges */
  std::map<SP::DynamicalSystem, bool> flag;

  unsigned int return_value = 0;

  SP::InteractionsGraph indexSet = _IG[inumber];

  InteractionsGraph::VIterator vi, viend;
  for(std11::tie(vi, viend) = indexSet->vertices();
      vi != viend; ++vi)
  {
    if(indexSet->properties(*vi).source)
    {
      if (flag.find(indexSet->properties(*vi).source) == flag.end())
      {
        flag[indexSet->properties(*vi).source] = true;
        return_value++;
      }
    }
    if(indexSet->properties(*vi).target)
    {
      if (flag.find(indexSet->properties(*vi).target) == flag.end())
      {
        flag[indexSet->properties(*vi).target] = true;
        return_value++;
      }
    }
  }

  InteractionsGraph::EIterator ei, eiend;

  for(std11::tie(ei, eiend) = indexSet->edges();
      ei != eiend; ++ei)
  {
    if (flag.find(indexSet->bundle(*ei)) == flag.end())
    {
      flag[indexSet->bundle(*ei)] = true;
      return_value++;
    }
  }

  return return_value;
}

SP::Interaction Topology::getInteraction(unsigned int requiredNumber)
{
  InteractionsGraph::VIterator vi, vdend;
  SP::Interaction inter;
  unsigned int currentNumber;
  for (std11::tie(vi, vdend) = _IG[0]->vertices(); vi != vdend; ++vi)
  {
    inter = _IG[0]->bundle(*vi);
    currentNumber = inter->number();
    if (currentNumber == requiredNumber)
      return inter;
  }

  RuntimeException::selfThrow("Topology::getInteractiob(n) inter not found.");

  return inter;
}
