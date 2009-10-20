/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include "Topology.h"
#include "NonSmoothLaw.h"
#include "NonSmoothDynamicalSystem.h"
#include "Interaction.h"
#include "UnitaryRelation.h"
#include "EqualityConditionNSL.h"

#include <boost/bind.hpp>
#include <algorithm>
#include <limits>

#define MAX_RELATIVE_DEGREE 999

// --- CONSTRUCTORS/DESTRUCTOR ---

// default
Topology::Topology(): isTopologyUpToDate(false), isTopologyTimeInvariant(true),
  numberOfConstraints(0)
{
  URG.resize(1);
  DSG.resize(1);

  URG[0].reset(new UnitaryRelationsGraph());
  DSG[0].reset(new DynamicalSystemsGraph());
  allInteractions.reset(new InteractionsSet());
  maxRelativeDegree = 0;
  minRelativeDegree = MAX_RELATIVE_DEGREE;
}

Topology::Topology(SP::InteractionsSet newInteractions) :
  isTopologyUpToDate(false), isTopologyTimeInvariant(true),
  numberOfConstraints(0)
{

  URG.resize(1);
  DSG.resize(1);

  URG[0].reset(new UnitaryRelationsGraph());
  DSG[0].reset(new DynamicalSystemsGraph());
  allInteractions.reset(new InteractionsSet());

  for (InteractionsIterator it = newInteractions->begin();
       it != newInteractions->end(); ++it)
  {
    addInteraction(*it);
  }

  maxRelativeDegree = 0;
  minRelativeDegree = MAX_RELATIVE_DEGREE;
  isTopologyUpToDate = false;
}


// a constructor with a DS set : when some DS may not be in interactions
Topology::Topology(SP::DynamicalSystemsSet newDSset, SP::InteractionsSet newInteractions) :
  isTopologyUpToDate(false), isTopologyTimeInvariant(true),
  numberOfConstraints(0)
{

  URG.resize(1);
  DSG.resize(1);

  URG[0].reset(new UnitaryRelationsGraph());
  DSG[0].reset(new DynamicalSystemsGraph());
  allInteractions.reset(new InteractionsSet());

  for (InteractionsIterator it = newInteractions->begin();
       it != newInteractions->end(); ++it)
  {
    addInteraction(*it);
  }

  for (DSIterator ids = newDSset->begin(); ids != newDSset->end() ; ++ids)
  {
    DSG[0]->add_vertex(*ids);
  }
  maxRelativeDegree = 0;
  minRelativeDegree = MAX_RELATIVE_DEGREE;
  isTopologyUpToDate = false;
}

// destructor
Topology::~Topology()
{
  clear();
}

void Topology::addInteractionInIndexSet(SP::Interaction inter)
{
  // Private function
  //
  // Creates UnitaryRelations corresponding to inter and add them into
  // URG

  // First, we get the number of relations in the interaction.  This
  // corresponds to inter->getNumberOfRelations but since Interaction
  // has not been initialized yet, this value is not set and we need
  // to get interaction size and nsLaw size.
  unsigned int nsLawSize = inter->nonSmoothLaw()->size();
  unsigned int m = inter->getSizeOfY() / nsLawSize;
  unsigned int pos; // relative position of the relation in the y
  // vector of the Interaction
  UnitaryRelationsGraph::EDescriptor ur_current_edge;
  UnitaryRelationsGraph::VDescriptor ui_current_vertex;
  UnitaryRelationsGraph::EDescriptor ds_current_edge;

  SP::DynamicalSystemsSet systems = inter->dynamicalSystems();

  numberOfConstraints += m * nsLawSize;

  // DSG is the hyper forest : (vertices : dynamical systems, edges :
  // unitary relations)
  //
  // URG is the hyper graph : (vertices : unitary relations, edges :
  // dynamical systems)

  // URG = L(DSG),  L is the line graph transformation


  // for all couples of ds in the interaction
  for (DSIterator i1ds = systems->begin(); i1ds != systems->end(); ++i1ds)
  {
    for (DSIterator i2ds = i1ds; i2ds != systems->end(); ++i2ds)
    {

      // build the unitary relations
      std::vector<SP::UnitaryRelation> current_urs;
      for (unsigned int i = 0, pos = 0; i < m; ++i, pos += nsLawSize)
      {
        SP::UnitaryRelation UR(new UnitaryRelation(inter, pos, i));
        current_urs.push_back(UR);
      };

      // one DS in the interaction is a special case : a self branch
      if ((i1ds == i2ds) && inter->dynamicalSystems()->size() == 1)
      {
        DynamicalSystemsGraph::VDescriptor dsgv;
        dsgv = DSG[0]->add_vertex(*i1ds);

        // this may be a multi edges graph
        for (std::vector<SP::UnitaryRelation>::iterator uri = current_urs.begin();
             uri != current_urs.end(); ++uri)
        {
          assert(!DSG[0]->is_edge(dsgv, dsgv, *uri));
          assert(!URG[0]->is_vertex(*uri));

          DSG[0]->add_edge(dsgv, dsgv, *uri, *URG[0]);

          assert(URG[0]->is_vertex(*uri));
          assert(DSG[0]->is_edge(dsgv, dsgv, *uri));

        }
      }
      else
        // multiples DS in the interaction : no self branch
        if (i1ds != i2ds)
        {
          DynamicalSystemsGraph::VDescriptor dsgv1, dsgv2;

          dsgv1 = DSG[0]->add_vertex(*i1ds);
          dsgv2 = DSG[0]->add_vertex(*i2ds);

          // this may be a multi edges graph
          for (std::vector<SP::UnitaryRelation>::iterator uri = current_urs.begin();
               uri != current_urs.end(); ++uri)
          {

            assert(!DSG[0]->is_edge(dsgv1, dsgv2, *uri));
            assert(!URG[0]->is_vertex(*uri));


            DSG[0]->add_edge(dsgv1, dsgv2, *uri, *URG[0]);

            assert(URG[0]->is_vertex(*uri));
            assert(DSG[0]->is_edge(dsgv1, dsgv2, *uri));
          }
        }
    }
  }
};

/* an edge is removed from DSG graph if the corresponding vertex is
   removed from the adjoint graph (URG)
*/
struct VertexIsRemoved
{
  VertexIsRemoved(SP::Interaction I,
                  SP::DynamicalSystemsGraph sg, SP::UnitaryRelationsGraph asg) :
    _I(I), _DSG(sg), _URG(asg) {};
  bool operator()(DynamicalSystemsGraph::EDescriptor ed)
  {

    if (_URG->is_vertex(_DSG->bundle(ed)))
    {
      UnitaryRelationsGraph::VDescriptor uvd = _URG->descriptor(_DSG->bundle(ed));

      if (_URG->bundle(uvd)->interaction() == _I)
      {
        _URG->remove_vertex(_DSG->bundle(ed));

        assert(_URG->size() == _DSG->edges_number() - 1);

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
  SP::DynamicalSystemsGraph _DSG;
  SP::UnitaryRelationsGraph _URG;
};


/* remove an interaction : remove edges (unitary relation) from DSG if
   corresponding vertices are removed from URG */
const bool Topology::removeInteractionFromIndexSet(SP::Interaction inter)
{

  for (DSIterator ids = inter->dynamicalSystems()->begin();
       ids != inter->dynamicalSystems()->end();
       ++ids)
  {
    DSG[0]->remove_out_edge_if
    (DSG[0]->descriptor(*ids),
     VertexIsRemoved(inter, DSG[0], URG[0]));
  };
};


void Topology::addInteraction(SP::Interaction inter)
{
  assert(allInteractions);
  assert(DSG[0]->edges_number() == URG[0]->size());

  allInteractions->insert(inter);
  addInteractionInIndexSet(inter);

  assert(DSG[0]->edges_number() == URG[0]->size());

};

void Topology::removeInteraction(SP::Interaction inter)
{
  assert(allInteractions);
  assert(DSG[0]->edges_number() == URG[0]->size());

  allInteractions->erase(inter);
  removeInteractionFromIndexSet(inter);

  assert(DSG[0]->edges_number() == URG[0]->size());
};

void Topology::removeDynamicalSystem(SP::DynamicalSystem ds)
{
  RuntimeException::selfThrow("remove dynamical system not implemented");
};


/* a visitor to set parameters depending on nslaw */

class ComplementarityConditionNSL;
class MixedComplementarityConditionNSL;
class NewtonImpactNSL;
class NewtonImpactFrictionNSL;

struct Topology::SetupFromNslaw : public SiconosVisitor
{
  SP::Topology parent;
  SP::Interaction interaction;
  SetupFromNslaw(SP::Topology p, SP::Interaction inter) :
    parent(p), interaction(inter) {};

  void visit(ComplementarityConditionNSL&)
  {
    parent->minRelativeDegree = std::min<int>(0, parent->minRelativeDegree);
    parent->maxRelativeDegree = std::max<int>(0, parent->maxRelativeDegree);
    interaction->setRelativeDegree(0);
  };
  void visit(EqualityConditionNSL&)
  {
    parent->minRelativeDegree = std::min<int>(0, parent->minRelativeDegree);
    parent->maxRelativeDegree = std::max<int>(0, parent->maxRelativeDegree);
    interaction->setRelativeDegree(0);
  };

  void visit(MixedComplementarityConditionNSL&)
  {
    parent->minRelativeDegree = std::min<int>(0, parent->minRelativeDegree);
    parent->maxRelativeDegree = std::max<int>(0, parent->maxRelativeDegree);
    interaction->setRelativeDegree(0);
  };

  void visit(NewtonImpactNSL&)
  {
    parent->minRelativeDegree = std::min<int>(2, parent->minRelativeDegree);
    parent->maxRelativeDegree = std::max<int>(2, parent->maxRelativeDegree);
    parent->isTopologyTimeInvariant = false;
    interaction->setRelativeDegree(2);
  };

  void visit(NewtonImpactFrictionNSL&)
  {
    parent->minRelativeDegree = std::min<int>(2, parent->minRelativeDegree);
    parent->maxRelativeDegree = std::max<int>(2, parent->maxRelativeDegree);
    parent->isTopologyTimeInvariant = false;
    interaction->setRelativeDegree(2);
  };

};


// Compute min & max of relative degrees
void Topology::computeRelativeDegrees()
{

  boost::shared_ptr<SetupFromNslaw> setupFromNslaw;
  UnitaryRelationsGraph::VIterator uv, uend;


  if (URG[0]->size() > 0)
  {
    minRelativeDegree = MAX_RELATIVE_DEGREE;
    maxRelativeDegree = 0;


    for (boost::tie(uv, uend) = URG[0]->vertices(); uv != uend; ++uv)
    {

      setupFromNslaw.reset(new SetupFromNslaw(shared_from_this(),
                                              URG[0]->bundle(*uv)->interaction()));

      URG[0]->bundle(*uv)->interaction()->nonSmoothLaw()->accept(*(setupFromNslaw.get()));

    }
  }
  else
  {
    // default values
    minRelativeDegree = 2;
    maxRelativeDegree = 2;
  }

}



const bool Topology::hasInteraction(SP::Interaction inter) const
{
  return allInteractions->isIn(inter);
}

const unsigned int Topology::getMaxRelativeDegree()
{
  return maxRelativeDegree;
}

const unsigned int Topology::getMinRelativeDegree()
{
  assert(minRelativeDegree != MAX_RELATIVE_DEGREE);
  return minRelativeDegree;
}

void Topology::initialize()
{

  computeRelativeDegrees();
  isTopologyUpToDate = true;
}

void Topology::clear()
{
  allInteractions->clear();

  URG.clear();
  DSG.clear();

  isTopologyUpToDate = false;
}
