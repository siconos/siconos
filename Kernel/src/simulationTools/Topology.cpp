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

#include <boost/bind.hpp>

const bool Topology::addInteractionInIndexSet(SP::Interaction inter)
{
  // Private function
  //
  // Creates UnitaryRelations corresponding to inter and add them into
  // URG

  // First, we get the number of relations in the interaction.  This
  // corresponds to inter->getNumberOfRelations but since Interaction
  // has not been initialized yet, this value is not set and we need
  // to get interaction size and nsLaw size.
  unsigned int nsLawSize = inter->getNonSmoothLawPtr()->getNsLawSize();
  unsigned int m = inter->getSizeOfY() / nsLawSize;
  unsigned int pos; // relative position of the relation in the y
  // vector of the Interaction
  UnitaryRelationsGraph::EDescriptor ur_current_edge;
  UnitaryRelationsGraph::VDescriptor ui_current_vertex;
  UnitaryRelationsGraph::EDescriptor ds_current_edge;

  SP::DynamicalSystemsSet systems = inter->getDynamicalSystemsPtr();

  bool inserted;

  bool res = true; // output value. False if insertion of one of the
  // relations fails.


  std::vector<SP::UnitaryRelation> current_urs;

  for (unsigned int i = 0, pos = 0; i < m; ++i, pos += nsLawSize)
  {
    SP::UnitaryRelation UR(new UnitaryRelation(inter, pos, i));
    current_urs.push_back(UR);
  };

  numberOfConstraints += m * nsLawSize;

  // the hypergraph (set of DS, Interaction) is made of a clique in
  // URG
  for (DSIterator i1ds = systems->begin(); i1ds != systems->end(); ++i1ds)
  {
    for (DSIterator i2ds = i1ds; i2ds != systems->end(); ++i2ds)
    {
      // one DS is a special case : a self branch
      if ((i1ds == i2ds) && inter->getDynamicalSystemsPtr()->size() == 1)
      {
        DynamicalSystemsGraph::VDescriptor dsgv;
        dsgv = DSG[0]->add_vertex(*i1ds);
        DSG[0]->color(dsgv) = boost::white_color;

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
      else if (i1ds != i2ds)
      {
        DynamicalSystemsGraph::VDescriptor dsgv1, dsgv2;

        dsgv1 = DSG[0]->add_vertex(*i1ds);
        dsgv2 = DSG[0]->add_vertex(*i2ds);

        DSG[0]->color(dsgv1) = boost::white_color;
        DSG[0]->color(dsgv2) = boost::white_color;


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
  return (res);
};

// an edge is removed if the corresponding vertex is removed in
// adjoint graph
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

      if (_URG->bundle(uvd)->getInteractionPtr() == _I)
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


const bool Topology::removeInteractionFromIndexSet(SP::Interaction inter)
{

  /* remove unitary relation in ur graph */



  for (DSIterator ids = inter->getDynamicalSystemsPtr()->begin();
       ids != inter->getDynamicalSystemsPtr()->end();
       ++ids)
  {
    DSG[0]->remove_out_edge_if
    (DSG[0]->descriptor(*ids),
     VertexIsRemoved(inter, DSG[0], URG[0]));
  };
};


void Topology::addInteraction(SP::Interaction inter)
{
  assert(DSG[0]->edges_number() == URG[0]->size());

  allInteractions->insert(inter);
  addInteractionInIndexSet(inter);

  assert(DSG[0]->edges_number() == URG[0]->size());

};

void Topology::removeInteraction(SP::Interaction inter)
{
  assert(DSG[0]->edges_number() == URG[0]->size());

  allInteractions->erase(inter);
  removeInteractionFromIndexSet(inter);

  assert(DSG[0]->edges_number() == URG[0]->size());
};


// Compute relative degrees map
void Topology::computeRelativeDegrees()
{
  // for each Unitary Relation relative degree vector depends on
  // NonSmooth Law and relations
  relativeDegrees.clear();
  std::string nslawType;

  // loop through URG
  UnitaryRelationsGraph::VIterator uv, uend;
  for (boost::tie(uv, uend) = URG[0]->vertices(); uv != uend; ++uv)
  {
    nslawType = URG[0]->bundle(*uv)->getNonSmoothLawType();
    if (nslawType == COMPLEMENTARITYCONDITIONNSLAW ||
        nslawType == MIXEDCOMPLEMENTARITYCONDITIONNSLAW)
      relativeDegrees[URG[0]->bundle(*uv)] = 0;

    else if (nslawType == NEWTONIMPACTNSLAW)
    {
      relativeDegrees[URG[0]->bundle(*uv)] = 2;
      isTopologyTimeInvariant = false;
    }
    else if (nslawType == NEWTONIMPACTFRICTIONNSLAW)
    {
      relativeDegrees[URG[0]->bundle(*uv)] = 2;
      isTopologyTimeInvariant = false;
    }
    else
      RuntimeException::selfThrow("Topology::computeRelativeDegree(...), not yet implemented for non smooth law of type" + nslawType);
  }

  // assert(!relativeDegrees.empty());
}

// --- CONSTRUCTORS/DESTRUCTOR ---

// default
Topology::Topology(): isTopologyUpToDate(false), isTopologyTimeInvariant(true),
  numberOfConstraints(0)
{}

Topology::Topology(SP::InteractionsSet newInteractions) :
  isTopologyUpToDate(false), isTopologyTimeInvariant(true),
  numberOfConstraints(0)
{
  if (URG.size() == 0)
  {
    URG.resize(1);
    URG[0].reset(new UnitaryRelationsGraph());
  }
  if (DSG.size() == 0)
  {
    DSG.resize(1);
    DSG[0].reset(new DynamicalSystemsGraph());
  }

  allInteractions = newInteractions;

  isTopologyUpToDate = false;
}

// destructor
Topology::~Topology()
{
  clear();
}

const bool Topology::hasInteraction(SP::Interaction inter) const
{
  return allInteractions->isIn(inter);
}

const unsigned int Topology::getMaxRelativeDegree()
{
  if (relativeDegrees.empty())
    return 2;
  //    RuntimeException::selfThrow("Topology::getMaxRelativeDegree, non-existent value, since the relative degrees map is empty.");

  ConstIteratorForRelativeDegrees it = max_element(relativeDegrees.begin(), relativeDegrees.end());
  return(it->second);
}

const unsigned int Topology::getMinRelativeDegree()
{
  if (relativeDegrees.empty())
    return 2;
  // RunimeException::selfThrow("Topology::getMinRelativeDegree, non-existent value, since the relative degrees map is empty.");

  ConstIteratorForRelativeDegrees it = min_element(relativeDegrees.begin(), relativeDegrees.end());
  return(it->second);
}

void Topology::initialize()
{

  assert(allInteractions && "Topology : allInteractions is NULL");
  //  assert(!allInteractions->isEmpty());

  // -- Creates Unitary Relations and put them in URG --- loop
  // through interactions list (from NSDS)
  for (InteractionsIterator it = allInteractions->begin();
       it != allInteractions->end(); ++it)
    addInteractionInIndexSet(*it);

  //-- Fills RelativeDegreesMaps in --
  computeRelativeDegrees();

  isTopologyUpToDate = true;
}

void Topology::clear()
{
  allInteractions->clear();

  URG.clear();
  DSG.clear();
  relativeDegrees.clear();

  isTopologyUpToDate = false;
}
