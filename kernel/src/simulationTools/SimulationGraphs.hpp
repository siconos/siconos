/* Siconos-Kernel, Copyright INRIA 2005-2015
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

/*! \file SiconosGraphs.hpp
 * \brief Definitions of the graphs used in Siconos
 */

#ifndef SimulationGraphs_H
#define SimulationGraphs_H

#include "SiconosGraph.hpp"
#include "SiconosProperties.hpp"
#include "SiconosPointers.hpp"
#include "SimulationTypeDef.hpp"

/** the graph structure :
 *
 * InteractionsGraph = L(DynamicalSystemsGraph)
 *
 * where L is the line graph
 * transformation */
typedef SiconosGraph < std11::shared_ptr<DynamicalSystem>, std11::shared_ptr<Interaction>,
        SystemProperties, InteractionProperties,
        GraphProperties > _DynamicalSystemsGraph;


typedef SiconosGraph < std11::shared_ptr<Interaction>, std11::shared_ptr<DynamicalSystem>,
        InteractionProperties, SystemProperties,
        GraphProperties > _InteractionsGraph;

struct DynamicalSystemsGraph : public _DynamicalSystemsGraph
{
  /** optional properties : memory is allocated only on first access */
  INSTALL_GRAPH_PROPERTIES(DynamicalSystems,
                           ((VertexSP, MatrixIntegrator, Ad)) // for ZOH Integration
                           ((VertexSP, MatrixIntegrator, AdInt)) // for ZOH Integration
                           ((VertexSP, MatrixIntegrator, Ld)) // For Observer (ZOH Integration)
                           ((VertexSP, MatrixIntegrator, Bd)) // For Controlled System (ZOH Integration)
                           ((VertexSP, SiconosMatrix, B)) // For Controlled System
                           ((VertexSP, SiconosMatrix, L)) // For Observer
                           ((VertexSP, PluggedObject, pluginB)) // For Controlled System
                           ((VertexSP, PluggedObject, pluginL)) // For Observer
                           ((VertexSP, SiconosVector, e)) // For Observer
                           ((VertexSP, SiconosVector, u)) // For Controlled System
                           ((VertexSP, PluggedObject, pluginU)) // For Controlled System (nonlinear w.r.t u)
                           ((VertexSP, PluggedObject, pluginJacgx)) // For Controlled System (nonlinear w.r.t u); compute nabla_x g(x, u)
                           ((VertexSP, SiconosVector, tmpXdot)) // For Controlled System (nonlinear w.r.t u); tmpXdot = g(x, u)
                           ((VertexSP, SimpleMatrix, jacgx)) // For Controlled System (nonlinear w.r.t u); jacgx = nabla_x g(x, u)
                           ((Vertex, std::string, name)) // a name for a dynamical system
                           ((Vertex, unsigned int, groupId))); // For group manipulations (example assign
                                                               // a material id for contact law
                                                               // determination
  // always needed -> SystemProperties

  /** serialization hooks */
  ACCEPT_SERIALIZATION(DynamicalSystemsGraph);

  // to be installed with INSTALL_GRAPH_PROPERTIES
  void eraseProperties(_DynamicalSystemsGraph::VDescriptor vd)
  {
    Ad._store->erase(vd);
    AdInt._store->erase(vd);
    Ld._store->erase(vd);
    Bd._store->erase(vd);
    B._store->erase(vd);
    L._store->erase(vd);
    pluginB._store->erase(vd);
    pluginL._store->erase(vd);
    e._store->erase(vd);
    u._store->erase(vd);
    pluginU._store->erase(vd);
    pluginJacgx._store->erase(vd);
    tmpXdot._store->erase(vd);
    jacgx._store->erase(vd);
    name._store->erase(vd);
    groupId._store->erase(vd);
  }
};

struct InteractionsGraph : public _InteractionsGraph
{
  /** optional properties : memory is allocated only on first access */
  INSTALL_GRAPH_PROPERTIES(Interactions,
                           ((Vertex, SP::SimpleMatrix, blockProj))        // ProjectOnConstraint
                           ((Edge, SP::SimpleMatrix, upper_blockProj))    // idem
                           ((Edge, SP::SimpleMatrix, lower_blockProj))  // idem
                           ((Vertex, std::string, name)));

  // to be installed with INSTALL_GRAPH_PROPERTIES
  void eraseProperties(_InteractionsGraph::VDescriptor vd)
  {
    blockProj._store->erase(vd);
    name._store->erase(vd);
  }

  // to be installed with INSTALL_GRAPH_PROPERTIES
  void eraseProperties(_InteractionsGraph::EDescriptor ed)
  {
    upper_blockProj._store->erase(ed);
    lower_blockProj._store->erase(ed);
  }

  /** serialization hooks */
  ACCEPT_SERIALIZATION(InteractionsGraph);
};

#endif
