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

/*! \file SiconosGraphs.hpp
 * \brief Definitions of the graphs used in Siconos
 */

#ifndef SimulationGraphs_H
#define SimulationGraphs_H

#include "SiconosGraph.hpp"
#include "SiconosProperties.hpp"

namespace siconos::algebra {

class SiconosMatrix;
class BlockVector;
class SimpleMatrix;
class SiconosVector;
}  // namespace siconos::algebra

namespace siconos::plugins {
class PluggedObject;
}

namespace siconos::modeling {

class DynamicalSystem;
class Interaction;
}  // namespace siconos::modeling

namespace siconos::simulation {

class MatrixIntegrator;
}  // namespace siconos::simulation

namespace siconos::integrators {

class OneStepIntegrator;
}  // namespace siconos::integrators

namespace siconos::graphs {

/** the graph structure :
 *
 * InteractionsGraph = L(DynamicalSystemsGraph)
 *
 * where L is the line graph
 * transformation
 *
 *
 * Properties on graph :
 * --------------------
 *
 * The properties on the graph enable to store the data that are specific to a simulation
 * strategy. It avoids to burden the modeling classes that should be as independent as possible
 * from the simulation choices.
 *
 * There are mainly two types of properties
 * <ul>
 * <li>  Mandatory properties DynamicalSystemProperties and InteractionProperties .
 *       These properties are always  instanciated for any kind of simulation.
 *       The accessors to the property are illustrated in the following example :
 *       For a given std::shared_ptr<DynamicalSystem> ds and a given graph
 * std::shared_ptr<DynamicalSystemsGraph> DSG
 *
 *       DynamicalSystemsGraph::VDescriptor dsv = DSG->descriptor(ds);
 *       auto osi = DSG->properties(dsv).osi;
 * </li>
 * <li> Optional Properties
 *      They are installed thanks to the macro INSTALL_GRAPH_PROPERTIES.
 *
 *      The accessors to the property are illustrated in the following example :
 *      For a given std::shared_ptr<DynamicalSystem> ds and a given graph
 * std::shared_ptr<DynamicalSystemsGraph> DSG
 *
 *      DynamicalSystemsGraph::VDescriptor dsv = DSG->descriptor(ds);
 *      DSG->name.insert(dsv, name); // insert the name in the property
 *      const std::string& name =  DSG[*dsv]->name;
 *
 *
 * </li>
 * </ul>
 */

/** \struct InteractionProperties mandatory properties for an Interaction  */
struct InteractionProperties {
  std::shared_ptr<siconos::algebra::SiconosMatrix> block{nullptr}; /**< diagonal block */
  std::shared_ptr<siconos::modeling::DynamicalSystem> source{nullptr};
  unsigned int source_pos;
  std::shared_ptr<siconos::modeling::DynamicalSystem> target{nullptr};
  unsigned int target_pos{0};
  unsigned int absolute_position{
      0}; /**< Absolute position of the interaction variables in the unknown vector in osnsp*/
  unsigned int absolute_position_proj{0}; /**< Absolute position of the interaction variables
                                          in the unknown vector in osnsp for projection*/
  bool forControl{false}; /**< true if the relation is used to add a control input to a DS */
  std::shared_ptr<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>
      workVectors; /**< set of SiconosVector, useful to ensure contiguous memory vectors, used
                      as buffers in OneStepIntegrator classes. */
  std::shared_ptr<std::vector<std::shared_ptr<siconos::algebra::BlockVector>>>
      workBlockVectors{
          nullptr}; /**< set of BlockVector, used as buffers in OneStepIntegrator classes. */
  std::shared_ptr<std::vector<std::shared_ptr<siconos::algebra::SimpleMatrix>>> workMatrices{
      nullptr}; /**< Internal buffers used on simulation size, to store
jacobians or other temporary matrices. */
  std::shared_ptr<siconos::integrators::OneStepIntegrator> osi1{
      nullptr}; /**< Integrator 1  used for the given Interaction */
  std::shared_ptr<siconos::integrators::OneStepIntegrator> osi2{
      nullptr}; /**< Integrator 2  used for the given Interaction */

  ACCEPT_SERIALIZATION(InteractionProperties);
};

/** \struct DynamicalSystemProperties mandatory properties for a DynamicalSystems */
struct DynamicalSystemProperties {
  std::shared_ptr<siconos::algebra::SiconosMatrix> upper_block{nullptr}; /**< i,j block i<j */
  std::shared_ptr<siconos::algebra::SiconosMatrix> lower_block{nullptr}; /**< i,j block i>j */
  std::shared_ptr<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>> workVectors{
      nullptr}; /**< Used for instance in Newton iteration */
  std::shared_ptr<std::vector<std::shared_ptr<siconos::algebra::SiconosMatrix>>> workMatrices{
      nullptr}; /**< Mostly for Lagrangian system.*/
  std::shared_ptr<siconos::integrators::OneStepIntegrator> osi{
      nullptr}; /**< Integrator used for the given DynamicalSystem */
  std::shared_ptr<siconos::algebra::SimpleMatrix> W{nullptr}; /**< Matrix for integration */
  std::shared_ptr<siconos::algebra::SimpleMatrix> WBoundaryConditions{
      nullptr}; /**< Matrix for integration of boundary conditions*/
  std::shared_ptr<siconos::algebra::SimpleMatrix> Winverse{
      nullptr};                      /**< Matrix for integration */
  unsigned int absolute_position{0}; /**< Absolute position of the ds variables in the unknown
                                     vector in osnsp*/
  //  std::shared_ptr<siconos::algebra::SiconosMemory> _xMemory            /**< old value of x,
  //  TBD */

  ACCEPT_SERIALIZATION(DynamicalSystemProperties);
};

// Note FP : workMatrices in DSProperties is used only in  NewmarkAlphaOSI => maybe it should
// be replaced with interprop workMat?

struct GraphProperties {
  bool symmetric{false};

  ACCEPT_SERIALIZATION(GraphProperties);
};

class _DynamicalSystemsGraph
    : public SiconosGraph<std::shared_ptr<siconos::modeling::DynamicalSystem>,
                          std::shared_ptr<siconos::modeling::Interaction>,
                          DynamicalSystemProperties, InteractionProperties, GraphProperties> {
  ACCEPT_SERIALIZATION(_DynamicalSystemsGraph);
};

class _InteractionsGraph
    : public SiconosGraph<std::shared_ptr<siconos::modeling::Interaction>,
                          std::shared_ptr<siconos::modeling::DynamicalSystem>,
                          InteractionProperties, DynamicalSystemProperties, GraphProperties> {
  ACCEPT_SERIALIZATION(_InteractionsGraph);
};

struct DynamicalSystemsGraph : public _DynamicalSystemsGraph {
  /** optional properties : memory is allocated only on first access */
  INSTALL_GRAPH_PROPERTIES(
      DynamicalSystems, ((siconos::graphs::VertexSP, siconos::simulation::MatrixIntegrator,
                          Ad))  // for ZOH Integration
      ((siconos::graphs::VertexSP, siconos::simulation::MatrixIntegrator,
        AdInt))  // for ZOH Integration
      ((siconos::graphs::VertexSP, siconos::simulation::MatrixIntegrator,
        Ld))  // For Observer (ZOH Integration)
      ((siconos::graphs::VertexSP, siconos::simulation::MatrixIntegrator,
        Bd))  // For Controlled System (ZOH Integration)
      ((siconos::graphs::VertexSP, siconos::algebra::SiconosMatrix,
        B))  // For Controlled System
      ((siconos::graphs::VertexSP, siconos::algebra::SiconosMatrix, L))  // For Observer
      ((siconos::graphs::VertexSP, siconos::plugins::PluggedObject,
        pluginB))  // For Controlled System
      ((siconos::graphs::VertexSP, siconos::plugins::PluggedObject,
        pluginL))                                                        // For Observer
      ((siconos::graphs::VertexSP, siconos::algebra::SiconosVector, e))  // For Observer
      ((siconos::graphs::VertexSP, siconos::algebra::SiconosVector,
        u))  // For Controlled System
      ((siconos::graphs::VertexSP, siconos::plugins::PluggedObject,
        pluginU))  // For Controlled System (nonlinear w.r.t u)
      ((siconos::graphs::VertexSP, siconos::plugins::PluggedObject,
        pluginJacgx))  // For Controlled System (nonlinear w.r.t u); compute nabla_x g(x, u)
      ((siconos::graphs::VertexSP, siconos::algebra::SiconosVector,
        tmpXdot))  // For Controlled System (nonlinear w.r.t u); tmpXdot = g(x, u)
      ((siconos::graphs::VertexSP, siconos::algebra::SimpleMatrix,
        jacgx))  // For Controlled System (nonlinear w.r.t u); jacgx = nabla_x g(x, u)
      ((Vertex, std::string, name))        // a name for a dynamical system
      ((Vertex, unsigned int, groupId)));  // For group manipulations (example assign
                                           // a material id for contact law
                                           // determination
  // always needed -> DynamicalSystemProperties

  ACCEPT_SERIALIZATION(DynamicalSystemsGraph);

  // to be installed with INSTALL_GRAPH_PROPERTIES
  void eraseProperties(_DynamicalSystemsGraph::VDescriptor vd) {
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

struct InteractionsGraph : public _InteractionsGraph {
  /** optional properties : memory is allocated only on first access */
  INSTALL_GRAPH_PROPERTIES(
      Interactions, ((siconos::graphs::Vertex, std::shared_ptr<siconos::algebra::SimpleMatrix>,
                      blockProj))  // ProjectOnConstraint
      ((siconos::graphs::Edge, std::shared_ptr<siconos::algebra::SimpleMatrix>,
        upper_blockProj))  // idem
      ((siconos::graphs::Edge, std::shared_ptr<siconos::algebra::SimpleMatrix>,
        lower_blockProj))  // idem
      ((siconos::graphs::Vertex, std::string, name)));

  // to be installed with INSTALL_GRAPH_PROPERTIES
  void eraseProperties(_InteractionsGraph::VDescriptor vd) {
    blockProj._store->erase(vd);
    name._store->erase(vd);
  }

  // to be installed with INSTALL_GRAPH_PROPERTIES
  void eraseProperties(_InteractionsGraph::EDescriptor ed) {
    upper_blockProj._store->erase(ed);
    lower_blockProj._store->erase(ed);
  }

  ACCEPT_SERIALIZATION(InteractionsGraph);
};
}  // namespace siconos::graphs
#endif
