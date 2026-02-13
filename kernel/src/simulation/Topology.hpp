/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*! \file Topology.hpp
  \brief Describes the topology of the DynamicalSystem and Interaction in the Simulation
*/
#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <memory>  // shared_from_this
#include <vector>

#include "SiconosSerialization.hpp"
#include "SimulationGraphs.hpp"

namespace siconos::simulation {

/**
   This class describes the topology of the non-smooth dynamical
   system. It holds all the "potential" Interactions".

   Topology is built in NSDS constructors but initialized in
   Simulation->initialize(), ie when all Interactions have been
   clearly defined.

   Note that indexSet0 holds all the possible relations (declared by
   user) not only those which are "actives".

   Construction consists in:
   - link with the NSDS that owns the topology.

   Initialization consists in:
   - scan of all the interactions of the NSDS
   - initialization of each interaction
   - insertion of the relations of all the Interaction into indexSet0

   Insertion of an Interaction into the set indexSet0:
   addInteractionInIndexSet0(std::shared_ptr<siconos::modeling::Interaction> inter) for each
   relation of the interaction, it creates a new Interaction and inserts it into indexSet0 It
   also counts the total number of "constraints" in the system.
*/
class Topology : public std::enable_shared_from_this<Topology> {
 private:
  ACCEPT_SERIALIZATION(Topology);

  /** dynamical systems graphs */
  std::vector<std::shared_ptr<siconos::graphs::DynamicalSystemsGraph>> _DSG;

  /** Interaction graphs (_IG[0]=L[DSG[0]], L is the line graph
      transformation) */
  std::vector<std::shared_ptr<siconos::graphs::InteractionsGraph>> _IG;

  /** check if topology is static or  not */
  bool _hasChanged = true;

  /** Total number of (scalar) constraints in the problem, ie sum of
      all nslaw sizes of Interactions of IndexSet0.*/
  unsigned int _numberOfConstraints = 0;

  /** symmetry in the blocks computation */
  bool _symmetric = false;

  /** initializations ( time invariance) from non
      smooth laws kind */
  struct SetupFromNslaw;
  friend struct Topology::SetupFromNslaw;

  // === PRIVATE FUNCTIONS ===

  /** schedules the relations of Interaction inter into IndexSet0 (ie
   *  creates the corresponding Interactions and add them into _DSG
   *  and _IG)
   *
   *  \param  inter the Interaction to add
   *  \param ds1 a std::shared_ptr<siconos::modeling::DynamicalSystem> first dynamical system
   * linked to the interaction \param ds2 a std::shared_ptr<siconos::modeling::DynamicalSystem>
   * second dynamical system linked to the interaction (default = None) \return a vertex
   * descriptor of the new vertex in IndexSet0
   */
  std::pair<siconos::graphs::DynamicalSystemsGraph::EDescriptor,
            siconos::graphs::InteractionsGraph::VDescriptor>
  __addInteractionInIndexSet0(std::shared_ptr<siconos::modeling::Interaction> inter,
                              std::shared_ptr<siconos::modeling::DynamicalSystem> ds1,
                              std::shared_ptr<siconos::modeling::DynamicalSystem> ds2 =
                                  std::shared_ptr<siconos::modeling::DynamicalSystem>());

  /** remove an Interaction from _IG and _DSG
   *
   *  \param inter a pointer to the Interaction to be removed
   */
  void __removeInteractionFromIndexSet(std::shared_ptr<siconos::modeling::Interaction> inter);

  /** remove a siconos::modeling::DynamicalSystem from _IG and _DSG
   *
   *  \param ds a pointer to the Dynamical System to be removed
   */
  void __removeDynamicalSystemFromIndexSet(
      std::shared_ptr<siconos::modeling::DynamicalSystem> ds);

  // Rules of five ...
  Topology(const Topology&) = delete;
  Topology(Topology&&) = delete;
  Topology& operator=(const Topology&) = delete;
  Topology& operator=(Topology&&) = delete;

 public:
  // --- CONSTRUCTORS/DESTRUCTOR ---

  /** default constructor */
  Topology();

  /** destructor */
  ~Topology() noexcept;

  // === GETTERS/SETTERS ===

  /** check if an dynamical system is already a vertex of the DSs graph.
   *
   *  \param ds the DS to test
   *  \return true if ds is in the graph
   */
  bool hasDynamicalSystem(std::shared_ptr<siconos::modeling::DynamicalSystem> ds) const;

  /** check if an interaction is already a vertex of the Interactions graph.
   *
   *  \param inter the Interaction to test
   *  \return true if inter is in the graph
   */
  bool hasInteraction(std::shared_ptr<siconos::modeling::Interaction> inter) const;

  /** remove an Interaction from the topology. The interaction is
   *  removed from Dynamical Systems graph and Interactions Graph.
   *  The interaction is not removed from actives subgraphs : see updateIndexSet
   *
   *  \param inter the interaction to remove
   */
  void removeInteraction(std::shared_ptr<siconos::modeling::Interaction> inter);

  /** add a dynamical system
   *
   *  \param ds the DynamicalSystem to add
   */
  void insertDynamicalSystem(std::shared_ptr<siconos::modeling::DynamicalSystem> ds);

  /** remove a Dynamical System from the topology. The dynamical
   *  system is removed from Dynamical Systems graph and Interactions
   *  Graph.  The dynamical system is not removed from actives
   *  subgraphs : see updateIndexSet
   *
   *  \param ds the dynamical system to remove
   */
  void removeDynamicalSystem(std::shared_ptr<siconos::modeling::DynamicalSystem> ds);

  /** set the name for this Dynamical System
   *
   *  \param ds a pointer to the system
   *  \param name the name of the DynamicalSystem
   */
  void setName(std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
               const std::string& name);

  /** get the name for this Dynamical System
   *
   *  \param ds a pointer to the system
   *  \return name the name of the DynamicalSystem, or empty string if not found.
   */
  std::string name(std::shared_ptr<siconos::modeling::DynamicalSystem> ds);

  /** set the name for an Interaction
   *
   *  \param inter a pointer to the Interaction
   *  \param name the name of the Interaction
   */
  void setName(std::shared_ptr<siconos::modeling::Interaction> inter, const std::string& name);

  /** get the name for this Interaction
   *
   *  \param inter a pointer to the Interaction
   *  \return name the name of the Interaction, or empty string if not found.
   */
  std::string name(std::shared_ptr<siconos::modeling::Interaction> inter);

  /** set the OSI for this DynamicalSystem
   *
   *  \param ds the DynamicalSystem
   *  \param OSI the integrator to use for this DS
   */
  void setOSI(std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
              std::shared_ptr<siconos::integrators::OneStepIntegrator> OSI);

  /** link two dynamical systems to a relation
   *
   *  \param inter a std::shared_ptr<siconos::modeling::Interaction>
   *  \param ds a std::shared_ptr<DynamicalSystem>
   *  \param ds2 a std::shared_ptr<DynamicalSystem> (optional)
   *  \return a vertex descriptor of the new vertex in IndexSet0
   */
  std::pair<siconos::graphs::DynamicalSystemsGraph::EDescriptor,
            siconos::graphs::InteractionsGraph::VDescriptor>
  link(std::shared_ptr<siconos::modeling::Interaction> inter,
       std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
       std::shared_ptr<siconos::modeling::DynamicalSystem> ds2 =
           std::shared_ptr<siconos::modeling::DynamicalSystem>());

  /** specify if the given Interaction is for controlling the DS
   *
   *  \param inter Interaction
   *  \param isControlInteraction true if the Interaction is used for
   *  control purposes
   **/
  void setControlProperty(std::shared_ptr<siconos::modeling::Interaction> inter,
                          const bool isControlInteraction);

  /** get a pointer to the graph of all Interactions.
   *
   *  \return a std::shared_ptr<siconos::graphs::InteractionsGraph>
   */
  inline std::shared_ptr<siconos::graphs::InteractionsGraph> indexSet0() const
  {
    return _IG[0];
  }

  std::shared_ptr<siconos::graphs::InteractionProperties> interaction_properties(
      unsigned int index, std::shared_ptr<siconos::modeling::Interaction> inter)
  {
    auto ui = indexSet(index)->descriptor(inter);
    return std::make_shared<siconos::graphs::InteractionProperties>(
        indexSet(index)->properties(ui));
  };

  /** get a pointer to the graph at level num of Interactions
   *
   *  \param num the number of indexSet
   *  \return a std::shared_ptr<siconos::graphs::InteractionsGraph>
   */
  std::shared_ptr<siconos::graphs::InteractionsGraph> indexSet(unsigned int num) const;

  /** get a pointer to the graph at level num of Interactions
   *
   *  \return a std::shared_ptr<siconos::graphs::InteractionsGraph>
   */
  inline size_t numberOfIndexSet() const { return _IG.size(); };

  /** reset graph at level num of Interactions
   *
   *  \param num the indexSet to reset
   */
  inline void resetIndexSetPtr(unsigned int num)
  {
    assert(num < _IG.size());

    // .. global properties may be defined here with
    // InteractionsSubGraphProperties(), see SiconosProperties.hpp
    // VertexSubProperties or EdgeSubProperties and the macros
    // INSTALL_GRAPH_PROPERTIES

    _IG[num] = std::make_shared<siconos::graphs::InteractionsGraph>();

    _IG[num]->properties().symmetric = _symmetric;
    _IG[num]->update_vertices_indices();
    _IG[num]->update_edges_indices();
  };

  /** get a pointer to the graph at level num of Dynamical System
   *
   *  \param num the level
   *  \return a shared ptr to the DynamicalSystemsGraph
   */
  inline std::shared_ptr<siconos::graphs::DynamicalSystemsGraph> dSG(unsigned int num) const
  {
    assert(num < _DSG.size());
    return _DSG[num];
  };

  /** get the number of Interactions Graphs
   *
   *  \return the number of Interactions Graphs
   */
  inline size_t indexSetsSize() const { return _IG.size(); };

  /** get the size of the  InteractionGraphs at a given level
   *
   *  \param level
   *  \return size of the  InteractionGraphs at a given level
   */
  inline size_t indexSetSize(unsigned int level) const { return _IG[level]->size(); };

  /** resize Interactions Graphs
   *
   *  \param newSize the new size
   */
  inline void indexSetsResize(unsigned int newSize) { _IG.resize(newSize); };

  // --- _hasChanged ---

  /** set _hasChanged to val
   *
   *  \param val a bool
   */
  inline void setHasChanged(const bool val) { _hasChanged = val; }

  /** check
   *
   *  \return a bool
   */
  inline bool hasChanged() const { return _hasChanged; }

  /** get the total number of scalar constraints
   *
   *  \return an unsigned int
   */
  inline unsigned int numberOfConstraints() const { return _numberOfConstraints; };

  void clear();

  /** set symmetry in the blocks computation
   *
   *  \param val a bool
   */
  void setSymmetric(bool val) { _symmetric = val; }

  /** initialize graphs properties */
  void setProperties();

  /** Get a dynamical system using its number
   *  \warning O(n) complexity
   *
   *  \param requiredNumber the required number
   *  \return a siconos::modeling::DynamicalSystem
   */
  std::shared_ptr<siconos::modeling::DynamicalSystem> getDynamicalSystem(
      unsigned int requiredNumber) const;

  /** list and display all dynamical systems
   */
  void displayDynamicalSystems() const;

  /** Get a dynamical system using its name
   *  \warning O(n) complexity
   *
   *  \param name the name of the dynamical system
   *  \return a siconos::modeling::DynamicalSystem
   */
  std::shared_ptr<siconos::modeling::DynamicalSystem> getDynamicalSystem(
      std::string name) const;

  /** Get a interaction using its number
   *  \warning O(n) complexity
   *
   *  \param requiredNumber the required number
   *  \return an Interaction
   */
  std::shared_ptr<siconos::modeling::Interaction> getInteraction(
      size_t requiredNumber) const;

  /** Get a interaction using its name
   *  \warning O(n) complexity
   *
   *  \param name the name of the Interaction
   *  \return an Interaction pointer
   */
  std::shared_ptr<siconos::modeling::Interaction> getInteraction(std::string name) const;

  /** get Interactions for a given DS
   *
   *  \return a vector of pointers to Interaction
   */
  std::vector<std::shared_ptr<siconos::modeling::Interaction>> interactionsForDS(
      std::shared_ptr<siconos::modeling::DynamicalSystem>) const;

  /** get Interactions for a given pair of DSs
   *
   *  \return a vector of pointers to Interaction
   */
  std::vector<std::shared_ptr<siconos::modeling::Interaction>> interactionsForPairOfDS(
      std::shared_ptr<siconos::modeling::DynamicalSystem> ds1,
      std::shared_ptr<siconos::modeling::DynamicalSystem> ds2 =
          std::shared_ptr<siconos::modeling::DynamicalSystem>()) const;

  /** get DynamicalSystems for a given Interaction
   *
   *  \return a vector of pointers to DynamicalSystem
   */
  std::vector<std::shared_ptr<siconos::modeling::DynamicalSystem>>
      dynamicalSystemsForInteraction(std::shared_ptr<siconos::modeling::Interaction>) const;

  /** Helper to get the descriptor in DSG0 from a DynamicalSystem
   *
   *  \param ds DynamicalSystem of which we want the descriptor
   *  \return the descriptor in DSG0 from a DynamicalSystem
   */
  siconos::graphs::DynamicalSystemsGraph::VDescriptor getDSG0Descriptor(
      std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
  {
    return _DSG[0]->descriptor(ds);
  }

  /** get the number of siconos::modeling::DynamicalSystem currently involved in an indexSet
   *
   *  \param inumber the indexSet number
   *  \return the number of siconos::modeling::DynamicalSystem involved
   */
  unsigned int numberOfInvolvedDS(unsigned int inumber);
};

}  // namespace siconos::simulation
#endif  // TOPOLOGY_H
