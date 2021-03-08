/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

#include "SiconosConst.hpp"
#include "SimulationTypeDef.hpp"
#include "SimulationGraphs.hpp"

/**  This class describes the topology of the non-smooth dynamical
 *  system. It holds all the "potential" Interactions".
 *
 *  Topology is built in NSDS constructors but initialized in
 *  Simulation->initialize(), ie when all Interactions have been
 *  clearly defined.
 *
 * Note that indexSet0 holds all the possible relations (declared by
 * user) not only those which are "actives".
 *
 * Construction consists in:
 *    - link with the NSDS that owns the topology.
 *
 * Initialization consists in:
 *    - scan of all the interactions of the NSDS
 *    - initialization of each interaction
 *    - insertion of the relations of all the Interaction into indexSet0
 *
 * Insertion of an Interaction into the set indexSet0:
 * addInteractionInIndexSet0(SP::Interaction inter) for each relation
 * of the interaction, it creates a new Interaction and inserts it
 * into indexSet0 It also counts the total number of "constraints" in
 * the system.
 *
 */
class Topology : public std::enable_shared_from_this<Topology>
{
public:

  typedef std::vector<SP::InteractionsGraph>::size_type IG_size_t;

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Topology);

  /** dynamical systems graphs */
  std::vector<SP::DynamicalSystemsGraph> _DSG;

  /** Interaction graphs (_IG[0]=L[DSG[0]], L is the line graph
      transformation) */
  std::vector<SP::InteractionsGraph> _IG;

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
  * creates the corresponding Interactions and add them into _DSG
  * and _IG)
  \param  inter the Interaction to add
  \param ds1 a SP::DynamicalSystem first dynamical system linked to the interaction
  \param ds2 a SP::DynamicalSystem second dynamical system linked to the interaction (default = None)
  \return a vertex descriptor of the new vertex in IndexSet0
  */
  std::pair<DynamicalSystemsGraph::EDescriptor, InteractionsGraph::VDescriptor>
  __addInteractionInIndexSet0(SP::Interaction inter, SP::DynamicalSystem ds1, SP::DynamicalSystem ds2 = SP::DynamicalSystem());

  /** remove an Interaction from _IG and _DSG
   * \param inter a pointer to the Interaction to be removed
   */
  void __removeInteractionFromIndexSet(SP::Interaction inter);

  /** remove a DynamicalSystem from _IG and _DSG
   * \param ds a pointer to the Dynamical System to be removed
   */
  void __removeDynamicalSystemFromIndexSet(SP::DynamicalSystem ds);

  /* forbid copy and assignment */
  Topology(const Topology&) = delete;
  Topology& operator=(const Topology&) = delete;
  
public:

  // --- CONSTRUCTORS/DESTRUCTOR ---

  /** default constructor */
  Topology();

  /** destructor */
  ~Topology();

  // === GETTERS/SETTERS ===

  /** check if an dynamical system is already a vertex of the DSs graph.
   *  \param ds the DS to test
   *  \return true if ds is in the graph
   */
  bool hasDynamicalSystem(SP::DynamicalSystem ds) const;

  /** check if an interaction is already a vertex of the Interactions graph.
   *  \param inter the Interaction to test
   *  \return true if inter is in the graph
   */
  bool hasInteraction(SP::Interaction inter) const;

  /** remove an Interaction from the topology. The interaction is
   *  removed from Dynamical Systems graph and Interactions Graph.
   *  The interaction is not removed from actives subgraphs : see updateIndexSet
   *  \param inter the interaction to remove
   */
  void removeInteraction(SP::Interaction inter);

  /** add a dynamical system
   * \param ds the DynamicalSystem to add
   */
  void insertDynamicalSystem(SP::DynamicalSystem ds);

  /** remove a Dynamical System from the topology. The dynamical
   *  system is removed from Dynamical Systems graph and Interactions
   *  Graph.  The dynamical system is not removed from actives
   *  subgraphs : see updateIndexSet
   *  \param ds the dynamical system to remove
   */
  void removeDynamicalSystem(SP::DynamicalSystem ds);

  /** set the name for this Dynamical System
   * \param ds a pointer to the system
   * \param name the name of the DynamicalSystem
   */
  void setName(SP::DynamicalSystem ds, const std::string& name);

  /** get the name for this Dynamical System
   * \param ds a pointer to the system
   * \return name the name of the DynamicalSystem, or empty string if not found.
   */
  std::string name(SP::DynamicalSystem ds);

  /** set the name for an Interaction
   * \param inter a pointer to the Interaction
   * \param name the name of the Interaction
   */
  void setName(SP::Interaction inter, const std::string& name);

  /** get the name for this Interaction
   * \param inter a pointer to the Interaction
   * \return name the name of the Interaction, or empty string if not found.
   */
  std::string name(SP::Interaction inter);

  /** set the OSI for this DynamicalSystem
   * \param ds the DynamicalSystem
   * \param OSI the integrator to use for this DS
   */
  void setOSI(SP::DynamicalSystem ds, SP::OneStepIntegrator OSI);

   /** link two dynamical systems to a relation
   * \param inter a SP::Interaction
   * \param ds a SP::DynamicalSystem
   * \param ds2 a SP::DynamicalSystem (optional)
   \return a vertex descriptor of the new vertex in IndexSet0
   */
  std::pair<DynamicalSystemsGraph::EDescriptor, InteractionsGraph::VDescriptor>
  link(SP::Interaction inter, SP::DynamicalSystem ds, SP::DynamicalSystem ds2 = SP::DynamicalSystem());

  /** specify if the given Interaction is for controlling the DS
   * \param inter Interaction
   * \param isControlInteraction true if the Interaction is used for
   * control purposes
   **/
  void setControlProperty(SP::Interaction inter,
                          const bool isControlInteraction);

  /** get a pointer to the graph of all Interactions.
   *  \return a SP::InteractionsGraph
   */
  inline SP::InteractionsGraph indexSet0() const
  {
    return _IG[0];
  }

  SP::InteractionProperties interaction_properties(unsigned int index, SP::Interaction inter)
  {
    InteractionsGraph::VDescriptor ui = indexSet(index)->descriptor(inter);
    SP::InteractionProperties inter_prop(new InteractionProperties(indexSet(index)->properties(ui)));
    return inter_prop;
  };

  /** get a pointer to the graph at level num of Interactions
   * \param num the number of indexSet
   * \return a SP::InteractionsGraph
   */
  SP::InteractionsGraph indexSet(unsigned int num) const;

  /** get a pointer to the graph at level num of Interactions
   *  \return a SP::InteractionsGraph
   */
  inline IG_size_t numberOfIndexSet() const
  {
    return _IG.size();
  };

  /** reset graph at level num of Interactions
   *  \param num the indexSet to reset
   */
  inline void resetIndexSetPtr(unsigned int num)
  {
    assert(num < _IG.size()) ;

    // .. global properties may be defined here with
    // InteractionsSubGraphProperties(), see SiconosProperties.hpp
    // VertexSubProperties or EdgeSubProperties and the macros
    // INSTALL_GRAPH_PROPERTIES

    _IG[num].reset(new InteractionsGraph());

    _IG[num]->properties().symmetric = _symmetric;
    _IG[num]->update_vertices_indices();
    _IG[num]->update_edges_indices();

  };

  /** get a pointer to the graph at level num of Dynamical System
   * \param num the level
   *\return a SP::DynamicalSystemsGraph
   */
  inline SP::DynamicalSystemsGraph dSG(unsigned int num) const
  {
    assert(num < _DSG.size()) ;
    return _DSG[num];
  };

  /** get the number of Interactions Graphs
   *  \return the number of Interactions Graphs
   */
  inline IG_size_t indexSetsSize() const
  {
    return _IG.size();
  };

  /** get the size of the  InteractionGraphs at a given level
   * \param level
   *  \return size of the  InteractionGraphs at a given level
   */
  inline IG_size_t indexSetSize(unsigned int level) const
  {
    return _IG[level]->size();
  };

  /** resize Interactions Graphs
   * \param newSize the new size
   */
  inline void indexSetsResize(unsigned int newSize)
  {
    _IG.resize(newSize);
  };


  // --- _hasChanged ---

  /** set _hasChanged to val
  *  \param val a bool
  */
  inline void setHasChanged(const bool val)
  {
    _hasChanged = val;
  }

  /** check
  *  \return a bool
  */
  inline bool hasChanged() const
  {
    return _hasChanged;
  }

  /** get the total number of scalar constraints
  *  \return an unsigned int
  */
  inline unsigned int numberOfConstraints() const
  {
    return _numberOfConstraints;
  };

  void clear();

  /** set symmetry in the blocks computation
   * \param val a bool
   */

  void setSymmetric(bool val)
  {
    _symmetric = val;
  }

  /** initialize graphs properties */
  void setProperties();

  /** Get a dynamical system using its number
   *   \warning O(n) complexity
   * \param requiredNumber the required number
   * \return a DynamicalSystem
   */
  SP::DynamicalSystem getDynamicalSystem(unsigned int requiredNumber) const;

  /** list and display all dynamical systems
   */
  void displayDynamicalSystems() const;

  /** Get a dynamical system using its name
   *  \warning O(n) complexity
   *  \param name the name of the dynamical system
   * \return a DynamicalSystem
   */
  SP::DynamicalSystem getDynamicalSystem(std::string name) const;

  /** Get a interaction using its number
   *  \warning O(n) complexity
   *  \param requiredNumber the required number
   *  \return an Interaction
   */
  SP::Interaction getInteraction(unsigned int requiredNumber) const;

  /** Get a interaction using its name
   *  \warning O(n) complexity
   *  \param name the name of the Interaction
   *  \return an Interaction pointer
   */
  SP::Interaction getInteraction(std::string name) const;

  /** get Interactions for a given DS
   * \return a vector of pointers to Interaction
   */
  std::vector<SP::Interaction> interactionsForDS(SP::DynamicalSystem) const;

  /** get Interactions for a given pair of DSs
   * \return a vector of pointers to Interaction
   */
  std::vector<SP::Interaction> interactionsForPairOfDS(
    SP::DynamicalSystem ds1,
    SP::DynamicalSystem ds2=SP::DynamicalSystem()) const;

  /** get DynamicalSystems for a given Interaction
   * \return a vector of pointers to DynamicalSystem
   */
  std::vector<SP::DynamicalSystem>
    dynamicalSystemsForInteraction(SP::Interaction) const;

  /** Helper to get the descriptor in DSG0 from a DynamicalSystem
   *  \param ds DynamicalSystem of which we want the descriptor
   *  \return the descriptor in DSG0 from a DynamicalSystem
   */
  DynamicalSystemsGraph::VDescriptor getDSG0Descriptor(SP::DynamicalSystem ds)
  {
    return _DSG[0]->descriptor(ds);
  }

  /** get the number of DynamicalSystem currently involved in an indexSet
   * \param inumber the indexSet number
   * \return the number of DynamicalSystem involved
   */
  unsigned int numberOfInvolvedDS(unsigned int inumber);


};

DEFINE_SPTR(Topology)

#endif // TOPOLOGY_H
