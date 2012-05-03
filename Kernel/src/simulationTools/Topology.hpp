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

/*! \file Topology.hpp

*/
#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "SiconosConst.hpp"
#include "InteractionsSet.hpp"
//#include "InteractionsSet.hpp"
#include "SimulationTypeDef.hpp"

class NonSmoothDynamicalSystem;
class Interaction;
class DynamicalSystem;
class SiconosMatrix;
class Interaction;


/**  This class describes the topology of the non-smooth dynamical
 *  system. It holds all the "potential" Interactions".
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) July 20, 2005
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
 * addInteractionInIndexSet(SP::Interaction inter) for each relation
 * of the interaction, it creates a new Interaction and inserts it
 * into indexSet0 It also counts the total number of "constraints" in
 * the system.
 *
 */
class Topology : public boost::enable_shared_from_this<Topology>
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Topology);


  // --- MEMBERS ---

  /** the set of all the interactions of the system */
  SP::InteractionsSet _allInteractions;

  /** dynamical systems graphs */
  std::vector<SP::DynamicalSystemsGraph> _DSG;

  /** Interaction graphs (_IG[0]=L[DSG[0]], L is the line graph
      transformation) */
  std::vector<SP::InteractionsGraph> _IG;

  /** check if topology has been updated since nsds modifications
      occur */
  bool _isTopologyUpToDate;

  /** check if topology is static or  not */
  bool _hasChanged;

  /** Total number of (scalar) constraints in the problem, ie sum of
      all nslaw sizes of Interactions of IndexSet0.*/
  unsigned int _numberOfConstraints;


  /** symmetry in the blocks computation */
  bool _symmetric;

  /** initializations ( time invariance) from non
      smooth laws kind */
  struct SetupFromNslaw;
  friend class Topology::SetupFromNslaw;


  // === PRIVATE FUNCTIONS ===

  /** schedules the relations of Interaction inter into IndexSet0 (ie
  * creates the corresponding Interactions and add them into _DSG
  * and _IG)
  \param  inter the Interaction to add
  \return a vertex descriptor of the new vertex in IndexSet0
  */
  InteractionsGraph::VDescriptor addInteractionInIndexSet(SP::Interaction inter);

  /** remove the Interactions of the interactions from _IG and
   * _DSG
   * \param inter the Interaction to remove
   */
  void removeInteractionFromIndexSet(SP::Interaction inter);

  /** default constructor
  */
  Topology();

public:

  // --- CONSTRUCTORS/DESTRUCTOR ---

  /** constructor from InteractionSet
  * \param a SP::InteractionSet
  */
  Topology(SP::InteractionsSet interSet);

  /** constructor from dynamical systems and interaction sets
  * \param newDSset a SP::DynamicalSystemsSet
  * \param newInteractions a SP::InteractionsSet
  */
  Topology(SP::DynamicalSystemsSet newDSset, SP::InteractionsSet newInteractions);


  /** destructor */
  ~Topology();


  // === GETTERS/SETTERS ===

  /** get all the Interactions of the Topology problem (saved in a set)
  *  \return an InteractionsSet
  */
  inline const SP::InteractionsSet interactions() const
  {
    return _allInteractions;
  }

  /** set the Interactions of the Topology problem (saved in a set)
   */
  inline void setInteractionsPtr(SP::InteractionsSet newInteractions)
  {
    _allInteractions->clear() ;
    _allInteractions = newInteractions;
  }


  /** check if Interaction inter is in the set
   *  \param inter an Interaction
   *  \return a bool
   */
  bool hasInteraction(SP::Interaction inter) const;

  /** add an Interaction in the topology. The interaction is both
   *  added in Dynamical Systems graph and Interactions Graph
   * \param a shared pointer to the interaction
   * \return a vertex descriptor to the new vertex in IndexSet0
   */
  InteractionsGraph::VDescriptor insertInteraction(SP::Interaction inter);

  /** remove an Interaction from the topology. The interaction is
   *  removed from Dynamical Systems graph and Interactions Graph.
   *  The interaction is not removed from actives subgraphs : see updateIndexSet
   *  \param inter the interaction to remove
   */
  void removeInteraction(SP::Interaction inter);

  /** add a dynamical system
   * \param ds the dynamical system to add
   */
  void insertDynamicalSystem(SP::DynamicalSystem ds);

  /** remove a dynamical system
   * \param ds the dynamical system to remove
   */
  void removeDynamicalSystem(SP::DynamicalSystem ds);

  /** link a dynamical system to a relation
   * \param inter a SP::Interaction
   * \param ds a SP::DynamicalSystem
   */
  void link(SP::Interaction inter, SP::DynamicalSystem ds);

  /** get a pointer to the graph of all Interactions.
   *  \return a SP::InteractionsGraph
   */
  inline SP::InteractionsGraph indexSet0()
  {
    return _IG[0];
  }

  /** get a pointer to the graph at level num of Interactions
   *  \return a SP::InteractionsGraph
   */
  inline SP::InteractionsGraph indexSet(unsigned int num)
  {
    assert(num < _IG.size()) ;
    return _IG[num];
  };

  /** get a pointer to the graph at level num of Interactions
   *  \return a SP::InteractionsGraph
   */
  inline unsigned int numberOfIndexSet()
  {
    return _IG.size();
  };

  /** reset graph at level num of Interactions
   *  \return a SP::InteractionsGraph
   */
  inline void resetIndexSetPtr(unsigned int num)
  {
    assert(num < _IG.size()) ;

    // .. global properties may be defined here with
    // InteractionsSubGraphProperties(), see SiconosProperties.hpp
    // VertexSubProperties or EdgeSubProperties and the macros
    // INSTALL_GRAPH_PROPERTIES

    _IG[num].reset(new InteractionsGraph());
    _IG[num]->properties().reset(new GraphProperties());

    _IG[num]->properties()->symmetric = _symmetric;

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

  /** get the number of Interactions Graphs */
  inline unsigned int indexSetsSize() const
  {
    return _IG.size();
  };

  /** resize Interactions Graphs
   * \param newSize the new size
   */
  inline void indexSetsResize(unsigned int newSize)
  {
    _IG.resize(newSize);
  };

  // --- isTopologyUpToDate ---

  /** set isTopologyUpToDate to val
  *  \param val a bool
  */
  inline void setUpToDate(const bool val)
  {
    _isTopologyUpToDate = val;
  }

  /** check if topology has been updated since modifications occurs on nsds
  *  \return a bool
  */
  inline bool isUpToDate() const
  {
    return _isTopologyUpToDate;
  }

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

  /** initializes the topology (called in Simulation->initialize)
  */
  void initialize();

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


};

DEFINE_SPTR(Topology);

#endif // TOPOLOGY_H
