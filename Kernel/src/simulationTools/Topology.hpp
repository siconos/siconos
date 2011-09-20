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

/*! \file Topology.h

*/
#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "SiconosConst.hpp"
#include "InteractionsSet.hpp"
#include "UnitaryRelationsSet.hpp"
#include "SimulationTypeDef.hpp"

class NonSmoothDynamicalSystem;
class Interaction;
class DynamicalSystem;
class SiconosMatrix;
class UnitaryRelation;


/**  This class describes the topology of the non-smooth dynamical
 *  system. It holds all the "potential" Unitary Relations".
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
 * of the interaction, it creates a new UnitaryRelation and inserts it
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

  /** unitary relations graphs (URG[0]=L[DSG[0]], L is the line graph
      transformation) */
  std::vector<SP::UnitaryRelationsGraph> _URG;

  /** check if topology has been updated since nsds modifications
      occur */
  bool _isTopologyUpToDate;

  /** check if topology is static or  not */
  bool _hasChanged;

  /** Total number of (scalar) constraints in the problem, ie sum of
      all nslaw sizes of Unitary Relations of IndexSet0.*/
  unsigned int _numberOfConstraints;


  /** symmetry in the blocks computation */
  bool _symmetric;

  /** initializations ( time invariance) from non
      smooth laws kind */
  struct SetupFromNslaw;
  friend class Topology::SetupFromNslaw;


  // === PRIVATE FUNCTIONS ===

  /** schedules the relations of Interaction inter into IndexSet0 (ie
  * creates the corresponding UnitaryRelations and add them into DSG
  * anr URG)
  \param: a pointer to Interaction
  */
  void addInteractionInIndexSet(SP::Interaction);

  /** remove the unitary relations of the interactions from URG and
   *   DSG */
  void removeInteractionFromIndexSet(SP::Interaction);

  /** default constructor
  */
  Topology();

public:

  // --- CONSTRUCTORS/DESTRUCTOR ---

  /** constructor from InteractionSet
  * \param: a SP::InteractionSet
  */
  Topology(SP::InteractionsSet);

  /** constructor from dynamical systems and interaction sets
  * \param: a SP::DynamicalSystemsSet
  * \param: a SP::InteractionsSet
  */
  Topology(SP::DynamicalSystemsSet, SP::InteractionsSet);


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
   *  \param a shared pointer to Interaction
   *  \return a bool
   */
  bool hasInteraction(SP::Interaction) const;

  /** add an Interaction in the topology. The interaction is both
   *  added in Dynamical Systems graph and Unitary Relations Graph
   * \param a shared pointer to the interaction
   */
  void insertInteraction(SP::Interaction);

  /** remove an Interaction from the topology. The interaction is
   *  removed from Dynamical Systems graph and Unitary Relations Graph.
   *  The interaction is not removed from actives subgraphs : see updateIndexSet
   *  \param a shared pointer to the interaction
   */
  void removeInteraction(SP::Interaction);

  /** add a dynamical system
   * \param a shared pointer to a dynamical system
   */
  void insertDynamicalSystem(SP::DynamicalSystem ds);

  /** remove a dynamical system
   * \param a shared pointer to a dynamical system
   */
  void removeDynamicalSystem(SP::DynamicalSystem ds);

  /** link a dynamical system to a relation
   * \param a SP::Interaction
   * \param a SP::DynamicalSystem
   */
  void link(SP::Interaction, SP::DynamicalSystem);

  /** get a pointer to the graph of all Unitary Relations.
   *  \return a SP::UnitaryRelationsGraph
   */
  inline SP::UnitaryRelationsGraph indexSet0()
  {
    return _URG[0];
  }

  /** get a pointer to the graph at level num of Unitary Relations
   *  \return a SP::UnitaryRelationsGraph
   */
  inline SP::UnitaryRelationsGraph indexSet(unsigned int num)
  {
    assert(num < _URG.size()) ;
    return _URG[num];
  };

  /** reset graph at level num of Unitary Relations
   *  \return a SP::UnitaryRelationsGraph
   */
  inline void resetIndexSetPtr(unsigned int num)
  {
    assert(num < _URG.size()) ;
    _URG[num].reset(new UnitaryRelationsGraph());
    _URG[num]->properties().symmetric = false;
  };

  /** get a pointer to the graph at level num of Dynamical System
   *  \return a SP::DynamicalSystemsGraph
   */
  inline SP::DynamicalSystemsGraph dSG(unsigned int num)
  {
    assert(num < _DSG.size()) ;
    return _DSG[num];
  };

  /** get the number of Unitary Relations Graphs */
  inline unsigned int indexSetsSize()
  {
    return _URG.size();
  };

  /** resize Unitary Relations Graphs */
  inline void indexSetsResize(unsigned int i)
  {
    return _URG.resize(i);
  };

  // --- isTopologyUpToDate ---

  /** set isTopologyUpToDate to val
  *  \param a bool
  */
  inline void setUpToDate(const bool val)
  {
    _isTopologyUpToDate = val;
  }

  /** check if topology has been updated since modifications occurs on nsds
  *  \return a bool
  */
  inline bool isUpToDate()
  {
    return _isTopologyUpToDate;
  }

  // --- _hasChanged ---

  /** set _hasChanged to val
  *  \param a bool
  */
  inline void setHasChanged(const bool val)
  {
    _hasChanged = val;
  }

  /** check
  *  \return a bool
  */
  inline bool hasChanged()
  {
    return _hasChanged;
  }

  /** get the total number of scalar constraints
  *  \return an unsigned int
  */
  inline unsigned int numberOfConstraints()
  {
    return _numberOfConstraints;
  };

  /** initializes the topology (called in Simulation->initialize)
  */
  void initialize();

  void clear();

  /** set symmetry in the blocks computation
   * \param a bool
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
