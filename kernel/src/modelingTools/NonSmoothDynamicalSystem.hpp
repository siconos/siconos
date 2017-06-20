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
/*! \file NonSmoothDynamicalSystem.hpp
 * \brief container for DynamicalSystem and Interaction
 */
#ifndef NSDS_H
#define NSDS_H

#include "SiconosPointers.hpp"
#include "Topology.hpp"
#include "DynamicalSystem.hpp"

/** the Non Smooth Dynamical System consists of DynamicalSystem
 *  and Interaction regrouped together in a Topology object,
 *  in the form of a graph of DynamicalSystem as nodes and Interaction as edges
 *  and its dual.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \date (Creation) Apr 23, 2004
 *
 */
class NonSmoothDynamicalSystem
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NonSmoothDynamicalSystem);


  /** TRUE if the NonSmoothDynamicalSystem is a boundary value problem*/
  bool _BVP;

  /** the topology of the system */
  SP::Topology _topology;

  NonSmoothDynamicalSystem(const NonSmoothDynamicalSystem& nsds);

  /** False is one of the interaction is non-linear.
   */
  bool _mIsLinear;

public:

  /** default constructor
   */
  NonSmoothDynamicalSystem();

  /** destructor
   */
  ~NonSmoothDynamicalSystem();

  // --- GETTERS/SETTERS ---

  /** get problem type (true if BVP)
   *  \return a bool
   */
  inline bool isBVP() const
  {
    return _BVP;
  }

  /** get problem type (true if IVP)
   *  \return a bool
   */
  inline bool isIVP() const
  {
    return !_BVP;
  }

  /** set the NonSmoothDynamicalSystem to BVP, else it is IVP
   *  \param newBvp true if BVP, false otherwise
   */
  inline void setBVP(const bool& newBvp)
  {
    _BVP = newBvp;
  }

  // === DynamicalSystems management ===

  /** get the number of Dynamical Systems present in the NSDS
      \return an unsigned int
   */
  inline unsigned int getNumberOfDS() const
  {
    return _topology->dSG(0)->size();
  }

  /** get all the dynamical systems declared in the NonSmoothDynamicalSystem.
   * \return a SP::DynamicalSystemsGraph
   */
  inline const SP::DynamicalSystemsGraph dynamicalSystems() const
  {
    return _topology->dSG(0);
  }

  /** add a dynamical system into the DS graph (as a vertex)
   * \param ds a pointer to the system to add
   */
  inline void insertDynamicalSystem(SP::DynamicalSystem ds)
  {
    _topology->insertDynamicalSystem(ds);
    _mIsLinear = ((ds)->isLinear() && _mIsLinear);
  };

  /** get Dynamical system number I
   * \param nb the identifier of the DynamicalSystem to get
   * \return a pointer on DynamicalSystem
   */
  inline SP::DynamicalSystem dynamicalSystem(int nb) const
  {
    return _topology->getDynamicalSystem(nb);
  }



  // === Interactions management ===

  /** get the number of Interactions present in the NSDS.
   *  \return an unsigned int
   */
  inline unsigned int getNumberOfInteractions() const
  {
    return _topology->indexSet0()->size();
  };

  /** return the graph of  Interactions present in the NSDS.
   *  \return SP::InteractionGraph
   */
  inline const SP::InteractionsGraph  interactions() const
  {
    return _topology->indexSet0();
  };


  /** remove an interaction to the system
   * \param inter a pointer to the interaction to remove
   */
  inline void removeInteraction(SP::Interaction inter)
  {
    _topology->removeInteraction(inter);
  };

  /** get Interaction number I
   * \param nb the identifier of the Interaction to get
   * \return a pointer to an Interaction
   */
  inline SP::Interaction interaction(int nb) const
  {
    return _topology->getInteraction(nb);
  }

  /** get Interaction named name
   * \param nb the name of the Interaction to get
   * \return a pointer to an Interaction
   */
  inline SP::Interaction interaction(std::string name) const
  {
    return _topology->getInteraction(name);
  }

  /** link an interaction to two dynamical systems
   * \param inter the interaction
   * \param ds1 a DynamicalSystem
   * \param ds2 a DynamicalSystem (optional)
   */
  void link(SP::Interaction inter, SP::DynamicalSystem ds1, SP::DynamicalSystem ds2 = SP::DynamicalSystem());

  /** set the name for this Dynamical System
   * \param ds a pointer to the system
   * \param name the name of the DynamicalSystem
   */
  inline void setName(SP::DynamicalSystem ds, const std::string& name)
  {
    _topology->setName(ds, name);
  };

  /** set the name for this Interaction
   * \param interaction a pointer to the Interaction
   * \param name the name of the Interaction
   */
  inline void setName(SP::Interaction interaction, const std::string& name)
  {
    _topology->setName(interaction, name);
  };


    /** specify id the given Interaction is for controlling the DS
   * \param inter the Interaction
   * \param isControlInteraction true if the Interaction is used for
   * control purposes
   **/
  void setControlProperty(SP::Interaction inter, const bool isControlInteraction)
  {
    _topology->setControlProperty(inter, isControlInteraction);
  }


  /** get the topology of the system
   *  \return a pointer on Topology
   */
  inline SP::Topology topology() const
  {
    return _topology;
  }

  /** display the data of the Non Smooth Dynamical System
   */
  void display() const;

  /** return false is one of the interations is not linear.  else
   *  return true.
   *  \return a bool
   */
  inline bool isLinear() const
  {
    return _mIsLinear;
  };

  void clear();

  /** set symmetry in the blocks computation
   * \param val a bool
   */
  void setSymmetric(bool val);

  /** Set all DS non-smooth part to zero.
   */
  void reset();

  /** Set all DS non-smooth part to zero for a given level.
   * \param level the level to will be zeroed
   */
  void reset(unsigned int level);

  /** save DynamicalSystems and Interactions states in Memories
   */
  void swapInMemory();

  /** save interaction states in memories. Applied to all interactions
   of the connected topology
  */
  void pushInteractionsInMemory();

  /** compute r thanks to lambda[level] for all Interactions
    * \param time
    * \param level lambda level
   */
  void updateInput(double time, unsigned int level);

  /** compute output for all the interactions
   * \param time
   *  \param level y min order to be computed
   */
  void updateOutput(double time, unsigned int level = 0);


  /** visit all dynamical systems in this system.
   * \param visitor an SP::SiconosVisitor that can visit classes derived from DS
   */
  void visitDynamicalSystems(SP::SiconosVisitor visitor);
};


#endif

