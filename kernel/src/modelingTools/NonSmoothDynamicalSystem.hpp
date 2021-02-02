/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

/** the NonSmoothDynamicalSystem consists in Dynamical Systems and Interactions
    structured into a graph defined in a Topology.
    In the DynamicalSystem graph, DynamicalSystem objects are nodes and Interaction objects
    are edges.

    To add a DynamicalSystem, use insertDynamicalSystem method.
    To add a new Interaction, use link method.

    A dual graph is also contructed, where Interactions are vertices and DynamicalSystems
    are edges.

*/
class NonSmoothDynamicalSystem
{
public:
  typedef enum
  {
    addDynamicalSystem, rmDynamicalSystem, addInteraction, rmInteraction, clearTopology
  } ChangeType;

  class Change
  {
  private:
    ACCEPT_SERIALIZATION(NonSmoothDynamicalSystem::Change);
    Change() = default;
  public:
    ChangeType typeOfChange;
    SP::DynamicalSystem ds;
    SP::Interaction i;

    Change(ChangeType t, SP::DynamicalSystem dsnew ):typeOfChange(t),ds(dsnew){};
    Change(ChangeType t, SP::Interaction inew):typeOfChange(t),i(inew){};
    Change(ChangeType t):typeOfChange(t){};
    void display() const;
  };

  typedef std::list<Change> ChangeLog;
  class ChangeLogIter
  {
    ACCEPT_SERIALIZATION(NonSmoothDynamicalSystem::Change);
  public:
    ChangeLogIter(){};
    ChangeLogIter(const ChangeLog& log,
                  const ChangeLog::const_iterator& i)
      : _log(&log), it(i) {};
    const ChangeLog *_log;
    ChangeLog::const_iterator it;
  };

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NonSmoothDynamicalSystem);

  /** initial time of the simulation */
  double _t0 = 0.;

  /** current time of the simulation
      Warning FP : it corresponds to the time
      at the end of the integration step.
      It means that _t corresponds to tkp1 of the
      simulation or nextTime().
   */
  double _t = _t0;

  /** final time of the simulation */
  double _T = 0.;

  /** information concerning the Model */
  std::string _title = "none",
    _author = "none",
    _description = "none",
    _date="unknown";

  /** TRUE if the NonSmoothDynamicalSystem is a boundary value problem*/
  bool _BVP = false;

  /** log list of the modifications of the nsds */
  std::list<Change> _changeLog;

  /** the topology of the system */
  SP::Topology _topology;

  /** False is one of the interaction is non-linear.
   */
  bool _mIsLinear = true;

  /* copy constructor, forbidden */
  NonSmoothDynamicalSystem(const NonSmoothDynamicalSystem& nsds) = delete;

  /* assignment, forbidden */
  OneStepNSProblem& operator=(const OneStepNSProblem& osnsp) = delete;

  /* Forbid default constructor except for serialization. No use to
     build a NSDS with t0=T.*/
protected:
  NonSmoothDynamicalSystem(){};

public:

  /** NSDS constructor.
   * \param t0 initial time
   * \param T final time
   */
  NonSmoothDynamicalSystem(double t0, double T);

  /** destructor
   */
  ~NonSmoothDynamicalSystem();

  // --- GETTERS/SETTERS ---
/** get the current time
   *  \return a double
   */
  inline double currentTime() const
  {
    return _t;
  }

  /** set the current time
   *  \param newValue the new time
   */
  inline void setCurrentTime(double newValue)
  {
    _t = newValue;
  }

  /** get initial time
   *  \return a double
   */
  inline double t0() const
  {
    return _t0;
  }

  /** set initial time of the time discretisation
   *  \param newT0
   */
  inline void sett0(double newT0)
  {
    _t0 = newT0;
  };

  /** get final time
   *  \return a double
   */
  inline double finalT() const
  {
    return _T;
  }

  /** set final time
   *  \param newValue the new final time for the Simulatiom
   */
  void setT(double newValue)
  {
    _T = newValue;
  };

  /** get the title of the simulation
   *  \return std::string : the title
   */
  inline const std::string  title() const
  {
    return _title;
  }

  /** set the title of the simulation
   *  \param s : the title
   */
  inline void setTitle(const std::string & s)
  {
    _title = s;
  }

  /** get the author of the simulation
   *  \return std::string : the author
   */
  inline const std::string  author() const
  {
    return _author;
  }

  /** set the author of the simulation
   *  \param s std::string : the author
   */
  inline void setAuthor(const std::string & s)
  {
    _author = s;
  }

  /** allows to get the description of the simulation
   *  \return std::string : the description
   */
  inline const std::string  description() const
  {
    return _description;
  }

  /** set the author of the simulation
   *  \param s std::string : the author
   */
  inline void setDescription(const std::string & s)
  {
    _description = s;
  }

  /** allows to get the date of the simulation
   *  \return std::string : the date
   */
  inline const std::string  date() const
  {
    return _date;
  }

  /** set the date of the simulation
   *  \param s std::string : the date
   */
  inline void setDate(const std::string & s)
  {
    _date = s;
  }

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


  /** get a reference to the changelog for an NSDS.
   * \return a reference to the changelog.
   */
  inline const ChangeLog& changeLog()
  {
    return _changeLog;
  };

  /** get an iterator to the last item in the changelog.
   * \return an iterator pointing at the last item in the changelog.
   */
  inline ChangeLogIter changeLogPosition()
  {
    ChangeLogIter it(_changeLog, _changeLog.end());
    // return iterator to last item, i.e. one less than end
    --it.it;
    return it;
  };

  /** get an iterator to the beginning of the changelog.
   * \return an iterator pointing at the beginning of the changelog.
   */
  inline ChangeLogIter changeLogBegin()
  {
    ChangeLogIter it(_changeLog, _changeLog.begin());
    return it;
  };

  /** clear the changelog up to a given position.
   *  \param it  This iterator must point to somewhere in the changelog
   *             for this NSDS.
   */
  void clearChangeLogTo(const ChangeLogIter& it);


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

  /** get all the dynamical systems declared in the NonSmoothDynamicalSystem.
   * into a std::vector<SP::DynamicalSystems>
   * Useful for iterates on DynamicalSystems in Python for instance
   * \return std::vector<SP::DynamicalSystems>
   */
  std::vector<SP::DynamicalSystem> dynamicalSystemsVector() const;


  /** add a dynamical system into the DS graph (as a vertex)
   * \param ds a pointer to the system to add
   */
  void insertDynamicalSystem(SP::DynamicalSystem ds);

  /** get Dynamical system number I
   * \param nb the identifier of the DynamicalSystem to get
   * \return a pointer on DynamicalSystem
   */
  inline SP::DynamicalSystem dynamicalSystem(int nb) const
  {
    return _topology->getDynamicalSystem(nb);
  }
  inline void displayDynamicalSystems() const
  {
    _topology->displayDynamicalSystems();
  }

  /** remove a dynamical system
   * \param ds a pointer to the dynamical system to remove
   */
  void removeDynamicalSystem(SP::DynamicalSystem ds);

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
  void removeInteraction(SP::Interaction inter);

  /** get Interaction number I
   * \param nb the identifier of the Interaction to get
   * \return a pointer to an Interaction
   */
  inline SP::Interaction interaction(int nb) const
  {
    return _topology->getInteraction(nb);
  }

  /** get Interaction named name
   * \param name of the Interaction to get
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

  /** get the name for this Dynamical System
   * \param ds a pointer to the system
   * \return name the name of the DynamicalSystem, or empty string if not found.
   */
  std::string name(SP::DynamicalSystem ds)
  {
    return _topology->name(ds);
  }

  /** set the name for this Interaction
   * \param interaction a pointer to the Interaction
   * \param name the name of the Interaction
   */
  inline void setName(SP::Interaction interaction, const std::string& name)
  {
    _topology->setName(interaction, name);
  };

  /** get the name for this Interaction
   * \param inter a pointer to the Interaction
   * \return name the name of the Interaction, or empty string if not found.
   */
  std::string name(SP::Interaction inter)
  {
    return _topology->name(inter);
  }

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

  /** compute output for all the interactions for a given level
   * \param time
   * \param level y order to be computed
   */
  void updateOutput(double time, unsigned int level = 0);

  /** compute output for all the interactions and for a level range
   * \param time
   * \param level_min y min order to be computed
   * \param level_max y max order to be computed
   */
  void updateOutput(double time, unsigned int level_min, unsigned int level_max);

  /** compute Jacobians for all the interactions (in indexSet0)
   * \param time
   */
  void computeInteractionJacobians(double time);

  /** compute Jacobians for all the interactions of a given index set.
   \param time
   \param indexSet InteractionsGraph of interest
   */
  void computeInteractionJacobians(double time, InteractionsGraph& indexSet);

  /** visit all dynamical systems in this system.
   * \param visitor an SP::SiconosVisitor that can visit classes derived from DS
   */
  void visitDynamicalSystems(SP::SiconosVisitor visitor);
};


#endif
