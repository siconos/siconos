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
/*! \file NonSmoothDynamicalSystem.hpp
 * container for DynamicalSystem and Interaction
 */
#ifndef NSDS_H
#define NSDS_H

#include <list>
#include <memory>
#include <string>
#include <vector>

#include "SiconosSerialization.hpp"

namespace siconos::internal {
struct SiconosVisitor;
}

namespace siconos::graphs {

struct DynamicalSystemsGraph;
struct InteractionsGraph;
}  // namespace siconos::graphs

namespace siconos::simulation {
class Topology;
}

namespace siconos::modeling {

class DynamicalSystem;
class Interaction;

/**
    the NonSmoothDynamicalSystem consists in Dynamical Systems and Interactions
    structured into a graph defined in a Topology.
    In the DynamicalSystem graph, DynamicalSystem objects are nodes and Interaction objects
    are edges.

    To add a DynamicalSystem, use insertDynamicalSystem method.
    To add a new Interaction, use link method.

    A dual graph is also contructed, where Interactions are vertices and DynamicalSystems
    are edges.

*/
class NonSmoothDynamicalSystem {
 public:
  enum class ChangeType {
    addDynamicalSystem,
    rmDynamicalSystem,
    addInteraction,
    rmInteraction,
    clearTopology
  };

  class Change {
   private:
    ACCEPT_SERIALIZATION(NonSmoothDynamicalSystem::Change);
    Change() = default;
    // Rule of five
    Change(const Change&) = delete;
    Change& operator=(Change&&) = delete;
    Change& operator=(const Change&) = delete;

   public:
    ChangeType typeOfChange;
    std::shared_ptr<DynamicalSystem> ds{nullptr};
    std::shared_ptr<Interaction> i{nullptr};

    Change(Change&&) = default;  // Required for push_back ...
    Change(ChangeType t, std::shared_ptr<DynamicalSystem> dsnew)
        : typeOfChange(t), ds(dsnew){};
    Change(ChangeType t, std::shared_ptr<Interaction> inew) : typeOfChange(t), i(inew){};
    Change(ChangeType t) : typeOfChange(t){};
    void display() const;
  };

 private:
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
  std::string _title{"none"};
  std::string _author = "none";
  std::string _description = "none";
  std::string _date = "unknown";

  /** TRUE if the NonSmoothDynamicalSystem is a boundary value problem*/
  bool _BVP = false;

  /** log list of the modifications of the nsds */
  std::list<Change> _changeLog = {};

  /** the topology of the system */
  std::shared_ptr<siconos::simulation::Topology> _topology{nullptr};

  /** False is one of the interaction is non-linear.
   */
  bool _mIsLinear = true;

  // Rule of five
  NonSmoothDynamicalSystem() = delete;
  NonSmoothDynamicalSystem(const NonSmoothDynamicalSystem&) = delete;
  NonSmoothDynamicalSystem(NonSmoothDynamicalSystem&&) = delete;
  NonSmoothDynamicalSystem& operator=(const NonSmoothDynamicalSystem&) = delete;
  NonSmoothDynamicalSystem& operator=(NonSmoothDynamicalSystem&&) = delete;

 public:
  /** NSDS constructor.
   *
   *  \param t0 initial time
   *  \param T final time
   */
  NonSmoothDynamicalSystem(double t0, double T);

  /** destructor
   */
  ~NonSmoothDynamicalSystem() noexcept;

  // --- GETTERS/SETTERS ---
  /** \return the current time value
   */
  inline double currentTime() const { return _t; }

  /** set the current time
   *
   *  \param newValue the new time
   */
  inline void setCurrentTime(double newValue) { _t = newValue; }

  /** \return initial time
   */
  inline double t0() const { return _t0; }

  /** set initial time of the time discretisation
   *
   *  \param newT0
   */
  inline void sett0(double newT0) { _t0 = newT0; };

  /** \return final time
   */
  inline double finalT() const { return _T; }

  /** set final time
   *
   *  \param newValue the new final time for the Simulatiom
   */
  void setT(double newValue) { _T = newValue; };

  /** get the title of the simulation
   *
   *  \return std::string : the title
   */
  inline const std::string title() const { return _title; }

  /** set the title of the simulation
   *
   *  \param s : the title
   */
  inline void setTitle(const std::string& s) { _title = s; }

  /** get the author of the simulation
   *
   *  \return std::string : the author
   */
  inline const std::string author() const { return _author; }

  /** set the author of the simulation
   *
   *  \param s std::string : the author
   */
  inline void setAuthor(const std::string& s) { _author = s; }

  /** allows to get the description of the simulation
   *
   *  \return std::string : the description
   */
  inline const std::string description() const { return _description; }

  /** set the author of the simulation
   *
   *  \param s std::string : the author
   */
  inline void setDescription(const std::string& s) { _description = s; }

  /** allows to get the date of the simulation
   *
   *  \return std::string : the date
   */
  inline const std::string date() const { return _date; }

  /** set the date of the simulation
   *
   *  \param s std::string : the date
   */
  inline void setDate(const std::string& s) { _date = s; }

  /** get problem type (true if BVP)
   *
   *  \return a bool
   */
  inline bool isBVP() const { return _BVP; }

  /** get problem type (true if IVP)
   *
   *  \return a bool
   */
  inline bool isIVP() const { return !_BVP; }

  /** set the NonSmoothDynamicalSystem to BVP, else it is IVP
   *
   *  \param newBvp true if BVP, false otherwise
   */
  inline void setBVP(const bool& newBvp) { _BVP = newBvp; }

  /** get a reference to the changelog for an NSDS.
   *
   *  \return a reference to the changelog.
   */
  inline const std::list<Change>& changeLog() { return _changeLog; };

  /** clear the changelog up to a given position.
   *
   *  \param it  This iterator must point to somewhere in the changelog
   *             for this NSDS.
   */
  void clearChangeLogTo(const std::list<Change>::const_iterator& it);

  // === DynamicalSystems management ===

  /** \return the number of Dynamical Systems present in the NSDS
   */
  size_t getNumberOfDS() const;

  /** get all the dynamical systems declared in the NonSmoothDynamicalSystem.
   *
   *  \return a std::shared_ptr<DynamicalSystemsGraph>
   */
  const std::shared_ptr<siconos::graphs::DynamicalSystemsGraph> dynamicalSystems() const;

  /** get all the dynamical systems declared in the NonSmoothDynamicalSystem.
   *  into a std::vector<std::shared_ptr<DynamicalSystem>>
   *  Useful for iterates on DynamicalSystems in Python for instance
   *
   *  \return std::vector<std::shared_ptr<DynamicalSystem>>
   */
  std::vector<std::shared_ptr<DynamicalSystem>> dynamicalSystemsVector() const;

  /** add a dynamical system into the DS graph (as a vertex)
   *
   *  \param ds a pointer to the system to add
   */
  void insertDynamicalSystem(std::shared_ptr<DynamicalSystem> ds);

  /** get Dynamical system number I
   *
   *  \param nb the identifier of the DynamicalSystem to get
   *  \return a pointer on DynamicalSystem
   */
  std::shared_ptr<siconos::modeling::DynamicalSystem> dynamicalSystem(unsigned int nb) const;

  void displayDynamicalSystems() const;

  /** remove a dynamical system
   *
   *  \param ds a pointer to the dynamical system to remove
   */
  void removeDynamicalSystem(std::shared_ptr<siconos::modeling::DynamicalSystem> ds);

  // === Interactions management ===

  /** get the number of Interactions present in the NSDS.
   *
   *  \return an unsigned int
   */
  size_t getNumberOfInteractions() const;

  /** return the graph of  Interactions present in the NSDS.
   *
   *  \return std::shared_ptr<InteractionsGraph>
   */
  const std::shared_ptr<siconos::graphs::InteractionsGraph> interactions() const;

  /** remove an interaction to the system
   *
   *  \param inter a pointer to the interaction to remove
   *
   */
  void removeInteraction(std::shared_ptr<siconos::modeling::Interaction> inter);

  /** get Interaction number I
   *
   *  \param nb the identifier of the Interaction to get
   *  \return a pointer to an Interaction
   */
  std::shared_ptr<siconos::modeling::Interaction> interaction(unsigned int nb) const;

  /** get Interaction named name
   *
   *  \param name of the Interaction to get
   *  \return a pointer to an Interaction
   */
  std::shared_ptr<siconos::modeling::Interaction> interaction(std::string name) const;

  /** get all the interactions declared in the NonSmoothDynamicalSystem.
   *  into a std::vector<std::shared_ptr<Interaction>>
   *  Useful for iterates on Interaction in Python for instance
   *
   *  \return std::vector<std::shared_ptr<Interaction>>
   */
  std::vector<std::shared_ptr<siconos::modeling::Interaction>> InteractionsVector() const;

  /** link an interaction to two dynamical systems
   *
   *  \param inter the interaction
   *  \param ds1 a DynamicalSystem
   *  \param ds2 a DynamicalSystem (optional)
   */
  void link(std::shared_ptr<siconos::modeling::Interaction> inter,
            std::shared_ptr<siconos::modeling::DynamicalSystem> ds1,
            std::shared_ptr<siconos::modeling::DynamicalSystem> ds2 =
                std::shared_ptr<siconos::modeling::DynamicalSystem>());

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

  /** set the name for this Interaction
   *
   *  \param interaction a pointer to the Interaction
   *  \param name the name of the Interaction
   */
  void setName(std::shared_ptr<siconos::modeling::Interaction> interaction,
               const std::string& name);

  /** get the name for this Interaction
   *
   *  \param inter a pointer to the Interaction
   *  \return name the name of the Interaction, or empty string if not found.
   */
  std::string name(std::shared_ptr<siconos::modeling::Interaction> inter);

  /** specify id the given Interaction is for controlling the DS
   *
   *  \param inter the Interaction
   *  \param isControlInteraction true if the Interaction is used for
   *  control purposes
   **/
  void setControlProperty(std::shared_ptr<siconos::modeling::Interaction> inter,
                          const bool isControlInteraction);

  /** get the topology of the system
   *
   *  \return a pointer on Topology
   */
  inline std::shared_ptr<siconos::simulation::Topology> topology() const { return _topology; }

  /** display the data of the Non Smooth Dynamical System
   */
  void display() const;

  /** return false is one of the interations is not linear.  else
   *  return true.
   *
   *  \return a bool
   */
  inline bool isLinear() const { return _mIsLinear; };

  void clear();

  /** set symmetry in the blocks computation
   *
   *  \param val a bool
   */
  void setSymmetric(bool val);

  /** Set all DS non-smooth part to zero for a given level.
   *
   *  \param level the level to will be zeroed
   */
  void resetNonSmoothPart(unsigned int level);

  /** save DynamicalSystems and Interactions states in Memories
   */
  void swapInMemory();

  /** save interaction states in memories. Applied to all interactions
      of the connected topology
  */
  void pushInteractionsInMemory();

  /** update the plugins of the DS
   *
   *  \param time to be used for plugins
   */
  void updateDSPlugins(double time);

  /** compute r thanks to lambda[level] for all Interactions
   *
   *  \param time
   *  \param level lambda level
   */
  void updateInput(double time, unsigned int level);

  /** compute output for all the interactions for a given level
   *
   *  \param time
   *  \param level y order to be computed
   */
  void updateOutput(double time, unsigned int level = 0);

  /** compute output for all the interactions and for a level range
   *
   *  \param time
   *  \param level_min y min order to be computed
   *  \param level_max y max order to be computed
   */
  void updateOutput(double time, unsigned int level_min, unsigned int level_max);

  /** compute Jacobians for all the interactions (in indexSet0)
   *
   *  \param time
   */
  void computeInteractionJacobians(double time);

  /** compute Jacobians for all the interactions of a given index set
   *
   *  \param time
   *  \param indexSet InteractionsGraph of interest
   */
  void computeInteractionJacobians(double time, siconos::graphs::InteractionsGraph& indexSet);

  /** visit all dynamical systems in this system
   *
   *  \param visitor a SiconosVisitor that can visit classes derived from DS
   */
  void visitDynamicalSystems(std::shared_ptr<siconos::internal::SiconosVisitor> visitor);
};

}  // namespace siconos::modeling
#endif
