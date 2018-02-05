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
public:
  enum
  {
    addDynamicalSystem, rmDynamicalSystem, addInteraction, rmInteraction, clearTopology
  };
  
  class Changes
  {
  public:
    int typeOfChange;
    SP::DynamicalSystem ds;
    SP::Interaction i;

    Changes(int t, SP::DynamicalSystem dsnew ):typeOfChange(t),ds(dsnew){};
    Changes(int t, SP::Interaction inew):typeOfChange(t),i(inew){};
    Changes(int t):typeOfChange(t){};
    void display() const
    {
      std::cout << "Changes display   " << this <<std::endl;
      if (typeOfChange == addDynamicalSystem)
      {
        std::cout << "typeOfChange : " << typeOfChange << " : addDynamicalSystem" << std::endl;
      }
      else if (typeOfChange == rmDynamicalSystem)
      {
        std::cout << "typeOfChange : " << typeOfChange << " : rmDynamicalSystem" << std::endl;
      }
      else if (typeOfChange == addInteraction)
      {
        std::cout << "typeOfChange : " << typeOfChange << " : addInteraction" << std::endl;
      }
      else if (typeOfChange == rmInteraction)
      {
        std::cout << "typeOfChange : " << typeOfChange << " : rmInteraction" << std::endl;
      }
      else if (typeOfChange == clearTopology)
      {
        std::cout << "typeOfChange : " << typeOfChange << " : clearTopology" << std::endl;
      }
    };
  };

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NonSmoothDynamicalSystem);


  /** current time of the simulation
      Warning FP : it corresponds to the time
      at the end of the integration step.
      It means that _t corresponds to tkp1 of the≈b
      simulation or nextTime().
   */
  double _t;

  /** initial time of the simulation */
  double _t0;

  /** final time of the simulation */
  double _T;

  /** information concerning the Model */
  std::string _title, _author, _description, _date;

  /** TRUE if the NonSmoothDynamicalSystem is a boundary value problem*/
  bool _BVP;

  /** log list of the modifications of the nsds */
  std::list<Changes> _changeLog;


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

  /** constructor with t0 and T
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


  inline const std::list<Changes>& changeLog()
  {
    return _changeLog;
  };

  inline std::list<Changes>::iterator changeLogPosition()
  {
    std::list<Changes>::iterator it = _changeLog.end();
    return --it;
  };

  




  
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
  void insertDynamicalSystem(SP::DynamicalSystem ds);

  /** get Dynamical system number I
   * \param nb the identifier of the DynamicalSystem to get
   * \return a pointer on DynamicalSystem
   */
  inline SP::DynamicalSystem dynamicalSystem(int nb) const
  {
    return _topology->getDynamicalSystem(nb);
  }

  /** remove a dynamical system
   * \param ds a pointer to the dynamical system to remove
   * \param removeInterations if true, all interactions connected to the ds will also be removed
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
   * \param ds a pointer to the system
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
