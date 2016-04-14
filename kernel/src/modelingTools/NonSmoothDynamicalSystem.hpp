/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
/*! \file NonSmoothDynamicalSystem.hpp
 * \brief container for DynamicalDystem and Interaction
 */
#ifndef NSDS_H
#define NSDS_H

#include "SiconosPointers.hpp"
#include "DynamicalSystemsSet.hpp"
#include "Topology.hpp"

/** the Non Smooth Dynamical System consists of DynamicalDystem
 *  and Interaction regrouped together in a Topology object,
 *  in the form of a graph of DynamicalDystem as nodes and Interaction as edges
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

  /** add a dynamical system
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
   * \return a pointer on Interaction
   */
  inline SP::Interaction interaction(int nb) const
  {
    return _topology->getInteraction(nb);
  }


  /** link an interaction to two dynamical systems
   * \param inter the interaction
   * \param ds1 a DynamicalSystem
   * \param ds2 a DynamicalSystem (optional)
   */
  void link(SP::Interaction inter, SP::DynamicalSystem ds1, SP::DynamicalSystem ds2 = SP::DynamicalSystem());

  

  

  inline void setOSI(SP::DynamicalSystem ds, SP::OneStepIntegrator OSI)
  {
    _topology->setOSI(ds, OSI);
    _mIsLinear = ((ds)->isLinear() && _mIsLinear);
  };

  /** set the name for this Dynamical System
   * \param ds a pointer to the system
   * \param name the name of the DynamicalSystem
   */
  inline void setName(SP::DynamicalSystem ds, const std::string& name)
  {
    _topology->setName(ds, name);
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

  /** calculate an indicator that gives convergence information for
   *  the DSs
   *  \return a double
   */
  double nsdsConvergenceIndicator();

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
};

#endif
