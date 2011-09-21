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
/*! \file NonSmoothDynamicalSystem.hpp
 */
#ifndef NSDS_H
#define NSDS_H

/** Available Dynamical Systems types*/
enum dynamicalsystem {LAGRANGIANNLDS, LAGRANGIANTIDS, LINEARTIDS};

#include "SiconosPointers.hpp"
#include "InteractionsSet.hpp"
#include "DynamicalSystemsSet.hpp"
#include "Topology.hpp"

class Interaction;
class DynamicalSystem;
class Topology;


/** the Non Smooth Dynamical System composed with dynamical systems
 *  that interact alltogether.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
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
  bool BVP;

  /** the topology of the system */
  SP::Topology _topology;

  /** the XML object linked to the NonSmoothDynamicalSystem to read XML data */
  SP::NonSmoothDynamicalSystemXML nsdsxml;

  /** default constructor
   */
  NonSmoothDynamicalSystem(): BVP(false) {};

private:

  NonSmoothDynamicalSystem(const NonSmoothDynamicalSystem&);

  /** False is one of the interaction is non-linear.
   */
  bool mIsLinear;

public:

  /** xml constructor
   *  \param: the XML object corresponding to the NonSmoothDynamicalSystem
   */
  NonSmoothDynamicalSystem(SP::NonSmoothDynamicalSystemXML);

  /** constructor from minimum data.
   *  \param: a pointer to DynamicalSystem
   *  \param: a pointer to Interaction
   *  \param: a bool
   */
  NonSmoothDynamicalSystem(SP::DynamicalSystem, SP::Interaction = SP::Interaction(), const bool& = false);

  /** constructor from data - Warning: DS and Interactions are not copied, but links are created
   *  between pointers of the two sets.
   *  \param: a set of DS
   *  \param: a set of Interactions
   *  \param: a bool
   */
  NonSmoothDynamicalSystem(DynamicalSystemsSet&, InteractionsSet&, const bool& = false);

  /** constructor from data (only DS, no Interactions)
   *  between pointers of the two sets.
   *  \param: a set of DS
   *  \param: a bool
   */
  NonSmoothDynamicalSystem(DynamicalSystemsSet&, const bool& = false);

  /** destructor
   */
  ~NonSmoothDynamicalSystem();

  // --- GETTERS/SETTERS ---

  /** get problem type (true if BVP)
   *  \return a bool
   */
  inline bool isBVP() const
  {
    return BVP;
  }

  /** get problem type (true if IVP)
   *  \return a bool
   */
  inline bool isIVP() const
  {
    return !BVP;
  }

  /** set the NonSmoothDynamicalSystem to BVP, else it is IVP
   *  \param bool : true if BVP, false otherwise
   */
  inline void setBVP(const bool& newBvp)
  {
    BVP = newBvp;
  }

  // === DynamicalSystems management ===

  /** get the number of Dynamical Systems present in the NSDS
      \return an unsigned int
   */
  inline unsigned int getNumberOfDS() const
  {
    return topology()->dSG(0)->size();
  };

  /** get all the DynamicalSystem of the NonSmoothDynamicalSystem
   *  problem
   * \return a DynamicalSystemsSet *
   */
  inline const SP::DynamicalSystemsGraph dynamicalSystems() const
  {
    return topology()->dSG(0);
  }

  // === Interactions management ===

  /** get the number of Interactions present in the NSDS (ie in allInteractions set)
   *  \return an unsigned int
   */
  inline unsigned int getNumberOfInteractions() const
  {
    return _topology->interactions()->size();
  };

  /** get all the Interactions of the NonSmoothDynamicalSystem problem (saved in a set)
   *  \return an InteractionsSet *
   */
  inline const SP::InteractionsSet interactions() const
  {
    return _topology->interactions();
  }

  /** add an interaction to the system
   * \param a shared pointer to the interaction
   */
  void insertInteraction(SP::Interaction inter)
  {
    _topology->insertInteraction(inter);
  };


  /** remove an interaction to the system
   * \param a shared pointer to the interaction
   */
  void removeInteraction(SP::Interaction inter)
  {
    _topology->removeInteraction(inter);
  };

  /** add a dynamical system
   * \param a shared pointer to a dynamical system
   */
  void insertDynamicalSystem(SP::DynamicalSystem ds)
  {
    _topology->insertDynamicalSystem(ds);
    mIsLinear = ((ds)->isLinear() && mIsLinear);
  };


  /** remove a dynamical system
   * \param a shared pointer to a dynamical system
   */
  void removeDynamicalSystem(SP::DynamicalSystem ds)
  {
    _topology->removeDynamicalSystem(ds);
  };

  /** link an interaction to a dynamical system
   * \param a SP::Interaction
   * \param a SP::DynamicalSystem
   */
  void link(SP::Interaction inter, SP::DynamicalSystem ds)
  {
    _topology->link(inter, ds);
    mIsLinear = ((inter)->relation()->isLinear() && mIsLinear);
  };


  /** get Dynamical system number I
      -   *  \param the identifier of the DynamicalSystem to get
      -   *  \return a pointer on DynamicalSystem
  */
  SP::DynamicalSystem dynamicalSystemNumber(int) const ;


  /** get the topology of the system
   *  \return a pointer on Topology
   */
  inline SP::Topology topology() const
  {
    return _topology;
  }

  /** get the xml linked object
   *  \return a pointer on NonSmoothDynamicalSystemXML
   */
  inline SP::NonSmoothDynamicalSystemXML nonSmoothDynamicalSystemXML()
  const
  {
    return nsdsxml;
  }

  /** set the xml linked object
   *  \param NonSmoothDynamicalSystemXML* : a pointer on
   *  NonSmoothDynamicalSystemXML* to link
   */
  inline void setNonSmoothDynamicalSystemXMLPtr(SP::NonSmoothDynamicalSystemXML newNsdsxml)
  {
    nsdsxml = newNsdsxml;
  }

  // --- OTHER FUNCTIONS ---

  /** copy the data of the NonSmoothDynamicalSystem to the XML tree
   *  \exception RuntimeException
   */
  void saveNSDSToXML();

  /** display the data of the Non Smooth Dynamical System
   */
  void display() const;

  /** calculate an indicator that gives convergence information for
   *  the DSs
   *  \return a double
   */
  double nsdsConvergenceIndicator() ;

  /** return false is one of the interations is not linear.  else
   *  return true.
   *  \return a bool
   */
  bool isLinear();

  void clear();

  /** set symmetry in the blocks computation
   * \param a bool
   */
  void setSymmetric(bool val)
  {
    topology()->setSymmetric(val);
  }

};

#endif
