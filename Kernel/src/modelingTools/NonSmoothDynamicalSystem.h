/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#ifndef NSDS_H
#define NSDS_H

#include "SiconosConst.h"

#include "DynamicalSystem.h"
#include "Interaction.h"
#include "EqualityConstraint.h"
#include "Topology.h"
#include "NonSmoothDynamicalSystemXML.h"
#include "DynamicalSystemsSet.h"
#include "InteractionsSet.h"
#include "check.h"
#include <iostream>
#include <map>
#include <vector>
#include <string>

enum dynamicalsystem {LAGRANGIANNLDS, LAGRANGIANTIDS, LINEARTIDS};

class Interaction;
class DynamicalSystem;
class EqualityConstraint;
class Topology;
class NonSmoothDynamicalSystemXML;

/** \class NonSmoothDynamicalSystem
 *  \brief This class describes the Non Smooth Dynamical System (NonSmoothDynamicalSystem)
 * composed with some dynamical systems, some algebraic constraints and some interactions between those systems
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
 *  \date (Creation) Apr 23, 2004
 *
 */

class NonSmoothDynamicalSystem
{
private:

  /** TRUE if the NonSmoothDynamicalSystem is a boundary value problem*/
  bool BVP;

  /** contains all the Dynamic Systems of the simulation */
  DynamicalSystemsSet allDS;

  /** inside-class allocation flags*/
  std::map<DynamicalSystem*, bool> isDSAllocatedIn;

  /** contains all the Interactions */
  InteractionsSet allInteractions;

  /** inside-class allocation flags*/
  std::map<Interaction*, bool> isInteractionAllocatedIn;

  /** contains the EqualityConstraints */
  std::vector<EqualityConstraint*> ecVector;

  /** the topology of the system */
  Topology * topology;

  /** the XML object linked to the NonSmoothDynamicalSystem to read XML data */
  NonSmoothDynamicalSystemXML *nsdsxml;

  /** Flags to check whether pointers were allocated in class constructors or not */
  std::deque<bool> isEcVectorAllocatedIn;
  bool isTopologyAllocatedIn;

  /** \fn NonSmoothDynamicalSystem(const bool&)
   *  \brief default constructor
   *  \param (optional) a bool which determines if the problem is BVP (true) or IVP (false)
   */
  NonSmoothDynamicalSystem(const bool& = false);

public:

  /** \fn NonSmoothDynamicalSystem(const NonSmoothDynamicalSystem&)
   *  \brief copy constructor
   *  \param  a ref to the NonSmoothDynamicalSystem to be copied
   */
  NonSmoothDynamicalSystem(const NonSmoothDynamicalSystem&);

  /** \fn NonSmoothDynamicalSystem(NonSmoothDynamicalSystemXML*)
   *  \brief xml constructor
   *  \param: the XML object corresponding to the NonSmoothDynamicalSystem
   */
  NonSmoothDynamicalSystem(NonSmoothDynamicalSystemXML*);

  /** \fn NonSmoothDynamicalSystem(DynamicalSystem* ds, Interaction* inter = NULL, const bool& isBvp = false)
   *  \brief constructor from minimum data.
   *  \param: a pointer to DynamicalSystem
   *  \param: a pointer to Interaction
   *  \param: a bool
   */
  NonSmoothDynamicalSystem(DynamicalSystem*, Interaction* = NULL, const bool& = false);

  /** \fn NonSmoothDynamicalSystem(DynamicalSystemsSet& ds, InteractionsSet& inter, const bool& isBvp = false)
   *  \brief constructor from data - Warning: DS and Interactions are not copied, but links are created
   *  between pointers of the two sets.
   *  \param: a set of DS
   *  \param: a set of Interactions
   *  \param: a bool
   */
  NonSmoothDynamicalSystem(DynamicalSystemsSet&, InteractionsSet&, const bool& = false);

  /** \fn NonSmoothDynamicalSystem(DynamicalSystemsSet& ds, const bool& isBvp = false)
   *  \brief constructor from data (only DS, no Interactions)
   *  between pointers of the two sets.
   *  \param: a set of DS
   *  \param: a bool
   */
  NonSmoothDynamicalSystem(DynamicalSystemsSet&, const bool& = false);

  /** \fn ~NonSmoothDynamicalSystem()
   *  \brief destructor
   */
  ~NonSmoothDynamicalSystem();

  // --- GETTERS/SETTERS ---

  /** \fn const bool isBVP(void)
   *  \brief get problem type (true if BVP)
   *  \return a bool
   */
  inline const bool isBVP() const
  {
    return BVP;
  }

  /** \fn const bool isIVP(void)
   *  \brief get problem type (true if IVP)
   *  \return a bool
   */
  inline const bool isIVP() const
  {
    return !BVP;
  }

  /** \fn void setBVP(const bool&)
   *  \brief set the NonSmoothDynamicalSystem to BVP, else it is IVP
   *  \param bool : true if BVP, false otherwise
   */
  inline void setBVP(const bool& newBvp)
  {
    BVP = newBvp;
  }

  // === DynamicalSystems management ===

  /** \fn inline const unsigned int getNumberOfDS() const
   *  \brief get the number of Dynamical Systems present in the NSDS (ie in allDS set)
   *  \return an unsigned int
   */
  inline const unsigned int getNumberOfDS() const
  {
    return allDS.size();
  };

  /** \fn const DynamicalSystemsSet getDynamicalSystems()
   *  \brief get all the DynamicalSystem of the NonSmoothDynamicalSystem problem (saved in a set)
   *  \return a DynamicalSystemsSet
   */
  inline const DynamicalSystemsSet getDynamicalSystems() const
  {
    return allDS;
  }

  /** \fn DynamicalSystem* getDynamicalSystemPtr(const int& position)
   *  \brief get DynamicalSystem at indix position in the set
   *  \param an int
   *  \return a pointer on DynamicalSystem
   */
  DynamicalSystem* getDynamicalSystemPtr(const int&) const ;

  /** \fn DynamicalSystem* getDynamicalSystemPtrNumber(const int& I)
   *  \brief get Dynamical system number I
   *  \param the identifier of the DynamicalSystem to get
   *  \return a pointer on DynamicalSystem
   */
  DynamicalSystem* getDynamicalSystemPtrNumber(const int&) const ;

  /** \fn void setDynamicalSystems(const DynamicalSystemsSet&)
   *  \brief to set allDS
   *  \param a DynamicalSystemsSet
   */
  void setDynamicalSystems(const DynamicalSystemsSet&) ;

  /** \fn const bool hasDynamicalSystemNumber(const int& N)  const
   *  \brief check if DynamicalSystem number N exists
   *  \param the identifier of the DynamicalSystem to get
   *  \return bool
   */
  const bool hasDynamicalSystemNumber(const int&) const ;

  /** \fn const bool hasDynamicalSystem(DynamicalSystem* ds)
   *  \brief check if DynamicalSystem ds is in the set
   *  \param a pointer to DynamicalSystem
   *  \return a bool
   */
  const bool hasDynamicalSystem(DynamicalSystem*) const;

  // === Interactions management ===

  /** \fn inline const unsigned int getNumberOfInteractions() const
   *  \brief get the number of Interactions present in the NSDS (ie in allInteractions set)
   *  \return an unsigned int
   */
  inline const unsigned int getNumberOfInteractions() const
  {
    return allInteractions.size();
  };

  /** \fn const InteractionsSet getInteractions()
   *  \brief get all the Interactions of the NonSmoothDynamicalSystem problem (saved in a set)
   *  \return an InteractionsSet
   */
  inline const InteractionsSet getInteractions() const
  {
    return allInteractions;
  }

  /** \fn Interaction* getInteractionPtr(const int& position)
   *  \brief get Interaction at indix position in the set
   *  \param an int
   *  \return a pointer on Interaction
   */
  Interaction* getInteractionPtr(const int&) const ;

  /** \fn Interaction* getInteractionPtrNumber(const int& I)
   *  \brief get Interaction number I
   *  \param the id-number of the Interaction to get
   *  \return a pointer on Interaction
   */
  Interaction* getInteractionPtrNumber(const int&) const ;

  /** \fn void setInteractions(const InteractionsSet&)
   *  \brief to set allInteractions
   *  \param an InteractionsSet
   */
  void setInteractions(const InteractionsSet&) ;

  /** \fn const bool hasInteractionNumber(const int& N)
   *  \brief check if Interaction number N exists
   *  \param the identifier of the Interaction to get
   *  \return bool
   */
  const bool hasInteractionNumber(const int&) const;

  /** \fn const bool hasInteraction(Interaction * inter)
   *  \brief check if Interaction inter is in the set
   *  \param a pointer to Interaction
   *  \return a bool
   */
  const bool hasInteraction(Interaction*) const;

  // === Equality constraints management ===

  /** \fn vector<EqualityConstraint*> getEqualityConstraints(void)
   *  \brief get the vector of algebraic constraints
   *  \return vector of EqualityConstraint
   */
  inline const std::vector<EqualityConstraint*> getEqualityConstraints(void) const
  {
    return ecVector;
  }

  /** \fn EqualityConstraint* getEqualityConstraintPtr(const int& N) const
   *  \brief get algebraic constraint at position N in vectorEc
   *  \param int : the position of the ec to get
   *  \return a pointer on EC
   */
  EqualityConstraint* getEqualityConstraintPtr(const int&) const;

  /** \fn void setEqualityConstraints(const vector<EqualityConstraint*>& )
   *  \brief set the vector of algebraic constraints
   *  \param vector<EqualityConstraint*> : new value for the vector
   */
  void setEqualityConstraints(const std::vector<EqualityConstraint*>& newEcVect) ;

  /** \fn Topology* getTopologyPtr() const
   *  \brief get the topology of the system
   *  \return a pointer on Topology
   */
  inline Topology* getTopologyPtr() const
  {
    return topology;
  }

  /** \fn inline NonSmoothDynamicalSystemXML* getNonSmoothDynamicalSystemXMLPtr()
   *  \brief get the xml linked object
   *  \return a pointer on NonSmoothDynamicalSystemXML
   */
  inline NonSmoothDynamicalSystemXML* getNonSmoothDynamicalSystemXMLPtr() const
  {
    return nsdsxml;
  }

  /** \fn inline void setNonSmoothDynamicalSystemXMLPtr( NonSmoothDynamicalSystemXML *nsdsxml )
   *  \brief set the xml linked object
   *  \param NonSmoothDynamicalSystemXML* : a pointer on NonSmoothDynamicalSystemXML* to link
   */
  inline void setNonSmoothDynamicalSystemXMLPtr(NonSmoothDynamicalSystemXML *newNsdsxml)
  {
    nsdsxml = newNsdsxml;
  }

  // --- OTHER FUNCTIONS ---

  /** \fn void saveNSDSToXML()
   *  \brief copy the data of the NonSmoothDynamicalSystem to the XML tree
   *  \exception RuntimeException
   */
  void saveNSDSToXML();

  /** \fn void display()
   *  \brief display the data of the Non Smooth Dynamical System
   */
  void display() const;

  /** \fn void addDynamicalSystemPtr(DynamicalSystem*)
   *  \brief add a DynamicalSystem into the NonSmoothDynamicalSystem (pointer link, no copy!)
   *  \param DynamicalSystem* : the DynamicalSystem to add
   */
  void addDynamicalSystemPtr(DynamicalSystem*);

  /** \fn void addInteractionPtr(Interaction*)
   *  \brief add an Interaction into the NonSmoothDynamicalSystem (pointer link, no copy!)
   *  \param Interaction : the Interaction to add
   */
  void addInteractionPtr(Interaction*);

  /** \fn void addEqualityConstraint(EqualityConstraint*)
   *  \brief add an EqualityConstraint to the NonSmoothDynamicalSystem
   *  \param EqualityConstraint* : the EqualityConstraint to add
   */
  void addEqualityConstraint(EqualityConstraint*);

  /** \fn double nsdsConvergenceIndicator() const
   *  \brief calculate an indicator that gives convergence information for the DSs
   *  \return a double
   */
  double nsdsConvergenceIndicator() ;
};

#endif
