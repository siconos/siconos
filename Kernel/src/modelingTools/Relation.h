/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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

/*! \file Relation.h
\brief General interface for relations.
 */

#ifndef RELATION_H
#define RELATION_H

#include "SimpleVector.h"
#include "SimpleMatrix.h"
#include "RuntimeException.h"
#include "Tools.h"
#include "SiconosSharedLibrary.h"

class Interaction;
class RelationXML;
class SimpleVector;

/** General Non Linear Relation (Virtual Base class for Relations).
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) Apr 27, 2004
 *
 *  A relation is a link between global variables of the Dynamical Systems and
 * some local ones, named y and lambda; belonging to one and only one Interaction.
 *
 * The present class is an interface to all relations provides tools to define and describe them.
 *
 * Each relation must have the two following functions:
 *
 *  - computeOutput(...) to compute y using DS global variables.
 *  - computeInput(...) to compute non-smooth DS part (r or p) using lambda.
 *
 * Depending on the DS class and the link type, various relations (ie derived classes) are available:
 *   - FirstOrder, for FirstOrderDS and derived classes.
 *   - Lagrangian, for LagrangianDS and derived classes.
 *
 *  The specific type (Linear, Scleronomous ...) is then given by the "subType". See Relation derived classes for details.
 *
 * The relation holds also:
 *  - a pointer to the Interaction that owns the present relation.
 *  - a pointer to an xml object
 *  - a VectorMap to handle links to DS variables (no copy!!!). Filled in during initialize.
 *
 */
class Relation
{

protected:

  /** type of the Relation: FirstOrder or Lagrangian */
  std::string  relationType;

  /** sub-type of the Relation (exple: LinearTIR or ScleronomousR ...) */
  std::string  subType;

  /** The Interaction linked to this Relation */
  Interaction * interaction;

  /** A map of vectors, used to save links (pointers) to DS objects of the interaction */
  VectorMap data;

  /** the object linked this Relation to read XML data */
  RelationXML *relationxml;

  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;

  /* contains the name of the plugin used to compute g function */
  /* contains the names of the various plug-in. Example: pluginNames["output"] is the function used to compute the output y.*/
  NamesList pluginNames;

  /** Flag to check if operators are plugged or not .*/
  BoolMap isPlugged;

  /** Flags to know if pointers have been allocated inside constructors or not */
  BoolMap isAllocatedIn;

  /** copy constructor => private, no copy nor pass-by-value.
   */
  Relation(const Relation&);

  /** default constructor
   *  \param a string that gives the type of the relation (optional)
   *  \param a string that gives the subtype of the relation (optional)
   */
  Relation(const std::string& = "Undefined", const std::string& = "Undefined");

  /** xml constructor
   *  \param RelationXML* : the XML object corresponding
   *  \param a string that gives the type of the relation (optional)
   *  \param a string that gives the subtype of the relation (optional)
   */
  Relation(RelationXML*, const std::string& = "Undefined", const std::string& = "Undefined");

public:

  /** destructor
   */
  virtual ~Relation();

  /** initialize the relation (check sizes, memory allocation ...)
   */
  virtual void initialize() = 0;

  /** To get the pointer to the Interaction linked to the present Relation
   *  \return a pointer to Interaction.
   */
  inline Interaction * getInteractionPtr()
  {
    return interaction;
  }

  /** To set the pointer to the Interaction linked to the present Relation
   *  \param a pointer to Interaction.
   */
  inline void setInteractionPtr(Interaction * newInter)
  {
    interaction = newInter;
  }

  /** To get the RelationXML* of the Relation
   *  \return a pointer on the RelationXML of the Relation
   */
  inline RelationXML* getRelationXML()
  {
    return relationxml;
  }

  /** To set the RelationXML* of the Relation
   *  \param RelationXML* : the pointer to set
   */
  inline void setRelationXML(RelationXML *rxml)
  {
    relationxml = rxml;
  }

  /** To get the type of the Relation (FirstOrder or Lagrangian)
   *  \return string : the type of the Relation
   */
  inline const std::string  getType() const
  {
    return relationType;
  }

  /** To get the subType of the Relation
   *  \return string : the sub-type of the Relation
   */
  inline const std::string  getSubType() const
  {
    return subType;
  }

  /** get the list of plug-in names
   *  \return a NamesList
   */
  inline NamesList getFunctionNames() const
  {
    return pluginNames;
  }

  /** get name of function that computes "name"
   *  \return a string
   */
  inline const std::string getFunctionName(const std::string& name) const
  {
    return (pluginNames.find(name))->second;
  }

  /** default function to compute y
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  virtual void computeOutput(double, unsigned int = 0) = 0;

  /** default function to compute y for the free state
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  virtual void computeFreeOutput(double, unsigned int = 0) = 0;

  /** default function to compute r
   *  \param double : current time
   *  \param unsigned int: "derivative" order of lambda used to compute input
   */
  virtual void computeInput(double, unsigned int) = 0;

  /** main relation members display
   */
  virtual void display() const;

  /** copy the data of the Relation to the XML tree
   */
  virtual void saveRelationToXML() const;

};

#endif // RELATION_H
