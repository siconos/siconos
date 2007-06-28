/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
/*! \file
*/

#ifndef __SIMULATIONXML__
#define __SIMULATIONXML__

#include "SiconosDOMTreeTools.h"

#include <set>

class Simulation;
class OneStepIntegratorXML;
class OneStepNSProblemXML;
class TimeDiscretisationXML;

/** set of OneStepIntegratorXML */
typedef std::set<OneStepIntegratorXML*> SetOfOSIXML;

/** iterator through SetOfInteractionXML */
typedef SetOfOSIXML::iterator SetOfOSIXMLIt;

/** XML management for Simulation
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.1.0.
 *   \date 05/17/2004
 *
 *
 *
 * SimulationXML allows to get OneStepIntegratorXMLs and OneStepNSProblemXMLs from a DOM tree.
 */
class SimulationXML
{

private:

  /* root node for simulation */
  xmlNodePtr rootNode;

  /* set of OneStepIntegratorXML* */
  SetOfOSIXML OSIXMLSet;

  /* OneStepNSProblemXML - Optional value, default = NULL */
  OneStepNSProblemXML *oneStepNSProblemXML;

  /* TimeDiscretisationXML */
  TimeDiscretisationXML *timeDiscretisationXML;

public:

  /** Default constructor
  */
  SimulationXML();

  /** Build a SimulationXML object from a DOM tree describing a Simulation
  *   \param the Simulation xml node
  */
  SimulationXML(xmlNodePtr rootSimulationNode);

  /** destructor
  */
  ~SimulationXML();

  // --- GETTERS/SETTERS ---

  /** Return the root node of the Simulation -> tag "Simulation"
  *  \return xmlNodePtr  : the root node
  */
  inline xmlNodePtr getRootNode()
  {
    return rootNode;
  };

  /** get the set of OneStepIntegratorXML
  *   \return a SetOfOSIXML
  */
  inline const SetOfOSIXML getOneStepIntegratorsXML() const
  {
    return OSIXMLSet;
  };

  /** true if set of OneStepIntegratorXML is not empty
  *   \return a bool
  */
  inline const bool hasOneStepIntegratorXML() const
  {
    return !(OSIXMLSet.empty());
  }

  /** Return the OneStepNSProblemXML pointer of the SimulationXML
  *   \return the OneStepNSProblemXML pointer of the SimulationXML ; NULL if SimulationXML does not have
  */
  inline OneStepNSProblemXML * getOneStepNSProblemXMLPtr()
  {
    return oneStepNSProblemXML;
  }

  /** Return the TimeDiscretisationXML of the SimulationXML
  *   \return the TimeDiscretisationXML of the SimulationXML
  */
  inline TimeDiscretisationXML* getTimeDiscretisationXMLPtr()
  {
    return timeDiscretisationXML;
  }

  /** Return the type of the Simulation
  *   \return the type of the SimulationXML
  */
  inline std::string  getSimulationXMLType()
  {
    return SiconosDOMTreeTools::getStringAttributeValue(rootNode, TYPE_ATTRIBUTE);
  }

  /** determines if the Simulation has a OneStepNSProblemXML
  *   \return bool :  false if the oneStepNSProblemXML* is NULL
  */
  inline const bool hasOneStepNSProblemXML() const
  {
    return (oneStepNSProblemXML != NULL);
  }

  /** save data of str into the DOM tree
  *   \param xmlNodePtr  : the root node of the SimulationXML
  *   \param Simulation* : the Simulation of this SimulationXML
  */
  void saveSimulation2XML(xmlNodePtr  , Simulation*);
};



#endif
