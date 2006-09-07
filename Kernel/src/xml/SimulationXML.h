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

/** \class SimulationXML
 *   \brief This class manages Simulation data part
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.3.0.
 *   \date 05/17/2004
 *
 *
 *
 * SimulationXML allows to get OneStepIntegratorXMLs and OneStepNSProblemXMLs from a DOM tree.
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

  /** \fn SimulationXML()
   *   \brief Default constructor
   */
  SimulationXML();

  /** \fn SimulationXML(xmlNodePtr rootSimulationNode)
   *   \brief Build a SimulationXML object from a DOM tree describing a Simulation
   *   \param the Simulation xml node
   */
  SimulationXML(xmlNodePtr rootSimulationNode);

  /** \fn ~SimulationXML()
   *  \brief destructor
   */
  ~SimulationXML();

  // --- GETTERS/SETTERS ---

  /** \fn xmlNodePtr getRootNode()
   *  \brief Return the root node of the Simulation -> tag "Simulation"
   *  \return xmlNodePtr  : the root node
   */
  inline xmlNodePtr getRootNode()
  {
    return rootNode;
  };

  /** \fn const SetOfOSIXML getOneStepIntegratorsXML()
   *   \brief get the set of OneStepIntegratorXML
   *   \return a SetOfOSIXML
   */
  inline const SetOfOSIXML getOneStepIntegratorsXML() const
  {
    return OSIXMLSet;
  };

  /** \fn const bool hasOneStepIntegratorXML() const
   *   \brief true if set of OneStepIntegratorXML is not empty
   *   \return a bool
   */
  inline const bool hasOneStepIntegratorXML() const
  {
    return !(OSIXMLSet.empty());
  }

  /** \fn * OneStepNSProblemXML getOneStepNSProblemXMLPtr()
   *   \brief Return the OneStepNSProblemXML pointer of the SimulationXML
   *   \return the OneStepNSProblemXML pointer of the SimulationXML ; NULL if SimulationXML does not have
   */
  inline OneStepNSProblemXML * getOneStepNSProblemXMLPtr()
  {
    return oneStepNSProblemXML;
  }

  /** \fn TimeDiscretisationXML* getTimeDiscretisationXMLPtr()
   *   \brief Return the TimeDiscretisationXML of the SimulationXML
   *   \return the TimeDiscretisationXML of the SimulationXML
   */
  inline TimeDiscretisationXML* getTimeDiscretisationXMLPtr()
  {
    return timeDiscretisationXML;
  }

  /** \fn inline string getSimulationXMLType()
   *   \brief Return the type of the Simulation
   *   \return the type of the SimulationXML
   */
  inline std::string  getSimulationXMLType()
  {
    return SiconosDOMTreeTools::getStringAttributeValue(rootNode, TYPE_ATTRIBUTE);
  }

  /** \fn bool hasOneStepNSProblemXML()
   *   \brief determines if the Simulation has a OneStepNSProblemXML
   *   \return bool :  false if the oneStepNSProblemXML* is NULL
   */
  inline const bool hasOneStepNSProblemXML() const
  {
    return (oneStepNSProblemXML != NULL);
  }

  /** \fn void saveSimulation2XML( xmlNodePtr  node, Simulation* str )
   *   \brief save data of str into the DOM tree
   *   \param xmlNodePtr  : the root node of the SimulationXML
   *   \param Simulation* : the Simulation of this SimulationXML
   */
  void saveSimulation2XML(xmlNodePtr  , Simulation*);
};



#endif
