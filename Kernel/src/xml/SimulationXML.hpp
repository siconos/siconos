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
/*! \file SimulationXML.hpp
  \brief XML management for simulation-type objects
 */

#ifndef __SIMULATIONXML__
#define __SIMULATIONXML__

#include "SiconosPointers.hpp"
#include "SiconosDOMTreeTools.hpp"
#include <set>

class Simulation;

/** set of OneStepIntegratorXML */
typedef std::set<SP::OneStepIntegratorXML> SetOfOSIXML;

/** iterator through set of OneStepIntegratorXML */
typedef SetOfOSIXML::iterator SetOfOSIXMLIt;

/** set of OneStepNSProblemXML */
typedef std::set<SP::OneStepNSProblemXML> SetOfOSNSPBXML;

/** iterator through  set of OneStepNSProblemXML*/
typedef SetOfOSNSPBXML::iterator SetOfOSNSPBXMLIt;

TYPEDEF_SPTR(SetOfOSIXML)
TYPEDEF_SPTR(SetOfOSNSPBXML)

/** XML management for Simulation
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 05/17/2004
 *
 *  XML interface for XML reading in classes TimeStepping and EventDriven.
 *
 * Attributes for Simulation tag:
 *  - type (required), simulation type, either TimeStepping or EventDriven.
 *
 * Tags for Simulation:
 *  - TimeDiscretisation (required): the time discretisation scheme
 *  - OneStepIntegrator_Definition (required): a list of one-step integrators
 *  - OneStepNSProblem (optional): a list of one-step non-smooth problems
 *
 * Template for XML input:
 \code
<Simulation type='EventDriven'>  <!-- Type = EventDriven or TimeStepping  -->
<TimeDiscretisation>
<h>0.05</h>
</TimeDiscretisation>
<OneStepIntegrator_Definition>
<!-- A list of OneStepIntegrators - See documentation of this class for details on XML input  -->
</OneStepIntegrator_Definition>
<OneStepNSProblem>
<!-- A list of OneStepNSProblem - See documentation of this class for details on XML input  -->
</OneStepNSProblem>
</Simulation>
\endcode
 *
 */
class SimulationXML
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SimulationXML);


  /* root node for simulation - Tag: Simulation */
  xmlNodePtr rootNode;

  /* set of OneStepIntegratorXML* */
  SetOfOSIXML OSIXMLSet;

  /* set of OneStepIntegratorXML* */
  SetOfOSNSPBXML OSNSPBXMLSet;

  /* TimeDiscretisationXML */
  SP::TimeDiscretisationXML _timeDiscretisationXML;

public:

  /** Default constructor */
  inline SimulationXML(): rootNode(NULL) {};

  /** Build a SimulationXML object from a DOM tree describing a Simulation
   *   \param the Simulation xml node
   */
  SimulationXML(xmlNodePtr);

  /** destructor
   */
  ~SimulationXML();

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
  inline bool hasOneStepIntegratorXML() const
  {
    return !(OSIXMLSet.empty());
  }

  /** get the set of OneStepNSProblemXML
   *  \return a SetOfOSNSPBXML
   */
  inline const SetOfOSNSPBXML getOneStepNSProblemsXML() const
  {
    return OSNSPBXMLSet;
  };

  /** true if set of OneStepNSProblemXML is not empty
   *  \return a bool
   */
  inline bool hasOneStepNSProblemXML() const
  {
    return !(OSNSPBXMLSet.empty());
  }

  /** Return the TimeDiscretisationXML of the SimulationXML
   *   \return the TimeDiscretisationXML of the SimulationXML
   */
  inline SP::TimeDiscretisationXML timeDiscretisationXML()
  {
    return _timeDiscretisationXML;
  }

  /** Return the type of the Simulation
   *   \return the type of the SimulationXML
   */
  inline std::string  getSimulationXMLType()
  {
    return SiconosDOMTreeTools::getStringAttributeValue(rootNode, TYPE_ATTRIBUTE);
  }

  /** save data of str into the DOM tree
   *   \param xmlNodePtr  : the root node of the SimulationXML
   *   \param SP::Simulation : the Simulation of this SimulationXML
   */
  void saveSimulation2XML(xmlNodePtr  , SP::Simulation);
};

#endif
