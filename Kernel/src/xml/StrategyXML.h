/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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

/** \class StrategyXML
 *   \brief This class manages Strategy data part
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.1.4.
 *   \date 05/17/2004
 *
 *
 *
 * StrategyXML allows to get OneStepIntegratorXMLs and OneStepNSProblemXMLs from a DOM tree.
 */
#ifndef __STRATEGYXML__
#define __STRATEGYXML__

#include "Strategy.h"
#include "SiconosDOMTreeTools.h"

#include "OneStepIntegratorXML.h"
#include "OneStepNSProblemXML.h"
#include "TimeDiscretisationXML.h"
#include <set>

class Strategy;
class OneStepIntegratorXML;
class OneStepNSProblemXML;
class TimeDiscretisationXML;

/** set of OneStepIntegratorXML */
typedef std::set<OneStepIntegratorXML*> SetOfOSIXML;

/** iterator through SetOfInteractionXML */
typedef SetOfOSIXML::iterator SetOfOSIXMLIt;



class StrategyXML
{

private:

  /* root node for strategy */
  xmlNodePtr rootNode;

  /* set of OneStepIntegratorXML* */
  SetOfOSIXML OSIXMLSet;

  /* OneStepNSProblemXML - Optional value, default = NULL */
  OneStepNSProblemXML *oneStepNSProblemXML;

  /* TimeDiscretisationXML */
  TimeDiscretisationXML *timeDiscretisationXML;

public:

  /** \fn StrategyXML()
   *   \brief Default constructor
   */
  StrategyXML();

  /** \fn StrategyXML(xmlNodePtr rootStrategyNode)
   *   \brief Build a StrategyXML object from a DOM tree describing a Strategy
   *   \param the Strategy xml node
   */
  StrategyXML(xmlNodePtr rootStrategyNode);

  /** \fn ~StrategyXML()
   *  \brief destructor
   */
  ~StrategyXML();

  // --- GETTERS/SETTERS ---

  /** \fn xmlNodePtr getRootNode()
   *  \brief Return the root node of the Strategy -> tag "Strategy"
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
   *   \brief Return the OneStepNSProblemXML pointer of the StrategyXML
   *   \return the OneStepNSProblemXML pointer of the StrategyXML ; NULL if StrategyXML does not have
   */
  inline OneStepNSProblemXML * getOneStepNSProblemXMLPtr()
  {
    return oneStepNSProblemXML;
  }

  /** \fn TimeDiscretisationXML* getTimeDiscretisationXMLPtr()
   *   \brief Return the TimeDiscretisationXML of the StrategyXML
   *   \return the TimeDiscretisationXML of the StrategyXML
   */
  inline TimeDiscretisationXML* getTimeDiscretisationXMLPtr()
  {
    return timeDiscretisationXML;
  }

  /** \fn inline string getStrategyXMLType()
   *   \brief Return the type of the Strategy
   *   \return the type of the StrategyXML
   */
  inline std::string  getStrategyXMLType()
  {
    return SiconosDOMTreeTools::getStringAttributeValue(rootNode, TYPE_ATTRIBUTE);
  }

  /** \fn bool hasOneStepNSProblemXML()
   *   \brief determines if the Strategy has a OneStepNSProblemXML
   *   \return bool :  false if the oneStepNSProblemXML* is NULL
   */
  inline const bool hasOneStepNSProblemXML() const
  {
    return (oneStepNSProblemXML != NULL);
  }

  /** \fn void saveStrategy2XML( xmlNodePtr  node, Strategy* str )
   *   \brief save data of str into the DOM tree
   *   \param xmlNodePtr  : the root node of the StrategyXML
   *   \param Strategy* : the Strategy of this StrategyXML
   */
  void saveStrategy2XML(xmlNodePtr  , Strategy*);
};



#endif
