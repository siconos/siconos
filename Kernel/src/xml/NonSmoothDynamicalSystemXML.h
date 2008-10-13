/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
/*! \file NonSmoothDynamicalSystemXML.h

*/


#ifndef __NonSmoothDynamicalSystemXML__
#define __NonSmoothDynamicalSystemXML__

#include "SiconosDOMTreeTools.h"
#include <set>

#include "SiconosPointers.h"

class NonSmoothDynamicalSystem;
class InteractionXML;
class DynamicalSystemXML;
const std::string NSDS_BVP = "bvp";

/** set of DSXML */
typedef std::set<SP::DynamicalSystemXML> SetOfDSXML;

/** set of InteractionXML */
typedef std::set<SP::InteractionXML> SetOfInteractionsXML;

/** iterator through SetOfInteractionXML */
typedef SetOfInteractionsXML::iterator SetOfInteractionsXMLIt;

/** XML management for NonSmoothDynamicalSystem
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date 04/04/2004
 *
 */
class NonSmoothDynamicalSystemXML
{
private:
  xmlNodePtr rootNode;

  /* Map of DSs */
  SetOfDSXML DSXMLSet;

  /* set of InteractionXML */
  SetOfInteractionsXML interactionsXMLSet;

  /** Builds DynamicalSystemXML objects from a DOM tree describing DSs
   *   \param xmlNodePtr  : the DSs DOM tree
   *   \exception XMLException : if a property of the NonSmoothDynamicalSystem lacks in the DOM tree
   */
  void loadDynamicalSystemXML(xmlNodePtr  rootDSNode);

  void loadNonSmoothDynamicalSystem();

public:

  /** Default constructor
  */
  NonSmoothDynamicalSystemXML();

  /** Build an NonSmoothDynamicalSystemXML object from a DOM tree describing an NonSmoothDynamicalSystem
  *   \param rootNSDSNode : the NSDS DOM tree
  */
  NonSmoothDynamicalSystemXML(xmlNodePtr  rootNSDSNode);

  /** Destructor
  */
  ~NonSmoothDynamicalSystemXML();

  /** Return the root node of the NonSmoothDynamicalSystemXML -> tag NSDS
  *   \return xmlNodePtr  : the root node
  */
  inline xmlNodePtr getRootNode()
  {
    return rootNode;
  };

  /** check if the NonSmoothDynamicalSystem is BVP or not
  *  \return a bool
  */
  inline const bool isBVP() const
  {
    return SiconosDOMTreeTools::hasAttributeValue(rootNode, NSDS_BVP);
  };

  /** set isbvp attribute of rootNode
  *   \param a bool
  */
  inline void setBVP(const bool& b)
  {
    if (!(SiconosDOMTreeTools::hasAttributeValue(rootNode, NSDS_BVP)))
    {
      if (b) xmlNewProp(rootNode, (xmlChar*)NSDS_BVP.c_str(), (xmlChar*)"true");
    }
    else
    {
      if (!b) xmlRemoveProp(xmlHasProp(rootNode, (xmlChar*)NSDS_BVP.c_str()));
    }
  }

  /** get the set of DynamicalSystemXML
  *   \return a SetOfDSXML
  */
  inline const SetOfDSXML getDynamicalSystemsXML() const
  {
    return DSXMLSet;
  };

  /** true if set of DynamicalSystemXML is not empty
  *   \return a bool
  */
  inline const bool hasDynamicalSystemXML() const
  {
    return !(DSXMLSet.empty());
  }

  /** get the set of InteractionsXML
  *   \return a SetOfInteractionsXML
  */
  inline const SetOfInteractionsXML getInteractionsXML() const
  {
    return interactionsXMLSet;
  };

  /** true if set of InteractionsXML is not empty
  *  \return a bool
  */
  inline const bool hasInteractionsXML() const
  {
    return !(interactionsXMLSet.empty());
  };

  /** makes the operations to add a NonSmoothDynamicalSystem to the SiconosModelXML
  *   \param xmlNodePtr  : the root node for the NonSmoothDynamicalSystemXML
  *   \param SP::NonSmoothDynamicalSystem : the NonSmoothDynamicalSystem of this NonSmoothDynamicalSystemXML
  */
  void updateNonSmoothDynamicalSystemXML(xmlNodePtr , SP::NonSmoothDynamicalSystem);

  /** loads the data of the NonSmoothDynamicalSystem into the NonSmoothDynamicalSystemXML (in the DOM tree)
  *   \param SP::NonSmoothDynamicalSystem : the NonSmoothDynamicalSystem of this NonSmoothDynamicalSystemXML
  */
  void loadNonSmoothDynamicalSystem(SP::NonSmoothDynamicalSystem);


};



#endif
