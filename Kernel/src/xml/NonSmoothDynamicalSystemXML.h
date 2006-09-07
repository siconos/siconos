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
/** \class NonSmoothDynamicalSystemXML
 *   \brief xml data management for NonSmoothDynamicalSystem
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
 *  \date 04/04/2004
 *
 * At the time, only DynamicalSystem and Interaction objects are handles by NonSmoothDynamicalSystem.
 * DSIO and EqualityConstraint are not. Thus all objects and functions related to those two classes are to be reviewed in the present class.
 * This should be done, if possible, in the same way as for DS and Interactions (with set in the place of map etc ...)
 *
 */

#ifndef __NonSmoothDynamicalSystemXML__
#define __NonSmoothDynamicalSystemXML__

#include "SiconosDOMTreeTools.h"

#include "EqualityConstraintXML.h"
#include "InteractionXML.h"
#include "DSInputOutputXML.h"
#include "DynamicalSystemXML.h"
#include <set>
#include <map>

class NonSmoothDynamicalSystem;
class InteractionXML;
class EqualityConstraintXML;
class DynamicalSystemXML;
class DSInputOutputXML;

const std::string NSDS_BVP = "bvp";

/** set of DSXML */
typedef std::set<DynamicalSystemXML*> SetOfDSXML;

/** iterator through setOfDSXML */
typedef SetOfDSXML::iterator SetOfDSXMLIt;

/** set of InteractionXML */
typedef std::set<InteractionXML*> SetOfInteractionsXML;

/** iterator through SetOfInteractionXML */
typedef SetOfInteractionsXML::iterator SetOfInteractionsXMLIt;

class NonSmoothDynamicalSystemXML
{
private:
  xmlNodePtr rootNode;

  /* Map of DSs */
  SetOfDSXML DSXMLSet;

  /* set of InteractionXML */
  SetOfInteractionsXML interactionsXMLSet;

  /* Map of EqualityConstraints */
  std::map<int, EqualityConstraintXML*> equalityConstraintXMLMap;

  /* Map of DSInputOutputs */
  std::map<int, DSInputOutputXML*> dsInputOutputXMLMap;

  /* vector of DSInputOutput numbers*/
  std::vector<int> definedDSInputOutputNumbers;

  /* vector of EqualityConstraint numbers*/
  std::vector<int> definedEqualityConstraintNumbers;

  /** \fn loadDynamicalSystemXML(xmlNodePtr  rootDSNode)
   *   \brief Builds DynamicalSystemXML objects from a DOM tree describing DSs
   *   \param xmlNodePtr  : the DSs DOM tree
   *   \exception XMLException : if a property of the NonSmoothDynamicalSystem lacks in the DOM tree
   */
  void loadDynamicalSystemXML(xmlNodePtr  rootDSNode);

  void loadNonSmoothDynamicalSystem();

  /** \fn void loadEqualityConstraintXML(xmlNodePtr  rootECNode)
   *   \brief Builds EqualityConstraintXML objects from a DOM tree describing EqualityConstraints
   *   \param xmlNodePtr  : the EqualityConstraints DOM tree
   *   \exception XMLException : if a number relating to an EqualityConstraint declares in the NonSmoothDynamicalSystem is already used
   */
  void loadEqualityConstraintXML(xmlNodePtr  rootECNode);

  /** \fn void loadDSInputOutputXML(xmlNodePtr  )
   *   \brief Builds DSInputOutputXML objects from a DOM tree describing DSInputOutputs
   *   \param xmlNodePtr  : the DSInputOutputs DOM tree
   *   \exception XMLException : if a number relating to an DSInputOutput declares in the NonSmoothDynamicalSystem is already used
   */
  void loadDSInputOutputXML(xmlNodePtr  rootdsioNode);

  /** \fn map<int, DSInputOutputXML*> getDSInputOutputXMLRelatingToDS( int number )
   *   \brief selects the DSInputOutputXML objects relating to a specific DynamicalSystem
   *   \param map<int, DSInputOutputXML*> : the map containing the DSInputOutputXML for a specific DynamicalSystem
   */
  std::map<int, DSInputOutputXML*> getDSInputOutputXMLRelatingToDS(int number);

public:

  /** \fn  NonSmoothDynamicalSystemXML();
   *   \brief Default constructor
   */
  NonSmoothDynamicalSystemXML();

  /** \fn NonSmoothDynamicalSystemXML(xmlNodePtr  rootNSDSNode)
   *   \brief Build an NonSmoothDynamicalSystemXML object from a DOM tree describing an NonSmoothDynamicalSystem
   *   \param rootNSDSNode : the NSDS DOM tree
   */
  NonSmoothDynamicalSystemXML(xmlNodePtr  rootNSDSNode);

  /** \fn ~NonSmoothDynamicalSystemXML();
   *   \brief Destructor
   */
  ~NonSmoothDynamicalSystemXML();

  /** \fn xmlNodePtr getRootNode()
   *   \brief Return the root node of the NonSmoothDynamicalSystemXML -> tag NSDS
   *   \return xmlNodePtr  : the root node
   */
  inline xmlNodePtr getRootNode()
  {
    return rootNode;
  };

  /** \fn const bool isBVP() const
   *  \brief check if the NonSmoothDynamicalSystem is BVP or not
   *  \return a bool
   */
  inline const bool isBVP() const
  {
    return SiconosDOMTreeTools::hasAttributeValue(rootNode, NSDS_BVP);
  };

  /** \fn void setBVP(const bool&)
   *   \brief set isbvp attribute of rootNode
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

  /** \fn const SetOfDSXML getDynamicalSystemsXML()
   *   \brief get the set of DynamicalSystemXML
   *   \return a SetOfDSXML
   */
  inline const SetOfDSXML getDynamicalSystemsXML() const
  {
    return DSXMLSet;
  };

  /** \fn const bool hasDynamicalSystemXML() const
   *   \brief true if set of DynamicalSystemXML is not empty
   *   \return a bool
   */
  inline const bool hasDynamicalSystemXML() const
  {
    return !(DSXMLSet.empty());
  }

  /** \fn const SetOfInteractionsXML getInteractionsXML()
   *   \brief get the set of InteractionsXML
   *   \return a SetOfInteractionsXML
   */
  inline const SetOfInteractionsXML getInteractionsXML() const
  {
    return interactionsXMLSet;
  };

  /** \fn const bool hasInteractionsXML() const
   *  \brief true if set of InteractionsXML is not empty
   *  \return a bool
   */
  inline const bool hasInteractionsXML() const
  {
    return !(interactionsXMLSet.empty());
  };

  /** \fn EqualityConstraintXML* getEqualityConstraintXML(int number)
   *   \brief Return the EqualityConstraintXML with id number
   *   \param number : int number : the number of the EqualityConstraintXML to return
   *  \exception XMLException
   *   \return the EqualityConstraintXML of number number, NULL if doesn't exist
   */
  EqualityConstraintXML* getEqualityConstraintXML(int number);

  /** \fn inline vector<int> getEqualityConstraintNumbers();
   *   \brief to get EqualityConstraintXML numbers
   *   \return vector EqualityConstraints integer numbers
   */
  inline std::vector<int> getEqualityConstraintNumbers()
  {
    return definedEqualityConstraintNumbers;
  }

  /** \fn void updateNonSmoothDynamicalSystemXML(xmlNodePtr , NonSmoothDynamicalSystem*)
   *   \brief makes the operations to add a NonSmoothDynamicalSystem to the SiconosModelXML
   *   \param xmlNodePtr  : the root node for the NonSmoothDynamicalSystemXML
   *   \param NonSmoothDynamicalSystem* : the NonSmoothDynamicalSystem of this NonSmoothDynamicalSystemXML
   */
  void updateNonSmoothDynamicalSystemXML(xmlNodePtr , NonSmoothDynamicalSystem*);

  /** \fn void loadNonSmoothDynamicalSystem( NonSmoothDynamicalSystem* )
   *   \brief loads the data of the NonSmoothDynamicalSystem into the NonSmoothDynamicalSystemXML (in the DOM tree)
   *   \param NonSmoothDynamicalSystem* : the NonSmoothDynamicalSystem of this NonSmoothDynamicalSystemXML
   */
  void loadNonSmoothDynamicalSystem(NonSmoothDynamicalSystem*);


};



#endif
