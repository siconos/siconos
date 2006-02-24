/* Siconos-Kernel version 1.1.2, Copyright INRIA 2005-2006.
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
 *   \brief This class manages NonSmoothDynamicalSystem data part
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.1.2.
 *   \date 04/04/2004
 *
 *
 *
 * NonSmoothDynamicalSystemXML allows to get DSXMLs and InteractionXMLs from a DOM tree.
 */

#ifndef __NonSmoothDynamicalSystemXML__
#define __NonSmoothDynamicalSystemXML__

#include "NonSmoothDynamicalSystem.h"

#include "SiconosDOMTreeTools.h"

#include "EqualityConstraintXML.h"
#include "InteractionXML.h"
#include "DSInputOutputXML.h"
#include "DynamicalSystemXML.h"

class NonSmoothDynamicalSystem;
class InteractionXML;
class EqualityConstraintXML;
class DynamicalSystemXML;
class DSInputOutputXML;

const std::string NSDS_BVP = "bvp";

class NonSmoothDynamicalSystemXML
{
public:
  NonSmoothDynamicalSystemXML();

  /** \fn NonSmoothDynamicalSystemXML(xmlNode * rootNSDSNode)
   *   \brief Build an NonSmoothDynamicalSystemXML object from a DOM tree describing an NonSmoothDynamicalSystem
   *   \param rootNSDSNode : the NSDS DOM tree
   */
  NonSmoothDynamicalSystemXML(xmlNode * rootNSDSNode);

  ~NonSmoothDynamicalSystemXML();

  ///* return a DSs Map */
  //inline map<int, DynamicalSystemXML> getDSMap() const;

  /** \fn xmlNode* getNonSmoothDynamicalSystemXMLNode()
   *   \brief Return the root node of the NonSmoothDynamicalSystemXML
   *   \return xmlNode* : the root node
   */
  inline xmlNode* getNonSmoothDynamicalSystemXMLNode()
  {
    return NonSmoothDynamicalSystemNode;
  }

  /** \fn DynamicalSystemXML* getDynamicalSystemXML(int number)
   *   \brief Return the DynamicalSystemXML with id number
   *   \param int number : the number of the DynamicalSystemXML to return
   *   \return the DynamicalSystemXML of number number NULL if doesn't exist
   */
  DynamicalSystemXML* getDynamicalSystemXML(int number);

  /** \fn InteractionXML* getInteractionXML(int number)
   *   \brief Return the InteracionXML with id number
   *   \param number : int number : the number of the InteractionXML to return
   *  \exception XMLException
   *   \return the InteractionXML of number number, NULL if doesn't exist
   */
  InteractionXML* getInteractionXML(int number);

  /** \fn EqualityConstraintXML* getEqualityConstraintXML(int number)
   *   \brief Return the EqualityConstraintXML with id number
   *   \param number : int number : the number of the EqualityConstraintXML to return
   *  \exception XMLException
   *   \return the EqualityConstraintXML of number number, NULL if doesn't exist
   */
  EqualityConstraintXML* getEqualityConstraintXML(int number);


  /** \fn inline vector<int> getDSNumbers();
   *   \brief Allows to know the defined DS
   *  \exception XMLException
   *   \return vector DS numbers
   */
  inline std::vector<int> getDSNumbers()
  {
    return definedDSNumbers;
  }

  /** \fn inline vector<int> getInteractionNumbers();
   *   \brief Allows to know the defined interactions
   *   \return vector Interactions integer numbers
   */
  inline std::vector<int> getInteractionNumbers()
  {
    return definedInteractionNumbers;
  }

  /** \fn inline vector<int> getEqualityConstraintNumbers();
   *   \brief Allows to know the defined EqualityConstraints
   *   \return vector EqualityConstraints integer numbers
   */
  inline std::vector<int> getEqualityConstraintNumbers()
  {
    return definedEqualityConstraintNumbers;
  }


  /** \fn bool isBVP()
   *   \brief Allows to know if the NonSmoothDynamicalSystem is BVP or not
   *   \return True if the NonSmoothDynamicalSystem is BVP false otherwise
   */
  inline bool isBVP()
  {
    if (SiconosDOMTreeTools::hasAttributeValue(this->NonSmoothDynamicalSystemNode, NSDS_BVP))
      return SiconosDOMTreeTools::getBooleanAttributeValue(NonSmoothDynamicalSystemNode, NSDS_BVP);
    else return false;
  }

  /** \fn void setBVP(bool)
   *   \brief Allows to define if the NonSmoothDynamicalSystem is BVP
   *   \param True if the NonSmoothDynamicalSystem is BVP false otherwise
   */
  inline void setBVP(bool b)
  {
    if (!(SiconosDOMTreeTools::hasAttributeValue(this->NonSmoothDynamicalSystemNode, NSDS_BVP)))
    {
      //if( b ) SiconosDOMTreeTools::setBooleanAttributeValue(this->NonSmoothDynamicalSystemNode, NSDS_BVP, true);
      if (b == true) xmlNewProp(NonSmoothDynamicalSystemNode, (xmlChar*)NSDS_BVP.c_str(), (xmlChar*)"true");
    }
    else
    {
      if (b == false) xmlRemoveProp(xmlHasProp(NonSmoothDynamicalSystemNode, (xmlChar*)NSDS_BVP.c_str()));
    }
  }

  /** \fn void updateNonSmoothDynamicalSystemXML(xmlNode*, NonSmoothDynamicalSystem*)
   *   \brief makes the operations to add a NonSmoothDynamicalSystem to the SiconosModelXML
   *   \param xmlNode* : the root node for the NonSmoothDynamicalSystemXML
   *   \param NonSmoothDynamicalSystem* : the NonSmoothDynamicalSystem of this NonSmoothDynamicalSystemXML
   */
  void updateNonSmoothDynamicalSystemXML(xmlNode*, NonSmoothDynamicalSystem*);

  /** \fn void loadNonSmoothDynamicalSystem( NonSmoothDynamicalSystem* )
   *   \brief loads the data of the NonSmoothDynamicalSystem into the NonSmoothDynamicalSystemXML (in the DOM tree)
   *   \param NonSmoothDynamicalSystem* : the NonSmoothDynamicalSystem of this NonSmoothDynamicalSystemXML
   */
  void loadNonSmoothDynamicalSystem(NonSmoothDynamicalSystem*);


private:
  xmlNode *NonSmoothDynamicalSystemNode;

  /* Map of DSs */
  std::map<int, DynamicalSystemXML*> DSXMLMap;

  /* Map of interactions */
  std::map<int, InteractionXML*> interactionXMLMap;

  /* Map of EqualityConstraints */
  std::map<int, EqualityConstraintXML*> equalityConstraintXMLMap;

  /* Map of DSInputOutputs */
  std::map<int, DSInputOutputXML*> dsInputOutputXMLMap;


  /* vector of DSInputOutput numbers*/
  std::vector<int> definedDSInputOutputNumbers;

  /* vector of DS numbers*/
  std::vector<int> definedDSNumbers;

  /* vector of Interaction numbers*/
  std::vector<int> definedInteractionNumbers;

  /* vector of EqualityConstraint numbers*/
  std::vector<int> definedEqualityConstraintNumbers;

  void loadNonSmoothDynamicalSystem();

  /** \fn loadDynamicalSystemXML(xmlNode * rootDSNode)
   *   \brief Builds DynamicalSystemXML objects from a DOM tree describing DSs
   *   \param xmlNode* : the DSs DOM tree
   *   \exception XMLException : if a property of the NonSmoothDynamicalSystem lacks in the DOM tree
   */
  void loadDynamicalSystemXML(xmlNode * rootDSNode);

  /** \fn loadInteractionXML(xmlNode * rootInteractionNode)
   *   \brief Builds InteractionXML objects from a DOM tree describing Interactions
   *   \param xmlNode* : the Interactions DOM tree
   *   \exception XMLException : if a number relating to an Interaction declares in the NonSmoothDynamicalSystem is already used
   */
  void loadInteractionXML(xmlNode * rootInteractionNode);

  /** \fn void loadEqualityConstraintXML(xmlNode * rootECNode)
   *   \brief Builds EqualityConstraintXML objects from a DOM tree describing EqualityConstraints
   *   \param xmlNode* : the EqualityConstraints DOM tree
   *   \exception XMLException : if a number relating to an EqualityConstraint declares in the NonSmoothDynamicalSystem is already used
   */
  void loadEqualityConstraintXML(xmlNode * rootECNode);

  /** \fn void loadDSInputOutputXML(xmlNode * )
   *   \brief Builds DSInputOutputXML objects from a DOM tree describing DSInputOutputs
   *   \param xmlNode* : the DSInputOutputs DOM tree
   *   \exception XMLException : if a number relating to an DSInputOutput declares in the NonSmoothDynamicalSystem is already used
   */
  void loadDSInputOutputXML(xmlNode * rootdsioNode);

  /** \fn map<int, DSInputOutputXML*> getDSInputOutputXMLRelatingToDS( int number )
   *   \brief selects the DSInputOutputXML objects relating to a specific DynamicalSystem
   *   \param map<int, DSInputOutputXML*> : the map containing the DSInputOutputXML for a specific DynamicalSystem
   */
  std::map<int, DSInputOutputXML*> getDSInputOutputXMLRelatingToDS(int number);

};



#endif
