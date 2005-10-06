/* Siconos version 1.0, Copyright INRIA 2005.
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
/** \class OneStepIntegratorXML
 *   \brief This class manages OneStepIntegrator data part
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 05/17/2004
 *
 *
 * OneStepIntegratorXML allows to manage data of a OneStepIntegrator DOM tree.
 */

#ifndef __OneStepIntegratorXML__
#define __OneStepIntegratorXML__

#include "SiconosDOMTreeTools.h"
#include "OneStepIntegrator.h"

class OneStepIntegrator;

// Tags
const std::string OSI_R = "r";
const std::string OSI_DS_CONCERNED = "DS_Concerned";

class OneStepIntegratorXML
{
public:

  OneStepIntegratorXML();

  /** \fn OneStepIntegratorXML(xmlNode * OneStepIntegratorNode, map<int, bool> definedDSNumbers)
   *   \brief Build a OneStepIntegratorXML object from a DOM tree describing a OneStepIntegrator
   *   \param OneStepIntegratorNode : the OneStepIntegrator DOM tree
   *   \param map<int, bool> definedDSNumbers : the DS numbers effectivly defined in the model
   *   \exception XMLException : if a property of the OneStepIntegrator lacks in the DOM tree
   */
  OneStepIntegratorXML(xmlNode * OneStepIntegratorNode, std::map<int, bool> definedDSNumbers);

  virtual ~OneStepIntegratorXML();


  /** \fn int getR()
   *   \brief Return r of the OneStepIntegrator
   *   \return The r integer of the OneStepIntegrator
   */
  inline int getR()
  {
    return SiconosDOMTreeTools::getIntegerContentValue(this->rNode);
  }

  /** \fn void setR(int r)
   *   \brief allows to save r of the OneStepIntegrator
   *   \param The r integer to save
   */
  inline void setR(int r)
  {
    if (this->hasR() == false)
    {
      this->rNode = SiconosDOMTreeTools::createIntegerNode(this->rootIntegratorXMLNode, OSI_R, r);
    }
    else SiconosDOMTreeTools::setIntegerContentValue(this->rNode, r);
  }

  /** \fn bool hasR()
   *   \brief determines if the r value is defined for this integrator
   *   \return bool : true if r is defined
   */
  inline bool hasR()
  {
    return (this->rNode != NULL);
  }

  /** \fn vector<int> getDSConcerned()
   *   \brief Return the DS numbers of the OneStepIntegrator
   *   \return The DS numbers vector of the OneStepIntegrator
   */
  inline std::vector<int> getDSConcerned()
  {
    return this->DSNumbersVector;
  }

  /** \fn vector<int> getDSConcerned()
   *   \brief Return the DS numbers of the OneStepIntegrator
   *   \return The DS numbers vector of the OneStepIntegrator
   */
  void setDSConcerned(std::vector<int>* ds);


  /** \fn string getType()
   *   \brief Return the type of the OneStepIntegratorXML
   *   \return The string type of the OneStepIntegratorXML
   */
  inline std::string getType()
  {
    //return SiconosDOMTreeTools::getStringAttributeValue(this->rootIntegratorXMLNode, OSI_TYPE);
    std::string type((char*)this->rootIntegratorXMLNode->name);
    return type;
  }

  /** \fn xmlNode* getNode()
   *   \brief allow to get the root node of the OneStepIntegratorXML
   *   \return xmlNode* : the root node of the OneStepIntegratorXML
   */
  inline xmlNode* getNode() const
  {
    return (xmlNode*)this->rootIntegratorXMLNode;
  }

  /** \fn void updateOneStepIntegratorXML( xmlNode* , OneStepIntegrator*  )
   *   \brief makes the operations to create a OneStepIntegratorXML to the StrategyXML
   *   \param xmlNode* : the root node of the OneStepIntegratorXML
   *   \param Strategy* : the OneStepIntegrator of this OneStepIntegratorXML
   */
  void updateOneStepIntegratorXML(xmlNode* node, OneStepIntegrator* osi);

  /** \fn bool hasAll()
   *  \brief All is an attribute of the DS_Concerned tag
   *  \return bool : true if attribute all is defined
   */
  inline bool hasAll()
  {
    if (SiconosDOMTreeTools::hasAttributeValue(this->DSConcernedNode, ALL_ATTRIBUTE))
      return SiconosDOMTreeTools::getBooleanAttributeValue(this->DSConcernedNode, ALL_ATTRIBUTE);
    else return false;
  }

  /** \fn void setAll(bool all)
   *   \brief Allows to modify the attribute "all" of the DS_concerned tag
   *   \param bool : the value to assign to the attribute
   */
  inline void setAll(bool all)
  {
    if (this->hasAll() == false)
    {
      if (all == true)
        xmlNewProp(this->DSConcernedNode, (xmlChar*)ALL_ATTRIBUTE.c_str(), (xmlChar*)"true");
    }
    else
    {
      if (all == false)
        xmlRemoveProp(xmlHasProp(this->DSConcernedNode, (xmlChar*)ALL_ATTRIBUTE.c_str()));
    }
  }

protected:
  //Nodes
  xmlNode * rNode;
  xmlNode * rootIntegratorXMLNode;
  xmlNode * DSConcernedNode;

private:

  //DSs (DS numbers)
  std::vector<int> DSNumbersVector;

  //Methods


  /** \fn loadOneStepIntegratonConcernedDS(xmlNode * DSConcernedNode, vector<int> definedDSNumbers)
   *   \brief load the DS numbers of the OneStepIntegrator
   *   \param xmlNode * DSConcernedNode : the DOM tree node of the concerned DS
   *   \param map<int, bool> definedDSNumbers : the DS numbers effectivly defined in the model
   *   \exception XMLException : if a DS of the OneStepIntegrator has not been defined in the model or is already used by another OneStepIntegrator
   */
  void loadOneStepIntegratorConcernedDS(xmlNode * DSConcernedNode, std::map<int, bool> definedDSNumbers);

};


#endif
