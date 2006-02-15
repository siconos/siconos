/* Siconos-Kernel version 1.1.1, Copyright INRIA 2005-2006.
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
#ifndef EQUALITYCONSTRAINTXML_H
#define EQUALITYCONSTRAINTXML_H

#include "SiconosDOMTreeTools.h"
#include "EqualityConstraint.h"

const std::string EQUALITYCONSTRAINT_G = "G";
const std::string EQUALITYCONSTRAINT_DSIO_CONCERNED = "DSInputOutput_Concerned";

class EqualityConstraint;

/** \class EqualityConstraintXML
 *  \brief object to manage XML data of an EqualityConstraint
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date 17/01/2005
 *
 *
 */
class EqualityConstraintXML
{
public:

  EqualityConstraintXML();

  /** \fn EqualityConstraintXML(xmlNode * , vector<int> )
  *   \brief Build a EqualityConstraintXML object from a DOM tree describing a EqualityConstraint
  *   \param xmlNode* : the EqualityConstraint DOM tree
  *   \param vector<int>  : vector of DSXML numbers to verify DS concerned by the EqualityConstraint (identified by number) exists
  */
  EqualityConstraintXML(xmlNode*, std::vector<int>);
  virtual ~EqualityConstraintXML();


  /** \fn int getNumber()
  *   \brief Return the number of the EqualityConstraintXML
  *   \return The integer number of the EqualityConstraintXML
  */
  inline int getNumber()
  {
    return SiconosDOMTreeTools::getIntegerAttributeValue(this->rootNode, NUMBER_ATTRIBUTE);
  }


  /** \fn int getType()
  *   \brief Return the type of the EqualityConstraintXML
  *   \return The string type of the EqualityConstraintXML
  */
  inline std::string getType()
  {
    std::string res((char*)this->rootNode->name);
    return res;
  }

  /** \fn SiconosMatrix getG()
  *   \brief Return G matrix of the EqualityConstraintXML
  *   \return SiconosMatrix : the G matrix of the EqualityConstraintXML
  */
  inline SiconosMatrix getG()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(this->GNode);
  }

  /** \fn void setG(SiconosMatrix *G)
  *   \brief allows to save the G matrix of the EqualityConstraintXML
  *   \param SiconosMatrix* : the G to save
  */
  inline void setG(SiconosMatrix *G)
  {
    if (this->GNode == NULL)
    {
      this->GNode = SiconosDOMTreeTools::createMatrixNode(this->rootNode, EQUALITYCONSTRAINT_G, *G);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(this->GNode, *G);
  }


  /** \fn vector<int> getDSIOConcerned()
  *   \brief Return the DSIOs concerned by the EqualityConstraintXML
  *   \return the integer vector who contains the DSs concerned by the EqualityConstraintXML
  */
  inline std::vector<int> getDSIOConcerned()
  {
    return this->definedDSIONumbers;
  }

  /** \fn void setDSIOConcerned( vector<int> )
  *   \brief allows to set the dynamical systems which are interacting together with this EqualityConstraintXML
  *   \param vector<int> : the dynamical system numbers
  */
  void setDSIOConcerned(std::vector<int>);

  /** \fn void updateEqualityConstraintXML( xmlNode* node, EqualityConstraint* ec );
  *   \brief makes the operations to create the EqualityConstraint of the Non Smooth Dynamical System
  *   \param xmlNode* : the root node of the EqualityConstraintXML
  *   \param EqualityConstraint* : the Relation of this EqualityConstraintXML
  */
  void updateEqualityConstraintXML(xmlNode* node, EqualityConstraint* ec);


  ///////////////////////////////
  /** \fn inline string getComputeInputPlugin()
    *   \brief Return the computeInput Plugin name of the EqualityConstraintXML
    *   \return The computeInput Plugin name of the EqualityConstraintXML
    *  \exception XMLException
    */
  inline std::string getComputeInputPlugin()
  {
    if (!isComputeInputPlugin())
      XMLException::selfThrow("EqualityConstraintXML - getComputeInputPlugin : computeInput is not calculated from a plugin ; Fint vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(this->computeInputNode, PLUGIN_ATTRIBUTE);
  }

  /** \fn inline string getComputeOutputPlugin()
  *   \brief Return the computeOutput Plugin name of the EqualityConstraintXML
  *   \return The computeOutput Plugin name of the EqualityConstraintXML
  *  \exception XMLException
  */
  inline std::string getComputeOutputPlugin()
  {
    if (!isComputeOutputPlugin())
      XMLException::selfThrow("EqualityConstraintXML - getComputeOutputPlugin : computeOutput is not calculated from a plugin ; Fint vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(this->computeOutputNode, PLUGIN_ATTRIBUTE);
  }

  /** \fn void setComputeInputPlugin(string plugin)
  *   \brief sets the computeInput Plugin name of the EqualityConstraintXML
  *   \param string :  The computeInput Plugin name of the EqualityConstraintXML
  *  \exception XMLException
  */
  inline void setComputeInputPlugin(std::string plugin)
  {
    if (computeInputNode == NULL)
    {
      this->computeInputNode = SiconosDOMTreeTools::createSingleNode(this->rootNode, COMPUTE_INPUT_TAG);
      xmlNewProp(this->computeInputNode, (xmlChar*)(PLUGIN_ATTRIBUTE.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(this->computeInputNode, PLUGIN_ATTRIBUTE, plugin);
  }

  /** \fn void setComputeOutputPlugin(string plugin)
  *   \brief sets the computeOutput Plugin name of the EqualityConstraintXML
  *   \param string :  The computeOutput Plugin name of the EqualityConstraintXML
  *  \exception XMLException
  */
  inline void setComputeOutputPlugin(std::string plugin)
  {
    if (this->computeOutputNode == NULL)
    {
      this->computeOutputNode = SiconosDOMTreeTools::createSingleNode(this->rootNode, COMPUTE_OUTPUT_TAG);
      xmlNewProp(this->computeOutputNode, (xmlChar*)(PLUGIN_ATTRIBUTE.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(this->computeOutputNode, PLUGIN_ATTRIBUTE, plugin);
  }

  /** \fn bool isComputeInputPlugin()
  *   \brief Return true if computeInput is calculated from a plugin
  *   \return True if computeInput is calculated from plugin
  */
  inline bool isComputeInputPlugin()
  {
    return xmlHasProp((xmlNodePtr)computeInputNode, (xmlChar *) PLUGIN_ATTRIBUTE.c_str());
  }

  /** \fn bool isComputeOutputPlugin()
  *   \brief Return true if computeOutput is calculated from a plugin
  *   \return True if computOutput is calculated from plugin
  */
  inline bool isComputeOutputPlugin()
  {
    return xmlHasProp((xmlNodePtr)computeOutputNode, (xmlChar *) PLUGIN_ATTRIBUTE.c_str());
  }

  /** \fn bool hasComputeInput()
   *  \brief return true if computeInputNode is defined
   *  \return true if computeInputNode is defined
   */
  inline bool hasComputeInput()
  {
    return (this->computeInputNode != NULL);
  }

  /** \fn bool hasComputeOutput()
   *  \brief return true if computeOutputNode is defined
   *  \return true if computeOutputNode is defined
   */
  inline bool hasComputeOutput()
  {
    return (this->computeOutputNode != NULL);
  }

protected :
  xmlNode * rootNode;
  xmlNode * GNode;
  xmlNode * dsioConcernedNode;

  xmlNode * computeInputNode;
  xmlNode * computeOutputNode;

  /* vector of DSIO numbers*/
  std::vector<int> definedDSIONumbers;

private :
  /** \fn loadECConcernedDSIO(xmlNode * , vector<int>)
  *   \brief load the DSs concerned by this interaction
  *   \param xmlNode * : the DOM tree node of DSIO concerned by the EqualityConstraint
  // *   \param vector<int> : vector of DSXML numbers to verify DS concerned by the interaction (identified by number) exists
  */
  void loadECConcernedDSIO(xmlNode * DSIOConcernedNode/*, std::vector<int> definedDSNumbers*/);
};

#endif // EQUALITYCONSTRAINTXML_H

