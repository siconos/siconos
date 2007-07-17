/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
/*! \file EqualityConstraintXML.h

*/

#ifndef EQUALITYCONSTRAINTXML_H
#define EQUALITYCONSTRAINTXML_H

#include "SiconosDOMTreeTools.h"
#include "SimpleMatrix.h"

const std::string EQUALITYCONSTRAINT_G = "G";
const std::string EQUALITYCONSTRAINT_DSIO_CONCERNED = "DSInputOutput_Concerned";
class SimpleMatrix;
class EqualityConstraint;

//! XML management for EqualityConstraint
/**  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date 17/01/2005
 *
 *
 */
class EqualityConstraintXML
{
public:

  EqualityConstraintXML();

  /** Build a EqualityConstraintXML object from a DOM tree describing a EqualityConstraint
  *   \param xmlNodePtr  : the EqualityConstraint DOM tree
  *   \param vector<int>  : vector of DSXML numbers to verify DS concerned by the EqualityConstraint (identified by number) exists
  */
  EqualityConstraintXML(xmlNodePtr , std::vector<int>);
  virtual ~EqualityConstraintXML();


  /** Return the number of the EqualityConstraintXML
  *   \return The integer number of the EqualityConstraintXML
  */
  inline int getNumber()
  {
    return SiconosDOMTreeTools::getAttributeValue<int>(rootNode, NUMBER_ATTRIBUTE);
  }


  /** Return the type of the EqualityConstraintXML
  *   \return The string type of the EqualityConstraintXML
  */
  inline std::string getType()
  {
    std::string res((char*)this->rootNode->name);
    return res;
  }

  /** Return G matrix of the EqualityConstraintXML
  *   \return SimpleMatrix : the G matrix of the EqualityConstraintXML
  */
  inline SimpleMatrix getG()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(this->GNode);
  }

  /** allows to save the G matrix of the EqualityConstraintXML
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


  /** Return the DSIOs concerned by the EqualityConstraintXML
  *   \return the integer vector who contains the DSs concerned by the EqualityConstraintXML
  */
  inline std::vector<int> getDSIOConcerned()
  {
    return this->definedDSIONumbers;
  }

  /** allows to set the dynamical systems which are interacting together with this EqualityConstraintXML
  *   \param vector<int> : the dynamical system numbers
  */
  void setDSIOConcerned(std::vector<int>);

  /** makes the operations to create the EqualityConstraint of the Non Smooth Dynamical System
  *   \param xmlNodePtr  : the root node of the EqualityConstraintXML
  *   \param EqualityConstraint* : the Relation of this EqualityConstraintXML
  */
  void updateEqualityConstraintXML(xmlNodePtr  node, EqualityConstraint* ec);


  ///////////////////////////////
  /** Return the computeInput Plugin name of the EqualityConstraintXML
  *   \return The computeInput Plugin name of the EqualityConstraintXML
  *  \exception XMLException
  */
  inline std::string getComputeInputPlugin()
  {
    if (!isComputeInputPlugin())
      XMLException::selfThrow("EqualityConstraintXML - getComputeInputPlugin : computeInput is not calculated from a plugin ; Fint vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(this->computeInputNode, PLUGIN_ATTRIBUTE);
  }

  /** Return the computeOutput Plugin name of the EqualityConstraintXML
  *   \return The computeOutput Plugin name of the EqualityConstraintXML
  *  \exception XMLException
  */
  inline std::string getComputeOutputPlugin()
  {
    if (!isComputeOutputPlugin())
      XMLException::selfThrow("EqualityConstraintXML - getComputeOutputPlugin : computeOutput is not calculated from a plugin ; Fint vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(this->computeOutputNode, PLUGIN_ATTRIBUTE);
  }

  /** sets the computeInput Plugin name of the EqualityConstraintXML
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

  /** sets the computeOutput Plugin name of the EqualityConstraintXML
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

  /** Return true if computeInput is calculated from a plugin
  *   \return True if computeInput is calculated from plugin
  */
  inline bool isComputeInputPlugin()
  {
    return xmlHasProp((xmlNodePtr)computeInputNode, (xmlChar *) PLUGIN_ATTRIBUTE.c_str());
  }

  /** Return true if computeOutput is calculated from a plugin
  *   \return True if computOutput is calculated from plugin
  */
  inline bool isComputeOutputPlugin()
  {
    return xmlHasProp((xmlNodePtr)computeOutputNode, (xmlChar *) PLUGIN_ATTRIBUTE.c_str());
  }

  /** return true if computeInputNode is defined
  *  \return true if computeInputNode is defined
  */
  inline bool hasComputeInput()
  {
    return (this->computeInputNode != NULL);
  }

  /** return true if computeOutputNode is defined
  *  \return true if computeOutputNode is defined
  */
  inline bool hasComputeOutput()
  {
    return (this->computeOutputNode != NULL);
  }

protected :
  xmlNodePtr  rootNode;
  xmlNodePtr  GNode;
  xmlNodePtr  dsioConcernedNode;

  xmlNodePtr  computeInputNode;
  xmlNodePtr  computeOutputNode;

  /* vector of DSIO numbers*/
  std::vector<int> definedDSIONumbers;

private :
  /** load the DSs concerned by this interaction
  *   \param xmlNodePtr  : the DOM tree node of DSIO concerned by the EqualityConstraint
  // *   \param vector<int> : vector of DSXML numbers to verify DS concerned by the interaction (identified by number) exists
  */
  void loadECConcernedDSIO(xmlNodePtr  DSIOConcernedNode/*, std::vector<int> definedDSNumbers*/);
};

#endif // EQUALITYCONSTRAINTXML_H

