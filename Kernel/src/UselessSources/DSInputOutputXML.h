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

/*! \file DSInputOutputXML.h
    \brief

*/
#ifndef __DSInputOutputXML__
#define __DSInputOutputXML__

#include "SiconosDOMTreeTools.h"

class DSInputOutput;

//! XML management for DSInputOutput
/**  \author SICONOS Development Team - copyright INRIA
 *   \version 2.1.0.
 *   \date 17/01/2005
 *
 * DSInputOutputXML allows to manage data of a DSInputOutput DOM tree.
 */
class DSInputOutputXML
{
public:
  DSInputOutputXML();

  /** Build a DSInputOutputXML object from a DOM tree describing a DSInputOutput
  *   \param xmlNodePtr  : the DSInputOutput DOM tree
  */
  DSInputOutputXML(xmlNodePtr);

  /** Destructor */
  virtual ~DSInputOutputXML();

  /** Return the computeInput Plugin name of the DSInputOutputXML
  *   \return The computeInput Plugin name of the DSInputOutputXML
  *  \exception XMLException
  */
  inline std::string getComputeInputPlugin()
  {
    if (!isComputeInputPlugin())
      XMLException::selfThrow("DSInputOutputXML - getComputeInputPlugin : computeInput is not calculated from a plugin ; Fint vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(this->computeInputNode, PLUGIN_ATTRIBUTE);
  }

  /** Return the computeOutput Plugin name of the DSInputOutputXML
  *   \return The computeOutput Plugin name of the DSInputOutputXML
  *  \exception XMLException
  */
  inline std::string getComputeOutputPlugin()
  {
    if (!isComputeOutputPlugin())
      XMLException::selfThrow("DSInputOutputXML - getComputeOutputPlugin : computeOutput is not calculated from a plugin ; Fint vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(this->computeOutputNode, PLUGIN_ATTRIBUTE);
  }

  /** sets the computeInput Plugin name of the DSInputOutputXML
  *   \param string :  The computeInput Plugin name of the DSInputOutputXML
  *  \exception XMLException
  */
  inline void setComputeInputPlugin(std::string plugin)
  {
    if (this->computeInputNode == NULL)
    {
      this->computeInputNode = SiconosDOMTreeTools::createSingleNode(this->rootDSIOXMLNode, COMPUTE_INPUT_TAG);
      xmlNewProp(this->computeInputNode, (xmlChar*)(PLUGIN_ATTRIBUTE.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(this->computeInputNode, PLUGIN_ATTRIBUTE, plugin);
  }

  /** sets the computeOutput Plugin name of the DSInputOutputXML
  *   \param string :  The computeOutput Plugin name of the DSInputOutputXML
  *  \exception XMLException
  */
  inline void setComputeOutputPlugin(std::string plugin)
  {
    if (this->computeOutputNode == NULL)
    {
      this->computeOutputNode = SiconosDOMTreeTools::createSingleNode(this->rootDSIOXMLNode, COMPUTE_OUTPUT_TAG);
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


  /** Return the number of the DSInputOutputXML
  *   \return The integer number of the DSInputOutputXML
  */
  inline int getNumber()
  {
    return SiconosDOMTreeTools::getAttributeValue<int>(rootDSIOXMLNode, NUMBER_ATTRIBUTE);
  }

  /** Return the type of the DSInputOutputXML
  *   \return The string type of the DSInputOutputXML
  */
  inline std::string getType()
  {
    std::string type((char*)this->rootDSIOXMLNode->name);
    return type;
  }

  /** Return the node of the DSInputOutputXML
  *   \return xmlNodePtr  : the node of the RelationXML in the DOM tree
  */
  inline xmlNodePtr  getNode()const
  {
    return this->rootDSIOXMLNode;
  }


  //    /** Return H matrix of the DSInputOutputXML
  //    *   \return SimpleMatrix : the H matrix of the DSInputOutputXML
  //    */
  //    inline SimpleMatrix getH()
  //    {
  //      return SiconosDOMTreeTools::getSiconosMatrixValue(this->HNode);
  //    }
  //
  //    /** allows to save the H matrix of the DSInputOutputXML
  //    *   \param SiconosMatrix* : the H to save
  //    */
  //    inline void setH(SiconosMatrix *H)
  //    {
  //      if( this->HNode == NULL )
  //      {
  //        this->HNode = SiconosDOMTreeTools::createMatrixNode(this->rootDSIOXMLNode, DSINPUTOUTPUT_H, H);
  //      }
  //      else SiconosDOMTreeTools::setSiconosMatrixValue(this->HNode, H);
  //    }

  /** makes the operations to create the DSInputOutput of the DynamicalSystem
  *   \param xmlNodePtr  : the root node of the DSInputOutputXML
  *   \param DSInputOutput* : the Relation of this DSInputOutputXML
  */
  void updateDSInputOutputXML(xmlNodePtr  node, DSInputOutput* dsio);

  /** Return the DSs concerned by the DSInputOutputXML
  *   \return the integer vector who contains the DSs concerned by the DSInputOutputXML
  */
  inline std::vector<int> getDSConcerned()
  {
    return this->definedDSNumbers;
  }

  /** allows to set the dynamical systems which are interacting together with this DSInputOutput
  *   \param vector<int> : the dynamical system numbers
  */
  void setDSConcerned(std::vector<int>);


protected:
  xmlNode * rootDSIOXMLNode;

  //    xmlNode * HNode;
  xmlNode * dsConcernedNode;
  xmlNode * computeInputNode;
  xmlNode * computeOutputNode;

  /* vector of DS numbers*/
  std::vector<int> definedDSNumbers;


private :

  /** load the DSs concerned by this interaction
  *   \param xmlNode * DSConcernedNode : the DOM tree node of DS concerned by the interaction
  *   \exception XMLException : if a DS number not exists
  */
  void loadDSIOConcernedDS(xmlNode * DSConcernedNode);
};


#endif
