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
/** \class DSInputOutputXML
 *   \brief This class manages Relation data part
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 17/01/2005
 *
 *
 *
 * DSInputOutputXML allows to manage data of a DSInputOutput DOM tree.
 */

#ifndef __DSInputOutputXML__
#define __DSInputOutputXML__

#include "SiconosDOMTreeTools.h"
#include "DSInputOutput.h"

class DSInputOutput;

class DSInputOutputXML
{
public:
  DSInputOutputXML();

  /** \fn DSInputOutputXML(xmlNode * , vector<int> )
   *   \brief Build a DSInputOutputXML object from a DOM tree describing a DSInputOutput
   *   \param xmlNode* : the DSInputOutput DOM tree
   //   *   \param vector<int>  : vector of DSXML numbers to verify DS concerned by the DSInputOutput (identified by number) exists
   */
  DSInputOutputXML(xmlNode*/*, vector<int>*/);
  virtual ~DSInputOutputXML();

  /** \fn inline string getComputeInputPlugin()
   *   \brief Return the computeInput Plugin name of the DSInputOutputXML
   *   \return The computeInput Plugin name of the DSInputOutputXML
   *  \exception XMLException
   */
  inline std::string getComputeInputPlugin()
  {
    if (!isComputeInputPlugin())
      XMLException::selfThrow("DSInputOutputXML - getComputeInputPlugin : computeInput is not calculated from a plugin ; Fint vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(this->computeInputNode, PLUGIN_ATTRIBUTE);
  }

  /** \fn inline string getComputeOutputPlugin()
   *   \brief Return the computeOutput Plugin name of the DSInputOutputXML
   *   \return The computeOutput Plugin name of the DSInputOutputXML
   *  \exception XMLException
   */
  inline std::string getComputeOutputPlugin()
  {
    if (!isComputeOutputPlugin())
      XMLException::selfThrow("DSInputOutputXML - getComputeOutputPlugin : computeOutput is not calculated from a plugin ; Fint vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(this->computeOutputNode, PLUGIN_ATTRIBUTE);
  }

  /** \fn void setComputeInputPlugin(string plugin)
   *   \brief sets the computeInput Plugin name of the DSInputOutputXML
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

  /** \fn void setComputeOutputPlugin(string plugin)
   *   \brief sets the computeOutput Plugin name of the DSInputOutputXML
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


  /** \fn int getNumber()
   *   \brief Return the number of the DSInputOutputXML
   *   \return The integer number of the DSInputOutputXML
   */
  inline int getNumber()
  {
    return SiconosDOMTreeTools::getIntegerAttributeValue(this->rootDSIOXMLNode, NUMBER_ATTRIBUTE);
  }

  /** \fn string getType()
   *   \brief Return the type of the DSInputOutputXML
   *   \return The string type of the DSInputOutputXML
   */
  inline std::string getType()
  {
    std::string type((char*)this->rootDSIOXMLNode->name);
    return type;
  }

  /** \fn xmlNode* getNode()
   *   \brief Return the node of the DSInputOutputXML
   *   \return xmlNode* : the node of the RelationXML in the DOM tree
   */
  inline xmlNode* getNode()const
  {
    return this->rootDSIOXMLNode;
  }


  //    /** \fn SiconosMatrix getH()
  //    *   \brief Return H matrix of the DSInputOutputXML
  //    *   \return SiconosMatrix : the H matrix of the DSInputOutputXML
  //    */
  //    inline SiconosMatrix getH()
  //    {
  //      return SiconosDOMTreeTools::getSiconosMatrixValue(this->HNode);
  //    }
  //
  //    /** \fn void setH(SiconosMatrix *H)
  //    *   \brief allows to save the H matrix of the DSInputOutputXML
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

  /** \fn void updateDSInputOutputXML( xmlNode* node, DSInputOutput* dsio );
   *   \brief makes the operations to create the DSInputOutput of the DynamicalSystem
   *   \param xmlNode* : the root node of the DSInputOutputXML
   *   \param DSInputOutput* : the Relation of this DSInputOutputXML
   */
  void updateDSInputOutputXML(xmlNode* node, DSInputOutput* dsio);

  /** \fn vector<int> getDSConcerned()
   *   \brief Return the DSs concerned by the DSInputOutputXML
   *   \return the integer vector who contains the DSs concerned by the DSInputOutputXML
   */
  inline std::vector<int> getDSConcerned()
  {
    return this->definedDSNumbers;
  }

  /** \fn void setDSConcerned( vector<int> )
   *   \brief allows to set the dynamical systems which are interacting together with this DSInputOutput
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

  /** \fn loadDSIOConcernedDS(xmlNode * , vector<int>)
   *   \brief load the DSs concerned by this interaction
   *   \param xmlNode * DSConcernedNode : the DOM tree node of DS concerned by the interaction
   //   *   \param vector<int> definedDSNumbers : vector of DSXML numbers to verify DS concerned by the interaction (identified by number) exists
   *   \exception XMLException : if a DS number not exists
   */
  void loadDSIOConcernedDS(xmlNode * DSConcernedNode/*, std::vector<int> definedDSNumbers*/);
};


#endif
