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

/** \class LagrangianDSIOXML
*   \brief This class manages LagrangianDSIO DSInputOutput data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.1.2.
*   \date 17/01/2005
*
*
*
* LagrangianDSIOXML allows to manage data of a LagrangianDSIO DOM tree.
*/


#ifndef __LNLRelationXML__
#define __LNLRelationXML__

#include "DSInputOutputXML.h"

class LagrangianDSIOXML : public DSInputOutputXML
{
public:

  LagrangianDSIOXML();

  /** \fn LagrangianDSIOXML(xmlNode * , vector<int> )
  *   \brief Build a DSInputOutputXML object from a DOM tree describing a DSInputOutput
  *   \param xmlNode* : the DSInputOutput DOM tree
  //    *   \param vector<int>  : vector of DSXML numbers to verify DS concerned by the DSInputOutput (identified by number) exists
  */
  LagrangianDSIOXML(xmlNode * dsioNode/*, vector<int> definedDSNumbers */);
  ~LagrangianDSIOXML();

  //    /** \fn int getComputeInputPlugin()
  //    *   \brief Returns the computeInput plugin of the LagrangianDSIOXML
  //    *   \return string which defines the plugin
  //    */
  //    inline string getComputeInputPlugin()
  //    {
  //      return  SiconosDOMTreeTools::getStringAttributeValue(this->computeInputNode, COMPUTE_INPUT_TAG);
  //    }
  //
  //    /** \fn inline void setComputeInputPlugin(string plugin)
  //    *   \brief saves the the computeInput plugin of the LagrangianDSIOXML
  //    *   \param The string corresponding to the plugin to save
  //    */
  //    inline void setComputeInputPlugin(string plugin)
  //    {
  //      if( this->computeInputNode == NULL )
  //      {
  //        this->computeInputNode = SiconosDOMTreeTools::createSingleNode(this->rootDSIOXMLNode, COMPUTE_INPUT_TAG);
  //        xmlNewProp(this->computeInputNode, (xmlChar*)(PLUGIN_ATTRIBUTE.c_str()), (xmlChar*)plugin.c_str() );
  //      }
  //      else SiconosDOMTreeTools::setStringAttributeValue(this->computeInputNode, PLUGIN_ATTRIBUTE, plugin);
  //    }
  //
  //    /** \fn int getComputeOutputPlugin()
  //    *   \brief Returns the computeOutput plugin of the LagrangianDSIOXML
  //    *   \return string which defines the plugin
  //    */
  //    inline string getComputeOutputPlugin()
  //    {
  //      return  SiconosDOMTreeTools::getStringAttributeValue(this->computeOutputNode, COMPUTE_OUTPUT_TAG);
  //    }
  //
  //    /** \fn inline void setComputeOutputPlugin(string plugin)
  //    *   \brief saves the the computeOutput plugin of the LagrangianDSIOXML
  //    *   \param The string corresponding to the plugin to save
  //    */
  //    inline void setComputeOutputPlugin(string plugin)
  //    {
  //      if( this->computeOutputNode == NULL )
  //      {
  //        this->computeOutputNode = SiconosDOMTreeTools::createSingleNode(this->rootDSIOXMLNode, COMPUTE_OUTPUT_TAG);
  //        xmlNewProp(this->computeOutputNode, (xmlChar*)(PLUGIN_ATTRIBUTE.c_str()), (xmlChar*)plugin.c_str() );
  //      }
  //      else SiconosDOMTreeTools::setStringAttributeValue(this->computeOutputNode, PLUGIN_ATTRIBUTE, plugin);
  //    }

private:
  //Nodes
};


#endif
