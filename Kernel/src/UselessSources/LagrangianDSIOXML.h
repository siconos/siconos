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

/*! \file LagrangianDSIOXML.h

*/


#ifndef __LNLRelationXML__
#define __LNLRelationXML__

#include "DSInputOutputXML.h"

//! XML management for LagrangianDSIO
/**  \author SICONOS Development Team - copyright INRIA
*   \version 2.1.0.
*   \date 17/01/2005
*
*
*
* LagrangianDSIOXML allows to manage data of a LagrangianDSIO DOM tree.
*/
class LagrangianDSIOXML : public DSInputOutputXML
{
public:

  LagrangianDSIOXML();

  /** Build a DSInputOutputXML object from a DOM tree describing a DSInputOutput
  *   \param xmlNode* : the DSInputOutput DOM tree
  //    *   \param vector<int>  : vector of DSXML numbers to verify DS concerned by the DSInputOutput (identified by number) exists
  */
  LagrangianDSIOXML(xmlNode * dsioNode/*, vector<int> definedDSNumbers */);
  ~LagrangianDSIOXML();

  //    /** Returns the computeInput plugin of the LagrangianDSIOXML
  //    *   \return string which defines the plugin
  //    */
  //    inline string getComputeInputPlugin()
  //    {
  //      return  SiconosDOMTreeTools::getStringAttributeValue(this->computeInputNode, COMPUTE_INPUT_TAG);
  //    }
  //
  //    /** saves the the computeInput plugin of the LagrangianDSIOXML
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
  //    /** Returns the computeOutput plugin of the LagrangianDSIOXML
  //    *   \return string which defines the plugin
  //    */
  //    inline string getComputeOutputPlugin()
  //    {
  //      return  SiconosDOMTreeTools::getStringAttributeValue(this->computeOutputNode, COMPUTE_OUTPUT_TAG);
  //    }
  //
  //    /** saves the the computeOutput plugin of the LagrangianDSIOXML
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
