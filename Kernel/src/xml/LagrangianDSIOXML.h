//$Id: LagrangianDSIOXML.h,v 1.5 2005/03/10 12:55:21 jbarbier Exp $

/** \class LagrangianDSIOXML
*   \brief This class manages LagrangianDSIO DSInputOutput data
*   \author J-M Barbier
*   \version 1.0
*   \date 17/01/2005
*
*
* $Date: 2005/03/10 12:55:21 $
* $Revision: 1.5 $
* $Author: jbarbier $
* $Source: /CVS/Siconos/SICONOS/src/xml/LagrangianDSIOXML.h,v $
*
* LagrangianDSIOXML allows to manage data of a LagrangianDSIO DOM tree.
*/


#ifndef __LNLRelationXML__
#define __LNLRelationXML__


#include <libxml/tree.h>

#include "DSInputOutputXML.h"
#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SiconosDOMTreeTools.h"


using namespace std;



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
//$Log: LagrangianDSIOXML.h,v $
//Revision 1.5  2005/03/10 12:55:21  jbarbier
//- implmentation of the EqualityConstraint and DSInputOutput classes in progress
//    attributes H (DSIO) et G (EC) added in XML and managed in XML objects
//
//Revision 1.4  2005/03/09 15:30:37  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.3  2005/03/07 13:17:21  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.2  2005/01/26 13:50:40  jbarbier
//
//- loading of an XML input file now loads EqualityConstraints and DSInputOutputs
//
//Revision 1.1  2005/01/17 10:56:27  jbarbier
//- classes EqualityConstraint and DSInputOutput added with inherited classes
//
//- classes EqualityConstraintXML and DSInputOutputXML added with inherited classes
//
