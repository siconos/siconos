//$Id: LagrangianNonLinearRXML.h,v 1.5 2005/03/07 13:17:21 jbarbier Exp $

/** \class LagrangianLinearRXML
*   \brief This class manages LagrangianLinear Relation data
*   \author J-M Barbier
*   \version 1.0
*   \date 08/12/2004
*
*
* $Date: 2005/03/07 13:17:21 $
* $Revision: 1.5 $
* $Author: jbarbier $
* $Source: /CVS/Siconos/SICONOS/src/xml/LagrangianNonLinearRXML.h,v $
*
* LagrangianNonLinearRXML allows to manage data of a LNLRelation DOM tree.
*/


#ifndef __LNLRelationXML__
#define __LNLRelationXML__


#include <libxml/tree.h>

#include "RelationXML.h"
#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SiconosDOMTreeTools.h"


using namespace std;



class LagrangianNonLinearRXML : public RelationXML
{
public:

  LagrangianNonLinearRXML();

  /** \fn LagrangianNonLinearRXML(xmlNode * LNLRelationNode)
  *   \brief Build a LagrangianNonLinearRXML object from a DOM tree describing a Relation with LNL type
  *   \param LagrangianNonLinearRXML : the LagrangianNonLinearR DOM tree
  *   \exception XMLException : if a property of the LagrangianLinear Relation lacks in the DOM tree
  */
  LagrangianNonLinearRXML(xmlNode * LNLRelationNode);

  ~LagrangianNonLinearRXML();

  /** \fn int getComputeInputPlugin()
  *   \brief Returns the computeInput plugin of the Relation
  *   \return string which defines the plugin
  */
  inline string getComputeInputPlugin()
  {
    return  SiconosDOMTreeTools::getStringAttributeValue(this->computeInputNode, COMPUTE_INPUT_TAG);
  }

  /** \fn inline void setComputeInputPlugin(string plugin)
  *   \brief saves the the computeInput plugin of the Relation
  *   \param The string corresponding to the plugin to save
  */
  inline void setComputeInputPlugin(string plugin)
  {
    if (this->computeInputNode == NULL)
    {
      this->computeInputNode = SiconosDOMTreeTools::createSingleNode(this->rootRelationXMLNode, COMPUTE_INPUT_TAG);
      xmlNewProp(this->computeInputNode, (xmlChar*)(PLUGIN_ATTRIBUTE.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(this->computeInputNode, PLUGIN_ATTRIBUTE, plugin);
  }

  /** \fn int getComputeOutputPlugin()
  *   \brief Returns the computeOutput plugin of the Relation
  *   \return string which defines the plugin
  */
  inline string getComputeOutputPlugin()
  {
    return  SiconosDOMTreeTools::getStringAttributeValue(this->computeOutputNode, COMPUTE_OUTPUT_TAG);
  }

  /** \fn inline void setComputeOutputPlugin(string plugin)
  *   \brief saves the the computeOutput plugin of the Relation
  *   \param The string corresponding to the plugin to save
  */
  inline void setComputeOutputPlugin(string plugin)
  {
    if (this->computeOutputNode == NULL)
    {
      this->computeOutputNode = SiconosDOMTreeTools::createSingleNode(this->rootRelationXMLNode, COMPUTE_OUTPUT_TAG);
      xmlNewProp(this->computeOutputNode, (xmlChar*)(PLUGIN_ATTRIBUTE.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(this->computeOutputNode, PLUGIN_ATTRIBUTE, plugin);
  }

private:
  //Nodes
};


#endif
//$Log: LagrangianNonLinearRXML.h,v $
//Revision 1.5  2005/03/07 13:17:21  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.4  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.3  2004/09/14 13:49:58  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.2  2004/09/10 11:26:27  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.1  2004/08/12 11:55:19  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
