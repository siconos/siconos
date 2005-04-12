
/** \class RelationXML
*   \brief This class manages Relation data part
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 04/13/2004
*
*
*
* RelationXML allows to manage data of a Relation DOM tree.
*/

#ifndef __RelationXML__
#define __RelationXML__


#include <string>
#include <libxml/tree.h>
#include "SiconosDOMTreeTools.h"

#include "Relation.h"



//using namespace std;

class Relation;


//const string RELATION_TYPE = "type";
//const string RELATION_INPUT = "computeInput";
//const string RELATION_OUTPUT = "computeOutput";
//const string RELATION_PLUGIN = "plugin";
//
//const string RELATION_LL = "LL";
//const string RELATION_LNL = "LNL";
//const string RELATION_LTI = "LTI";
#include "XMLTagsName.h"

class RelationXML
{
public:
  RelationXML();
  RelationXML(xmlNode*);
  ~RelationXML();

  /** \fn inline string getComputeInputPlugin()
  *   \brief Return the computeInput Plugin name of the RelationXML
  *   \return The computeInput Plugin name of the RelationXML
  *  \exception XMLException
  */
  inline string getComputeInputPlugin()
  {
    if (this->isComputeInputPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->computeInputNode, PLUGIN_ATTRIBUTE);
    XMLException::selfThrow("RelationXML - getComputeInputPlugin : computeInput is not calculated from a plugin ; Fint vector is given");
  }

  /** \fn inline string getComputeOutputPlugin()
  *   \brief Return the computeOutput Plugin name of the RelationXML
  *   \return The computeOutput Plugin name of the RelationXML
  *  \exception XMLException
  */
  inline string getComputeOutputPlugin()
  {
    if (this->isComputeOutputPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->computeOutputNode, PLUGIN_ATTRIBUTE);
    XMLException::selfThrow("RelationXML - getComputeOutputPlugin : computeOutput is not calculated from a plugin ; Fint vector is given");
  }

  /** \fn void setComputeInputPlugin(string plugin)
  *   \brief sets the computeInput Plugin name of the RelationXML
  *   \param string :  The computeInput Plugin name of the RelationXML
  *  \exception XMLException
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

  /** \fn void setComputeOutputPlugin(string plugin)
  *   \brief sets the computeOutput Plugin name of the RelationXML
  *   \param string :  The computeOutput Plugin name of the RelationXML
  *  \exception XMLException
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


  /** \fn string getType()
  *   \brief Return the type of the RelationXML
  *   \return The string type of the RelationXML
  */
  inline string getType()
  {
    //return SiconosDOMTreeTools::getStringAttributeValue(this->rootRelationXMLNode, RELATION_TYPE);
    string type((char*)this->rootRelationXMLNode->name);
    return type;
  }

  /** \fn xmlNode* getNode()
  *   \brief Return the node of the RelationXML
  *   \return xmlNode* : the node of the RelationXML in the DOM tree
  */
  inline xmlNode* getNode()const
  {
    return this->rootRelationXMLNode;
  }

  /** \fn void updateRelationXML( xmlNode* node, Relation* rel );
  *   \brief makes the operations to create the Relation of the Interaction
  *   \param xmlNode* : the root node of the RelationXML
  *   \param Relation* : the Relation of this RelationXML
  */
  void updateRelationXML(xmlNode* node, Relation* rel);


protected:
  xmlNode * rootRelationXMLNode;

  xmlNode * computeInputNode;
  xmlNode * computeOutputNode;
};


#endif
//$Log: RelationXML.h,v $
//Revision 1.16  2005/03/07 13:17:21  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.15  2004/12/08 12:49:39  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.14  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.13  2004/09/14 13:49:59  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.12  2004/07/29 14:25:45  jbarbier
