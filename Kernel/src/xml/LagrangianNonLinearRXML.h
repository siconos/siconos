
/** \class LagrangianLinearRXML
*   \brief This class manages LagrangianLinear Relation data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 08/12/2004
*
*
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
