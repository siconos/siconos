//$Id: EqualityConstraintXML.h,v 1.10 2005/03/21 16:48:04 jbarbier Exp $
#ifndef EQUALITYCONSTRAINTXML_H
#define EQUALITYCONSTRAINTXML_H


#include "SiconosDOMTreeTools.h"
#include "EqualityConstraint.h"
#include "XMLTagsName.h"

const string EQUALITYCONSTRAINT_G = "G";
const string EQUALITYCONSTRAINT_DSIO_CONCERNED = "DSInputOutput_Concerned";


class EqualityConstraint;

/** \class EqualityConstraintXML
 *  \brief object to manage XML data of an EqualityConstraint
 *  \author Jean-Michel Barbier
 *  \version 0.1
 *  \date 17/01/2005
 *
 * $Date: 2005/03/21 16:48:04 $
 * $Revision: 1.10 $
 * $Author: jbarbier $
 * $Source: /CVS/Siconos/SICONOS/src/xml/EqualityConstraintXML.h,v $
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
  EqualityConstraintXML(xmlNode*, vector<int>);
  ~EqualityConstraintXML();


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
  inline string getType()
  {
    string res((char*)this->rootNode->name);
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
      this->GNode = SiconosDOMTreeTools::createMatrixNode(this->rootNode, EQUALITYCONSTRAINT_G, G);
    }
    else SiconosDOMTreeTools::setSiconosMatrixValue(this->GNode, G);
  }


  /** \fn vector<int> getDSIOConcerned()
  *   \brief Return the DSIOs concerned by the EqualityConstraintXML
  *   \return the integer vector who contains the DSs concerned by the EqualityConstraintXML
  */
  inline vector<int> getDSIOConcerned()
  {
    return this->definedDSIONumbers;
  }

  /** \fn void setDSIOConcerned( vector<int> )
  *   \brief allows to set the dynamical systems which are interacting together with this EqualityConstraintXML
  *   \param vector<int> : the dynamical system numbers
  */
  void setDSIOConcerned(vector<int>);

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
  inline string getComputeInputPlugin()
  {
    if (this->isComputeInputPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->computeInputNode, PLUGIN_ATTRIBUTE);
    XMLException::selfThrow("EqualityConstraintXML - getComputeInputPlugin : computeInput is not calculated from a plugin ; Fint vector is given");
  }

  /** \fn inline string getComputeOutputPlugin()
  *   \brief Return the computeOutput Plugin name of the EqualityConstraintXML
  *   \return The computeOutput Plugin name of the EqualityConstraintXML
  *  \exception XMLException
  */
  inline string getComputeOutputPlugin()
  {
    if (this->isComputeOutputPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->computeOutputNode, PLUGIN_ATTRIBUTE);
    XMLException::selfThrow("EqualityConstraintXML - getComputeOutputPlugin : computeOutput is not calculated from a plugin ; Fint vector is given");
  }

  /** \fn void setComputeInputPlugin(string plugin)
  *   \brief sets the computeInput Plugin name of the EqualityConstraintXML
  *   \param string :  The computeInput Plugin name of the EqualityConstraintXML
  *  \exception XMLException
  */
  inline void setComputeInputPlugin(string plugin)
  {
    if (this->computeInputNode == NULL)
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
  inline void setComputeOutputPlugin(string plugin)
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
  vector<int> definedDSIONumbers;

private :
  /** \fn loadECConcernedDSIO(xmlNode * , vector<int>)
  *   \brief load the DSs concerned by this interaction
  *   \param xmlNode * : the DOM tree node of DSIO concerned by the EqualityConstraint
  // *   \param vector<int> : vector of DSXML numbers to verify DS concerned by the interaction (identified by number) exists
  */
  void loadECConcernedDSIO(xmlNode * DSIOConcernedNode/*, vector<int> definedDSNumbers*/);
};

#endif // EQUALITYCONSTRAINTXML_H

//$Log: EqualityConstraintXML.h,v $
//Revision 1.10  2005/03/21 16:48:04  jbarbier
//- EqualityConstraint : computeInput and computeOutput functions added (plugin funcitons)
//
//- link OneStepNSProblem - EqualityConstraint established
//
//- modification of OneStepNSProblem save according to change to normType[64] in SiconosNumerics.h
//
//Revision 1.9  2005/03/15 09:57:48  jbarbier
//- EqualityConstraint save OK
//
//Revision 1.8  2005/03/14 16:05:27  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.7  2005/03/10 12:55:21  jbarbier
//- implmentation of the EqualityConstraint and DSInputOutput classes in progress
//    attributes H (DSIO) et G (EC) added in XML and managed in XML objects
//
//Revision 1.6  2005/03/09 15:30:36  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.5  2005/03/08 12:41:38  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.4  2005/01/26 13:50:40  jbarbier
//
//- loading of an XML input file now loads EqualityConstraints and DSInputOutputs
//
//Revision 1.3  2005/01/25 14:51:47  jbarbier
//- attributes id, type and XML object added to EqualityConstraint
//
//Revision 1.2  2005/01/18 10:35:17  jbarbier
//- attribute "r" no longer used for Moreau integrator
//
//- modificatoin in the tests for Moreau integrator
//
//- file XMLTagsName.h for further use to regroup all xml tags name...
//
//Revision 1.1  2005/01/17 10:56:27  jbarbier
//- classes EqualityConstraint and DSInputOutput added with inherited classes
//
//- classes EqualityConstraintXML and DSInputOutputXML added with inherited classes
//