//$Id: DSInputOutputXML.h,v 1.10 2005/03/15 09:57:48 jbarbier Exp $

/** \class DSInputOutputXML
*   \brief This class manages Relation data part
*   \author Jean-Michel Barbier
*   \version 1.0
*   \date 17/01/2005
*
*
* $Date: 2005/03/15 09:57:48 $
* $Revision: 1.10 $
* $Author: jbarbier $
* $Source: /CVS/Siconos/SICONOS/src/xml/DSInputOutputXML.h,v $
*
* DSInputOutputXML allows to manage data of a DSInputOutput DOM tree.
*/

#ifndef __DSInputOutputXML__
#define __DSInputOutputXML__


#include <string>
#include <libxml/tree.h>
#include "SiconosDOMTreeTools.h"
#include "DSInputOutput.h"

#include "XMLTagsName.h"


using namespace std;

const string DSINPUTOUTPUT_H = "H";
//const string DSINPUTOUTPUT_DS_CONCERNED = "DS_Concerned";


class DSInputOutput;


class DSInputOutputXML
{
public:
  DSInputOutputXML();

  /** \fn DSInputOutputXML(xmlNode * , vector<int> )
  *   \brief Build a DSInputOutputXML object from a DOM tree describing a DSInputOutput
  *   \param xmlNode* : the DSInputOutput DOM tree
  //    *   \param vector<int>  : vector of DSXML numbers to verify DS concerned by the DSInputOutput (identified by number) exists
  */
  DSInputOutputXML(xmlNode*/*, vector<int>*/);
  ~DSInputOutputXML();

  //    /** \fn inline string getComputeInputPlugin()
  //    *   \brief Return the computeInput Plugin name of the DSInputOutputXML
  //    *   \return The computeInput Plugin name of the DSInputOutputXML
  //    *  \exception XMLException
  //    */
  //    inline string getComputeInputPlugin()
  //    {
  //      if (this->isComputeInputPlugin())
  //        return  SiconosDOMTreeTools::getStringAttributeValue(this->computeInputNode, PLUGIN_ATTRIBUTE);
  //      XMLException::selfThrow("DSInputOutputXML - getComputeInputPlugin : computeInput is not calculated from a plugin ; Fint vector is given");
  //    }
  //
  //    /** \fn inline string getComputeOutputPlugin()
  //    *   \brief Return the computeOutput Plugin name of the DSInputOutputXML
  //    *   \return The computeOutput Plugin name of the DSInputOutputXML
  //    *  \exception XMLException
  //    */
  //    inline string getComputeOutputPlugin()
  //    {
  //      if (this->isComputeOutputPlugin())
  //        return  SiconosDOMTreeTools::getStringAttributeValue(this->computeOutputNode, PLUGIN_ATTRIBUTE);
  //      XMLException::selfThrow("DSInputOutputXML - getComputeOutputPlugin : computeOutput is not calculated from a plugin ; Fint vector is given");
  //    }
  //
  //    /** \fn void setComputeInputPlugin(string plugin)
  //    *   \brief sets the computeInput Plugin name of the DSInputOutputXML
  //    *   \param string :  The computeInput Plugin name of the DSInputOutputXML
  //    *  \exception XMLException
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
  //    /** \fn void setComputeOutputPlugin(string plugin)
  //    *   \brief sets the computeOutput Plugin name of the DSInputOutputXML
  //    *   \param string :  The computeOutput Plugin name of the DSInputOutputXML
  //    *  \exception XMLException
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
  //
  //    /** \fn bool isComputeInputPlugin()
  //    *   \brief Return true if computeInput is calculated from a plugin
  //    *   \return True if computeInput is calculated from plugin
  //    */
  //    inline bool isComputeInputPlugin()
  //    {
  //      return xmlHasProp((xmlNodePtr)computeInputNode,(xmlChar *) PLUGIN_ATTRIBUTE.c_str());
  //    }
  //
  //    /** \fn bool isComputeOutputPlugin()
  //    *   \brief Return true if computeOutput is calculated from a plugin
  //    *   \return True if computOutput is calculated from plugin
  //    */
  //    inline bool isComputeOutputPlugin()
  //    {
  //      return xmlHasProp((xmlNodePtr)computeOutputNode,(xmlChar *) PLUGIN_ATTRIBUTE.c_str());
  //    }
  //
  //
  //    /** \fn bool hasComputeInput()
  //     *  \brief return true if computeInputNode is defined
  //     *  \return true if computeInputNode is defined
  //     */
  //    inline bool hasComputeInput()
  //    {
  //      return (this->computeInputNode!=NULL);
  //    }
  //
  //    /** \fn bool hasComputeOutput()
  //     *  \brief return true if computeOutputNode is defined
  //     *  \return true if computeOutputNode is defined
  //     */
  //    inline bool hasComputeOutput()
  //    {
  //      return (this->computeOutputNode!=NULL);
  //    }


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
  inline string getType()
  {
    string type((char*)this->rootDSIOXMLNode->name);
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


  /** \fn SiconosMatrix getH()
  *   \brief Return H matrix of the DSInputOutputXML
  *   \return SiconosMatrix : the H matrix of the DSInputOutputXML
  */
  inline SiconosMatrix getH()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(this->HNode);
  }

  /** \fn void setH(SiconosMatrix *H)
  *   \brief allows to save the H matrix of the DSInputOutputXML
  *   \param SiconosMatrix* : the H to save
  */
  inline void setH(SiconosMatrix *H)
  {
    if (this->HNode == NULL)
    {
      this->HNode = SiconosDOMTreeTools::createMatrixNode(this->rootDSIOXMLNode, DSINPUTOUTPUT_H, H);
    }
    else SiconosDOMTreeTools::setSiconosMatrixValue(this->HNode, H);
  }

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
  inline vector<int> getDSConcerned()
  {
    return this->definedDSNumbers;
  }

  /** \fn void setDSConcerned( vector<int> )
  *   \brief allows to set the dynamical systems which are interacting together with this DSInputOutput
  *   \param vector<int> : the dynamical system numbers
  */
  void setDSConcerned(vector<int>);


protected:
  xmlNode * rootDSIOXMLNode;

  xmlNode * HNode;
  xmlNode * dsConcernedNode;
  //    xmlNode * computeInputNode;
  //    xmlNode * computeOutputNode;

  /* vector of DS numbers*/
  vector<int> definedDSNumbers;


private :

  /** \fn loadDSIOConcernedDS(xmlNode * , vector<int>)
  *   \brief load the DSs concerned by this interaction
  *   \param xmlNode * DSConcernedNode : the DOM tree node of DS concerned by the interaction
  //    *   \param vector<int> definedDSNumbers : vector of DSXML numbers to verify DS concerned by the interaction (identified by number) exists
  *   \exception XMLException : if a DS number not exists
  */
  void loadDSIOConcernedDS(xmlNode * DSConcernedNode/*, vector<int> definedDSNumbers*/);
};


#endif
//$Log: DSInputOutputXML.h,v $
//Revision 1.10  2005/03/15 09:57:48  jbarbier
//- EqualityConstraint save OK
//
//Revision 1.9  2005/03/14 16:05:27  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.8  2005/03/10 12:55:21  jbarbier
//- implmentation of the EqualityConstraint and DSInputOutput classes in progress
//    attributes H (DSIO) et G (EC) added in XML and managed in XML objects
//
//Revision 1.7  2005/03/09 15:30:35  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.6  2005/03/08 12:41:37  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.5  2005/03/07 13:17:21  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.4  2005/02/24 15:50:20  jbarbier
//- LCP prepared to changes needed for several interactions
//
//- new function for the SiconosMatrices to copy a block matrix into another matrix
//
//- tag names of BoundaryConditionXML, DSInputOutputXML, DSXML, InteractionXML, LagrangianLinearRXML, LagrangianNLDSXML put in XMLTagNames.h
//
//Revision 1.3  2005/01/26 13:50:40  jbarbier
//
//- loading of an XML input file now loads EqualityConstraints and DSInputOutputs
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
