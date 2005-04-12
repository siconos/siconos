
/** \class OneStepIntegratorXML
*   \brief This class manages OneStepIntegrator data part
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/17/2004
*
*
* OneStepIntegratorXML allows to manage data of a OneStepIntegrator DOM tree.
*/


#ifndef __OneStepIntegratorXML__
#define __OneStepIntegratorXML__


#include <vector>
#include <string>
#include <map>
#include <libxml/tree.h>

//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SiconosMatrix.h"
#include "SiconosDOMTreeTools.h"

#include "OneStepIntegrator.h"


//using namespace std;

class OneStepIntegrator;


//Tags
const string OSI_R = "r";
const string OSI_DS_CONCERNED = "DS_Concerned";


#include "XMLTagsName.h"

class OneStepIntegratorXML
{
public:

  OneStepIntegratorXML();

  /** \fn OneStepIntegratorXML(xmlNode * OneStepIntegratorNode, map<int, bool> definedDSNumbers)
  *   \brief Build a OneStepIntegratorXML object from a DOM tree describing a OneStepIntegrator
  *   \param OneStepIntegratorNode : the OneStepIntegrator DOM tree
  *   \param map<int, bool> definedDSNumbers : the DS numbers effectivly defined in the model
  *   \exception XMLException : if a property of the OneStepIntegrator lacks in the DOM tree
  */
  OneStepIntegratorXML(xmlNode * OneStepIntegratorNode, map<int, bool> definedDSNumbers);

  ~OneStepIntegratorXML();


  /** \fn int getR()
  *   \brief Return r of the OneStepIntegrator
  *   \return The r integer of the OneStepIntegrator
  */
  inline int getR()
  {
    return SiconosDOMTreeTools::getIntegerContentValue(this->rNode);
  }

  /** \fn void setR(int r)
  *   \brief allows to save r of the OneStepIntegrator
  *   \param The r integer to save
  */
  inline void setR(int r)
  {
    if (this->hasR() == false)
    {
      this->rNode = SiconosDOMTreeTools::createIntegerNode(this->rootIntegratorXMLNode, OSI_R, r);
    }
    else SiconosDOMTreeTools::setIntegerContentValue(this->rNode, r);
  }

  /** \fn bool hasR()
  *   \brief determines if the r value is defined for this integrator
  *   \return bool : true if r is defined
  */
  inline bool hasR()
  {
    return (this->rNode != NULL);
  }

  /** \fn vector<int> getDSConcerned()
  *   \brief Return the DS numbers of the OneStepIntegrator
  *   \return The DS numbers vector of the OneStepIntegrator
  */
  inline vector<int> getDSConcerned()
  {
    return this->DSNumbersVector;
  }

  /** \fn vector<int> getDSConcerned()
  *   \brief Return the DS numbers of the OneStepIntegrator
  *   \return The DS numbers vector of the OneStepIntegrator
  */
  void setDSConcerned(vector<int>* ds);


  /** \fn string getType()
  *   \brief Return the type of the OneStepIntegratorXML
  *   \return The string type of the OneStepIntegratorXML
  */
  inline string getType()
  {
    //return SiconosDOMTreeTools::getStringAttributeValue(this->rootIntegratorXMLNode, OSI_TYPE);
    string type((char*)this->rootIntegratorXMLNode->name);
    return type;
  }

  /** \fn xmlNode* getNode()
  *   \brief allow to get the root node of the OneStepIntegratorXML
  *   \return xmlNode* : the root node of the OneStepIntegratorXML
  */
  inline xmlNode* getNode() const
  {
    return (xmlNode*)this->rootIntegratorXMLNode;
  }

  /** \fn void updateOneStepIntegratorXML( xmlNode* , OneStepIntegrator*  )
  *   \brief makes the operations to create a OneStepIntegratorXML to the StrategyXML
  *   \param xmlNode* : the root node of the OneStepIntegratorXML
  *   \param Strategy* : the OneStepIntegrator of this OneStepIntegratorXML
  */
  void updateOneStepIntegratorXML(xmlNode* node, OneStepIntegrator* osi);

  /** \fn bool hasAll()
   *  \brief All is an attribute of the DS_Concerned tag
   *  \return bool : true if attribute all is defined
   */
  inline bool hasAll()
  {
    if (SiconosDOMTreeTools::hasAttributeValue(this->DSConcernedNode, ALL_ATTRIBUTE))
      return SiconosDOMTreeTools::getBooleanAttributeValue(this->DSConcernedNode, ALL_ATTRIBUTE);
    else return false;
  }

  /** \fn void setAll(bool all)
  *   \brief Allows to modify the attribute "all" of the DS_concerned tag
  *   \param bool : the value to assign to the attribute
  */
  inline void setAll(bool all)
  {
    if (this->hasAll() == false)
    {
      if (all == true)
        xmlNewProp(this->DSConcernedNode, (xmlChar*)ALL_ATTRIBUTE.c_str(), (xmlChar*)"true");
    }
    else
    {
      if (all == false)
        xmlRemoveProp(xmlHasProp(this->DSConcernedNode, (xmlChar*)ALL_ATTRIBUTE.c_str()));
    }
  }

protected:
  //Nodes
  xmlNode * rNode;
  xmlNode * rootIntegratorXMLNode;
  xmlNode * DSConcernedNode;

private:



  //DSs (DS numbers)
  vector<int> DSNumbersVector;

  //Methods


  /** \fn loadOneStepIntegratonConcernedDS(xmlNode * DSConcernedNode, vector<int> definedDSNumbers)
  *   \brief load the DS numbers of the OneStepIntegrator
  *   \param xmlNode * DSConcernedNode : the DOM tree node of the concerned DS
  *   \param map<int, bool> definedDSNumbers : the DS numbers effectivly defined in the model
  *   \exception XMLException : if a DS of the OneStepIntegrator has not been defined in the model or is already used by another OneStepIntegrator
  */
  void loadOneStepIntegratorConcernedDS(xmlNode * DSConcernedNode, map<int, bool> definedDSNumbers);

};


#endif
//$Log: OneStepIntegratorXML.h,v $
//Revision 1.20  2005/03/08 14:23:45  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.19  2005/01/14 11:51:25  jbarbier
//- attribute "all" of the OneStepNSProblem terminated
//in OneStepIntegrator and Interaction the attribute is available
//
//Revision 1.18  2004/12/08 12:49:39  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.17  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.16  2004/09/21 11:49:10  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.15  2004/09/14 13:49:58  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.14  2004/09/10 11:26:28  charlety
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
//Revision 1.13  2004/07/29 14:25:44  jbarbier
