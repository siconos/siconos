/** \class InteractionXML
*   \brief This class manages Interaction data part
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date (Creation) 04/12/2004
*

*
* InteractionXML allows to manage data of a Interaction DOM tree.
*/


#ifndef __INTERACTIONXML__
#define __INTERACTIONXML__


#include <vector>
#include <string>
#include <libxml/tree.h>
#include <stdlib.h>

//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"
#include "RelationXML.h"
#include "NonSmoothLawXML.h"

#include "SiconosDOMTreeTools.h"

#include "Interaction.h"

using namespace std;


class Interaction;
class NonSmoothLawXML;
class RelationXML;


const string INTERACTION_STATUS = "Status";
const string INTERACTION_Y = "y";
const string INTERACTION_LAMBDA = "lambda";
const string INTERACTION_NINTER = "nInter";
const string INTERACTION_DS_CONCERNED = "DS_Concerned";
const string INTERACTION_INTERACTWITHDS_NUMBER = "interactsWithDS_Number";
#include "XMLTagsName.h"

class InteractionXML
{
public:
  InteractionXML();

  /** \fn InteractionXML(xmlNode * rootNodeInteraction, int number)
  *   \brief Build a InteractionXML object from a DOM tree describing a Interaction
  *   \param xmlNode * interactionNode : the Interaction DOM tree
  *   \param vector<int> definedDSNumbers : vector of DSXML numbers to verify DS concerned by the interaction (identified by number) exists
  */
  InteractionXML(xmlNode * interactionNode, vector<int> definedDSNumbers);

  ~InteractionXML();


  /** \fn xmlNode* getInteractionXMLNode()
  *   \brief Return the root node of the InteractionXML
  *   \return xmlNode* : the root node
  */
  inline xmlNode* getInteractionXMLNode()
  {
    return this->rootInteractionXMLNode;
  }

  /** \fn int getNumber()
  *   \brief Return the number of the InteractionXML
  *   \return The number of the InteractionXML
  */
  inline int getNumber()
  {
    return SiconosDOMTreeTools::getIntegerAttributeValue(this->rootInteractionXMLNode, NUMBER_ATTRIBUTE);
  }

  /** \fn void setNumber(int i)
  *   \brief allows to save the number of the InteractionXML
  *   \return The number to save
  */
  inline void setNumber(int i)
  {
    SiconosDOMTreeTools::setIntegerAttributeValue(this->rootInteractionXMLNode, NUMBER_ATTRIBUTE, i);
  }

  /** \fn bool hasId()
   *  \brief return true if idNode is defined
   *  \return true if idNode is defined
   */
  inline bool hasId()
  {
    return (this->idNode != NULL);
  }

  /** \fn string getId()
  *   \brief Return the id of the InteractionXML
  *   \return The string id of the InteractionXML
  */
  inline string getId()
  {
    return SiconosDOMTreeTools::getStringContentValue(this->idNode);
  }

  /** \fn void setId(string s)
  *   \brief allows to save the id of the InteractionXML
  *   \return The string id to save
  */
  inline void setId(string s)
  {
    if (hasId() == false)
    {
      this->idNode = SiconosDOMTreeTools::createStringNode(this->rootInteractionXMLNode, ID_ATTRIBUTE, s);
    }
    else SiconosDOMTreeTools::setStringContentValue(this->idNode, s);
  }

  /** \fn int getNumber()
  *   \brief Return the number of the InteractionXML
  *   \return The number of the InteractionXML
  */
  inline int getNInter()
  {
    return SiconosDOMTreeTools::getIntegerContentValue(this->nInterNode);
  }

  /** \fn void setNumber(int i)
  *   \brief allows to save the number of the InteractionXML
  *   \return The number to save
  */
  inline void setNInter(int nInter)
  {
    if (this->nInterNode == false)
    {
      this->nInterNode = SiconosDOMTreeTools::createIntegerNode(this->rootInteractionXMLNode, INTERACTION_NINTER, nInter);
    }
    else SiconosDOMTreeTools::setIntegerContentValue(this->nInterNode, nInter);
  }

  /** \fn vector<int> getStatus()
  *   \brief Return the status of the InteractionXML
  *   \return vector<int> status : the status of the InteractionXML
  */
  inline vector<int> getStatus()
  {
    return SiconosDOMTreeTools::getVectorIntContentValue(this->statusNode);
  }

  /** \fn void setStatus(vector<int> i)
  *   \brief allows to save the status of the InteractionXML
  *   \return vector<int> status : the status to save
  */
  inline void setStatus(vector<int> i)
  {
    if (this->statusNode == false)
    {
      this->statusNode = SiconosDOMTreeTools::createVectorIntNode(this->rootInteractionXMLNode, INTERACTION_STATUS, i);
    }
    else SiconosDOMTreeTools::setVectorIntContentValue(this->statusNode, i);
  }

  /** \fn bool hasY()
   *  \brief return true if yNode is defined
   *  \return true if yNode is defined
   */
  inline bool hasY()
  {
    return (this->yNode != NULL);
  }

  /** \fn SimpleVector getY()
  *   \brief Return y vector of the InteractionXML
  *   \return SimpleVector : the y of the InteractionXML
  */
  inline /*SiconosVector*/SimpleVector getY()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(this->yNode);
  }

  /** \fn void setY(SiconosVector *v)
  *   \brief allows to save the y of the InteractionXML
  *   \param SiconosVector* : the y to save
  */
  inline void setY(SiconosVector *v)
  {
    if (hasY() == false)
    {
      this->yNode = SiconosDOMTreeTools::createVectorNode(this->rootInteractionXMLNode, INTERACTION_Y, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->yNode, v);
  }

  /** \fn bool hasLambda()
   *  \brief return true if lambdaNode is defined
   *  \return true if lambdaNode is defined
   */
  inline bool hasLambda()
  {
    return (this->lambdaNode != NULL);
  }

  /** \fn SimpleVector getLambda()
  *   \brief Return the lambda vector of the InteractionXML
  *   \return SimpleVector : the lambda of the InteractionXML
  */
  inline /*SiconosVector*/ SimpleVector getLambda()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(this->lambdaNode);
  }

  /** \fn void setLambda(SiconosVector *v)
  *   \brief allows to save the lambda of the InteractionXML
  *   \return SiconosVector* : the lambda to save
  */
  inline void setLambda(SiconosVector *v)
  {
    if (this->hasLambda() == false)
    {
      this->lambdaNode = SiconosDOMTreeTools::createVectorNode(this->rootInteractionXMLNode, INTERACTION_LAMBDA, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->lambdaNode, v);

    /* affichage pour les tests d'integration */
  }


  /** \fn vector< vector<int> > getDSConcerned()
  *   \brief Return the DSs concerned by the InteractionXML
  *   \return the 2D integer vector who contains the DSs concerned coulpes by the InteractionXML
  */
  inline vector< vector<int> > getDSConcerned()
  {
    return this->DSCouples;
  }

  /** \fn void setDSConcerned( vector<DynamicalSystem*> )
  *   \brief allows to set the dynamical systems which are interacting together with this interaction
  *   \param vector<DynamicalSystem*> : the dynamical systems in interaction
  */
  void setDSConcerned(vector<DynamicalSystem*>);


  /** \fn RelationXML* getRelationXML()
  *   \brief Return the relationXML of the InteractionXML
  *   \return The relationXML of the InteractionXML
  */
  inline RelationXML* getRelationXML()
  {
    return this->relationXML;
  }


  /** \fn NonSmoothLawXML* getNonSmoothLawXML()
  *   \brief Return the NonSmoothLawXML of the InteractionXML
  *   \return The NonSmoothLawXML of the InteractionXML
  */
  inline NonSmoothLawXML* getNonSmoothLawXML()
  {
    return this->nSLawXML;
  }

  /** \fn void updateInteractionXML( xmlNode* node, Interaction* inter );
  *   \brief makes the operations to add an Interaction to the NSDS
  *   \param xmlNode* : the root node of the InteractionXML
  *   \param Interaction* : the Interaction of this InteractionXML
  */
  void updateInteractionXML(xmlNode* node, Interaction* inter);

  /** \fn void loadInteraction( Interaction* )
  *   \brief loads the data of the Interaction into the InteractionXML (in the DOM tree)
  *   \param NSDS* : the Interaction of this InteractionXML
  */
  void loadInteraction(Interaction*);

  /** \fn bool hasAll()
   *  \brief All is an attribute of the DS_Concerned tag
   *  \return bool : true if attribute all is defined
   */
  inline bool hasAll()
  {
    if (SiconosDOMTreeTools::hasAttributeValue(this->dsConcernedNode, ALL_ATTRIBUTE))
      return SiconosDOMTreeTools::getBooleanAttributeValue(this->dsConcernedNode, ALL_ATTRIBUTE);
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
        xmlNewProp(this->dsConcernedNode, (xmlChar*)ALL_ATTRIBUTE.c_str(), (xmlChar*)"true");
    }
    else
    {
      if (all == false)
        xmlRemoveProp(xmlHasProp(this->dsConcernedNode, (xmlChar*)INTERACTION_DS_CONCERNED.c_str()));
    }
  }


private:

  //Nodes
  xmlNode * rootInteractionXMLNode;

  xmlNode * idNode;
  xmlNode * nInterNode;
  xmlNode * statusNode;
  xmlNode * yNode;
  xmlNode * lambdaNode;
  xmlNode * isActiveNode;
  xmlNode * dsConcernedNode;

  //Couples of DSs (DS numbers)
  vector< vector<int> > DSCouples;

  //Relation
  RelationXML *relationXML;

  //NSLAW
  NonSmoothLawXML *nSLawXML;

  //Methods

  /** \fn loadInteractionProperties(xmlNode * , vector<int>)
  *   \brief load the different properties of a Interaction
  *   \param xmlNode * nodeInteraction : the DOM tree node of the concern Interaction
  *   \param vector<int> definedDSNumbers : vector of DSXML numbers to verify DS concerned by the interaction (identified by number) exists
  *   \exception XMLException : if a property of the Interaction lacks in the DOM tree
  */
  void loadInteractionProperties(xmlNode * interactionNode, vector<int> definedDSNumbers);

  /** \fn loadInteractionConcernedDS(xmlNode * , vector<int>)
  *   \brief load the DSs concerned by this interaction
  *   \param xmlNode * DSConcernedNode : the DOM tree node of DS concerned by the interaction
  *   \param vector<int> definedDSNumbers : vector of DSXML numbers to verify DS concerned by the interaction (identified by number) exists
  *   \exception XMLException : if a DS number not exists
  */
  void loadInteractionConcernedDS(xmlNode * DSConcernedNode, vector<int> definedDSNumbers);


};

#endif
//$Log: InteractionXML.h,v $
//Revision 1.39  2005/03/10 12:55:21  jbarbier
//- implmentation of the EqualityConstraint and DSInputOutput classes in progress
//    attributes H (DSIO) et G (EC) added in XML and managed in XML objects
//
//Revision 1.38  2005/03/08 12:41:38  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.37  2005/03/04 15:35:27  jbarbier
//- README files added for some samples
//
//- beginning of the refactoring of XML module constants
//
//Revision 1.36  2005/02/24 15:50:21  jbarbier
//- LCP prepared to changes needed for several interactions
//
//- new function for the SiconosMatrices to copy a block matrix into another matrix
//
//- tag names of BoundaryConditionXML, DSInputOutputXML, DSXML, InteractionXML, LagrangianLinearRXML, LagrangianDSXML put in XMLTagNames.h
//
//Revision 1.35  2005/01/26 13:50:40  jbarbier
//
//- loading of an XML input file now loads EqualityConstraints and DSInputOutputs
//
//Revision 1.34  2005/01/18 10:35:17  jbarbier
//- attribute "r" no longer used for Moreau integrator
//
//- modificatoin in the tests for Moreau integrator
//
//- file XMLTagsName.h for further use to regroup all xml tags name...
//
//Revision 1.33  2005/01/14 11:51:25  jbarbier
//- attribute "all" of the OneStepNSProblem terminated
//in OneStepIntegrator and Interaction the attribute is available
//
//Revision 1.32  2005/01/13 14:14:40  jbarbier
//- correction in the XML output about size attribute in tags DS_Concerned and Interactoin _Concerned
//
//- modifications in creation of XML objects when saving data with partial XML input file
//
//Revision 1.31  2004/12/08 12:49:39  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.30  2004/09/27 08:24:26  charlety
//
//_ Modifications in doxygen comments.
//
//Revision 1.29  2004/09/21 11:49:10  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.28  2004/09/14 13:49:56  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.27  2004/09/10 11:26:26  charlety
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
//Revision 1.26  2004/09/10 08:04:49  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.25  2004/07/29 14:04:00  jbarbier
//- new test on SiconosMemoryXML
//
//- last functions hasAttribute() in the XML part added
//
//Revision 1.24  2004/07/08 12:51:07  jbarbier
//-attribut status of the Interaction changed to vector<int> in the InteractionXML
//-modifications of the xml schema for the vector<int> status of an interaction
//-modifications of the xml schema on the DynamicalSystems to determine for each
//attribut if it is optional or required.
//-attribut BVP of the NSDS moved to optional
//
//Revision 1.23  2004/07/07 08:14:54  jbarbier
//-modifications on the test after renamming
//
//-modification of the XML schema, attributs row, col and size of Matrices and
//Vector changed from 'positiveInteger' to 'nonNegativeInteger'
//
//-function setSiconosVector/Matrix which take a SiconosVector/Matrix* in parameter to avoid
//unnecessary vector and matrix copies
//
