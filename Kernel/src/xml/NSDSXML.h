//$Id: NSDSXML.h,v 1.42 2005/03/14 16:05:27 jbarbier Exp $

/** \class NSDSXML
*   \brief This class manages NSDS data part
*   \author J. Blanc-Tranchant
*   \version 1.0
*   \date 04/04/2004
*
*
* $Date: 2005/03/14 16:05:27 $
* $Revision: 1.42 $
* $Author: jbarbier $
* $Source: /CVS/Siconos/SICONOS/src/xml/NSDSXML.h,v $
*
* NSDSXML allows to get DSXMLs and InteractionXMLs from a DOM tree.
*/


#ifndef __NSDSXML__
#define __NSDSXML__

#include <string>
#include <vector>
#include <map>
#include <libxml/tree.h>


#include "NonSmoothDynamicalSystem.h"

#include "DSXML.h"
#include "InteractionXML.h"
#include "EqualityConstraintXML.h"
//#include "DSInputOutputXML.h"
#include "SiconosDOMTreeTools.h"

#include "XMLTagsName.h"


using namespace std;


class NonSmoothDynamicalSystem;
class SiconosModelXML;
class DSXML;
class InteractionXML;
class EqualityConstraintXML;


const string NSDS_BVP = "bvp";


class NSDSXML
{
public:
  NSDSXML();

  /** \fn NSDSXML(xmlNode * rootNSDSNode)
  *   \brief Build an NSDSXML object from a DOM tree describing an NSDS
  *   \param rootNSDSNode : the NSDS DOM tree
  */
  NSDSXML(xmlNode * rootNSDSNode);

  ~NSDSXML();

  ///* return a DSs Map */
  //inline map<int, DSXML> getDSMap() const;

  /** \fn xmlNode* getNSDSXMLNode()
  *   \brief Return the root node of the NSDSXML
  *   \return xmlNode* : the root node
  */
  inline xmlNode* getNSDSXMLNode()
  {
    return this->NSDSNode;
  }

  /** \fn DSXML* getDSXML(int number)
  *   \brief Return the DSXML with id number
  *   \param int number : the number of the DSXML to return
  *   \return the DSXML of number number NULL if doesn't exist
  */
  DSXML* getDSXML(int number);

  /** \fn InteractionXML* getInteractionXML(int number)
  *   \brief Return the InteracionXML with id number
  *   \param number : int number : the number of the InteractionXML to return
  *   \exception XMLException
  *   \return the InteractionXML of number number, NULL if doesn't exist
  */
  InteractionXML* getInteractionXML(int number);

  /** \fn EqualityConstraintXML* getEqualityConstraintXML(int number)
  *   \brief Return the EqualityConstraintXML with id number
  *   \param number : int number : the number of the EqualityConstraintXML to return
  *   \exception XMLException
  *   \return the EqualityConstraintXML of number number, NULL if doesn't exist
  */
  EqualityConstraintXML* getEqualityConstraintXML(int number);


  /** \fn inline vector<int> getDSNumbers();
  *   \brief Allows to know the defined DS
  *   \exception XMLException
  *   \return vector DS numbers
  */
  inline vector<int> getDSNumbers()
  {
    return this->definedDSNumbers;
  }

  /** \fn inline vector<int> getInteractionNumbers();
  *   \brief Allows to know the defined interactions
  *   \return vector Interactions integer numbers
  */
  inline vector<int> getInteractionNumbers()
  {
    return this->definedInteractionNumbers;
  }

  /** \fn inline vector<int> getEqualityConstraintNumbers();
  *   \brief Allows to know the defined EqualityConstraints
  *   \return vector EqualityConstraints integer numbers
  */
  inline vector<int> getEqualityConstraintNumbers()
  {
    return this->definedEqualityConstraintNumbers;
  }


  /** \fn bool isBVP()
  *   \brief Allows to know if the NSDS is BVP or not
  *   \return True if the NSDS is BVP false otherwise
  */
  inline bool isBVP()
  {
    if (SiconosDOMTreeTools::hasAttributeValue(this->NSDSNode, NSDS_BVP))
      return SiconosDOMTreeTools::getBooleanAttributeValue(this->NSDSNode, NSDS_BVP);
    else return false;
  }

  /** \fn void setBVP(bool)
  *   \brief Allows to define if the NSDS is BVP
  *   \param True if the NSDS is BVP false otherwise
  */
  inline void setBVP(bool b)
  {
    if (!(SiconosDOMTreeTools::hasAttributeValue(this->NSDSNode, NSDS_BVP)))
    {
      //if( b ) SiconosDOMTreeTools::setBooleanAttributeValue(this->NSDSNode, NSDS_BVP, true);
      if (b == true) xmlNewProp(this->NSDSNode, (xmlChar*)NSDS_BVP.c_str(), (xmlChar*)"true");
    }
    else
    {
      if (b == false) xmlRemoveProp(xmlHasProp(this->NSDSNode, (xmlChar*)NSDS_BVP.c_str()));
    }
  }

  /** \fn void updateNSDSXML(xmlNode*, NonSmoothDynamicalSystem*)
  *   \brief makes the operations to add a NSDS to the SiconosModelXML
  *   \param xmlNode* : the root node for the NSDSXML
  *   \param NonSmoothDynamicalSystem* : the NonSmoothDynamicalSystem of this NSDSXML
  */
  void updateNSDSXML(xmlNode*, NonSmoothDynamicalSystem*);

  /** \fn void loadNonSmoothDynamicalSystem( NonSmoothDynamicalSystem* )
  *   \brief loads the data of the NSDS into the NSDSXML (in the DOM tree)
  *   \param NonSmoothDynamicalSystem* : the NSDS of this NSDSXML
  */
  void loadNonSmoothDynamicalSystem(NonSmoothDynamicalSystem*);


private:
  xmlNode *NSDSNode;

  /* Map of DSs */
  map<int, DSXML*> DSXMLMap;

  /* Map of interactions */
  map<int, InteractionXML*> interactionXMLMap;

  /* Map of EqualityConstraints */
  map<int, EqualityConstraintXML*> equalityConstraintXMLMap;


  /* vector of DS numbers*/
  vector<int> definedDSNumbers;

  /* vector of Interaction numbers*/
  vector<int> definedInteractionNumbers;

  /* vector of EqualityConstraint numbers*/
  vector<int> definedEqualityConstraintNumbers;


  /** \fn loadNonSmoothDynamicalSystem()
  *   \brief Load the NonSmoothDynamicalSystem : Interactions and DSs
  *   \exception XMLException
  */
  void loadNonSmoothDynamicalSystem();

  /** \fn loadDSXML(xmlNode * rootDSNode)
  *   \brief Builds DSXML objects from a DOM tree describing DSs
  *   \param xmlNode* : the DSs DOM tree
  *   \exception XMLException : if a property of the NSDS lacks in the DOM tree
  */
  void loadDSXML(xmlNode * rootDSNode);

  /** \fn loadInteractionXML(xmlNode * rootInteractionNode)
  *   \brief Builds InteractionXML objects from a DOM tree describing Interactions
  *   \param xmlNode* : the Interactions DOM tree
  *   \exception XMLException : if a number relating to an Interaction declares in the NSDS is already used
  */
  void loadInteractionXML(xmlNode * rootInteractionNode);

  /** \fn void loadEqualityConstraintXML(xmlNode * rootECNode)
  *   \brief Builds EqualityConstraintXML objects from a DOM tree describing EqualityConstraints
  *   \param xmlNode* : the EqualityConstraints DOM tree
  *   \exception XMLException : if a number relating to an EqualityConstraint declares in the NSDS is already used
  */
  void loadEqualityConstraintXML(xmlNode * rootECNode);

};



#endif
//$Log: NSDSXML.h,v $
//Revision 1.42  2005/03/14 16:05:27  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.41  2005/03/09 15:30:38  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.40  2005/03/08 12:41:39  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.39  2005/03/04 15:35:27  jbarbier
//- README files added for some samples
//
//- beginning of the refactoring of XML module constants
//
//Revision 1.38  2005/01/26 13:50:40  jbarbier
//
//- loading of an XML input file now loads EqualityConstraints and DSInputOutputs
//
//Revision 1.37  2005/01/20 14:44:49  jbarbier
//- NSDS class renamed NonSmoothDynamicalSystem
//
//- code reduce, some comments remove
//
//Revision 1.36  2004/09/21 11:49:10  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.35  2004/09/16 11:35:25  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.34  2004/09/10 08:04:51  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.33  2004/08/23 14:30:03  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NSDS. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.32  2004/08/20 15:26:46  jbarbier
//- creation of a Model and save in the XML is ok
//- creation of a NSDS and save in the XML is ok
//- creation of a NonLinearSystemDS and save in the XML is OK
//
//Revision 1.31  2004/08/20 07:34:23  jbarbier
//- creation of Model, NSDS in comand program succeed in creating SiconosModelXML,
//NSDSXML
//
//Revision 1.30  2004/08/05 12:44:44  jbarbier
//- loading XML file with no OneStepNSProblem succesfull
//
//- NonLinearSystemDS is now available
//
//Revision 1.29  2004/07/29 14:25:44  jbarbier
//- $Log: NSDSXML.h,v $
//- Revision 1.42  2005/03/14 16:05:27  jbarbier
//- - manual creation of DSInputOutput saving OK
//-
//- - in progress for EqualityConstraint
//-
//- Revision 1.41  2005/03/09 15:30:38  jbarbier
//- - add of LagrangianEC class
//-
//- - in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//-
//- Revision 1.40  2005/03/08 12:41:39  jbarbier
//- - constant variables files modified :
//- Some constants added in SiconosConst
//-
//- all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//-
//- Revision 1.39  2005/03/04 15:35:27  jbarbier
//- - README files added for some samples
//-
//- - beginning of the refactoring of XML module constants
//-
//- Revision 1.38  2005/01/26 13:50:40  jbarbier
//-
//- - loading of an XML input file now loads EqualityConstraints and DSInputOutputs
//-
//- Revision 1.37  2005/01/20 14:44:49  jbarbier
//- - NSDS class renamed NonSmoothDynamicalSystem
//-
//- - code reduce, some comments remove
//-
//- Revision 1.36  2004/09/21 11:49:10  jbarbier
//- - correction in the XML save for a manual construction of the platform :
//-     DS_Concerned of the Interaction
//-     DS_Concerned of the Integrator
//-
//- - test updated for these changes
//-
//- Revision 1.35  2004/09/16 11:35:25  jbarbier
//- - save of the TimeDiscretisation in a XML file in manual creation of the
//- platform which was forgotten is now available.
//-
//- - the save of the platform's data can be done when the platform is created with
//- an XML input file and completed with dynmical systems, interactions, one-step
//- non smooth problem and one-step integrator.
//-
//- Revision 1.34  2004/09/10 08:04:51  jbarbier
//- - XML save available for BoundaryCondition and Interaction
//-
//- Revision 1.33  2004/08/23 14:30:03  jbarbier
//- - All the dynamical systems can be created in a comand program and added to a
//- NSDS. The save is OK, but the creation of the boundary conditions is not yet
//- finished.
//-
//- Revision 1.32  2004/08/20 15:26:46  jbarbier
//- - creation of a Model and save in the XML is ok
//- - creation of a NSDS and save in the XML is ok
//- - creation of a NonLinearSystemDS and save in the XML is OK
//-
//- Revision 1.31  2004/08/20 07:34:23  jbarbier
//- - creation of Model, NSDS in comand program succeed in creating SiconosModelXML,
//- NSDSXML
//-
//- Revision 1.30  2004/08/05 12:44:44  jbarbier
//- - loading XML file with no OneStepNSProblem succesfull
//-
//- - NonLinearSystemDS is now available
//- and $Id: NSDSXML.h,v 1.42 2005/03/14 16:05:27 jbarbier Exp $ added
//
