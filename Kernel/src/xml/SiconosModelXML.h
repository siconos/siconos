//$Id: SiconosModelXML.h,v 1.42 2005/03/08 14:23:46 jbarbier Exp $

/** \class SiconosModelXML
*   \brief This class allows to verify and download a Siconos XML data file, manage it data and save in a new Siconos XML data file
*   \author J. Blanc-Tranchant
*   \version 1.0
*   \date 04/04/2004
*
*
* $Date: 2005/03/08 14:23:46 $
* $Revision: 1.42 $
* $Author: jbarbier $
* $Source: /CVS/Siconos/SICONOS/src/xml/SiconosModelXML.h,v $
*
* SiconosModelXML allows the verification of a Siconos XML data file thanks to an XML Schema.
* All verifications can't be done with schema : others are done in the code, as the download of the XML file (download is stopped if error).
* This class calls NSDSXML and StrategyXML classes to manage Siconos NSDS and Strategy data of a simulation.
*/


#ifndef __MODELXML__
#define __MODELXML__

#include "Model.h"

#include "SiconosDOMTreeTools.h"
#include "NSDSXML.h"
#include "StrategyXML.h"
#include "RuntimeException.h"

//#include "KernelDefaultConfig.h"

#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xmlschemas.h>

#include <stdio.h>
#include <iostream>
#include <string>
#include <map>


using namespace std;

extern int MATRIX_MAX_SIZE;
extern int VECTOR_MAX_SIZE;
extern string FILE_STORAGE;
extern string XML_SCHEMA;

class Model;
class NSDSXML;
class StrategyXML;

const string ISO = "ISO-8859-1";
const string SM_T_CURRENT = "t";
const string SM_T = "T";
const string SM_T0 = "t0";
const string SM_TIME = "Time";

const string SM_TITLE = "Title";
const string SM_AUTHOR = "Author";
const string SM_DESCRIPTION = "Description";
const string SM_DATE = "Date";
const string SM_XMLSCHEMA = "SchemaXML";

#include "XMLTagsName.h"

//const string DEFAULT_XMLSCHEMA="/config/xmlschema/SiconosModelSchema-V1.2.xsd";



class SiconosModelXML
{
public:
  SiconosModelXML();

  /** \fn SiconosModelXML(char * siconosModelXMLFilePath)
  *   \brief Build an SiconosModelXML object from a Siconos XML data file
  *   \param siconosModelXMLFilePath : the path string of the Siconos XML data file
  *   \exception XMLException : exception may be the XML file does not exist ; XML schema has wrong syntax  or the XML file does not respect it ; etc.
  */
  SiconosModelXML(char * siconosModelXMLFilePath);

  /** \fn ~SiconosModelXML()
  *   \brief Destroy an SiconosModelXML object
  */
  ~SiconosModelXML();


  /** \fn xmlNode* getRootNode()
  *   \brief Gets the rootNode of the SiconosModelXML
  *   \return xmlNode* : the rootNode
  */
  inline xmlNode* getRootNode()
  {
    return this->rootNode;
  }

  /** \fn saveSiconosModelInXMLFile(char *siconosModelXMLFilePath)
  *   \brief Saves the Siconos Model creating a Siconos XML data file
  *   \param siconosModelXMLFilePath the path string of the Siconos XML data file to save
  */
  void saveSiconosModelInXMLFile(char *siconosModelXMLFilePath);

  /** \fn getNSDSXML()
  *   \brief This function allows to get the NSDSXML
  *   \return The NSDSXML of the SiconosModelXML
  */
  inline NSDSXML* getNSDSXML()
  {
    return this->nsdsXML;
  }

  /** \fn getStrategyXML()
  *   \brief This function allows to get the StrategyXML
  *   \return The StrategyXML of the SiconosModelXML
  */
  inline StrategyXML* getStrategyXML()
  {
    return this->strategyXML;
  }

  /** \fn bool hasT()
  *   \brief determines if the current time T is in the DOM tree
  *   \return bool :  true if T is in the DOM tree
  */
  inline bool hasT()
  {
    return (this->TNode != NULL);
  }

  /** \fn bool hasTCurrent()
  *   \brief determines if the current time t is in the DOM tree
  *   \return bool :  true if t is in the DOM tree
  */
  inline bool hasTCurrent()
  {
    return (this->tNode != NULL);
  }

  /** \fn void getTCurrent()
  *   \brief Gets the value of t
  *   \return t value
  */
  inline double getTCurrent()
  {
    return SiconosDOMTreeTools::getDoubleContentValue(this->tNode);
  }

  /** \fn double getT0()
  *   \brief Gets the value of t0
  *   \return t0 value
  */
  inline double getT0()
  {
    return  SiconosDOMTreeTools::getDoubleContentValue(this->t0Node);
  }

  /** \fn double getT()
  *   \brief Gets the value of T
  *   \return T value
  */
  inline double getT()
  {
    return  SiconosDOMTreeTools::getDoubleContentValue(this->TNode);
  }

  /** \fn void setT0(double)
  *   \brief Sets the value of t0
  *   \param The new t0 value
  */
  inline void setT0(double t0)
  {
    if (this->t0Node == NULL)
    {
      this->t0Node = SiconosDOMTreeTools::createDoubleNode(this->timeNode, SM_T0, t0);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(this->t0Node, t0);
  }

  /** \fn void setT(double)
  *   \brief Sets the value of T
  *   \param The new T value
  */
  inline void setT(double T)
  {
    if (this->hasT() == false)
    {
      this->TNode = SiconosDOMTreeTools::createDoubleNode(this->timeNode, SM_T, T);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(this->TNode, T);
  }

  /** \fn void setTCurrent(double)
  *   \brief Sets the value of t
  *   \param The new t value
  */
  inline void setTCurrent(double t)
  {
    if (this->hasTCurrent() == false)
    {
      this->tNode = SiconosDOMTreeTools::createDoubleNode(this->timeNode, SM_T_CURRENT, t);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(this->tNode, t);
  }

  /** \fn bool checkSiconosDOMTree()
  *   \brief checks the content of the DOM tree associated to the data of the SiconosModelXML
  * by using the XML schema of the platform
  *   \return bool : true if the XML schema is respected, else false
  *   \exception SiconosException
  */
  bool checkSiconosDOMTree();

  /** \fn bool checkSiconosDOMTreeCoherency()
  *   \brief checks the coherency of the content of the DOM tree associated to the data of the SiconosModelXML
  * by using the XML schema of the platform
  *   \return bool : true if the data are coherent, else false
  *   \exception SiconosException
  */
  bool checkSiconosDOMTreeCoherency();

  /** \fn bool hasStrategy()
  *   \brief determines if the Strategy is defined
  *   \return bool :  false if the strategyXML* is NULL
  */
  inline bool hasStrategy()
  {
    return (this->strategyXML != NULL);
  }


  /** \fn string getTitle()
  *   \brief Gets the title of the model
  *   \return string
  */
  inline string getTitle()
  {
    return  SiconosDOMTreeTools::getStringContentValue(this->titleNode);
  }

  /** \fn string getAuthor()
  *   \brief Gets the author of the model
  *   \return string
  */
  inline string getAuthor()
  {
    return  SiconosDOMTreeTools::getStringContentValue(this->authorNode);
  }

  /** \fn string getDescription()
  *   \brief Gets the Description of the model
  *   \return string
  */
  inline string getDescription()
  {
    return  SiconosDOMTreeTools::getStringContentValue(this->descriptionNode);
  }

  /** \fn string getDate()
  *   \brief Gets the date of the model
  *   \return string
  */
  inline string getDate()
  {
    return  SiconosDOMTreeTools::getStringContentValue(this->dateNode);
  }

  /** \fn string getXMLSchema()
  *   \brief Gets the XML Schema of the model
  *   \return string
  */
  inline string getXMLSchema()
  {
    return  SiconosDOMTreeTools::getStringContentValue(this->xmlSchemaNode);
  }

  /** \fn void setTitle(string s)
  *   \brief allows to save the title of the DSXML
  *   \param string : The string s of the DSXML
  */
  inline void setTitle(string s)
  {
    if (this->hasTitle() == false)
    {
      this->titleNode = SiconosDOMTreeTools::createStringNode(this->rootNode, SM_TITLE, s);
    }
    else SiconosDOMTreeTools::setStringContentValue(this->titleNode, s);
  }

  /** \fn void setAuthor(string s)
  *   \brief allows to save the author of the DSXML
  *   \param string : The string s of the DSXML
  */
  inline void setAuthor(string s)
  {
    if (this->hasAuthor() == false)
    {
      this->authorNode = SiconosDOMTreeTools::createStringNode(this->rootNode, SM_AUTHOR, s);
    }
    else SiconosDOMTreeTools::setStringContentValue(this->authorNode, s);
  }

  /** \fn void setDescription(string s)
  *   \brief allows to save the Description of the DSXML
  *   \param string : The string s of the DSXML
  */
  inline void setDescription(string s)
  {
    if (this->hasDescription() == false)
    {
      this->descriptionNode = SiconosDOMTreeTools::createStringNode(this->rootNode, SM_DESCRIPTION, s);
    }
    else SiconosDOMTreeTools::setStringContentValue(this->descriptionNode, s);
  }

  /** \fn void setDate(string s)
  *   \brief allows to save the date of the DSXML
  *   \param string : The string s of the DSXML
  */
  inline void setDate(string s)
  {
    if (this->hasDate() == false)
    {
      this->dateNode = SiconosDOMTreeTools::createStringNode(this->rootNode, SM_DATE, s);
    }
    else SiconosDOMTreeTools::setStringContentValue(this->dateNode, s);
  }

  /** \fn void setXMLSchema(string s)
  *   \brief allows to save the xml schema of the DSXML
  *   \param string : The string s of the DSXML
  */
  inline void setXMLSchema(string s)
  {
    if (this->hasXMLSchema() == false)
    {
      this->xmlSchemaNode = SiconosDOMTreeTools::createStringNode(this->rootNode, SM_XMLSCHEMA, s);
    }
    else SiconosDOMTreeTools::setStringContentValue(this->xmlSchemaNode, s);
  }

  /** \fn bool hasTitle()
  *   \brief determines if the title of the model is in the DOM tree
  *   \return bool :  true if the author of the model is in the DOM tree
  */
  inline bool hasTitle()
  {
    return (this->titleNode != NULL);
  }

  /** \fn bool hasAuthor()
  *   \brief determines if the author of the model is in the DOM tree
  *   \return bool :  true if the title of the model is in the DOM tree
  */
  inline bool hasAuthor()
  {
    //return SiconosDOMTreeTools::hasAttributeValue(this->rootNode, SM_AUTHOR);
    return (this->authorNode != NULL);
  }

  /** \fn bool hasDescription()
  *   \brief determines if the description of the model is in the DOM tree
  *   \return bool :  true if T is in the DOM tree
  */
  inline bool hasDescription()
  {
    //return SiconosDOMTreeTools::hasAttributeValue(this->rootNode, SM_DESCRIPTION);
    return (this->descriptionNode != NULL);
  }
  /** \fn bool hasDate()
  *   \brief determines if the date of the model is in the DOM tree
  *   \return bool :  true if date of the model is in the DOM tree
  */
  inline bool hasDate()
  {
    //return SiconosDOMTreeTools::hasAttributeValue(this->rootNode, SM_DATE);
    return (this->dateNode != NULL);
  }

  /** \fn bool hasXMLSchema()
  *   \brief determines if the xml schema of the model is in the DOM tree
  *   \return bool :  true if the xml schema of the model is in the DOM tree
  */
  inline bool hasXMLSchema()
  {
    return (this->xmlSchemaNode != NULL);
  }


  /** \fn loadModel(Model*)
  *   \brief Loads the model with all the data required to construct it
  *  \param Model* : the Model which contains all the data to build the SiconosModelXML
  */
  void loadModel(Model*);

  /** \fn int validateXmlFile(string xmlFile, string xmlSchema)
  *   \brief
  *  \param string :
  *  \param string :
  *  \return int :
  */
  int validateXmlFile(string xmlFile, string xmlSchema);



private:

  //string XML_SCHEMA_FILE = "/config/xmlschema/SiconosModelSchema-V1.0.xsd";
  string xmlSchemaFile;// = "/config/xmlschema/SiconosModelSchema-V1.0.xsd";

  NSDSXML *nsdsXML;
  StrategyXML *strategyXML;

  /*Time*/
  xmlNode * tNode;
  xmlNode * t0Node;
  xmlNode * TNode;

  xmlNode * titleNode;
  xmlNode * authorNode;
  xmlNode * descriptionNode;
  xmlNode * dateNode;
  xmlNode * xmlSchemaNode;

  xmlNode * rootNode;
  xmlNode * timeNode;

  xmlDoc * doc;


  /** \fn loadModel(xmlNode *)
  *   \brief Loads the model : Time (to and T) and the NSDSXML and StrategyXML of the model in order to manage NSDS and Strategy data of the simulation
  *   \param xmlNode * : the root node of the Model
  *   \exception XMLException : if a property of the SiconosModel (NSDS or Strategy or Time) lacks in the DOM tree
  */
  void loadModel(xmlNode *rootNode);

  /** \fn loadTime(xmlNode *)
  *   \brief Loads Time (to and T)
  *   \param xmlNode * : the root node of the Time element
  *   \exception XMLException : if a property (t, T0) of the Time lacks in the DOM tree
  */
  void loadTime(xmlNode *timeNode);

  /** \fn loadNSDSXML(xmlNode *)
  *   \brief Loads NSDSXML
  *   \param xmlNode * : the root node of the
  */
  void loadNSDSXML(xmlNode *NSDSNode);

  /** \fn loadStrategyXML(xmlNode *)
  *   \brief Loads StrategyXML
  *   \param xmlNode * : the root node of the
  */
  void loadStrategyXML(xmlNode *strategyNode);


  const string SICONOSPATH;
};

#endif

//$Log: SiconosModelXML.h,v $
//Revision 1.42  2005/03/08 14:23:46  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.41  2005/01/25 10:33:15  jbarbier
//- modifications for test purpose
//
//Revision 1.40  2005/01/20 14:44:49  jbarbier
//- NSDS class renamed NonSmoothDynamicalSystem
//
//- code reduce, some comments remove
//
//Revision 1.39  2005/01/20 09:05:35  jbarbier
//- configuration file available and usable
//
//- save of vectors and matrices into external files (version 0.1)
//
//Revision 1.38  2005/01/11 17:08:31  jbarbier
//- last modification about the BoundaryCondition
//<BoundaryCondition>
//  <[type]>
//  <[type]/>
//<BoundaryCondition/>
//
//- modification of the xml files for this modification
//
//- version 1.2 of the xml schema
//
//Revision 1.37  2005/01/10 14:07:42  jbarbier
//- file ChangeLog added for the autotools
//
//- xml schema corrected : BoundaryCondition modified without "choice"
//
//- new function of the Model to check at every moment the validity of an xml file according to the xml Schema
//
//Revision 1.36  2004/12/20 15:01:26  jbarbier
//- schema XML renamed V1.1
//
//- schema XML corrected about OneStepNSProblem:
//  tag OneStepNSProblem contains tags LCP, QP, ... and other tags to add
//  further
//
//Revision 1.35  2004/11/26 14:10:58  jbarbier
//- the xml schema used to check the xml file is the one given in the input xml
//file
//
//Revision 1.34  2004/09/27 13:27:14  jbarbier
//
//- Siconos schema renamed : SiconosModelSchema-V1.0.xsd
//
//- new required tags of the model : title, author, description, date, xmlSchema.
//They replace previous attributes author, description and date of the Model.
//
//Revision 1.33  2004/09/14 13:49:59  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.32  2004/08/20 07:34:23  jbarbier
//- creation of Model, NSDS in comand program succeed in creating SiconosModelXML,
//NSDSXML
//
//Revision 1.31  2004/08/12 14:28:38  jbarbier
//- createTimeDiscretisation in progress
//
//Revision 1.30  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.29  2004/08/10 14:51:49  jbarbier
//- functions initialize() of the Lsodar and Adams OneStepIntegrator completed to
//call the function initialize() of the base class
//
//Revision 1.28  2004/08/03 12:07:12  jbarbier
//- all test on th eModel are successfull
//
//- new tests on the Model with the opening of XML file
//
//- link TimeDiscretisation -> Strategy
//
//- attribute T of the Model is now optional
//
//Revision 1.27  2004/07/29 14:25:45  jbarbier
//- $Log: SiconosModelXML.h,v $
//- Revision 1.42  2005/03/08 14:23:46  jbarbier
//- - modification of constant variables :
//- in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//-
//- in simualtion tools, some constants have been moved to SiconosConst.h
//-
//- Revision 1.41  2005/01/25 10:33:15  jbarbier
//- - modifications for test purpose
//-
//- Revision 1.40  2005/01/20 14:44:49  jbarbier
//- - NSDS class renamed NonSmoothDynamicalSystem
//-
//- - code reduce, some comments remove
//-
//- Revision 1.39  2005/01/20 09:05:35  jbarbier
//- - configuration file available and usable
//-
//- - save of vectors and matrices into external files (version 0.1)
//-
//- Revision 1.38  2005/01/11 17:08:31  jbarbier
//- - last modification about the BoundaryCondition
//- <BoundaryCondition>
//-   <[type]>
//-   <[type]/>
//- <BoundaryCondition/>
//-
//- - modification of the xml files for this modification
//-
//- - version 1.2 of the xml schema
//-
//- Revision 1.37  2005/01/10 14:07:42  jbarbier
//- - file ChangeLog added for the autotools
//-
//- - xml schema corrected : BoundaryCondition modified without "choice"
//-
//- - new function of the Model to check at every moment the validity of an xml file according to the xml Schema
//-
//- Revision 1.36  2004/12/20 15:01:26  jbarbier
//- - schema XML renamed V1.1
//-
//- - schema XML corrected about OneStepNSProblem:
//-   tag OneStepNSProblem contains tags LCP, QP, ... and other tags to add
//-   further
//-
//- Revision 1.35  2004/11/26 14:10:58  jbarbier
//- - the xml schema used to check the xml file is the one given in the input xml
//- file
//-
//- Revision 1.34  2004/09/27 13:27:14  jbarbier
//-
//- - Siconos schema renamed : SiconosModelSchema-V1.0.xsd
//-
//- - new required tags of the model : title, author, description, date, xmlSchema.
//- They replace previous attributes author, description and date of the Model.
//-
//- Revision 1.33  2004/09/14 13:49:59  jbarbier
//- - files added in sample/ to run run the main_siconos test program
//-
//- - all the platform can now be saved in an XML file when it is created manually
//-
//- Revision 1.32  2004/08/20 07:34:23  jbarbier
//- - creation of Model, NSDS in comand program succeed in creating SiconosModelXML,
//- NSDSXML
//-
//- Revision 1.31  2004/08/12 14:28:38  jbarbier
//- - createTimeDiscretisation in progress
//-
//- Revision 1.30  2004/08/11 14:43:45  jbarbier
//- - beginning of the mechanism of creation without XML input file of the objects of the platform with the
//- creatObjtect methods
//-
//- - function saveWToXML for Moreau integrator, and same specific functions to save
//- M,q and Q,p for LCP and QP
//-
//- - function to check coherency of the Model
//-
//- Revision 1.29  2004/08/10 14:51:49  jbarbier
//- - functions initialize() of the Lsodar and Adams OneStepIntegrator completed to
//- call the function initialize() of the base class
//-
//- Revision 1.28  2004/08/03 12:07:12  jbarbier
//- - all test on th eModel are successfull
//-
//- - new tests on the Model with the opening of XML file
//-
//- - link TimeDiscretisation -> Strategy
//-
//- - attribute T of the Model is now optional
//- and $Id: SiconosModelXML.h,v 1.42 2005/03/08 14:23:46 jbarbier Exp $ added
//
