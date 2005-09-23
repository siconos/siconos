
/** \class SiconosModelXML
 *   \brief This class allows to verify and download a Siconos XML data file, manage it data and save in a new Siconos XML data file
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 04/04/2004
 *
 *
 *
 * SiconosModelXML allows the verification of a Siconos XML data file thanks to an XML Schema.
 * All verifications can't be done with schema : others are done in the code, as the download of the XML file (download is stopped if error).
 * This class calls NonSmoothDynamicalSystemXML and StrategyXML classes to manage Siconos NSDS and Strategy data of a simulation.
 */
#ifndef __MODELXML__
#define __MODELXML__

#include "Model.h"

#include "SiconosDOMTreeTools.h"

#include "NonSmoothDynamicalSystemXML.h"
#include "StrategyXML.h"

#include <libxml/parser.h>
#include <libxml/xmlschemas.h>


extern unsigned int MATRIX_MAX_SIZE;
extern unsigned int VECTOR_MAX_SIZE;
extern std::string  FILE_STORAGE;
extern std::string  XML_SCHEMA;

class Model;
class NonSmoothDynamicalSystemXML;
class StrategyXML;

const std::string  ISO = "ISO-8859-1";
const std::string  SM_T_CURRENT = "t";
const std::string  SM_T = "T";
const std::string  SM_T0 = "t0";
const std::string  SM_TIME = "Time";

const std::string  SM_TITLE = "Title";
const std::string  SM_AUTHOR = "Author";
const std::string  SM_DESCRIPTION = "Description";
const std::string  SM_DATE = "Date";
const std::string  SM_XMLSCHEMA = "SchemaXML";

//const std::string  DEFAULT_XMLSCHEMA="/config/xmlschema/SiconosModelSchema-V1.2.xsd";

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

  /** \fn getNonSmoothDynamicalSystemXML()
   *   \brief This function allows to get the NonSmoothDynamicalSystemXML
   *   \return The NonSmoothDynamicalSystemXML of the SiconosModelXML
   */
  inline NonSmoothDynamicalSystemXML* getNonSmoothDynamicalSystemXML()
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
    return SiconosDOMTreeTools::getDoubleContentValue(this->TNode);
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
  inline std::string  getTitle()
  {
    return  SiconosDOMTreeTools::getStringContentValue(this->titleNode);
  }

  /** \fn string getAuthor()
   *   \brief Gets the author of the model
   *   \return string
   */
  inline std::string  getAuthor()
  {
    return  SiconosDOMTreeTools::getStringContentValue(this->authorNode);
  }

  /** \fn string getDescription()
   *   \brief Gets the Description of the model
   *   \return string
   */
  inline std::string  getDescription()
  {
    return  SiconosDOMTreeTools::getStringContentValue(this->descriptionNode);
  }

  /** \fn string getDate()
   *   \brief Gets the date of the model
   *   \return string
   */
  inline std::string  getDate()
  {
    return  SiconosDOMTreeTools::getStringContentValue(this->dateNode);
  }

  /** \fn string getXMLSchema()
   *   \brief Gets the XML Schema of the model
   *   \return string
   */
  inline std::string  getXMLSchema()
  {
    return  SiconosDOMTreeTools::getStringContentValue(this->xmlSchemaNode);
  }

  /** \fn void setTitle(string s)
   *   \brief allows to save the title of the DynamicalSystemXML
   *   \param string : The string s of the DynamicalSystemXML
   */
  inline void setTitle(std::string  s)
  {
    if (this->hasTitle() == false)
    {
      this->titleNode = SiconosDOMTreeTools::createStringNode(this->rootNode, SM_TITLE, s);
    }
    else SiconosDOMTreeTools::setStringContentValue(this->titleNode, s);
  }

  /** \fn void setAuthor(string s)
   *   \brief allows to save the author of the DynamicalSystemXML
   *   \param string : The string s of the DynamicalSystemXML
   */
  inline void setAuthor(std::string  s)
  {
    if (this->hasAuthor() == false)
    {
      this->authorNode = SiconosDOMTreeTools::createStringNode(this->rootNode, SM_AUTHOR, s);
    }
    else SiconosDOMTreeTools::setStringContentValue(this->authorNode, s);
  }

  /** \fn void setDescription(string s)
   *   \brief allows to save the Description of the DynamicalSystemXML
   *   \param string : The string s of the DynamicalSystemXML
   */
  inline void setDescription(std::string  s)
  {
    if (this->hasDescription() == false)
    {
      this->descriptionNode = SiconosDOMTreeTools::createStringNode(this->rootNode, SM_DESCRIPTION, s);
    }
    else SiconosDOMTreeTools::setStringContentValue(this->descriptionNode, s);
  }

  /** \fn void setDate(string s)
   *   \brief allows to save the date of the DynamicalSystemXML
   *   \param string : The string s of the DynamicalSystemXML
   */
  inline void setDate(std::string  s)
  {
    if (this->hasDate() == false)
    {
      this->dateNode = SiconosDOMTreeTools::createStringNode(this->rootNode, SM_DATE, s);
    }
    else SiconosDOMTreeTools::setStringContentValue(this->dateNode, s);
  }

  /** \fn void setXMLSchema(string s)
   *   \brief allows to save the xml schema of the DynamicalSystemXML
   *   \param string : The string s of the DynamicalSystemXML
   */
  inline void setXMLSchema(std::string  s)
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
  int validateXmlFile(std::string  xmlFile, std::string  xmlSchema);



private:

  std::string  xmlSchemaFile;

  NonSmoothDynamicalSystemXML *nsdsXML;
  StrategyXML *strategyXML;

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
   *   \brief Loads the model : Time (to and T) and the NonSmoothDynamicalSystemXML and StrategyXML of the model in order to manage NSDS and Strategy data of the simulation
   *  \param xmlNode * : the root node of the Model
   *   \exception XMLException : if a property of the SiconosModel (NSDS or Strategy or Time) lacks in the DOM tree
   */
  void loadModel(xmlNode *rootNode);

  /** \fn loadTime(xmlNode *)
   *   \brief Loads Time (to and T)
   *    \param xmlNode * : the root node of the Time element
   *   \exception XMLException : if a property (t, T0) of the Time lacks in the DOM tree
   */
  void loadTime(xmlNode *timeNode);

  /** \fn loadNonSmoothDynamicalSystemXML(xmlNode *)
   *   \brief Loads NonSmoothDynamicalSystemXML
   *    \param xmlNode * : the root node of the
   */
  void loadNonSmoothDynamicalSystemXML(xmlNode *NSDSNode);

  /** \fn loadStrategyXML(xmlNode *)
   *   \brief Loads StrategyXML
   *    \param xmlNode * : the root node of the
   */
  void loadStrategyXML(xmlNode *strategyNode);


  const std::string  SICONOSPATH;
};

#endif

