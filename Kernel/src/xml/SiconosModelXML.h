/* Siconos-Kernel version 1.1.2, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/

/** \class SiconosModelXML
 *   \brief This class allows to verify and download a Siconos XML data file, manage it data and save in a new Siconos XML data file
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.1.2.
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

class Model;
class NonSmoothDynamicalSystemXML;
class StrategyXML;

// warning: the xml_schema location must corresponds to the package name
// provided in configure.ac (AC_INIT(package name ...)).
// Here : siconos-kernel
const std::string XML_SCHEMA = "/share/siconos-kernel/SiconosModelSchema-V1.2.xsd";
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

private:

  /**  Root node => named "SiconosModel". */
  xmlNodePtr rootNode;
  /** Node named "Time". */
  xmlNodePtr timeNode;
  /** Main xml document (ie node corresponding to the parsed input file). */
  xmlDocPtr doc;

  /** Xml Schema file name and location. */
  std::string  xmlSchemaFile;
  /** Model title node. */
  xmlNodePtr titleNode;
  /** Author node. */
  xmlNodePtr authorNode;
  /** Model description node. */
  xmlNodePtr descriptionNode;
  /** Model date. */
  xmlNodePtr dateNode;
  /** schema name and location node. */
  xmlNodePtr xmlSchemaNode;
  /**  */
  xmlNodePtr tNode;
  /** Initial time node. */
  xmlNodePtr t0Node;
  /** Final time node. */
  xmlNodePtr TNode;
  /** Non smooth dynamical system node. */
  NonSmoothDynamicalSystemXML *nsdsXML;
  /** Strategy node. */
  StrategyXML *strategyXML;

  const std::string  SICONOSPATH;

  /** \fn loadModel(xmlNode *)
   *  \brief Load model data from xml file
   *  \param xmlNode * : the root node of the Model
   *  \exception XMLException : if a property of the SiconosModel (NSDS or Strategy or Time) lacks in the DOM tree
   */
  void loadModel(xmlNode *rootNode);

  /** \fn loadTime(xmlNode *)
   *  \brief Loads Time boundary values
   *  \param xmlNode * : the root node of the Time element
   *  \exception XMLException : if a property (t, T0) of the Time lacks in the DOM tree
   */
  void loadTime(xmlNode *timeNode);

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
    return rootNode;
  }

  /** \fn string getTitle()
   *   \brief Gets the title of the model
   *   \return string
   */
  inline const std::string  getTitle()
  {
    return  SiconosDOMTreeTools::getStringContentValue(titleNode);
  }

  /** \fn string getAuthor()
   *   \brief Gets the author of the model
   *   \return string
   */
  inline const std::string  getAuthor()
  {
    return  SiconosDOMTreeTools::getStringContentValue(authorNode);
  }

  /** \fn string getDescription()
   *   \brief Gets the Description of the model
   *   \return string
   */
  inline const std::string  getDescription()
  {
    return  SiconosDOMTreeTools::getStringContentValue(descriptionNode);
  }

  /** \fn string getDate()
   *   \brief Gets the date of the model
   *   \return string
   */
  inline const std::string  getDate()
  {
    return  SiconosDOMTreeTools::getStringContentValue(dateNode);
  }

  /** \fn string getXMLSchema()
   *   \brief Gets the XML Schema of the model
   *   \return string
   */
  inline const std::string  getXMLSchema()
  {
    return  SiconosDOMTreeTools::getStringContentValue(xmlSchemaNode);
  }

  /** \fn void setTitle(const string& s)
   *   \brief allows to save the title of the DynamicalSystemXML
   *   \param string : The string s of the DynamicalSystemXML
   */
  inline void setTitle(const std::string&  s)
  {
    if (!hasTitle())
      titleNode = SiconosDOMTreeTools::createStringNode(rootNode, SM_TITLE, s);
    else SiconosDOMTreeTools::setStringContentValue(titleNode, s);
  }

  /** \fn void setAuthor(const string& s)
   *   \brief allows to save the author of the DynamicalSystemXML
   *   \param string : The string s of the DynamicalSystemXML
   */
  inline void setAuthor(const std::string&  s)
  {
    if (!hasAuthor())
      authorNode = SiconosDOMTreeTools::createStringNode(rootNode, SM_AUTHOR, s);
    else SiconosDOMTreeTools::setStringContentValue(authorNode, s);
  }

  /** \fn void setDescription(string s)
   *   \brief allows to save the Description of the DynamicalSystemXML
   *   \param string : The string s of the DynamicalSystemXML
   */
  inline void setDescription(const std::string&  s)
  {
    if (! hasDescription())
      descriptionNode = SiconosDOMTreeTools::createStringNode(rootNode, SM_DESCRIPTION, s);
    else SiconosDOMTreeTools::setStringContentValue(descriptionNode, s);
  }

  /** \fn void setDate(string s)
   *   \brief allows to save the date of the DynamicalSystemXML
   *   \param string : The string s of the DynamicalSystemXML
   */
  inline void setDate(const std::string&  s)
  {
    if (!hasDate())
      dateNode = SiconosDOMTreeTools::createStringNode(rootNode, SM_DATE, s);
    else SiconosDOMTreeTools::setStringContentValue(dateNode, s);
  }

  /** \fn void setXMLSchema(string s)
   *   \brief allows to save the xml schema of the DynamicalSystemXML
   *   \param string : The string s of the DynamicalSystemXML
   */
  inline void setXMLSchema(const std::string&  s)
  {
    if (!hasXMLSchema())
      xmlSchemaNode = SiconosDOMTreeTools::createStringNode(rootNode, SM_XMLSCHEMA, s);
    else SiconosDOMTreeTools::setStringContentValue(xmlSchemaNode, s);
  }

  /** \fn bool hasTitle()
   *   \brief determines if the title of the model is in the DOM tree
   *   \return bool :  true if the author of the model is in the DOM tree
   */
  inline const bool hasTitle() const
  {
    return (titleNode != NULL);
  }

  /** \fn bool hasAuthor()
   *   \brief determines if the author of the model is in the DOM tree
   *   \return bool :  true if the title of the model is in the DOM tree
   */
  inline const bool hasAuthor() const
  {
    return (authorNode != NULL);
  }

  /** \fn bool hasDescription()
   *   \brief determines if the description of the model is in the DOM tree
   *   \return bool :  true if T is in the DOM tree
   */
  inline const bool hasDescription() const
  {
    return (descriptionNode != NULL);
  }
  /** \fn bool hasDate()
   *   \brief determines if the date of the model is in the DOM tree
   *   \return bool :  true if date of the model is in the DOM tree
   */
  inline const bool hasDate() const
  {
    return (dateNode != NULL);
  }

  /** \fn bool hasXMLSchema()
   *   \brief determines if the xml schema of the model is in the DOM tree
   *   \return bool :  true if the xml schema of the model is in the DOM tree
   */
  inline const bool hasXMLSchema() const
  {
    return (xmlSchemaNode != NULL);
  }

  /** \fn bool hasT()
   *   \brief determines if the current time T is in the DOM tree
   *   \return bool :  true if T is in the DOM tree
   */
  inline const bool hasT() const
  {
    return (TNode != NULL);
  }

  /** \fn bool hasTCurrent()
   *   \brief determines if the current time t is in the DOM tree
   *   \return bool :  true if t is in the DOM tree
   */
  inline const bool hasTCurrent() const
  {
    return (tNode != NULL);
  }

  /** \fn void getTCurrent()
   *   \brief Gets the value of t
   *   \return t value
   */
  inline const double getTCurrent() const
  {
    return SiconosDOMTreeTools::getDoubleContentValue(tNode);
  }

  /** \fn double getT0()
   *   \brief Gets the value of t0
   *   \return t0 value
   */
  inline const double getT0() const
  {
    return  SiconosDOMTreeTools::getDoubleContentValue(t0Node);
  }

  /** \fn double getT()
   *   \brief Gets the value of T
   *   \return T value
   */
  inline const double getT() const
  {
    return SiconosDOMTreeTools::getDoubleContentValue(TNode);
  }

  /** \fn void setT0(double)
   *   \brief Sets the value of t0
   *   \param The new t0 value
   */
  inline void setT0(const double& t0)
  {
    if (t0Node == NULL)
      t0Node = SiconosDOMTreeTools::createDoubleNode(timeNode, SM_T0, t0);
    else SiconosDOMTreeTools::setDoubleContentValue(t0Node, t0);
  }

  /** \fn void setT(double)
   *   \brief Sets the value of T
   *   \param The new T value
   */
  inline void setT(const double& T)
  {
    if (!hasT())
      TNode = SiconosDOMTreeTools::createDoubleNode(timeNode, SM_T, T);
    else SiconosDOMTreeTools::setDoubleContentValue(TNode, T);
  }

  /** \fn void setTCurrent(double)
   *   \brief Sets the value of t
   *   \param The new t value
   */
  inline void setTCurrent(const double& t)
  {
    if (!hasTCurrent())
      tNode = SiconosDOMTreeTools::createDoubleNode(timeNode, SM_T_CURRENT, t);
    else SiconosDOMTreeTools::setDoubleContentValue(tNode, t);
  }

  /** \fn getNonSmoothDynamicalSystemXML()
   *   \brief This function allows to get the NonSmoothDynamicalSystemXML
   *   \return The NonSmoothDynamicalSystemXML of the SiconosModelXML
   */
  inline NonSmoothDynamicalSystemXML* getNonSmoothDynamicalSystemXML() const
  {
    return nsdsXML;
  }

  /** \fn getStrategyXML()
   *   \brief This function allows to get the StrategyXML
   *   \return The StrategyXML of the SiconosModelXML
   */
  inline StrategyXML* getStrategyXML() const
  {
    return strategyXML;
  }

  /** \fn bool hasStrategy()
   *   \brief determines if the Strategy is defined
   *   \return bool :  false if the strategyXML* is NULL
   */
  inline const bool hasStrategy() const
  {
    return (strategyXML != NULL);
  }

  /** \fn saveSiconosModelInXMLFile(char *siconosModelXMLFilePath)
   *   \brief Saves the Siconos Model creating a Siconos XML data file
   *   \param siconosModelXMLFilePath the path string of the Siconos XML data file to save
   */
  void saveSiconosModelInXMLFile(char *siconosModelXMLFilePath);

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
  int validateXmlFile(const std::string&  xmlFile, const std::string&  xmlSchema);

};

#endif

