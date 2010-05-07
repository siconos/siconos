/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/

/*! \file
*/
#ifndef __MODELXML__
#define __MODELXML__

#include "SiconosPointers.hpp"
#include "SiconosDOMTreeTools.hpp"
#include "KernelConfig.h"
#include <libxml/parser.h>
#include <libxml/xmlschemas.h>

// warning: the xml_schema location must corresponds to the package name
// provided in configure.ac (AC_INIT(package name ...)).
// Here : siconos-kernel
//const std::string XML_SCHEMA = "/share/siconos-kernel/SiconosModelSchema-V1.2.xsd";
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

/** XML management for Model
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 04/04/2004
 *
 *
 *
 * SiconosModelXML allows the verification of a Siconos XML data file thanks to an XML Schema.
 * All verifications can't be done with schema : others are done in the code, as the download of the XML file (download is stopped if error).
 * This class calls NonSmoothDynamicalSystemXML and SimulationXML classes to manage Siconos NSDS and Simulation data of a simulation.
 */
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
  SP::NonSmoothDynamicalSystemXML nsdsXML;
  /** Simulation node. */
  SP::SimulationXML simulationXML;

  const std::string  SICONOSPATH;

  /** Load model data from xml file
  *  \param xmlNode * : the root node of the Model
  *  \exception XMLException : if a property of the SiconosModel (NSDS or Simulation or Time) lacks in the DOM tree
  */
  void loadModel(xmlNode *rootNode);

  /** Loads Time boundary values
  *  \param xmlNode * : the root node of the Time element
  *  \exception XMLException : if a property (t, T0) of the Time lacks in the DOM tree
  */
  void loadTime(xmlNode *timeNode);

public:
  SiconosModelXML();

  /** Build an SiconosModelXML object from a Siconos XML data file
  *   \param siconosModelXMLFilePath : the path string of the Siconos XML data file
  *   \exception XMLException : exception may be the XML file does not exist ; XML schema has wrong syntax  or the XML file does not respect it ; etc.
  */
  SiconosModelXML(const std::string& siconosModelXMLFilePath);

  /** Destroy an SiconosModelXML object
  */
  ~SiconosModelXML();


  /** Gets the rootNode of the SiconosModelXML
  *   \return xmlNode* : the rootNode
  */
  inline xmlNode* getRootNode()
  {
    return rootNode;
  }

  /** Gets the title of the model
  *   \return string
  */
  inline const std::string  title()
  {
    return  SiconosDOMTreeTools::getStringContentValue(titleNode);
  }

  /** Gets the author of the model
  *   \return string
  */
  inline const std::string  author()
  {
    return  SiconosDOMTreeTools::getStringContentValue(authorNode);
  }

  /** Gets the Description of the model
  *   \return string
  */
  inline const std::string  description()
  {
    return  SiconosDOMTreeTools::getStringContentValue(descriptionNode);
  }

  /** Gets the date of the model
  *   \return string
  */
  inline const std::string  date()
  {
    return  SiconosDOMTreeTools::getStringContentValue(dateNode);
  }

  /** Gets the XML Schema of the model
  *   \return string
  */
  inline const std::string  getXMLSchema()
  {
    return  SiconosDOMTreeTools::getStringContentValue(xmlSchemaNode);
  }

  /** allows to save the title of the DynamicalSystemXML
  *   \param string : The string s of the DynamicalSystemXML
  */
  inline void setTitle(const std::string&  s)
  {
    if (!hasTitle())
      titleNode = SiconosDOMTreeTools::createStringNode(rootNode, SM_TITLE, s);
    else SiconosDOMTreeTools::setStringContentValue(titleNode, s);
  }

  /** allows to save the author of the DynamicalSystemXML
  *   \param string : The string s of the DynamicalSystemXML
  */
  inline void setAuthor(const std::string&  s)
  {
    if (!hasAuthor())
      authorNode = SiconosDOMTreeTools::createStringNode(rootNode, SM_AUTHOR, s);
    else SiconosDOMTreeTools::setStringContentValue(authorNode, s);
  }

  /** allows to save the Description of the DynamicalSystemXML
  *   \param string : The string s of the DynamicalSystemXML
  */
  inline void setDescription(const std::string&  s)
  {
    if (! hasDescription())
      descriptionNode = SiconosDOMTreeTools::createStringNode(rootNode, SM_DESCRIPTION, s);
    else SiconosDOMTreeTools::setStringContentValue(descriptionNode, s);
  }

  /** allows to save the date of the DynamicalSystemXML
  *   \param string : The string s of the DynamicalSystemXML
  */
  inline void setDate(const std::string&  s)
  {
    if (!hasDate())
      dateNode = SiconosDOMTreeTools::createStringNode(rootNode, SM_DATE, s);
    else SiconosDOMTreeTools::setStringContentValue(dateNode, s);
  }

  /** allows to save the xml schema of the DynamicalSystemXML
  *   \param string : The string s of the DynamicalSystemXML
  */
  inline void setXMLSchema(const std::string&  s)
  {
    if (!hasXMLSchema())
      xmlSchemaNode = SiconosDOMTreeTools::createStringNode(rootNode, SM_XMLSCHEMA, s);
    else SiconosDOMTreeTools::setStringContentValue(xmlSchemaNode, s);
  }

  /** determines if the title of the model is in the DOM tree
  *   \return bool :  true if the author of the model is in the DOM tree
  */
  inline bool hasTitle() const
  {
    return titleNode ;
  }

  /** determines if the author of the model is in the DOM tree
  *   \return bool :  true if the title of the model is in the DOM tree
  */
  inline bool hasAuthor() const
  {
    return authorNode;
  }

  /** determines if the description of the model is in the DOM tree
  *   \return bool :  true if T is in the DOM tree
  */
  inline bool hasDescription() const
  {
    return descriptionNode;
  }
  /** determines if the date of the model is in the DOM tree
  *   \return bool :  true if date of the model is in the DOM tree
  */
  inline bool hasDate() const
  {
    return dateNode;
  }

  /** determines if the xml schema of the model is in the DOM tree
  *   \return bool :  true if the xml schema of the model is in the DOM tree
  */
  inline bool hasXMLSchema() const
  {
    return xmlSchemaNode;
  }

  /** determines if the current time T is in the DOM tree
  *   \return bool :  true if T is in the DOM tree
  */
  inline bool hasT() const
  {
    return (TNode);
  }

  /** determines if the current time t is in the DOM tree
  *   \return bool :  true if t is in the DOM tree
  */
  inline bool hasTCurrent() const
  {
    return (tNode);
  }

  /** Gets the value of t
  *   \return t value
  */
  inline double getTCurrent() const
  {
    return SiconosDOMTreeTools::getContentValue<double>(tNode);
  }

  /** Gets the value of t0
  *   \return t0 value
  */
  inline double t0() const
  {
    return  SiconosDOMTreeTools::getContentValue<double>(t0Node);
  }

  /** Gets the value of T
  *   \return T value
  */
  inline double getT() const
  {
    return SiconosDOMTreeTools::getContentValue<double>(TNode);
  }

  /** Sets the value of t0
  *   \param The new t0 value
  */
  inline void sett0(const double& t0)
  {
    if (!t0Node)
      t0Node = SiconosDOMTreeTools::createDoubleNode(timeNode, SM_T0, t0);
    else SiconosDOMTreeTools::setDoubleContentValue(t0Node, t0);
  }

  /** Sets the value of T
  *   \param The new T value
  */
  inline void setT(const double& T)
  {
    if (!hasT())
      TNode = SiconosDOMTreeTools::createDoubleNode(timeNode, SM_T, T);
    else SiconosDOMTreeTools::setDoubleContentValue(TNode, T);
  }

  /** Sets the value of t
  *   \param The new t value
  */
  inline void setTCurrent(const double& t)
  {
    if (!hasTCurrent())
      tNode = SiconosDOMTreeTools::createDoubleNode(timeNode, SM_T_CURRENT, t);
    else SiconosDOMTreeTools::setDoubleContentValue(tNode, t);
  }

  /** This function allows to get the NonSmoothDynamicalSystemXML
  *   \return The NonSmoothDynamicalSystemXML of the SiconosModelXML
  */
  inline SP::NonSmoothDynamicalSystemXML getNonSmoothDynamicalSystemXML() const
  {
    return nsdsXML;
  }

  /** This function allows to get the SimulationXML
  *   \return The SimulationXML of the SiconosModelXML
  */
  inline SP::SimulationXML getSimulationXML() const
  {
    return simulationXML;
  }

  /** determines if the Simulation is defined
  *   \return bool :  false if the simulationXML* is NULL
  */
  inline bool hasSimulation() const
  {
    return simulationXML ;
  }

  /** Saves the Siconos Model creating a Siconos XML data file
  *   \param siconosModelXMLFilePath the path string of the Siconos XML data file to save
  */
  void saveSiconosModelInXMLFile(const char *siconosModelXMLFilePath);

  /** checks the content of the DOM tree associated to the data of the SiconosModelXML
  * by using the XML schema of the platform
  *   \return bool : true if the XML schema is respected, else false
  *   \exception SiconosException
  */
  bool checkSiconosDOMTree();

  /** checks the coherency of the content of the DOM tree associated to the data of the SiconosModelXML
  * by using the XML schema of the platform
  *   \return bool : true if the data are coherent, else false
  *   \exception SiconosException
  */
  bool checkSiconosDOMTreeCoherency();

  /** Loads the model with all the data required to construct it
  *  \param SP::Model : the Model which contains all the data to build the SiconosModelXML
  */
  void loadModel(SP::Model);

  /**
  *  \param string :
  *  \param string :
  *  \return int :
  */
  int validateXmlFile(const std::string&  xmlFile, const std::string&  xmlSchema);

};

#endif

