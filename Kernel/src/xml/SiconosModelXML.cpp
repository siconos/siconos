/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include "SiconosModelXML.h"
#include "Model.h"
#include "Simulation.h"
#include "NonSmoothDynamicalSystem.h"
#include "SimulationXML.h"
#include "NonSmoothDynamicalSystemXML.h"

using namespace std;

SiconosModelXML::SiconosModelXML():
  rootNode(NULL), timeNode(NULL), doc(NULL), xmlSchemaFile(XML_SCHEMA),
  titleNode(NULL), authorNode(NULL), descriptionNode(NULL), dateNode(NULL), xmlSchemaNode(NULL),
  tNode(NULL), t0Node(NULL), TNode(NULL)
{
  doc = xmlNewDoc((xmlChar*)"1.0");
  if (doc == NULL)
    XMLException::selfThrow("SiconosModelXML - Creation of the document aborted");

  //=====Create and set the root element node=====//
  rootNode = xmlNewNode(NULL, (xmlChar*)MODEL_TAG.c_str());
  xmlDocSetRootElement(doc, rootNode);

  setTitle("no title defined");
  setAuthor("no author defined");
  setDescription("no description");
  setDate("no date defined");
  setXMLSchema("no specific XML Schema defined");
}


SiconosModelXML::SiconosModelXML(char * siconosModelXMLFilePath):
  rootNode(NULL), timeNode(NULL), doc(NULL), xmlSchemaFile(XML_SCHEMA),
  titleNode(NULL), authorNode(NULL), descriptionNode(NULL), dateNode(NULL), xmlSchemaNode(NULL),
  tNode(NULL), t0Node(NULL), TNode(NULL)
{
  if (siconosModelXMLFilePath != NULL)
  {
    //===== XML input file loading =====
    doc = xmlParseFile(siconosModelXMLFilePath);

    if (doc == NULL) // check if parsing is ok
      XMLException::selfThrow("SiconosModelXML - Document " + (string)siconosModelXMLFilePath + " not parsed successfully.");

    //===== Retrieve the document's root element =====
    xmlNode *newRootNode = xmlDocGetRootElement(doc);

    if (newRootNode == NULL)
    {
      xmlFreeDoc(doc);
      XMLException::selfThrow("SiconosModelXML - Empty xml document");
    }
    else rootNode = newRootNode;

    //=====- Check if document is of the right type (SiconosModel) - //
    if (xmlStrcmp(rootNode->name, (const xmlChar *) "SiconosModel"))
    {
      XMLException::selfThrow("SiconosModelXML - Wrong xml document type, root node !=SiconosModel.");
      xmlFreeDoc(doc);
    }

    //===== Load the XML schema =====//
    // get schema location and name from the corresponding node
    xmlNodePtr xmlSchemaGiven = SiconosDOMTreeTools::findNodeChild(rootNode, SM_XMLSCHEMA);
    if (xmlSchemaGiven != NULL)
    {
      // First case: the full-path for xmlSchema is given by user in input file
      xmlSchemaFile = SiconosDOMTreeTools::getStringContentValue(xmlSchemaGiven);
    }
    else
    {
      // Second case: no schema is given by user -> default one is used: SICONOSPATH/XML_SCHEMA

      xmlSchemaFile = XML_SCHEMA;
    }

    string schemaFile = xmlSchemaFile; //getenv("SICONOSPATH") + (string)"/share/SICONOS/";

    cout << " **** The xml schema used is: " << schemaFile << " **** " << endl;

    //=====Load Schema=====
    std::ifstream givenXMLSchema(schemaFile.c_str());
    if (givenXMLSchema == NULL)
      XMLException::selfThrow("SiconosModelXML constructor. File \"" + schemaFile + "\" does not exist.");

    xmlSchemaParserCtxtPtr ctxt = xmlSchemaNewParserCtxt(schemaFile.c_str());
    xmlSchemaSetParserErrors(ctxt, (xmlSchemaValidityErrorFunc) fprintf, (xmlSchemaValidityWarningFunc) fprintf, stderr);

    //=====Verifify the Schema validity=====
    xmlSchemaPtr schema = xmlSchemaParse(ctxt);
    xmlSchemaFreeParserCtxt(ctxt);

    if (schema == NULL)
      XMLException::selfThrow("SiconosModelXML - Unvalid schema: " + schemaFile + ".");

    xmlSchemaValidCtxtPtr validctxt = xmlSchemaNewValidCtxt(schema);

    xmlSchemaSetValidErrors(validctxt, (xmlSchemaValidityErrorFunc) fprintf, (xmlSchemaValidityWarningFunc) fprintf, stderr);

    //=====Verifify that the XML file respects the schema=====

    int xmlValid = 0;

    xmlSchemaFreeValidCtxt(validctxt);
    xmlSchemaFree(schema);

    if (xmlValid == 0)
      cout << "SiconosModelXML - Your XML input file, " << siconosModelXMLFilePath << ", is valid with respect to the schema of reference." << endl;
    else if (xmlValid == -1)
      XMLException::selfThrow("SiconosModelXML - Internal or API error to verify your XML model file : " + (string)siconosModelXMLFilePath + ".");
    else //positive error code number returned
      XMLException::selfThrow("SiconosModelXML - Your XML inout file " + (string)siconosModelXMLFilePath + " does not respect the Siconos schema.");

    //===== read model-nodes values in xml files ====
    loadModel(rootNode);
  }
  else // if siconosModelXMLFilePath == NULL
  {
    doc = xmlNewDoc((xmlChar*)"3.0.0");
    if (doc == NULL)
      XMLException::selfThrow("SiconosModelXML - Creation of the document aborted");

    //===== Create and set the root element node =====//
    rootNode = xmlNewNode(NULL, (xmlChar*)MODEL_TAG.c_str());
    xmlDocSetRootElement(doc, rootNode);
  }
}

SiconosModelXML::~SiconosModelXML()
{
  /* Free up the resulting document */
  if (doc != NULL) xmlFreeDoc(doc);

  /* Free the global variables that may
   *  have been allocated by the parser.
   */
  xmlCleanupParser();

}

void SiconosModelXML::saveSiconosModelInXMLFile(const char * siconosModelXMLFilePath)
{
  xmlIndentTreeOutput = 1;//to have the indentation

  //Save in XML file if the tree changes :
  xmlSaveFormatFileEnc(siconosModelXMLFilePath, doc, ISO.c_str(), xmlIndentTreeOutput); //last parameter 1 to formatting spaces and have indentation
}


void SiconosModelXML::loadModel(xmlNode *rootNode)
{
  xmlNode *node;

  // title, author, description and date tags (required)
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, SM_TITLE)) != NULL)
    titleNode = node;
  else XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_TITLE + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, SM_AUTHOR)) != NULL)
    authorNode = node;
  else XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_AUTHOR + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, SM_DESCRIPTION)) != NULL)
    descriptionNode = node;
  else XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_DESCRIPTION + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, SM_DATE)) != NULL)
    dateNode = node;
  else XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_DATE + " not found.");

  // xml schema location and name tag (optional)
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, SM_XMLSCHEMA)) != NULL)
    xmlSchemaNode = node;

  // Non smooth dynamical system tag (required)
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, NSDS_TAG)) != NULL)
    nsdsXML.reset(new NonSmoothDynamicalSystemXML(node));
  else
    XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + NSDS_TAG + " not found.");
  // Time interval (required)
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, SM_TIME)) != NULL)
  {
    loadTime(node);
    timeNode = node;
  }
  else XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_TIME + " not found.");

  // Simulation tag (optional)
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, SIMULATION_TAG)) != NULL)
    simulationXML.reset(new SimulationXML(node));
  else
  {
    cout << " /!\\ Warning: SiconosModelXML - loadModel: no tag " << SIMULATION_TAG << " found. This may not be a problem since this tag is optional ./!\\" << endl;
    timeNode = node;
    simulationXML.reset();
  }
}

void SiconosModelXML::loadModel(SP::Model model)
{
  if (model)
  {
    xmlNode* node;
    /*
     * creation of the time node
     */
    timeNode = xmlNewChild(rootNode, NULL, (xmlChar*)SM_TIME.c_str(), NULL);
    setT0(model->getT0());
    setTCurrent(model->getCurrentTime());
    setT(model->getFinalT());
    setTitle(model->getTitle());
    setAuthor(model->getAuthor());
    setDescription(model->getDescription());
    setDate(model->getDate());
    setXMLSchema(model->getXmlSchema());

    /*
     * creation of the NSDS node
     */
    if (model->getNonSmoothDynamicalSystemPtr() != NULL)
    {
      //nsdsXML = new NonSmoothDynamicalSystemXML( xmlNewChild(rootNode, NULL, (xmlChar*)SM_SIMULATION.c_str(), NULL) );
      nsdsXML.reset(new NonSmoothDynamicalSystemXML());

      // linkage between the NSDS and his NonSmoothDynamicalSystemXML
      model->getNonSmoothDynamicalSystemPtr()->setNonSmoothDynamicalSystemXMLPtr(nsdsXML);

      // creation of the nodes of the NSDS with the right data
      node = xmlNewChild(rootNode, NULL, (xmlChar*)NSDS_TAG.c_str(), NULL);
      nsdsXML->updateNonSmoothDynamicalSystemXML(node, model->getNonSmoothDynamicalSystemPtr());
    }
    else nsdsXML.reset();

    /*
     * creation of the Simulation node
     */
    if (model->getSimulationPtr() != NULL)
    {
      simulationXML.reset(new SimulationXML());
      // linkage between the Simulation and his SimulationXML
      model->getSimulationPtr()->setSimulationXMLPtr(simulationXML);

      // creation of the nodes of the Simulation with the right data
      node = xmlNewChild(rootNode, NULL, (xmlChar*)SIMULATION_TAG.c_str(), NULL);
      xmlNewProp(node, (xmlChar*)TYPE_ATTRIBUTE.c_str(), (xmlChar*)model->getSimulationPtr()->getType().c_str());
      simulationXML->saveSimulation2XML(node, model->getSimulationPtr());
    }
  }
  else XMLException::selfThrow("SiconosModelXML::loadModel(Model * model) : no Model has been given.");
}


void SiconosModelXML::loadTime(xmlNode *timeNode)
{
  xmlNode * node;

  if ((node = SiconosDOMTreeTools::findNodeChild(timeNode, SM_T0)) != NULL) t0Node = node;
  else XMLException::selfThrow("SiconosModelXML - loadTime error : tag " + SM_T0 + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(timeNode, SM_T)) != NULL) TNode = node;
  else TNode = NULL;

  if ((node = SiconosDOMTreeTools::findNodeChild(timeNode, SM_T_CURRENT)) != NULL) tNode = node;
  else tNode = NULL;
}

bool SiconosModelXML::checkSiconosDOMTree()
{
  bool res = false;
  string schemaFile = getenv("SICONOSPATH") + xmlSchemaFile;

  if (doc == NULL)
    cout << "Warning: no DOM tree has been initialized yet." << endl;
  else
  {
    //=====Load Schema=====//
    xmlSchemaParserCtxtPtr ctxt = xmlSchemaNewParserCtxt(schemaFile.c_str());
    xmlSchemaSetParserErrors(ctxt, (xmlSchemaValidityErrorFunc) fprintf, (xmlSchemaValidityWarningFunc) fprintf, stderr);

    //=====Verifify the Schema validity=====//
    xmlSchemaPtr schema = xmlSchemaParse(ctxt);
    xmlSchemaFreeParserCtxt(ctxt);

    if (schema == NULL)
      XMLException::selfThrow("SiconosModelXML - please correct the xml schema : " + schemaFile + ".");
    xmlSchemaValidCtxtPtr validctxt = xmlSchemaNewValidCtxt(schema);
    xmlSchemaSetValidErrors(validctxt, (xmlSchemaValidityErrorFunc) fprintf, (xmlSchemaValidityWarningFunc) fprintf, stderr);

    //=====Verifify the XML file respects the schema=====//
    int xmlValid = 0;
    // xmlSchemaValidateDoc(validctxt, doc);

    xmlSchemaFreeValidCtxt(validctxt);
    xmlSchemaFree(schema);

    if (xmlValid == 0)
    {
      printf("SiconosModelXML - The DOM tree in memory respects the XML schema.\n");
      res = true;
    }
    else if (xmlValid == -1)
      XMLException::selfThrow("SiconosModelXML - Internal or API error to verify your DOM tree.");
    else //positive error code number returned
    {
      saveSiconosModelInXMLFile("invalidDOMtree.xml");
      cout << "Invalid DOM tree saved in \"invalidDOMtree.xml\" file." << endl;
      XMLException::selfThrow("SiconosModelXML - The DOM tree in memory doesn't respect the XML schema.");
    }
  }
  cout << "<< SiconosModelXML::checkSiconosDOMTree()" << endl;
  return res;
}

bool SiconosModelXML::checkSiconosDOMTreeCoherency()
{
  bool res = true;

  // checks if the NSDS is BVP or not
  // get the set of DSXML of the NonSmoothDynamicalSystemXML.
  SetOfDSXML dsSet = getNonSmoothDynamicalSystemXML()->getDynamicalSystemsXML();
  SetOfDSXMLIt it;

  // checks the Matrix and Vector type
  // \todo : verification of the type of all the matrices and vectors (from XML file, from external file, from a plugin)
  return res;
}

int SiconosModelXML::validateXmlFile(const string& xmlFile, const string& xmlSchema)
{
  int res = 0;
  string schemaXML;

  if (xmlSchema == "")
  {
    schemaXML = XML_SCHEMA;//DEFAULT_XMLSCHEMA;DEFAULT_XMLSCHEMA;
    cout << "/!\\No XML Schema file has been given !\n# Generic schema will be used" << endl;
  }
  else schemaXML = xmlSchema;

  xmlSchemaParserCtxtPtr ctxt = xmlSchemaNewParserCtxt(schemaXML.c_str());
  xmlSchemaSetParserErrors(ctxt, (xmlSchemaValidityErrorFunc) fprintf, (xmlSchemaValidityWarningFunc) fprintf, stderr);

  //=====Verifify the Schema validity=====//
  xmlSchemaPtr schema = xmlSchemaParse(ctxt);
  xmlSchemaFreeParserCtxt(ctxt);

  if (schema == NULL)
    XMLException::selfThrow("SiconosModelXML - please correct the xml schema : " + schemaXML + ".");

  xmlSchemaValidCtxtPtr validctxt = xmlSchemaNewValidCtxt(schema);
  xmlSchemaSetValidErrors(validctxt, (xmlSchemaValidityErrorFunc) fprintf, (xmlSchemaValidityWarningFunc) fprintf, stderr);

  //=====Loads the XML input file=====//
  xmlDoc *doc = xmlParseFile(xmlFile.c_str());
  if (doc == NULL)
    XMLException::selfThrow("SiconosModelXML - Could not find or error(s) in your model XML file : " + xmlFile + ".");

  //=====Verififys the XML file respects the schema=====//
  int xmlValid = 0;
  //xmlSchemaValidateDoc(validctxt, doc);

  xmlSchemaFreeValidCtxt(validctxt);
  xmlSchemaFree(schema);

  if (xmlValid == 0)
  {
    cout << "SiconosModelXML - Your XML model file : " << xmlFile << " is valid with respect to the schema" << endl;
    res = 1;
  }
  else if (xmlValid == -1)
    XMLException::selfThrow("SiconosModelXML - Internal or API error to verify your XML model file : " + xmlFile + ".");
  else //positive error code number returned
    XMLException::selfThrow("SiconosModelXML - Your XML model file " + xmlFile + " doesn't respect the siconos schema.");

  return res;
}

