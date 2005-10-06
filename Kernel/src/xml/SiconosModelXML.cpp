/* Siconos version 1.0, Copyright INRIA 2005.
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
using namespace std;

SiconosModelXML::SiconosModelXML():
  xmlSchemaFile(XML_SCHEMA), nsdsXML(NULL), strategyXML(NULL), tNode(NULL),
  t0Node(NULL), TNode(NULL), titleNode(NULL), authorNode(NULL), descriptionNode(NULL),
  dateNode(NULL), xmlSchemaNode(NULL), rootNode(NULL), timeNode(NULL), doc(NULL)
{
  IN("SiconosModelXML::SiconosModelXML()\n");
  if (getenv("SICONOSPATH") == NULL)
    RuntimeException::selfThrow("Environment variable SICONOSPATH is not defined.");

  doc = xmlNewDoc((xmlChar*)"1.0");
  if (doc == NULL)
    XMLException::selfThrow("SiconosModelXML - Creation of the document aborted");

  //------------Create and set the root element node----------//
  rootNode = xmlNewNode(NULL, (xmlChar*)MODEL_TAG.c_str());
  xmlDocSetRootElement(doc, rootNode);

  setTitle("no title defined");
  setAuthor("no author defined");
  setDescription("no description");
  setDate("no date defined");
  setXMLSchema("no specific XML Schema defined");
  OUT("SiconosModelXML::SiconosModelXML()\n");
}


SiconosModelXML::SiconosModelXML(char * siconosModelXMLFilePath):
  xmlSchemaFile("none"), nsdsXML(NULL), strategyXML(NULL), tNode(NULL),
  t0Node(NULL), TNode(NULL), titleNode(NULL), authorNode(NULL),   descriptionNode(NULL),
  dateNode(NULL), xmlSchemaNode(NULL), rootNode(NULL), timeNode(NULL), doc(NULL)
{
  /*
   * this initialize the library and check potential API mismatches
   * between the version it was compiled for and the actual shared
   * library used.
   */
  //xmlKeepBlanksDefault(0);

  if (getenv("SICONOSPATH") == NULL)
    RuntimeException::selfThrow("Environment variable SICONOSPATH is not defined.");

  if (siconosModelXMLFilePath != NULL)
  {
    //      LIBXML_TEST_VERSION

    //---------- XML input file loading ------------//
    doc = xmlParseFile(siconosModelXMLFilePath);

    if (doc == NULL)
      XMLException::selfThrow("SiconosModelXML - Document " + (string)siconosModelXMLFilePath + " not parsed successfully.");

    //------------ Retrieve the document's root element ----------//
    xmlNode *newRootNode = xmlDocGetRootElement(doc);

    if (newRootNode == NULL)
    {
      xmlFreeDoc(doc);
      XMLException::selfThrow("SiconosModelXML - Empty xml document");
    }
    else rootNode = newRootNode;

    //------------- Check document is of the right type (SiconosModel) ----- //
    if (xmlStrcmp(rootNode->name, (const xmlChar *) "SiconosModel"))
    {
      XMLException::selfThrow("SiconosModelXML - Wrong xml document type, root node !=SiconosModel.");
      xmlFreeDoc(doc);
    }

    //----------Loads the XML schema------------//
    /* 1- gets the value of the node which contains the tag "XMLSchema"
     * 2- gets the name and location of the XML shcema file
     * 3- loads the XML schema
     */
    xmlNode *xmlSchemaGiven = SiconosDOMTreeTools::findNodeChild(rootNode, SM_XMLSCHEMA);
    if (xmlSchemaGiven != NULL)
    {
      xmlSchemaFile = SiconosDOMTreeTools::getStringContentValue(xmlSchemaGiven);
      // checks if the path of the file is well formed
      unsigned int cpt = 0;
      bool goodFile = true;

      /*
       * the spaces are deleted at the end and at the begining of the string value
       */
      while ((xmlSchemaFile[cpt] == ' ') && (cpt < xmlSchemaFile.size()))
        cpt++;
      // delete spaces at the beginning
      if (cpt != 0)
        xmlSchemaFile = xmlSchemaFile.substr(cpt, xmlSchemaFile.size() - cpt);

      cpt = xmlSchemaFile.size();
      if (cpt > 0)
      {
        while ((xmlSchemaFile[cpt - 1] == ' ') && (cpt > 0))
          cpt--;
        // delete spaces at the end
        if (cpt != 0)
          xmlSchemaFile = xmlSchemaFile.substr(0, cpt);
      }
      else goodFile = false;

      if (xmlSchemaFile[0] != '/')
        xmlSchemaFile = "/" + xmlSchemaFile;

      /*
       * if there's no schema given (empty field)
       */
      if (!goodFile)
      {
        //xmlSchemaFile = DEFAULT_XMLSCHEMA;
        xmlSchemaFile = XML_SCHEMA;//DEFAULT_XMLSCHEMA;
        cout << "/!\\No XML Schema file has been given !\n# Generic schema will be used" << endl;
      }
    }
    else
    {
      //xmlSchemaFile = "/config/xmlschema/SiconosModelSchema-V1.0.xsd";
      XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_XMLSCHEMA + " not found.");
    }

    string schemaFile = getenv("SICONOSPATH") + xmlSchemaFile;

    //----------Load Schema--------//
    std::ifstream givenXMLSchema(schemaFile.c_str());
    if (givenXMLSchema == NULL)
      XMLException::selfThrow("SiconosModelXML - File \"" + schemaFile + "\" doesn't exist.");

    xmlSchemaParserCtxtPtr ctxt = xmlSchemaNewParserCtxt(schemaFile.c_str());
    xmlSchemaSetParserErrors(ctxt, (xmlSchemaValidityErrorFunc) fprintf, (xmlSchemaValidityWarningFunc) fprintf, stderr);


    //----------Verifify the Schema validity------------//
    xmlSchemaPtr schema = xmlSchemaParse(ctxt);
    xmlSchemaFreeParserCtxt(ctxt);

    if (schema == NULL)
      XMLException::selfThrow("SiconosModelXML - please correct the xml schema : " + schemaFile + ".");

    xmlSchemaValidCtxtPtr validctxt = xmlSchemaNewValidCtxt(schema);

    xmlSchemaSetValidErrors(validctxt, (xmlSchemaValidityErrorFunc) fprintf, (xmlSchemaValidityWarningFunc) fprintf, stderr);

    //----------Verififys the XML file respects the schema------------//

    int xmlValid = 0;
    //  xmlSchemaValidateDoc(validctxt, doc);


    xmlSchemaFreeValidCtxt(validctxt);
    xmlSchemaFree(schema);

    if (xmlValid == 0) cout << "SiconosModelXML - Your XML model file : " << siconosModelXMLFilePath << " is valid with respect to the schema" << endl;
    else if (xmlValid == -1)
      XMLException::selfThrow("SiconosModelXML - Internal or API error to verify your XML model file : " + (string)siconosModelXMLFilePath + ".");
    else //positive error code number returned
      XMLException::selfThrow("SiconosModelXML - Your XML model file " + (string)siconosModelXMLFilePath + " doesn't respect the siconos schema.");

    //-------------------Load NonSmoothDynamicalSystemXML, StrategyXML and to and T ---------------------//
    loadModel(rootNode);
  }
  else
  {
    // \to do, when "siconosModelXMLFilePath == NULL"

    doc = xmlNewDoc((xmlChar*)"1.0");

    if (doc == NULL)
      XMLException::selfThrow("SiconosModelXML - Creation of the document aborted");

    //------------Create and set the root element node----------//
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

  delete nsdsXML;
  nsdsXML = NULL;
  delete strategyXML;
  strategyXML = NULL;
}

void SiconosModelXML::saveSiconosModelInXMLFile(char * siconosModelXMLFilePath)
{
  xmlIndentTreeOutput = 1;//to have the indentation

  //Save in XML file if the tree changes :
  xmlSaveFormatFileEnc(siconosModelXMLFilePath, doc, ISO.c_str(), xmlIndentTreeOutput); //last parameter 1 to formatting spaces and have indentation
}


void SiconosModelXML::loadModel(xmlNode *rootNode)
{
  xmlNode *node;

  /*
   * load the information about the Model
   */
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

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, SM_XMLSCHEMA)) != NULL)
    xmlSchemaNode = node;
  else XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_XMLSCHEMA + " not found.");
  /*
   * loadModel function is called, so an xml file has been given to load the data of the system
   */
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, NSDS_TAG)) != NULL)
  {
    loadNonSmoothDynamicalSystemXML(node);
  }
  else
    XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + NSDS_TAG + " not found.");

  /*
   * loadTime function is called, so an xml file has been given to load the data of the system
   */
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, SM_TIME)) != NULL)
  {
    loadTime(node);
    timeNode = node;
  }
  else XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_TIME + " not found.");
  /*
   * Warning, in an xml input data file, the strategy is only optional
   */
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, STRATEGY_TAG)) != NULL)
  {
    loadStrategyXML(node);
  }
  else
  {
    //XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_STRATEGY + " not found.");

    /*
     * the strategy is optional so if there's no strategy, we only print to the screen a warning
     */
    cout << "SiconosModelXML - loadModel Warning : tag " << STRATEGY_TAG << " not found. The Strategy is optional" << endl;
    timeNode = node;
    strategyXML = NULL;
  }
}

void SiconosModelXML::loadModel(Model * model)
{
  IN("SiconosModelXML::loadModel(Model * model)\n");
  if (model != NULL)
  {
    xmlNode* node;
    /*
     * creation of the time node
     */
    timeNode = xmlNewChild(rootNode, NULL, (xmlChar*)SM_TIME.c_str(), NULL);
    setT0(model->getT0());
    setTCurrent(model->getCurrentT());
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
      //nsdsXML = new NonSmoothDynamicalSystemXML( xmlNewChild(rootNode, NULL, (xmlChar*)SM_STRATEGY.c_str(), NULL) );
      nsdsXML = new NonSmoothDynamicalSystemXML();

      // linkage between the NSDS and his NonSmoothDynamicalSystemXML
      model->getNonSmoothDynamicalSystemPtr()->setNonSmoothDynamicalSystemXMLPtr(nsdsXML);

      // creation of the nodes of the NSDS with the right data
      node = xmlNewChild(rootNode, NULL, (xmlChar*)NSDS_TAG.c_str(), NULL);
      nsdsXML->updateNonSmoothDynamicalSystemXML(node, model->getNonSmoothDynamicalSystemPtr());
    }
    else nsdsXML = NULL;

    /*
     * creation of the Strategy node
     */
    if (model->getStrategyPtr() != NULL)
    {
      strategyXML = new StrategyXML();
      // linkage between the Strategy and his StrategyXML
      model->getStrategyPtr()->setStrategyXMLPtr(strategyXML);

      // creation of the nodes of the Strategy with the right data
      node = xmlNewChild(rootNode, NULL, (xmlChar*)STRATEGY_TAG.c_str(), NULL);
      xmlNewProp(node, (xmlChar*)TYPE_ATTRIBUTE.c_str(), (xmlChar*)model->getStrategyPtr()->getType().c_str());
      strategyXML->updateStrategyXML(node, model->getStrategyPtr());
    }
  }
  else XMLException::selfThrow("SiconosModelXML::loadModel(Model * model) : no Model has been given.");
  OUT("SiconosModelXML::loadModel(Model * model)\n");
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


void SiconosModelXML::loadNonSmoothDynamicalSystemXML(xmlNode *NSDSNode)
{
  nsdsXML = new NonSmoothDynamicalSystemXML(NSDSNode);
}


void SiconosModelXML::loadStrategyXML(xmlNode *strategyNode)
{
  strategyXML = new StrategyXML(strategyNode, nsdsXML->getDSNumbers(), nsdsXML->getInteractionNumbers());
}

bool SiconosModelXML::checkSiconosDOMTree()
{
  //  cout<<">> SiconosModelXML::checkSiconosDOMTree()"<<endl;
  bool res = false;
  string schemaFile = getenv("SICONOSPATH") + xmlSchemaFile;

  if (doc == NULL)
    cout << "Warning: no DOM tree has been initialized yet." << endl;
  else
  {
    //----------Load Schema--------//
    xmlSchemaParserCtxtPtr ctxt = xmlSchemaNewParserCtxt(schemaFile.c_str());
    xmlSchemaSetParserErrors(ctxt, (xmlSchemaValidityErrorFunc) fprintf, (xmlSchemaValidityWarningFunc) fprintf, stderr);

    //----------Verifify the Schema validity------------//
    xmlSchemaPtr schema = xmlSchemaParse(ctxt);
    xmlSchemaFreeParserCtxt(ctxt);

    if (schema == NULL)
      XMLException::selfThrow("SiconosModelXML - please correct the xml schema : " + schemaFile + ".");
    xmlSchemaValidCtxtPtr validctxt = xmlSchemaNewValidCtxt(schema);
    xmlSchemaSetValidErrors(validctxt, (xmlSchemaValidityErrorFunc) fprintf, (xmlSchemaValidityWarningFunc) fprintf, stderr);

    //----------Verifify the XML file respects the schema------------//
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
  //  string errormsg;
  char errormsg[256];
  unsigned int i;
  IN("SiconosModelXML::checkSiconosDOMTreeCoherency\n");
  /*
   * Matrices and Vector can come from the XML file, an external file or a plugin,
   * but can only come from one of these way of storage
   */

  // checks if the NSDS is BVP or not
  bool bvp = getNonSmoothDynamicalSystemXML()->isBVP();
  if (bvp)
  {
    // the NSDS is BVP so the DynamicalSystems must have a BoundaryCondition
    vector<int> definedDS = getNonSmoothDynamicalSystemXML()->getDSNumbers();
    for (i = 0; i < definedDS.size(); i++)
    {
      if (getNonSmoothDynamicalSystemXML()->getDynamicalSystemXML(definedDS[i])->getBoundaryConditionXML() == NULL)
      {
        sprintf(errormsg, "SiconosModelXML::checkSiconosDOMTreeCoherency - the DynamicalSystem which 'id' is %i has no BoundariCondition whereas the NSDS is BVP", definedDS[i]);
        XMLException::selfThrow(errormsg);
      }
    }
  }
  else
  {
    // the NSDS is not BVP so the DynamicalSystems must have no BoundaryCondition
    vector<int> definedDS = getNonSmoothDynamicalSystemXML()->getDSNumbers();
    for (i = 0; i < definedDS.size(); i++)
    {
      if (getNonSmoothDynamicalSystemXML()->getDynamicalSystemXML(definedDS[i])->getBoundaryConditionXML() != NULL)
      {
        sprintf(errormsg, "SiconosModelXML::checkSiconosDOMTreeCoherency - Warning : the DynamicalSystem which 'id' is %i has BoundariCondition whereas the NSDS is not BVP", definedDS[i]);
        cout << errormsg << endl;
      }
    }
  }

  // checks the Matrix and Vector type
  // \todo : verification of the type of all the matrices and vectors (from XML file, from external file, from a plugin)

  OUT("SiconosModelXML::checkSiconosDOMTreeCoherency\n");
  return res;
}

int SiconosModelXML::validateXmlFile(string xmlFile, string xmlSchema)
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

  //----------Verifify the Schema validity------------//
  xmlSchemaPtr schema = xmlSchemaParse(ctxt);
  xmlSchemaFreeParserCtxt(ctxt);

  if (schema == NULL)
    XMLException::selfThrow("SiconosModelXML - please correct the xml schema : " + schemaXML + ".");

  xmlSchemaValidCtxtPtr validctxt = xmlSchemaNewValidCtxt(schema);
  xmlSchemaSetValidErrors(validctxt, (xmlSchemaValidityErrorFunc) fprintf, (xmlSchemaValidityWarningFunc) fprintf, stderr);

  //----------Loads the XML input file------------//
  xmlDoc *doc = xmlParseFile(xmlFile.c_str());
  if (doc == NULL)
    XMLException::selfThrow("SiconosModelXML - Could not find or error(s) in your model XML file : " + xmlFile + ".");

  //----------Verififys the XML file respects the schema------------//
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

