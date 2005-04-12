
#include "SiconosModelXML.h"
#include <stdio.h>
#include <stdlib.h>
#include "check.h"

//#include "KernelDefaultConfig.h"

//#ifdef LIBXML_TREE_ENABLED


SiconosModelXML::SiconosModelXML()
{
  IN("SiconosModelXML::SiconosModelXML()\n");
  //xmlKeepBlanksDefault(0);
  if (getenv("SICONOSPATH") == NULL)
    RuntimeException::selfThrow("Environment variable SICONOSPATH is not defined.");

  this->titleNode = NULL;
  this->authorNode = NULL;
  this->descriptionNode = NULL;
  this->dateNode = NULL;
  this->xmlSchemaNode = NULL;

  this->t0Node = NULL;
  this->tNode = NULL;
  this->TNode = NULL;

  this->nsdsXML = NULL;
  this->strategyXML = NULL;

  this->doc = xmlNewDoc((xmlChar*)"1.0");
  if (this->doc == NULL)
    XMLException::selfThrow("SiconosModelXML - Creation of the document aborted");

  //------------Create and set the root element node----------//
  this->rootNode = xmlNewNode(NULL, (xmlChar*)MODEL_TAG.c_str());
  xmlDocSetRootElement(this->doc, this->rootNode);

  this->setTitle("no title defined");
  this->setAuthor("no author defined");
  this->setDescription("no description");
  this->setDate("no date defined");
  this->setXMLSchema("no specific XML Schema defined");

  // the default schema will be used
  xmlSchemaFile = XML_SCHEMA;//DEFAULT_XMLSCHEMA;

  //  // creation of the time node
  //  this->timeNode = xmlNewChild(this->rootNode, NULL, (xmlChar*)SM_TIME.c_str(), NULL);
  //
  //  // we must construct the NSDSXML which is required
  //  this->nsdsXML = new NSDSXML( xmlNewChild(this->rootNode, NULL, (xmlChar*)SM_STRATEGY.c_str(), NULL) );

  OUT("SiconosModelXML::SiconosModelXML()\n");
}


SiconosModelXML::SiconosModelXML(char * siconosModelXMLFilePath)
{
  /*
   * this initialize the library and check potential API mismatches
   * between the version it was compiled for and the actual shared
   * library used.
   */
  //xmlKeepBlanksDefault(0);
  this->nsdsXML = NULL;
  this->strategyXML = NULL;

  this->titleNode = NULL;
  this->authorNode = NULL;
  this->descriptionNode = NULL;
  this->dateNode = NULL;
  this->xmlSchemaNode = NULL;

  if (getenv("SICONOSPATH") == NULL)
    RuntimeException::selfThrow("Environment variable SICONOSPATH is not defined.");

  if (siconosModelXMLFilePath != NULL)
  {
    //    string schemaFile = getenv("SICONOSPATH") + XML_SCHEMA_FILE;
    //      LIBXML_TEST_VERSION

    this->doc = NULL;

    //----------Load Schema--------//
    /*
    xmlSchemaParserCtxtPtr ctxt = xmlSchemaNewParserCtxt(schemaFile.c_str());
    xmlSchemaSetParserErrors(ctxt,(xmlSchemaValidityErrorFunc) fprintf,(xmlSchemaValidityWarningFunc) fprintf, stderr);
    */

    //----------Verifify the Schema validity------------//
    /*
     xmlSchemaPtr schema = xmlSchemaParse(ctxt);
     xmlSchemaFreeParserCtxt(ctxt);

     if (schema == NULL)
     XMLException::selfThrow("SiconosModelXML - please correct the xml schema : " + schemaFile + ".");

     xmlSchemaValidCtxtPtr validctxt = xmlSchemaNewValidCtxt(schema);

     xmlSchemaSetValidErrors(validctxt, (xmlSchemaValidityErrorFunc) fprintf, (xmlSchemaValidityWarningFunc) fprintf, stderr);
    */


    //----------Loads the XML input file------------//
    this->doc = xmlParseFile(siconosModelXMLFilePath);

    if (this->doc == NULL)
      XMLException::selfThrow("SiconosModelXML - Could not find or error(s) in your model XML file : " + (string)siconosModelXMLFilePath + ".");

    //------------Verify if the root element node exists----------//
    xmlNode *rootNode = xmlDocGetRootElement(doc);

    if (rootNode == NULL)
    {
      xmlFreeDoc(doc);
      XMLException::selfThrow("SiconosModelXML - Internal error : get the root node.");
    }
    else this->rootNode = rootNode;

    //----------Loads the XML schema------------//
    /* 1- gets the value of the node which contains the tag "XMLSchema"
     * 2- gets the name and location of the XML shcema file
     * 3- loads the XML schema
     */
    xmlNode *xmlSchemaGiven = SiconosDOMTreeTools::findNodeChild(rootNode, SM_XMLSCHEMA);
    if (xmlSchemaGiven != NULL)
    {
      this->xmlSchemaFile = SiconosDOMTreeTools::getStringContentValue(xmlSchemaGiven);
      // checks if the path of the file is well formed
      int cpt = 0;
      bool wellFormed = false;
      bool goodFile = true;

      /*
       * the spaces are deleted at the end and at the begining of the string value
       */
      while ((this->xmlSchemaFile[cpt] == ' ') && (cpt < this->xmlSchemaFile.size()))
        cpt++;
      // delete spaces at the beginning
      if (cpt != 0)
        this->xmlSchemaFile = this->xmlSchemaFile.substr(cpt, this->xmlSchemaFile.size() - cpt);

      cpt = this->xmlSchemaFile.size();
      if (cpt > 0)
      {
        while ((this->xmlSchemaFile[cpt - 1] == ' ') && (cpt > 0))
          cpt--;
        // delete spaces at the end
        if (cpt != 0)
          this->xmlSchemaFile = this->xmlSchemaFile.substr(0, cpt);
      }
      else goodFile = false;

      if (this->xmlSchemaFile[0] != '/')
        this->xmlSchemaFile = "/" + this->xmlSchemaFile;

      /*
       * if there's no schema given (empty field)
       */
      if (!goodFile)
      {
        //this->xmlSchemaFile = DEFAULT_XMLSCHEMA;
        this->xmlSchemaFile = XML_SCHEMA;//DEFAULT_XMLSCHEMA;
        cout << "/!\\No XML Schema file has been given !\n# Generic schema will be used" << endl;
      }
    }
    else
    {
      //this->xmlSchemaFile = "/config/xmlschema/SiconosModelSchema-V1.0.xsd";
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

    // CORRECTION POUR EVITER seg fault sur fedora core 2 --- 06/01/2005
    //int xmlValid = 0;
    int xmlValid = xmlSchemaValidateDoc(validctxt, doc);


    xmlSchemaFreeValidCtxt(validctxt);
    xmlSchemaFree(schema);

    if (xmlValid == 0) cout << "SiconosModelXML - Your XML model file : " << siconosModelXMLFilePath << " is valid with respect to the schema" << endl;
    //     if (xmlValid == 0) printf("SiconosModelXML - Your XML model file : %s is valid with respect to the schema %s \n", siconosModelXMLFilePath,schemaFile);
    else if (xmlValid == -1)
      XMLException::selfThrow("SiconosModelXML - Internal or API error to verify your XML model file : " + (string)siconosModelXMLFilePath + ".");
    else //positive error code number returned
      XMLException::selfThrow("SiconosModelXML - Your XML model file " + (string)siconosModelXMLFilePath + " doesn't respect the siconos schema.");

    //------------Verify if the root element node exists----------//
    /*      xmlNode *rootNode = xmlDocGetRootElement(doc);

    if (rootNode==NULL)
    {
    xmlFreeDoc(doc);
    XMLException::selfThrow("SiconosModelXML - Internal error : get the root node.");
    }
    else this->rootNode = rootNode;
    */
    //-------------------Load NSDSXML, StrategyXML and to and T ---------------------//
    this->loadModel(rootNode);
  }
  else
  {
    // \to do, when "siconosModelXMLFilePath == NULL"

    this->doc = xmlNewDoc((xmlChar*)"1.0");

    if (this->doc == NULL)
      XMLException::selfThrow("SiconosModelXML - Creation of the document aborted");

    //------------Create and set the root element node----------//
    this->rootNode = xmlNewNode(NULL, (xmlChar*)MODEL_TAG.c_str());
    xmlDocSetRootElement(this->doc, this->rootNode);
  }
}

SiconosModelXML::~SiconosModelXML()
{
  /* Free up the resulting document */
  if (this->doc != NULL) xmlFreeDoc(doc);

  /* Free the global variables that may
   *  have been allocated by the parser.
   */
  xmlCleanupParser();

  if (this->nsdsXML != NULL) delete this->nsdsXML;
  if (this->strategyXML != NULL) delete this->strategyXML;
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
    this->titleNode = node;
  else XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_TITLE + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, SM_AUTHOR)) != NULL)
    this->authorNode = node;
  else XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_AUTHOR + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, SM_DESCRIPTION)) != NULL)
    this->descriptionNode = node;
  else XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_DESCRIPTION + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, SM_DATE)) != NULL)
    this->dateNode = node;
  else XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_DATE + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, SM_XMLSCHEMA)) != NULL)
    this->xmlSchemaNode = node;
  else XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_XMLSCHEMA + " not found.");


  /*
   * loadModel function is called, so an xml file has been given to load the data of the system
   */
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, NSDS_TAG)) != NULL)
  {
    this->loadNSDSXML(node);
  }
  else
    XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + NSDS_TAG + " not found.");

  /*
   * loadTime function is called, so an xml file has been given to load the data of the system
   */
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, SM_TIME)) != NULL)
  {
    this->loadTime(node);
    this->timeNode = node;
  }
  else XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_TIME + " not found.");

  /*
   * Warning, in an xml input data file, the strategy is only optional
   */
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, STRATEGY_TAG)) != NULL)
  {
    this->loadStrategyXML(node);
  }
  else
  {
    //XMLException::selfThrow("SiconosModelXML - loadModel error : tag " + SM_STRATEGY + " not found.");

    /*
     * the strategy is optional so if there's no strategy, we only print to the screen a warning
     */
    cout << "SiconosModelXML - loadModel Warning : tag " << STRATEGY_TAG << " not found. The Strategy is optional" << endl;
    this->timeNode = node;
    this->strategyXML = NULL;
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
    this->timeNode = xmlNewChild(this->rootNode, NULL, (xmlChar*)SM_TIME.c_str(), NULL);
    this->setT0(model->getT0());
    this->setTCurrent(model->getCurrentT());
    this->setT(model->getFinalT());
    this->setTitle(model->getTitle());
    this->setAuthor(model->getAuthor());
    this->setDescription(model->getDescription());
    this->setDate(model->getDate());
    this->setXMLSchema(model->getXMLSchema());

    /*
     * creation of the NSDS node
     */
    if (model->getNonSmoothDynamicalSystem() != NULL)
    {
      //this->nsdsXML = new NSDSXML( xmlNewChild(this->rootNode, NULL, (xmlChar*)SM_STRATEGY.c_str(), NULL) );
      this->nsdsXML = new NSDSXML();

      // linkage between the NSDS and his NSDSXML
      model->getNonSmoothDynamicalSystem()->setNSDSXML(this->nsdsXML);

      // creation of the nodes of the NSDS with the right data
      //this->nsdsXML->loadNSDS( model->getNSDS() );
      node = xmlNewChild(this->rootNode, NULL, (xmlChar*)NSDS_TAG.c_str(), NULL);
      this->nsdsXML->updateNSDSXML(node, model->getNonSmoothDynamicalSystem());
    }
    else this->nsdsXML = NULL;

    /*
     * creation of the Strategy node
     */
    if (model->getStrategy() != NULL)
    {
      this->strategyXML = new StrategyXML();
      // linkage between the Strategy and his StrategyXML
      model->getStrategy()->setStrategyXML(this->strategyXML);

      // creation of the nodes of the Strategy with the right data
      node = xmlNewChild(this->rootNode, NULL, (xmlChar*)STRATEGY_TAG.c_str(), NULL);
      xmlNewProp(node, (xmlChar*)TYPE_ATTRIBUTE.c_str(), (xmlChar*)model->getStrategy()->getType().c_str());
      this->strategyXML->updateStrategyXML(node, model->getStrategy());
    }
  }
  else XMLException::selfThrow("SiconosModelXML::loadModel(Model * model) : no Model has been given.");
  OUT("SiconosModelXML::loadModel(Model * model)\n");
}


void SiconosModelXML::loadTime(xmlNode *timeNode)
{
  xmlNode * node;

  if ((node = SiconosDOMTreeTools::findNodeChild(timeNode, SM_T0)) != NULL) this->t0Node = node;
  else XMLException::selfThrow("SiconosModelXML - loadTime error : tag " + SM_T0 + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(timeNode, SM_T)) != NULL) this->TNode = node;
  else TNode = NULL;

  if ((node = SiconosDOMTreeTools::findNodeChild(timeNode, SM_T_CURRENT)) != NULL) this->tNode = node;
  else tNode = NULL;
}


void SiconosModelXML::loadNSDSXML(xmlNode *NSDSNode)
{
  this->nsdsXML = new NSDSXML(NSDSNode);
}


void SiconosModelXML::loadStrategyXML(xmlNode *strategyNode)
{
  this->strategyXML = new StrategyXML(strategyNode, this->nsdsXML->getDSNumbers(), this->nsdsXML->getInteractionNumbers());
}

bool SiconosModelXML::checkSiconosDOMTree()
{
  //  cout<<">> SiconosModelXML::checkSiconosDOMTree()"<<endl;
  bool res = false;
  string schemaFile = getenv("SICONOSPATH") + xmlSchemaFile;

  if (this->doc == NULL)
    cout << "/!\ No DOM tree has been initialized yet." << endl;
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
    int xmlValid = xmlSchemaValidateDoc(validctxt, doc);

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
      this->saveSiconosModelInXMLFile("invalidDOMtree.xml");
      cout << "Invalid DOM tree saved in \"invalidDOMtree.xml\" file." << endl;
      XMLException::selfThrow("SiconosModelXML - The DOM tree in memory doesn't respect the XML schema.");
      //cout<<"WARNING : SiconosModelXML - The DOM tree in memory doesn't respect the XML schema."<<endl;
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
  int i;
  IN("SiconosModelXML::checkSiconosDOMTreeCoherency\n");
  /*
   * Matrices and Vector can come from the XML file, an external file or a plugin,
   * but can only come from one of these way of storage
   */

  // checks if the NSDS is BVP or not
  bool bvp = this->getNSDSXML()->isBVP();
  if (bvp)
  {
    // the NSDS is BVP so the DynamicalSystems must have a BoundaryCondition
    vector<int> definedDS = this->getNSDSXML()->getDSNumbers();
    for (i = 0; i < definedDS.size(); i++)
    {
      if (this->getNSDSXML()->getDSXML(definedDS[i])->getBoundaryConditionXML() == NULL)
      {
        sprintf(errormsg, "SiconosModelXML::checkSiconosDOMTreeCoherency - the DynamicalSystem which 'id' is %i has no BoundariCondition whereas the NSDS is BVP", definedDS[i]);
        XMLException::selfThrow(errormsg);
      }
    }
  }
  else
  {
    // the NSDS is not BVP so the DynamicalSystems must have no BoundaryCondition
    vector<int> definedDS = this->getNSDSXML()->getDSNumbers();
    for (i = 0; i < definedDS.size(); i++)
    {
      if (this->getNSDSXML()->getDSXML(definedDS[i])->getBoundaryConditionXML() != NULL)
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
  int xmlValid = xmlSchemaValidateDoc(validctxt, doc);

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

