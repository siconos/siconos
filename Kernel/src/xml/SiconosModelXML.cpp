
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

//$Log: SiconosModelXML.cpp,v $
//Revision 1.74  2005/03/17 16:01:10  jbarbier
//- bug in xmlSchema attribute saving of the Model fixed
//
//- bug in overwriting solving algorithm fixed
//
//- nice indentation for new nodes added into input xml file removed because of limitation of libxml. The way to have the nice indentation was creating phantom tags.
//
//Revision 1.73  2005/03/17 10:58:34  jbarbier
//- "xmlKeepBlanksDefault (0);" added to SiconosModelXML to have a nice indented xml output when an xml input file
//
//Revision 1.72  2005/03/15 14:44:04  jbarbier
//- pySiconos.i edited to remove local paths
//
//- checkCoherency checks whether the DSInputOutputs and EqualityConstraints have unique numbers
//
//Revision 1.71  2005/03/15 09:57:48  jbarbier
//- EqualityConstraint save OK
//
//Revision 1.70  2005/03/08 14:23:45  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.69  2005/01/20 14:44:49  jbarbier
//- NSDS class renamed NonSmoothDynamicalSystem
//
//- code reduce, some comments remove
//
//Revision 1.68  2005/01/20 09:05:35  jbarbier
//- configuration file available and usable
//
//- save of vectors and matrices into external files (version 0.1)
//
//Revision 1.67  2005/01/13 14:14:40  jbarbier
//- correction in the XML output about size attribute in tags DS_Concerned and Interactoin _Concerned
//
//- modifications in creation of XML objects when saving data with partial XML input file
//
//Revision 1.66  2005/01/12 09:01:30  jbarbier
//- xmlSchemaValidateDoc now available and used in the platform
//
//Revision 1.65  2005/01/11 17:08:31  jbarbier
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
//Revision 1.64  2005/01/10 14:07:42  jbarbier
//- file ChangeLog added for the autotools
//
//- xml schema corrected : BoundaryCondition modified without "choice"
//
//- new function of the Model to check at every moment the validity of an xml file according to the xml Schema
//
//Revision 1.63  2005/01/06 10:18:09  charlety
//
//_ modifications to be autotools compliant
//
//Revision 1.62  2004/12/20 15:01:26  jbarbier
//- schema XML renamed V1.1
//
//- schema XML corrected about OneStepNSProblem:
//  tag OneStepNSProblem contains tags LCP, QP, ... and other tags to add
//  further
//
//Revision 1.61  2004/12/08 12:49:39  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.60  2004/12/06 14:10:59  jbarbier
//- more test on the path given for the xml schema in XML input files
//  (empty field is replaced by default xml schema file, else invalid files
//  launch XMLExcpetion)
//
//Revision 1.59  2004/11/26 14:10:58  jbarbier
//- the xml schema used to check the xml file is the one given in the input xml
//file
//
//Revision 1.58  2004/09/27 13:27:14  jbarbier
//
//- Siconos schema renamed : SiconosModelSchema-V1.0.xsd
//
//- new required tags of the model : title, author, description, date, xmlSchema.
//They replace previous attributes author, description and date of the Model.
//
//Revision 1.57  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.56  2004/09/21 11:49:10  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.55  2004/09/16 11:35:25  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.54  2004/09/14 13:49:59  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.53  2004/09/10 08:05:24  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.52  2004/08/23 14:30:03  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NSDS. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.51  2004/08/20 15:26:46  jbarbier
//- creation of a Model and save in the XML is ok
//- creation of a NSDS and save in the XML is ok
//- creation of a NonLinearSystemDS and save in the XML is OK
//
//Revision 1.50  2004/08/20 07:34:23  jbarbier
//- creation of Model, NSDS in comand program succeed in creating SiconosModelXML,
//NSDSXML
//
//Revision 1.49  2004/08/13 11:26:58  jbarbier
//- function createNSDS complete
//
//- function createDynamicalSystem and createLinearSystemDS complete
//
//- function  createLagrangianNLDS in progress
//
//Revision 1.48  2004/08/12 14:28:37  jbarbier
//- createTimeDiscretisation in progress
//
//Revision 1.47  2004/08/12 11:55:19  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.46  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.45  2004/08/10 14:51:49  jbarbier
//- functions initialize() of the Lsodar and Adams OneStepIntegrator completed to
//call the function initialize() of the base class
//
//Revision 1.44  2004/08/06 11:27:53  jbarbier
//- new tests with the XML and the optional attributes
//
//- tests on the save of the XML data
//
//Revision 1.43  2004/08/03 12:07:12  jbarbier
//- all test on th eModel are successfull
//
//- new tests on the Model with the opening of XML file
//
//- link TimeDiscretisation -> Strategy
//
//- attribute T of the Model is now optional
//
//Revision 1.42  2004/07/30 14:37:16  jbarbier
//- saving methods for DynamicalSystemXML and LagrangianNLDSXML
//
//Revision 1.41  2004/07/29 14:25:45  jbarbier
