#include "Model.h"

#include "DynamicalSystem.h"
#include "LagrangianDS.h"
#include "LagrangianLinearTIDS.h"
#include "LinearSystemDS.h"

#include "TimeStepping.h"
#include "EventDriven.h"

#include "OneStepIntegrator.h"
#include "Moreau.h"
#include "Lsodar.h"
#include "Adams.h"

#include "OneStepNSProblem.h"
#include "LCP.h"
#include "QP.h"
#include "Relay.h"

#include "check.h"


Model::Model()
{
  this->strategy = NULL;
  this->modelxml = NULL;
  this->nsds = NULL;

  // initialisation to pass through READ_UNINIT_MEM
  this->t = 0.0;
  this->t0 = 0.0;
  this->T = 0.0;

}

Model::Model(char *xmlFile, float t, float t0, float T, NonSmoothDynamicalSystem* nsds, Strategy* strategy)
{
  this->strategy = NULL;
  this->modelxml = NULL;
  this->nsds = NULL;

  // initialisation to pass through READ_UNINIT_MEM
  this->t = 0.0;
  this->t0 = 0.0;
  this->T = 0.0;

  if (xmlFile != NULL)
  {
    this->readModel(xmlFile);
    this->linkModelXML();
    this->fillModelWithModelXML();
  }
  else
  {
    /*
     * no xml file in input
     * the DOM tree must be created
     * no "linkModelXML" to do
     *
     * the call to the readModel function which must create the SiconosModelXML
     */
    this->readModel(xmlFile);

    /*
     * construction of the NonSmoothDynamicalSystem and the Strategy
     */
  }

  if (T  != -1) this->T  = T;
  if (t0 != -1) this->t0 = t0;
  if (t  != -1) this->t  = t;
}

Model::~Model()
{
  IN("Model::~Model()\n");
  if (this->strategy != NULL) delete this->strategy;
  if (this->modelxml != NULL) delete this->modelxml;
  if (this->nsds != NULL) delete this->nsds;
  OUT("Model::~Model()\n");
}


/////////////////////////

bool Model::isModelComplete(void)
{
  // different things to check, to be defined
  /**
   *\WARNING incomplete
   */
  return true;
}

void Model::readModel(char* xmlFile)
{
  IN("Model::readModel\n");
  /** if there's no xml file the constructor of SiconosModelXML must create the DOM tree */
  this->modelxml = new SiconosModelXML(xmlFile);

  OUT("Model::readModel\n");
}

void Model::saveToXMLFile(char* xmlFile)
{
  IN("Model::saveToXMLFile\n");

  cout << "## Model->checkXMLPlatform()" << endl;
  /*   * the first operation to do is to check the XML objects   */
  this->checkXMLPlatform();

  cout << "## Model->savePlatformToXML()" << endl;
  /*   * copy the values of the platform to the DOM tree   */
  this->savePlatformToXML();

  cout << "## Model->checkXMLDOMTree()" << endl;
  /*   * verifies that the DOM tree respects the XML schema   */
  this->checkXMLDOMTree();

  cout << "## Model->saveSiconosModelInXMLFile()" << endl;
  /*   * saves in a file the DOM tree   */
  this->modelxml->saveSiconosModelInXMLFile(xmlFile);
  OUT("Model::saveToXMLFile\n");
}

void Model::saveToDOMTree()
{
  IN("Model::saveToDOMTree\n");
  this->checkXMLPlatform();
  this->savePlatformToXML();
  this->checkXMLDOMTree();
  OUT("Model::saveToDOMTree\n");
}

void Model::savePlatformToXML()
{
  int size, i;

  IN("Model::savePlatformToXML\n");
  /*
   * update of the data of the Model
   */
  this->modelxml->setT0(this->t0);
  this->modelxml->setT(this->T);
  this->modelxml->setTCurrent(this->t);

  this->modelxml->setTitle(this->title);
  this->modelxml->setAuthor(this->author);
  this->modelxml->setDescription(this->description);
  this->modelxml->setDate(this->date);
  this->modelxml->setXMLSchema(this->xmlSchema);


  /*
   * save of the NonSmoothDynamicalSystem
   */
  this->nsds->saveNSDSToXML();

  //  size = this->nsds->getDSVectorSize();
  //  for(i = 0; i<size; i++)
  //    {
  //      if( this->nsds->getDynamicalSystem(i)->getType() == LNLDS )
  //  (static_cast<LagrangianDS*>(this->nsds->getDynamicalSystem(i)))->saveDSToXML();
  //      else if( this->nsds->getDynamicalSystem(i)->getType() == LTIDS )
  //  (static_cast<LagrangianLinearTIDS*>(this->nsds->getDynamicalSystem(i)))->saveDSToXML();
  //      else if( this->nsds->getDynamicalSystem(i)->getType() == LSDS )
  //  (static_cast<LinearSystemDS*>(this->nsds->getDynamicalSystem(i)))->saveDSToXML();
  //      else if( this->nsds->getDynamicalSystem(i)->getType() == NLSDS )
  //  this->nsds->getDynamicalSystem(i)->saveDSToXML();
  //      else RuntimeException::selfThrow("Model::saveToXML - bad kind of DS");
  //    }
  //
  //  size = this->nsds->getInteractionVectorSize();
  //  for(i = 0; i<size; i++)
  //    {
  //      this->nsds->getInteraction(i)->saveInteractionToXML();
  //    }

  /*
   * save of the Strategy
   */

  if (this->strategy != NULL)
  {
    this->strategy->getTimeDiscretisation()->saveTimeDiscretisationToXML();

    if (this->strategy->getType() == TIMESTEPPING_STRATEGY)
      (static_cast<TimeStepping*>(this->strategy))->saveStrategyToXML();
    else if (this->strategy->getType() == EVENTDRIVEN_STRATEGY)
      (static_cast<EventDriven*>(this->strategy))->saveStrategyToXML();
    else RuntimeException::selfThrow("Model::saveToXML - bad kind of Strategy");

    //    size = this->strategy->getOneStepIntegratorVectorSize();
    //    for(i = 0; i<size; i++)
    //    {
    //      if( this->strategy->getOneStepIntegrator(i)->getType() == MOREAU_INTEGRATOR )
    //        (static_cast<Moreau*>(this->strategy->getOneStepIntegrator(i)))->saveIntegratorToXML();
    //      else if( this->strategy->getOneStepIntegrator(i)->getType() == ADAMS_INTEGRATOR )
    //        (static_cast<Adams*>(this->strategy->getOneStepIntegrator(i)))->saveIntegratorToXML();
    //      else if( this->strategy->getOneStepIntegrator(i)->getType() == LSODAR_INTEGRATOR )
    //        (static_cast<Lsodar*>(this->strategy->getOneStepIntegrator(i)))->saveIntegratorToXML();
    //      else RuntimeException::selfThrow("Model::saveToXML - bad kind of OneStepIntegrator");
    //    }
    //
    //    if( this->strategy->getStrategyXML()->hasOneStepNSProblemXML() )
    //    {
    //      if( this->strategy->getOneStepNSProblem()->getType() == LCP_OSNSP )
    //        (static_cast<LCP*>(this->strategy->getOneStepNSProblem()))->saveNSProblemToXML();
    //      else if( this->strategy->getOneStepNSProblem()->getType() == QP_OSNSP )
    //        (static_cast<QP*>(this->strategy->getOneStepNSProblem()))->saveNSProblemToXML();
    //      else if( this->strategy->getOneStepNSProblem()->getType() == RELAY_OSNSP )
    //        (static_cast<Relay*>(this->strategy->getOneStepNSProblem()))->saveNSProblemToXML();
    //      else RuntimeException::selfThrow("Model::saveToXML - bad kind of OneStepNSProblem");
    //    }
  }
  else //RuntimeException::selfThrow("Model::saveToXML - object StrategyXML does not exist");
    cout << "Model::saveToXML - Warnig : No Strategy is defined" << endl;

  OUT("Model::savePlatformToXML\n");
}

void Model::linkModelXML(void)
{
  IN("Model::linkModelXML\n");
  if (this->modelxml != NULL)
  {
    this->nsds = new NonSmoothDynamicalSystem();
    this->nsds->createNonSmoothDynamicalSystem(this->modelxml->getNSDSXML());

    if (this->modelxml->hasStrategy())
    {
      if (this->modelxml->getStrategyXML()->getStrategyXMLType() == TIMESTEPPING_TAG)
      {
        this->strategy = new TimeStepping();
        static_cast<TimeStepping*>(this->strategy)->createStrategy(this->modelxml->getStrategyXML(), this);
      }
      else if (this->modelxml->getStrategyXML()->getStrategyXMLType() == EVENTDRIVEN_TAG)
      {
        this->strategy = new EventDriven();
        static_cast<EventDriven*>(this->strategy)->createStrategy(this->modelxml->getStrategyXML(), this);
      }
    }
    else cout << "Warning - No Strategy is defined." << endl;
  }
  else RuntimeException::selfThrow("Model::linkModelXML - modelxml == NULL");

  OUT("Model::linkModelXML\n");

}

void Model::fillModelWithModelXML()
{
  OUT("Model::fillModelWithModelXML\n");
  if (this->modelxml != NULL)
  {
    /*
     * the value defined for T when it is not defined is -1
     */
    if (this->modelxml->hasT()) this->T = this->modelxml->getT();
    else this->T = -1;

    this->t0 = this->modelxml->getT0();

    if (this->modelxml->hasTCurrent()) this->t = this->modelxml->getTCurrent();
    else this->t = this->t0;

    this->title = this->modelxml->getTitle();
    this->author = this->modelxml->getAuthor();
    this->description = this->modelxml->getDescription();
    this->date = this->modelxml->getDate();
    this->xmlSchema = this->modelxml->getXMLSchema();
  }
  else RuntimeException::selfThrow("Model::fillModelWithModelXML - object ModelXML does not exist");
  IN("Model::fillModelWithModelXML\n");

}

void Model::runSimulation(void)
{
  IN("Model::runSimulation\n");
  OUT("Model::runSimulation\n");
}

void Model::doOneStep(void)
{
  IN("Model::doOneStep\n");
  OUT("Model::doOneStep\n");
}

void Model::createModel(char *xmlFile, /*float t,*/ float t0, float T, string title, string author, string description, string date, string schema)
{
  IN("Model::createModel\n");
  cout << "Model::createModel" << endl;
  if (xmlFile != NULL)
  {
    /*
     * load of the data from a xml file
     * if values are given for t, t0 and T ( != -1 ), these values will be saved
     */
    this->readModel(xmlFile);
    cout << "Model read" << endl;
    this->fillModelWithModelXML();
    cout << "Model filled" << endl;
    this->linkModelXML();
    cout << "Model linked" << endl;

    if (title  != "title") this->title = title;
    if (author  != "author") this->author  = author;
    if (description != "description") this->description = description;
    if (date  != "date") this->date  = date;
    if ((schema != "none") && (schema  != XML_SCHEMA)) this->xmlSchema = schema;

    //this->getStrategy()->getTimeDiscretisation()->checkTimeDiscretisation();
  }
  else
  {
    /*
     * in this case, needed data are given in parameters
     */

    /*
     * no xml file in input
     * the DOM tree must be created
     * no "linkModelXML" to do
     */

    /** T final */
    if (T  != -1)
    {
      this->T  = T;
    }
    else this->T = -1;

    /** t0 initial time */
    if (t0 != -1) this->t0 = t0;
    else RuntimeException::selfThrow("Model::createModel - a value for t0 must be given");

    /** t current id optionnal, the default value is t0 which is required */
    this->t = t0;

    this->title = title;
    this->author  = author;
    this->description = description;
    this->date  = date;
    if (schema == "none") this->xmlSchema = XML_SCHEMA; //xmlSchema;
    else this->xmlSchema = schema;
  }
  OUT("Model::createModel\n");
}

void Model::createModel(float t0, float T, string title, string author, string description, string date, string xmlSchema)
{
  IN("Model::createModel\n");
  cout << "Model::createModel" << endl;
  /*
   * in this case, needed data are given in parameters
   */

  /*
   * no xml file in input
   * the DOM tree must be created
   * no "linkModelXML" to do
   */

  /** T final */
  if (T  != -1)
  {
    this->T  = T;
  }
  else this->T = -1;

  /** t0 initial time */
  if (t0 != -1) this->t0 = t0;
  else RuntimeException::selfThrow("Model::createModel - a value for t0 must be given");

  /** t current id optionnal, the default value is t0 which is required */
  this->t = t0;

  this->title = title;
  this->author  = author;
  this->description = description;
  this->date  = date;
  if (xmlSchema == "none") this->xmlSchema = XML_SCHEMA; //xmlSchema;
  else this->xmlSchema = xmlSchema;

  OUT("Model::createModel\n");
}

void Model::checkModelCoherency()
{
  int number;
  int i, j, k, cpt;
  char num[32];
  string error;

  /*
   * at first, checking the XML
   * if matrix and vector are well defined by example
   * if DynamicalSystems have BoundaryConditions when the NonSmoothDynamicalSystem is BVP by example
   * ...
   */
  if (this->modelxml->checkSiconosDOMTreeCoherency() == true) cout << "Data of the XML DOM tree are coherent." << endl;
  else cout << "Warning : Data of the XML DOM tree are not coherent." << endl;

  /*
   * we can check here other properties that the platform must have
   */
  // the number of each EqualityConstraint must be unique
  for (i = 0; i < this->nsds->getEqualityConstraints().size(); i++)
  {
    for (j = i + 1; j < this->nsds->getEqualityConstraints().size(); j++)
    {
      if (this->nsds->getEqualityConstraint(i)->getNumber() == this->nsds->getEqualityConstraint(j)->getNumber())
      {
        number = this->nsds->getEqualityConstraint(j)->getNumber();
        sprintf(num, "%d", number);
        error = "/!\\ Error, 2 EqualityConstraints have the same number : ";
        error += num;
        error += " \n";
        RuntimeException::selfThrow(error);
      }
    }
  }

  // the number of each DSInputOutput must be unique
  /*
   * we get all the DSInputOutput numbers from all the DynamicalSystems
   * and we can check if there are redundant numbers
   */
  cpt = 0; // cpt corresponds to the size of 'vec'
  vector<int> vec;
  for (i = 0; i < this->nsds->getDynamicalSystems().size(); i++)
  {
    for (j = 0; j < this->nsds->getDynamicalSystem(i)->getDSInputOutputs().size(); j++)
    {
      if (cpt == 0)
      {
        vec.push_back(this->nsds->getDynamicalSystem(i)->getDSInputOutput(0)->getNumber());
        cpt++;
      }
      else
      {
        for (k = 0; k < cpt; k++)
          if (vec[k] != this->nsds->getDynamicalSystem(i)->getDSInputOutput(j)->getNumber())
          {
            vec.push_back(this->nsds->getDynamicalSystem(i)->getDSInputOutput(j)->getNumber());
            cpt++;
            break;
          }
          else
          {
            number = vec[k];
            sprintf(num, "%d", number);
            error = "/!\\ Error, 2 DSInputOutputs have the same number : ";
            error += num;
            error += " \n";
            RuntimeException::selfThrow(error);
          }
      }
    }
  }
}

void Model::display() const
{
  IN("Model::display\n");

  cout << "| t = " << t << endl;
  cout << "| t0= " << t0 << endl;
  cout << "| T = " << T << endl;
  cout << "| &strategy = " << this->strategy << endl;
  cout << "| &nsds = " << this->nsds << endl;
  cout << "| &modelxml = " << this->modelxml << endl;
  cout << "| author = " << this->author << endl;
  cout << "| description = " << this->description << endl;
  cout << "| date = " << this->date << endl;
  cout << "| title = " << this->title << endl;
  cout << "| xmlSchema = " << this->xmlSchema << endl;
  cout << "|===========================" << endl;

  OUT("Model::display\n");
}

bool Model::checkXMLDOMTree()
{
  bool res = false;
  if (this->modelxml != NULL)
    res = this->modelxml->checkSiconosDOMTree();

  cout << " # checkModelCoherency()" << endl;
  this->checkModelCoherency();
  return res;
}

void Model::checkXMLPlatform()
{
  IN("Model::checkXMLPlatform\n");

  int i;
  if (this->modelxml != NULL)
  {
    if (this->modelxml->getNSDSXML() != NULL)
    {
      /*
       * we must check if each DynamicalSystem has an DynamicalSystemXML
       */
      //    vector<DynamicalSystem*> vectDS = this->nsds->getDynamicalSystems();
      //    for(i=0; i<vectDS.size(); i++)
      //      {
      //        //if( this->modelxml->getNSDSXML()->getDSXML( vectDS[i]->getNumber() ) == NULL )
      //        if( vectDS[i]->getDynamicalSystemXML() == NULL )
      //    {
      /*
       * we must create/update the DSXMLs
       */
      this->nsds->getNSDSXML()->updateNSDSXML(this->modelxml->getNSDSXML()->getNSDSXMLNode(), this->nsds);
      //    }
      //      }

      /*
       * we must check if each Interaction has an InteractionXML
       */
      //    vector<Interaction*> vectInter = this->nsds->getInteractions();
      //    for(i=0; i<vectInter.size(); i++)
      //      {
      //        //if( this->modelxml->getNSDSXML()->getInteractionXML( vectInter[i]->getNumber() ) == NULL )
      //        if( vectInter[i]->getInteractionXML() == NULL )
      //    {
      /*
       * we must create/update the InteractionXMLs
       */
      //      this->nsds->getNSDSXML()->updateNSDSXML( this->modelxml->getNSDSXML()->getNSDSXMLNode(), this->nsds );
      //    }
      //        else
      //    {
      //      /*
      //       * we must check if the Interaction contains
      //       * a Relation and his RelationXML
      //       * and and NonSmoothLaw with his NonSmoothLawXML
      //       */
      //      if( vectInter[i]->getRelation() != NULL )
      //        {
      //          if( vectInter[i]->getRelation()->getRelationXML() == NULL )
      //      {
      //        /*
      //         * creation of the RelationXML for this Relation
      //         */
      //        vectInter[i]->getInteractionXML()->updateInteractionXML( this->modelxml->getNSDSXML()->getInteractionXML(vectInter[i]->getNumber())->getInteractionXMLNode(), vectInter[i] );
      //      }
      //        }
      //      else RuntimeException::selfThrow("Model::checkXMLPlatform - There's no Relation in this Interaction, the XML platform can't be built");
      //
      //      if( vectInter[i]->getNonSmoothLaw() != NULL )
      //        {
      //          if( vectInter[i]->getNonSmoothLaw()->getNonSmoothLawXML() == NULL )
      //      {
      //        /*
      //         * creation of the NonSmoothLawXML for this NonSmoothLaw
      //         */
      //        vectInter[i]->getInteractionXML()->updateInteractionXML( this->modelxml->getNSDSXML()->getInteractionXML(vectInter[i]->getNumber())->getInteractionXMLNode(), vectInter[i] );
      //      }
      //        }
      //      else RuntimeException::selfThrow("Model::checkXMLPlatform - There's no NonSmoothLaw in this Interaction, the XML platform can't be built");
      //    }
      //      }
    }
    else if (this->nsds != NULL)
    {
      /*
       * creation of the NSDSXML and of all the DynamicalSystemXML and InteractionXML
       */
      this->modelxml->loadModel(this);
      // \todo to be tested !!
    }
    else RuntimeException::selfThrow("Model::checkXMLPlatform - There's no NonSmoothDynamicalSystem in the Platform, the XML platform can't be built");

    if ((this->strategy != NULL)) //&& (this->modelxml->getStrategyXML() != NULL) )
    {
      if (this->modelxml->getStrategyXML() == NULL)
      {
        /*
         * no StrategyXML already exists, so no TimeDiscretisationXML, OneStepIntegratorXML and OneStepNSProblemXML are existing
         * because these objects are required when a Strategy is defined in the XML input file
         */

        /*
         * we must update all the Model to do
         * the creation of the StrategyXML and of all the OneStepIntegratorXML and OneStepNSProblemXML
         */
        this->modelxml->loadModel(this);
        // \todo to be tested !!
      }
      else
      {
        this->strategy->getStrategyXML()->updateStrategyXML(this->modelxml->getStrategyXML()->getNode(), this->strategy);
      }
    }
  }
  else
  {
    /*
     * in this case, we must create all the XML objects
     * SiconosModelXML, NSDSXML, StrategyXML, ...
     *
     * to build all the XML objects, we must fold all the objects of the platform
     *
     */
    this->modelxml = new SiconosModelXML();
    this->modelxml->loadModel(this);
  }

  OUT("Model::checkXMLPlatform\n");
}


/*******************************************************
 *
 * function to create the platform from a C++ programm
 *
 *//////////////////////////////////////////////////////

NonSmoothDynamicalSystem* Model::createNonSmoothDynamicalSystem(bool bvp)
{
  this->nsds = new NonSmoothDynamicalSystem();
  //  this->nsds->setBVP( bvp );
  this->nsds->createNonSmoothDynamicalSystem(NULL, bvp);//, NULL);

  return this->nsds;
}


Strategy* Model::createStrategy(string type)
{
  if (type == TIMESTEPPING_STRATEGY || type == "timestepping" || type == "Timestepping")
  {
    this->strategy = new TimeStepping();
    static_cast<TimeStepping*>(this->strategy)->createStrategy(NULL, this);
    return this->strategy;
  }
  else if (type == EVENTDRIVEN_STRATEGY || type == "Eventdriven" || type == "eventdriven")
  {
    this->strategy = new EventDriven();
    static_cast<EventDriven*>(this->strategy)->createStrategy(NULL, this);
    return this->strategy;
  }
  else
  {
    cout << "Warning : The Strategy type '" << type << "'doesn't exist !" << endl;
    return NULL;
  }
}

Strategy* Model::createTimeStepping()
{
  this->strategy = new TimeStepping();
  static_cast<TimeStepping*>(this->strategy)->createStrategy(NULL, this);
  return this->strategy;
}

Strategy* Model::createTimeEventDriven()
{
  this->strategy = new EventDriven();
  static_cast<EventDriven*>(this->strategy)->createStrategy(NULL, this);
  return this->strategy;
}

void Model::setXMLSchema(const string str)
{
  this->xmlSchema = str;
}

int Model::xmlSchemaValidated(string xmlFile, string xmlSchema)
{
  int res;
  cout << "int Model::xmlSchemaValidated(string xmlFile, string xmlSchema)" << endl;
  res = this->modelxml->validateXmlFile(xmlFile, xmlSchema);
  return res;
}


void Model::setMatrixMaxSize(const int max)
{
  MATRIX_MAX_SIZE = max;
}

void Model::setVectorMaxSize(const int max)
{
  VECTOR_MAX_SIZE = max;
}

void Model::setFileStorage(const string fs)
{
  if ((fs == N_ASCII) || (fs == N_BINARY)) FILE_STORAGE = fs;
  else cout << "/!\\ file storage not changed, new file storage method doesn't exist!" << endl;
}

