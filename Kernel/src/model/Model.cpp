//$Id: Model.cpp,v 1.83 2005/03/17 16:01:10 jbarbier Exp $
#include "Model.h"

#include "DynamicalSystem.h"
#include "LagrangianNLDS.h"
#include "LagrangianTIDS.h"
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
  //  (static_cast<LagrangianNLDS*>(this->nsds->getDynamicalSystem(i)))->saveDSToXML();
  //      else if( this->nsds->getDynamicalSystem(i)->getType() == LTIDS )
  //  (static_cast<LagrangianTIDS*>(this->nsds->getDynamicalSystem(i)))->saveDSToXML();
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

//$Log: Model.cpp,v $
//Revision 1.83  2005/03/17 16:01:10  jbarbier
//- bug in xmlSchema attribute saving of the Model fixed
//
//- bug in overwriting solving algorithm fixed
//
//- nice indentation for new nodes added into input xml file removed because of limitation of libxml. The way to have the nice indentation was creating phantom tags.
//
//Revision 1.82  2005/03/16 10:49:17  jbarbier
//- update of the pySiconos.i with the new header files of the platform
//
//- new createModel function without paramater for an xml input file
//
//Revision 1.81  2005/03/15 14:44:03  jbarbier
//- pySiconos.i edited to remove local paths
//
//- checkCoherency checks whether the DSInputOutputs and EqualityConstraints have unique numbers
//
//Revision 1.80  2005/03/14 16:05:26  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.79  2005/03/11 15:06:20  jbarbier
//- save to XML methods of EqualityConstraint and DSInputOutput added
//
//- XML loading process modified : Model loads NSDS, then NSDS loads the DynamicalSystems, EqualityConstraints, Interactions; Modle loads Strategy, then Strategy loads TimeDiscretisation, then the Integrators, then the OneStepNSProblem
//
//Revision 1.78  2005/03/09 15:30:20  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.77  2005/03/08 14:23:41  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.76  2005/02/14 09:52:25  charlety
//_ getters / setters put inline
//
//Revision 1.75  2005/02/01 11:08:41  charlety
//
//_ some displays of values during computations suppressed.
//
//Revision 1.74  2005/01/20 14:44:47  jbarbier
//- NSDS class renamed NonSmoothDynamicalSystem
//
//- code reduce, some comments remove
//
//Revision 1.73  2005/01/20 09:05:34  jbarbier
//- configuration file available and usable
//
//- save of vectors and matrices into external files (version 0.1)
//
//Revision 1.72  2005/01/18 17:07:36  charlety
//
//_ added autotools makefiles for sample directory
//
//Revision 1.71  2005/01/14 08:46:38  jbarbier
//- in class SiconosModelXML, hasOneStepNSProblem renamed hasOneStepNSProblemXML => name change in the call in class Model
//
//Revision 1.70  2005/01/13 14:14:38  jbarbier
//- correction in the XML output about size attribute in tags DS_Concerned and Interactoin _Concerned
//
//- modifications in creation of XML objects when saving data with partial XML input file
//
//Revision 1.69  2005/01/11 17:08:29  jbarbier
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
//Revision 1.68  2005/01/10 14:07:42  jbarbier
//- file ChangeLog added for the autotools
//
//- xml schema corrected : BoundaryCondition modified without "choice"
//
//- new function of the Model to check at every moment the validity of an xml file according to the xml Schema
//
//Revision 1.67  2004/12/20 15:01:25  jbarbier
//- schema XML renamed V1.1
//
//- schema XML corrected about OneStepNSProblem:
//  tag OneStepNSProblem contains tags LCP, QP, ... and other tags to add
//  further
//
//Revision 1.66  2004/12/08 13:53:54  jbarbier
//- Numerics : name of the methode in structures of the SiconosNumerics.h changed
//from char* to char[64] to simplify the use
//
//- Kernel : t, T, and t0 initialised to pass through READ_UNINIT_MEM
//
//Revision 1.65  2004/12/08 12:49:33  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.64  2004/09/28 14:19:39  jbarbier
//
//- new test of the model : manual creation of a Model, a NonSmoothDynamicalSystem, some dynamical
//systems, a Strategy and some integrators.
//
//- new function : Model::saveToDOMTree(), to create/update the XML objects
//and to save the data of the platform in the DOM tree, without saving to a file.
//
//Revision 1.63  2004/09/27 13:27:12  jbarbier
//
//- Siconos schema renamed : SiconosModelSchema-V1.0.xsd
//
//- new required tags of the model : title, author, description, date, xmlSchema.
//They replace previous attributes author, description and date of the Model.
//
//Revision 1.62  2004/09/23 14:09:23  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.61  2004/09/22 10:54:43  jbarbier
//- light modification according to the attribute mass of the lagrangian dynamical
//systems. The lagrangianNLDS take always an function from a plugin to compute the
//mass, whereas the lagrangianTIDS needs only a matrix.
//
//- xml input files have been modified in consequence
//
//Revision 1.60  2004/09/16 11:35:24  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.59  2004/09/14 13:49:53  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.58  2004/09/10 11:26:05  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.57  2004/09/09 08:57:44  jbarbier
//- functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//createTimeDiscretisation of the Strategy done.
//
//=> all functions to create manually the objects of the platform are done
//
//Revision 1.56  2004/09/07 07:31:10  jbarbier
//- create strategies methods of the Model done
//
//Revision 1.55  2004/08/26 11:15:55  jbarbier
//- functions of the Model : saveToXML renamed savePlatformToXML, and saveModel
//renamed saveToXMLFile
//
//Revision 1.54  2004/08/20 15:26:43  jbarbier
//- creation of a Model and save in the XML is ok
//- creation of a NonSmoothDynamicalSystem and save in the XML is ok
//- creation of a NonLinearSystemDS and save in the XML is OK
//
//Revision 1.53  2004/08/20 07:34:22  jbarbier
//- creation of Model, NonSmoothDynamicalSystem in comand program succeed in creating SiconosModelXML,
//NSDSXML
//
//Revision 1.52  2004/08/18 14:37:15  jbarbier
//- creation of Model, NonSmoothDynamicalSystem, Strategy(TimeStepping and EventDriven) and
//DynamicalSystem available when the creation is in a command program
//
//Revision 1.51  2004/08/13 11:26:57  jbarbier
//- function createNSDS complete
//
//- function createDynamicalSystem and createLinearSystemDS complete
//
//- function  createLagrangianNLDS in progress
//
//Revision 1.50  2004/08/12 11:55:14  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.49  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.48  2004/08/10 14:51:48  jbarbier
//- functions initialize() of the Lsodar and Adams OneStepIntegrator completed to
//call the function initialize() of the base class
//
//Revision 1.47  2004/08/09 15:00:49  jbarbier
//- changes in the cardinality of some attributes of the DynamicalSystem,
//OneStepIntegrator
//
//- modifications in classes Moreau, Lsodar, Adams for these new cardinalities
//
//- corrections in the test xml files
//
//Revision 1.46  2004/08/06 11:27:53  jbarbier
//- new tests with the XML and the optional attributes
//
//- tests on the save of the XML data
//
//Revision 1.45  2004/08/06 10:46:30  charlety
//
//_ example Oscillator in progress
//_ corrected a bug : theXML  save of a system without OneStepIntegrator was not OK.
//
//Revision 1.44  2004/08/06 08:38:00  jbarbier
//- saveToXML available for NonLinearSystemDS
//
//Revision 1.43  2004/08/05 15:13:24  charlety
//
//_ LSODAR --> Lsodar (end)
//
//Revision 1.42  2004/08/05 09:31:37  jbarbier
//- test successfull on the TimeDiscretisation for the triplet (T, t0, N),
//...
//
//- cjecking of this triplet in the createTimeDiscretisation method
//
//Revision 1.41  2004/08/04 14:51:01  jbarbier
//- new test using xml_uncomplete7.xml, test with no interaction defined in the
//XML input file
//
//- for the TimeDiscretisation, the triplet (t0,T,h), (t0,T,N) or (t0,h,N) ids
//required, the missing element is now computed
//
//Revision 1.40  2004/08/04 11:03:22  jbarbier
//- about the SiconosMemory : when a SiconosMemory has a maxSize greater than the
//number of steps in memory required by an integrator, the oldest SiconosVector
//are deleted
//
//- the way to initialize the SiconosMemory by the integrator has been updated to
//match with these changes
//
//Revision 1.39  2004/08/03 12:07:11  jbarbier
//- all test on th eModel are successfull
//
//- new tests on the Model with the opening of XML file
//
//- link TimeDiscretisation -> Strategy
//
//- attribute T of the Model is now optional
//
//Revision 1.38  2004/07/30 14:37:13  jbarbier
//- saving methods for DynamicalSystemXML and LagrangianNLDSXML
//
//Revision 1.37  2004/07/28 14:40:18  jbarbier
//- tests on the platform
//
//Revision 1.36  2004/07/27 14:56:01  jbarbier
//- functions createStrategy, createTimeDiscretisation and createIntegrator done
//
//Revision 1.35  2004/07/23 14:39:25  jbarbier
//- createModel, createNSDS, createDynamicalSystem, createBoundaryCondition OK
//it's the first step, it do the same thing that before, but the process is
//unified and it must simply add the traitment for the creation of the nodes in
//the DOM tree
//
//Revision 1.34  2004/07/12 11:22:53  charlety
//
//_ the state vector x of the dynamical system is now plugged to q and the velocity when this system is a Lagrangian one.
//
//_ A function to force a SiconosVector to be composite has been written. This is a temporary solution. We should use the operator = and the constructor by copy instead. This problem will be fixed later in the summer.
//
//Revision 1.33  2004/07/09 14:18:55  jbarbier
//-management of the optional attributes of the DynamicalSystem
//-new node t for current time in the time of the NonSmoothDynamicalSystem
//
//Revision 1.32  2004/07/05 12:45:13  jbarbier
//optionnal Matrix and Vector attributes save ... ok
//XML save  ok even without boundarycondition
//
//Revision 1.31  2004/07/02 14:40:20  jbarbier
//some IN and OUT added
//BoundaryConditon saveToXML ok
//problem after, during the call of the destructors
//
//Revision 1.30  2004/06/30 09:43:32  acary
//Correction on the Log messages
//
//Revision 1.29  2004/06/29 10:41:08  acary
//Ajout des Macros IN et OUT
//Ajout des Tag CVS Id et Log
//
