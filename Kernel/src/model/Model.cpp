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

#include "Model.h"

// includes to be deleted thanks to factories ...
#include "TimeStepping.h"
#include "EventDriven.h"

using namespace std;

// --- CONSTRUCTORS ---

// -> xml
Model::Model(char *xmlFile):
  t(0.0), t0(0.0), T(-1.0), strat(NULL), nsds(NULL),
  modelxml(NULL), title("none"), author("nobody"), description("none"),
  date("none"), xmlSchema(XML_SCHEMA),
  isNsdsAllocatedIn(true), isStrategyAllocatedIn(false), isModelXmlAllocatedIn(true)
{
  IN("Model::Model xml\n");
  if (xmlFile != NULL)
  {
    // Built DOMtree
    modelxml = new SiconosModelXML(xmlFile);
    // Load data (default value for T is -1)
    if (modelxml->hasT()) T = modelxml->getT();
    t0 = modelxml->getT0();
    if (modelxml->hasTCurrent()) t = modelxml->getTCurrent();
    else t = t0;
    title = modelxml->getTitle();
    author = modelxml->getAuthor();
    description = modelxml->getDescription();
    date = modelxml->getDate();
    if (modelxml->hasXMLSchema())  xmlSchema = modelxml->getXMLSchema();
    else xmlSchema = XML_SCHEMA;
    // Memory allocation for nsds and strategy
    nsds = new NonSmoothDynamicalSystem(modelxml->getNonSmoothDynamicalSystemXML());
    if (modelxml->hasStrategy())
    {
      isStrategyAllocatedIn = true;
      if (modelxml->getStrategyXML()->getStrategyXMLType() == TIMESTEPPING_TAG)
        strat = new TimeStepping(modelxml->getStrategyXML(), this);
      else if (modelxml->getStrategyXML()->getStrategyXMLType() == EVENTDRIVEN_TAG)
        strat = new EventDriven(modelxml->getStrategyXML(), this);
      else RuntimeException::selfThrow("Model: xml constructor, wrong type of strategy" + (modelxml->getStrategyXML()->getStrategyXMLType()));
    }
  }
  else RuntimeException::selfThrow("Model: xml constructor, xmlfile = NULL");
  OUT("Model::Model xml\n");
}

// --- From a minimum set of data ---
// -> without precision on strategy ... (usefull ?)
Model::Model(const double& newT0, const double& newT, const string& newTitle, const string& newAuthor,
             const string& newDescription, const string& newDate, const string& newSchema):
  t(newT0), t0(newT0), T(-1), strat(NULL), nsds(NULL), modelxml(NULL), title(newTitle),
  author(newAuthor), description(newDescription), date(newDate), xmlSchema(newSchema),
  isNsdsAllocatedIn(false), isStrategyAllocatedIn(false), isModelXmlAllocatedIn(false)
{
  IN("Model::constructor from min data\n");
  if (newT > t0) T = newT;
  else  RuntimeException::selfThrow("Model::constructor from min data: Warning, final T lower than t0");
  OUT("Model::constructor from min data\n");
}

// -> with a clear strategy defined from a string
/*Model::Model(double newT0, string strategy_name, double newT, string newTitle, string newAuthor,
       string newDescription, string newDate, string newSchema):
  t(newT0),t0(newT0),T(-1),title(newTitle), strat(NULL), nsds(NULL), modelxml(NULL) ,
  author(newAuthor), description(newDescription),date(newDate),xmlSchema(newSchema)
{
  IN("Model::constructor from min data\n");*/
// check time boundaries ...
//if( newT > t0 ) T = newT;
//  else cout << "Model::constructor from min data: Warning, wrong new value for final T -> ignored"<< endl;
// Memory allocation for strategy object
//  if(strategy_name == TIMESTEPPING_TAG )
//  strat = new TimeStepping(this);
//else if(strategy_name == EVENTDRIVEN_TAG )
//  strat = new EventDriven(this);
//else RuntimeException::selfThrow("Model: constructor from min data, wrong strategy type");
//OUT("Model::constructor from min data\n");
//}


// -> with a clear strategy given with a pointer

//Model::Model(double newT0, Strategy * newStrategy, double newT, string newTitle, string newAuthor,
//       string newDescription, string newDate, string newSchema):
//  t(newT0),t0(newT0),T(-1),title(newTitle), strat(newStrategy), nsds(NULL), modelxml(NULL) ,
//  author(newAuthor), description(newDescription),date(newDate),xmlSchema(newSchema)
//{
//  IN("Model::constructor from min data\n");
// check time boundaries ...
//  if( newT > t0 ) T = newT;
//  else cout << "Model::constructor from min data: Warning, wrong new value for final T -> ignored"<< endl;
// link with model object in strategy
//  newStrategy->setModelPtr(this);
//  OUT("Model::constructor from min data\n");
//}

Model::~Model()
{
  IN("Model::~Model()\n");
  if (isNsdsAllocatedIn)
  {
    delete nsds;
    nsds = NULL;
  }
  if (isStrategyAllocatedIn)
  {
    delete strat;
    strat = NULL;
  }
  if (isModelXmlAllocatedIn)
  {
    delete modelxml;
    modelxml = NULL;
  }
  OUT("Model::~Model()\n");
}

// --- SETTERS ---
void Model::setStrategyPtr(Strategy *newPtr)
{
  if (isStrategyAllocatedIn) delete strat;
  strat = newPtr;
  isStrategyAllocatedIn = false;
  newPtr->setModelPtr(this);
}

void Model::setNonSmoothDynamicalSystemPtr(NonSmoothDynamicalSystem *newPtr)
{
  if (isNsdsAllocatedIn) delete nsds;
  nsds = newPtr;
  isNsdsAllocatedIn = false;
}

void Model::setSiconosModelXMLPtr(SiconosModelXML *newPtr)
{
  if (isModelXmlAllocatedIn) delete modelxml;
  modelxml = newPtr;
  isModelXmlAllocatedIn = false;
}

// --- OTHER FUNCTIONS ---


// --- XML RELATED FUNCTIONS ---
void Model::saveToXMLFile(char* xmlFile)
{
  IN("Model::saveToXMLFile\n");

  cout << "## Model->checkXMLPlatform()" << endl;
  //   the first operation to do is to check the XML objects
  checkXMLPlatform();

  cout << "## Model->savePlatformToXML()" << endl;
  //   copy the values of the platform to the DOM tree
  savePlatformToXML();

  cout << "## Model->checkXMLDOMTree()" << endl;
  //   verifies that the DOM tree respects the XML schema
  checkXMLDOMTree();

  cout << "## Model->saveSiconosModelInXMLFile()" << endl;
  //   saves in a file the DOM tree
  modelxml->saveSiconosModelInXMLFile(xmlFile);
  OUT("Model::saveToXMLFile\n");
}

void Model::saveToDOMTree()
{
  IN("Model::saveToDOMTree\n");
  checkXMLPlatform();
  savePlatformToXML();
  checkXMLDOMTree();
  OUT("Model::saveToDOMTree\n");
}

void Model::savePlatformToXML()
{
  IN("Model::savePlatformToXML\n");

  // update of the data of the Model
  modelxml->setT0(t0);
  modelxml->setT(T);
  modelxml->setTCurrent(t);
  modelxml->setTitle(title);
  modelxml->setAuthor(author);
  modelxml->setDescription(description);
  modelxml->setDate(date);
  modelxml->setXMLSchema(xmlSchema);

  // save of the NonSmoothDynamicalSystem
  nsds->saveNSDSToXML();

  // save of the Strategy

  if (strat != NULL)
  {
    strat->getTimeDiscretisationPtr()->saveTimeDiscretisationToXML();

    if (strat->getType() == TIMESTEPPING_STRATEGY)
      (static_cast<TimeStepping*>(strat))->saveStrategyToXML();
    else if (strat->getType() == EVENTDRIVEN_STRATEGY)
      (static_cast<EventDriven*>(strat))->saveStrategyToXML();
    else RuntimeException::selfThrow("Model::savePlatformToXML - bad kind of Strategy");
  }
  else //RuntimeException::selfThrow("Model::saveToXML - object StrategyXML does not exist");
    cout << "Model::saveToXML - Warning : No Strategy is defined" << endl;

  OUT("Model::savePlatformToXML\n");
}

bool Model::checkXMLDOMTree()
{
  bool res = false;
  if (modelxml != NULL)
    res = modelxml->checkSiconosDOMTree();

  cout << " # checkModelCoherency()" << endl;
  checkModelCoherency();
  return res;
}

void Model::checkXMLPlatform()
{
  IN("Model::checkXMLPlatform\n");

  if (modelxml != NULL)
  {
    if (modelxml->getNonSmoothDynamicalSystemXML() != NULL)
    {
      // we must create/update the DynamicalSystemXMLs
      nsds->getNonSmoothDynamicalSystemXMLPtr()->updateNonSmoothDynamicalSystemXML(modelxml->getNonSmoothDynamicalSystemXML()->getNonSmoothDynamicalSystemXMLNode(), nsds);
    }
    else if (nsds != NULL)
    {
      // creation of the NonSmoothDynamicalSystemXML and of all the DynamicalSystemXML and InteractionXML
      modelxml->loadModel(this);
      // \todo to be tested !!
    }
    else RuntimeException::selfThrow("Model::checkXMLPlatform - There's no NonSmoothDynamicalSystem in the Platform, the XML platform can't be built");

    if ((strat != NULL))
    {
      if (modelxml->getStrategyXML() == NULL)
      {
        //
        // no StrategyXML already exists, so no TimeDiscretisationXML, OneStepIntegratorXML and OneStepNSProblemXML are existing
        // because these objects are required when a Strategy is defined in the XML input file

        // we must update all the Model to do
        // the creation of the StrategyXML and of all the OneStepIntegratorXML and OneStepNSProblemXML
        //
        modelxml->loadModel(this);
        // \todo to be tested !!
      }
      else
      {
        strat->getStrategyXMLPtr()->saveStrategy2XML(modelxml->getStrategyXML()->getNode(), strat);
      }
    }
  }
  else
  {
    // in this case, we must create all the XML objects
    // SiconosModelXML, NonSmoothDynamicalSystemXML, StrategyXML, ...

    // to build all the XML objects, we must fold all the objects of the platform

    modelxml = new SiconosModelXML();
    modelxml->loadModel(this);
  }

  OUT("Model::checkXMLPlatform\n");
}


void Model::checkModelCoherency()
{
  int number;
  unsigned int i, j, k, cpt;
  char num[32];
  string error;


  // at first, checking the XML
  // if matrix and vector are well defined by example
  // if DynamicalSystems have BoundaryConditions when the NonSmoothDynamicalSystem is BVP for example

  if (modelxml->checkSiconosDOMTreeCoherency() == true) cout << "Data of the XML DOM tree are coherent." << endl;
  else cout << "Warning : Data of the XML DOM tree are not coherent." << endl;

  // we can check here other properties that the platform must have
  // the number of each EqualityConstraint must be unique
  for (i = 0; i < nsds->getEqualityConstraints().size(); i++)
  {
    for (j = i + 1; j < nsds->getEqualityConstraints().size(); j++)
    {
      if (nsds->getEqualityConstraintPtr(i)->getNumber() == nsds->getEqualityConstraintPtr(j)->getNumber())
      {
        number = nsds->getEqualityConstraintPtr(j)->getNumber();
        sprintf(num, "%d", number);
        error = "/!\\ Error, 2 EqualityConstraints have the same number : ";
        error += num;
        error += " \n";
        RuntimeException::selfThrow(error);
      }
    }
  }

  // the number of each DSInputOutput must be unique

  // we get all the DSInputOutput numbers from all the DynamicalSystems
  // and we can check if there are redundant numbers

  cpt = 0; // cpt corresponds to the size of 'vec'
  vector<int> vec;
  for (i = 0; i < nsds->getDynamicalSystems().size(); i++)
  {
    for (j = 0; j < nsds->getDynamicalSystemPtr(i)->getDSInputOutputs().size(); j++)
    {
      if (cpt == 0)
      {
        vec.push_back(nsds->getDynamicalSystemPtr(i)->getDSInputOutput(0)->getNumber());
        cpt++;
      }
      else
      {
        for (k = 0; k < cpt; k++)
          if (vec[k] != nsds->getDynamicalSystemPtr(i)->getDSInputOutput(j)->getNumber())
          {
            vec.push_back(nsds->getDynamicalSystemPtr(i)->getDSInputOutput(j)->getNumber());
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

int Model::xmlSchemaValidated(string xmlFile, string xmlSchema)
{
  int res;
  cout << "int Model::xmlSchemaValidated(string xmlFile, string xmlSchema)" << endl;
  res = modelxml->validateXmlFile(xmlFile, xmlSchema);
  return res;
}


// --- OTHER FUNCTIONS ---

void Model::display() const
{
  IN("Model::display\n");
  cout << " ===== Model display =====" << endl;
  cout << "| current time = " << t << endl;
  cout << "| t0 (initial time) = " << t0 << endl;
  cout << "| T (final time) = " << T << endl;
  cout << "| &strategy = " << endl;
  if (strat != NULL) cout << strat << endl;
  else cout << "-> NULL" << endl;
  cout << "| &nsds = " << endl;
  if (nsds != NULL) cout << nsds << endl;
  else cout << "-> NULL" << endl;
  cout << "| &modelxml = " << endl;
  if (modelxml != NULL) cout << modelxml << endl;
  else cout << "-> NULL" << endl;
  cout << "| author = " << author << endl;
  cout << "| description = " << description << endl;
  cout << "| date = " << date << endl;
  cout << "| title = " << title << endl;
  cout << "| xmlSchema = " << xmlSchema << endl;
  cout << "============================" << endl;
  OUT("Model::display\n");
}

//=======================================================
//
// function to create the platform from a C++ programm
//
//=======================================================

Strategy* Model::createStrategy(string type)
{
  if (isStrategyAllocatedIn) delete strat;
  isStrategyAllocatedIn = false;
  strat = NULL;
  if (type == TIMESTEPPING_STRATEGY || type == "timestepping" || type == "Timestepping")
  {
    strat = new TimeStepping(NULL, this);
    isStrategyAllocatedIn = true ;
  }
  else if (type == EVENTDRIVEN_STRATEGY || type == "Eventdriven" || type == "eventdriven")
  {
    strat = new EventDriven(NULL, this);
    isStrategyAllocatedIn = true ;
  }
  else RuntimeException::selfThrow("Model::create Strategy:wrong type of strategy:" + type);
  return strat;
}

Strategy* Model::createTimeStepping()
{
  if (isStrategyAllocatedIn) delete strat;
  strat = new TimeStepping(NULL, this);
  isStrategyAllocatedIn = true ;
  return strat;
}

Strategy* Model::createTimeEventDriven()
{
  if (isStrategyAllocatedIn) delete strat;
  strat = new EventDriven(NULL, this);
  isStrategyAllocatedIn = true ;
  return strat;
}

// --- Default (private) constructor ---
Model::Model(): t(0.0), t0(0.0), T(0.0), strat(NULL), nsds(NULL),
  modelxml(NULL), title("none"), author("nobody"), description("none"),
  date("none"), xmlSchema("none")
{}
