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

#include "Model.h"
#include "SiconosModelXML.h"
#include "config.h"
#include "SimulationXML.h"
#include "NonSmoothDynamicalSystem.h"
#include "NonSmoothDynamicalSystemXML.h"
#include "TimeDiscretisation.h"
#include "TimeStepping.h"
#include "EventDriven.h"

using namespace std;

// --- CONSTRUCTORS ---

// --- Default (private) constructor ---
Model::Model(): t(0.0), t0(0.0), T(0.0), strat(NULL), nsds(NULL),
  modelxml(NULL), title("none"), author("nobody"), description("none"),
  date("none"), xmlSchema("none")
{}

// -> xml
Model::Model(char *xmlFile):
  t(0.0), t0(0.0), T(-1.0), strat(NULL), nsds(NULL),
  modelxml(NULL), title("none"), author("nobody"), description("none"),
  date("none"), xmlSchema(XML_SCHEMA)
{
  if (xmlFile == NULL)
    RuntimeException::selfThrow("Model: xml constructor, xmlfile = NULL");

  // Built DOMtree
  modelxml = new SiconosModelXML(xmlFile);
  isAllocatedIn["ModelXml"] = true;

  // Load data (default value for T is -1)
  if (modelxml->hasT()) T = modelxml->getT();
  t0 = modelxml->getT0();
  if (modelxml->hasTCurrent()) t = modelxml->getTCurrent();
  else t = t0;
  title = modelxml->getTitle();
  author = modelxml->getAuthor();
  description = modelxml->getDescription();
  date = modelxml->getDate();
  if (modelxml->hasXMLSchema())
    xmlSchema = modelxml->getXMLSchema();

  // Memory allocation for nsds and simulation
  nsds = new NonSmoothDynamicalSystem(modelxml->getNonSmoothDynamicalSystemXML());
  isAllocatedIn["Nsds"] = true;
  if (modelxml->hasSimulation())
  {
    isAllocatedIn["Simulation"] = true;
    if (modelxml->getSimulationXML()->getSimulationXMLType() == TIMESTEPPING_TAG)
      strat = new TimeStepping(modelxml->getSimulationXML(), this);
    else if (modelxml->getSimulationXML()->getSimulationXMLType() == EVENTDRIVEN_TAG)
      strat = new EventDriven(modelxml->getSimulationXML(), this);
    else RuntimeException::selfThrow("Model: xml constructor, wrong type of simulation" + (modelxml->getSimulationXML()->getSimulationXMLType()));
  }
}

// --- From a minimum set of data ---
Model::Model(double newT0, double newT, const string& newTitle, const string& newAuthor,
             const string& newDescription, const string& newDate, const string& newSchema):
  t(newT0), t0(newT0), T(-1), strat(NULL), nsds(NULL), modelxml(NULL), title(newTitle),
  author(newAuthor), description(newDescription), date(newDate), xmlSchema(newSchema)
{
  if (newT > t0) T = newT;
  else if (newT > 0 && newT <= t0)
    RuntimeException::selfThrow("Model::constructor from min data: Warning, final T lower than t0");
  // else no T in the model!
}

Model::~Model()
{
  if (isAllocatedIn["Nsds"]) delete nsds;
  nsds = NULL;
  if (isAllocatedIn["Simulation"]) delete strat;
  strat = NULL;
  if (isAllocatedIn["ModelXml"]) delete modelxml;
  modelxml = NULL;
}

void Model::setSimulationPtr(Simulation *newPtr)
{
  // Warning: this function may be used carefully because of the links between Model and TimeDiscretisation
  // The model of the simulation input MUST be the current model.
  if (isAllocatedIn["Simulation"]) delete strat;
  strat = newPtr;
  isAllocatedIn["Simulation"] = false;
}

void Model::setNonSmoothDynamicalSystemPtr(NonSmoothDynamicalSystem *newPtr)
{
  if (isAllocatedIn["Nsds"]) delete nsds;
  nsds = newPtr;
  isAllocatedIn["Nsds"] = false;
}

void Model::setSiconosModelXMLPtr(SiconosModelXML *newPtr)
{
  if (isAllocatedIn["ModelXml"]) delete modelxml;
  modelxml = newPtr;
  isAllocatedIn["ModelXml"] = false;
}

// --- OTHER FUNCTIONS ---


// --- XML RELATED FUNCTIONS ---
void Model::saveToXMLFile(char* xmlFile)
{
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
}

void Model::saveToDOMTree()
{
  checkXMLPlatform();
  savePlatformToXML();
  checkXMLDOMTree();
}

void Model::savePlatformToXML()
{
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

  // save of the Simulation

  if (strat != NULL)
  {
    strat->getTimeDiscretisationPtr()->saveTimeDiscretisationToXML();

    if (strat->getType() == "TimeStepping")
      (static_cast<TimeStepping*>(strat))->saveSimulationToXML();
    else if (strat->getType() == "EventDriven")
      (static_cast<EventDriven*>(strat))->saveSimulationToXML();
    else RuntimeException::selfThrow("Model::savePlatformToXML - bad kind of Simulation");
  }
  else //RuntimeException::selfThrow("Model::saveToXML - object SimulationXML does not exist");
    cout << "Model::saveToXML - Warning : No Simulation is defined" << endl;
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
  if (modelxml != NULL)
  {
    if (modelxml->getNonSmoothDynamicalSystemXML() != NULL)
    {
      // we must create/update the DynamicalSystemXMLs
      nsds->getNonSmoothDynamicalSystemXMLPtr()->updateNonSmoothDynamicalSystemXML(modelxml->getNonSmoothDynamicalSystemXML()->getRootNode(), nsds);
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
      if (modelxml->getSimulationXML() == NULL)
      {
        //
        // no SimulationXML already exists, so no TimeDiscretisationXML, OneStepIntegratorXML and OneStepNSProblemXML are existing
        // because these objects are required when a Simulation is defined in the XML input file

        // we must update all the Model to do
        // the creation of the SimulationXML and of all the OneStepIntegratorXML and OneStepNSProblemXML
        //
        modelxml->loadModel(this);
        // \todo to be tested !!
      }
      else
      {
        strat->getSimulationXMLPtr()->saveSimulation2XML(modelxml->getSimulationXML()->getRootNode(), strat);
      }
    }
  }
  else
  {
    // in this case, we must create all the XML objects
    // SiconosModelXML, NonSmoothDynamicalSystemXML, SimulationXML, ...

    // to build all the XML objects, we must fold all the objects of the platform

    modelxml = new SiconosModelXML();
    isAllocatedIn["ModelXml"] = true;
    modelxml->loadModel(this);
  }
}


void Model::checkModelCoherency()
{
  // at first, checking the XML
  // if matrix and vector are well defined by example
  // if DynamicalSystems have BoundaryConditions when the NonSmoothDynamicalSystem is BVP for example

  if (modelxml->checkSiconosDOMTreeCoherency() == true) cout << "Data of the XML DOM tree are coherent." << endl;
  else cout << "Warning : Data of the XML DOM tree are not coherent." << endl;
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
  cout << " =========> Model named " << title << ", written by " << author << " (" << date << ")." << endl;
  cout << " ----- Description: " << description << endl;
  cout << " ----- xml schema: " << xmlSchema << endl;
  cout << endl;
  cout << " Time runs from " << t0 << " to " << T << endl;
  cout << " Current time is " << t << endl;
  cout << endl;
  if (nsds == NULL) cout << "No NSDS linked to the Model" << endl;
  if (strat != NULL) cout << "The simulation (name: " << strat->getName() << ") is a " << strat->getType() << "." << endl;
  else cout << "No simulation attached to this model." << endl;
  cout << endl;
  cout << " ============================" << endl;
}

