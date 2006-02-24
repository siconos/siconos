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
#include "OneStepNSProblem.h"
using namespace std;

// --- CONSTRUCTORS/DESTRUCTOR ---
// Default constructor
OneStepNSProblem::OneStepNSProblem():
  nspbType("undefined"), dim(0), solver(NULL), isSolverAllocatedIn(false),  strategy(NULL), onestepnspbxml(NULL)
{}

// xml constructor
OneStepNSProblem::OneStepNSProblem(OneStepNSProblemXML* osnspbxml, Strategy* newStrat):
  nspbType("undefined"), dim(0), solver(NULL), isSolverAllocatedIn(false), strategy(newStrat), onestepnspbxml(osnspbxml)
{
  if (onestepnspbxml != NULL) // get dimension of the problem ...
  {
    if (onestepnspbxml->hasDim()) dim = onestepnspbxml->getDimNSProblem();
  }
  else RuntimeException::selfThrow("OneStepNSProblem::xml constructor, xml file=NULL");

  // read list of interactions and equality constraints concerned
  if (strategy != NULL)
  {
    interactionVector = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractions();
    ecVector = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getEqualityConstraints();
  }
  else cout << "OneStepNSPb xml-constructor - Warning: no strategy linked to OneStepPb" << endl;
}

// Constructor with given strategy and a pointer on Solver
OneStepNSProblem::OneStepNSProblem(Strategy * newStrat, Solver* newSolver):
  nspbType("undefined"), dim(0), solver(newSolver), isSolverAllocatedIn(false), strategy(newStrat), onestepnspbxml(NULL)
{
  if (strategy != NULL)
  {
    interactionVector = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractions();
    ecVector = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getEqualityConstraints();
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem:: constructor from strategy, given strategy == NULL");
}

OneStepNSProblem::~OneStepNSProblem()
{
  map< Interaction* , SiconosMatrix*>::iterator it;
  for (it = diagonalBlocksMap.begin(); it != diagonalBlocksMap.end(); it++)
  {
    SiconosMatrix * tmp = (*it).second;
    if (tmp != NULL)  delete tmp ;
    tmp = NULL;
  }

  map< Interaction* , map<Interaction *, SiconosMatrix*> >::iterator it2;
  map<Interaction *, SiconosMatrix*>::iterator it3;
  for (it2 = extraDiagonalBlocksMap.begin(); it2 != extraDiagonalBlocksMap.end(); it2++)
  {

    for (it3 = ((*it2).second).begin(); it3 != ((*it2).second).end(); it3++)
    {
      SiconosMatrix * tmp = (*it3).second;
      if (tmp != NULL)  delete tmp ;
      tmp = NULL;
    }
  }
  if (isSolverAllocatedIn) delete solver;
  solver = NULL;
  strategy = NULL;
  onestepnspbxml = NULL;
}

Interaction* OneStepNSProblem::getInteractionPtr(const unsigned int& nb)
{
  if (nb >= interactionVector.size())
    RuntimeException::selfThrow("OneStepNSProblem::getInteractionPtr(const int& nb) - number greater than size of interaction vector");
  return interactionVector[nb];
}

void OneStepNSProblem::setSolverPtr(Solver * newSolv)
{
  if (isSolverAllocatedIn) delete solver;
  solver = newSolv;
  isSolverAllocatedIn = false;
}

void OneStepNSProblem::addInteraction(Interaction *interaction)
{
  interactionVector.push_back(interaction);
}

void OneStepNSProblem::initialize()
{
  // update topology if necessary (ie take into account modifications in the NonSmoothDynamicalSystem)
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();
  if (!(topology->isUpToDate()))
    topology->updateTopology();

  updateOutput();
  updateInput();
}

void OneStepNSProblem::computeEffectiveOutput()
{

  // 3 steps to update the effective output, this for each interaction:
  //  - compute prediction for y ( -> yp), this for the r-1 first derivatives, r being
  //    the relative degree
  //  - compute indexMax using this prediction
  //  - compute effectiveIndexes, a list of the indexes for which constraints will be applied
  //
  //

  // get topology of the NonSmooth Dynamical System
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();
  // get time step
  double pasH = strategy->getTimeDiscretisationPtr()->getH();

  unsigned int i ; // index of derivation
  unsigned int j ; // relation number
  unsigned int sizeOutput; // effective size of vector y for a specific interaction
  unsigned int globalSizeOutput = 0; // effective size of global vector y (ie including all interactions) = sum sizeOutput over all interactions
  unsigned int k;
  // === loop over the interactions ===
  vector<Interaction*>::iterator it;
  for (it = interactionVector.begin(); it != interactionVector.end(); it++)
  {
    // get the output vector (values for previous time step)
    vector<SimpleVector *> yOld = (*it)->getYOld();

    // get relative degrees vector of this interaction (one relative degree for each relation!)
    vector<unsigned int> relativeDegree = topology->getRelativeDegrees(*it);
    unsigned int numberOfRelations = relativeDegree.size(); // this size corresponds to the interaction size, ie the number of relations

    // --- prediction vector ---

    // we compute yp[i], i =0..r-2. r is equal to the maximum value of all the relative degrees.
    // For the moment we consider that the interaction is homegeneous, ie all the relations have the same degree.
    // if r<2, no prediction, all relations are effective.
    unsigned int sizeYp;

    if (relativeDegree[0] != 0)  sizeYp = relativeDegree[0] - 1;
    else sizeYp = 0;

    if (sizeYp > 0)
    {
      // --- prediction vector ---

      vector<SimpleVector *> yp;
      yp.resize(sizeYp, NULL);
      // allocate and initialize yp with yOld.
      for (i = 0; i < sizeYp ; i++)
        yp[i] = new SimpleVector(*yOld[i]);
      // \todo the way prediction is calculated should be defined by user elsewhere
      *(yp[0]) = *(yOld[0]) +  0.5 * pasH * *(yOld[1]) ;

      // --- indexMax ---

      // loop from 0 to relative degree to find the first yp>0
      vector<unsigned int> indexMax;
      indexMax.resize(numberOfRelations, 0);

      NonSmoothLaw * nslaw = (*it) ->getNonSmoothLawPtr();
      string nslawType = nslaw -> getType();

      // loop over various non smooth law types
      if (nslawType == COMPLEMENTARITYCONDITIONNSLAW || nslawType == NEWTONIMPACTNSLAW)
      {
        for (j = 0; j < numberOfRelations; j++)
        {
          for (i = 0; i < sizeYp; i++)
          {
            if ((*(yp[i]))(j) <= 0)
              indexMax[j]++;
            else
              break;
          }
        }
      }
      else if (nslawType == NEWTONIMPACTFRICTIONNSLAW)
      {
        if (nspbType == "FrictionContact2D")
        {
          for (j = 0; j < numberOfRelations / 2; j++)
          {
            for (i = 0; i < sizeYp; i++)
            {
              if ((*(yp[i]))(2 * j) <= 0)
                indexMax[2 * j]++;
              else
                break;
            }
            indexMax[2 * j + 1] = indexMax[2 * j];
          }
        }
        else if (nspbType == "FrictionContact3D")
        {
          for (j = 0; j < numberOfRelations; j = j + 3)
          {
            for (i = 0; i < sizeYp; i++)
            {
              if ((*(yp[i]))(j) <= 0)
                indexMax[j]++;
              else
                break;
            }
            indexMax[j + 1] = indexMax[j];
            indexMax[j + 2] = indexMax[j];
          }
        }
        else
          RuntimeException::selfThrow("OneStepNSProblem::computeEffectiveOutput unknown Non smooth problem type: " + nslawType);
      }
      else
        RuntimeException::selfThrow("OneStepNSProblem::computeEffectiveOutput not yet implemented for non smooth law of type " + nslawType);

      topology->setIndexMax(*it, indexMax);

      for (i = 0; i < sizeYp ; i++)
        delete yp[i];

      // --- effective indexes ---

      // compute sizeOutput for the current interaction
      sizeOutput = topology->computeEffectiveSizeOutput(*it);

      vector<unsigned int> effectiveIndexes, blockIndexes, indexMin;
      indexMin = topology->getIndexMin(*it);

      effectiveIndexes.resize(sizeOutput, 0);
      blockIndexes.resize(sizeOutput, 0);

      k = 0;

      for (j = 0; j < numberOfRelations; j++)
      {
        for (i = indexMin[j]; i < indexMax[j] + 1; i++)
        {
          effectiveIndexes[k] = i + j * (relativeDegree[j]);
          blockIndexes[k] = i - indexMin[j] + j * (relativeDegree[j] - indexMin[j]);
          k++;
        }
      }
      topology->setEffectiveIndexes(*it, effectiveIndexes);
      blockIndexesMap[*it] = blockIndexes;

    }
    else
    {
      // compute sizeOutput for the current interaction
      sizeOutput = topology->computeEffectiveSizeOutput(*it);
    }
    globalSizeOutput   += sizeOutput;


  }// == end of interactions loop ==

  topology->setEffectiveSizeOutput(globalSizeOutput);

  // compute effective positions map
  topology->computeInteractionEffectivePositionMap();
}
void OneStepNSProblem::nextStep()
{
  vector<Interaction*>::iterator it;
  for (it = interactionVector.begin(); it != interactionVector.end(); it++)
    (*it)->swapInMemory();
}

void OneStepNSProblem::updateInput()
{
  vector<Interaction*>::iterator it;
  double currentTime = strategy->getModelPtr()->getCurrentT();

  for (it = interactionVector.begin(); it != interactionVector.end(); it++)
    (*it)->getRelationPtr() -> computeInput(currentTime);
}

void OneStepNSProblem::updateOutput()
{
  vector<Interaction*>::iterator it;
  double currentTime = strategy->getModelPtr()->getCurrentT();
  for (it = interactionVector.begin(); it != interactionVector.end(); it++)
    (*it)->getRelationPtr()->computeOutput(currentTime);
}

void OneStepNSProblem::compute(const double& time)
{
  RuntimeException::selfThrow("OneStepNSProblem::compute - not yet implemented for problem type =" + getType());
}

void OneStepNSProblem::saveNSProblemToXML()
{
  IN("OneStepNSProblem::saveNSProblemToXML\n");
  if (onestepnspbxml != NULL)
  {
    onestepnspbxml->setDimNSProblem(dim);
    vector<int> v;
    for (unsigned int i = 0; i < interactionVector.size(); i++)
      v.push_back(interactionVector[i]->getNumber());
    //onestepnspbxml->setInteractionConcerned( v, allInteractionConcerned() );

    /*
     * save of the solving method to XML
     */

    //    onestepnspbxml->setSolver(solvingFormalisation, methodName, normType, tolerance, maxIter, searchDirection );
  }
  else RuntimeException::selfThrow("OneStepNSProblem::saveNSProblemToXML - OneStepNSProblemXML object not exists");
  OUT("OneStepNSProblem::saveNSProblemToXML\n");
}

bool OneStepNSProblem::isOneStepNsProblemComplete() const
{
  bool isComplete = true;

  if (nspbType != "LCP" || nspbType != "FrictionContact2D" || nspbType != "FrictionContact3D" || nspbType != "QP" || nspbType != "Relay")
  {
    cout << "OneStepNSProblem is not complete: unknown problem type " << nspbType << endl;
    isComplete = false;
  }

  if (dim == 0)
  {
    cout << "OneStepNSProblem warning: problem size == 0" << endl;
    isComplete = false;
  }

  if (!(interactionVector.size() > 0))
  {
    cout << "OneStepNSProblem warning: interaction vector is empty" << endl;
    isComplete = false;
  }
  else
  {
    vector< Interaction* >::const_iterator it;
    for (it = interactionVector.begin(); it != interactionVector.end(); it++)
      if (*it == NULL)
        cout << "OneStepNSProblem warning: an interaction points to NULL" << endl;
  }

  if (!(ecVector.size() > 0))
  {
    cout << "OneStepNSProblem warning: equality constraints vector is empty" << endl;
    isComplete = false;
  }
  else
  {
    vector< EqualityConstraint* >::const_iterator it;
    for (it = ecVector.begin(); it != ecVector.end(); it++)
      if (*it == NULL)
        cout << "OneStepNSProblem warning: an equalityConstraint of the problem points to NULL" << endl;
  }

  if (strategy == NULL)
  {
    cout << "OneStepNSProblem warning: no strategy linked with the problem" << endl;
    isComplete = false;
  }

  if (solver == NULL)
  {
    cout << "OneStepNSProblem warning: no solver defined in the problem" << endl;
    isComplete = false;
  }

  return isComplete;
}

void OneStepNSProblem::check_solver(const int& info) const
{
  string solverName = solver->getSolverAlgorithmName();
  // info = 0 => ok
  // else: depend on solver
  // \todo: chose a standard output for solver parameter "info", that do not depend on the solver -> Numerics part
  if (info != 0)
  {
    cout << "OneStepNS computeOutput warning: output message from solver is equal to " << info << " => may have failed?" << endl;
    RuntimeException::selfThrow(" Non smooth problem, solver convergence failed");
    /*      if(info == 1)
    cout <<" reach max iterations number with solver " << solverName << endl;
    else if (info == 2)
    {
    if (solverName == "LexicoLemke" || solverName == "CPG" || solverName == "NLGS")
    RuntimeException::selfThrow(" negative diagonal term with solver "+solverName);
    else if (solverName == "QP" || solverName == "NSQP" )
    RuntimeException::selfThrow(" can not satisfy convergence criteria for solver "+solverName);
    else if (solverName == "Latin")
    RuntimeException::selfThrow(" Choleski factorisation failed with solver Latin");
    }
    else if (info == 3 && solverName == "CPG")
       cout << "pWp null in solver CPG" << endl;
    else if (info == 3 && solverName == "Latin")
    RuntimeException::selfThrow("Null diagonal term with solver Latin");
    else if (info == 5 && (solverName == "QP" || solverName == "NSQP"))
    RuntimeException::selfThrow("Length of working array insufficient in solver "+solverName);
    else
    RuntimeException::selfThrow("Unknown error type in solver "+ solverName);
    */
  }
}
