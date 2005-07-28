#include "OneStepNSProblem.h"
using namespace std;

// --- CONSTRUCTORS/DESTRUCTOR ---

// Default constructor
OneStepNSProblem::OneStepNSProblem():
  nspbType("none"), n(0), solver(""), strategy(NULL), onestepnspbxml(NULL)
{}

// xml constructor
OneStepNSProblem::OneStepNSProblem(OneStepNSProblemXML* osnspbxml, Strategy* newStrat):
  nspbType("none"), n(0), solver(""), strategy(newStrat), onestepnspbxml(osnspbxml)
{
  if (onestepnspbxml != NULL)
  {
    if (onestepnspbxml->hasN()) n = onestepnspbxml->getN();
    if (onestepnspbxml->hasSolver())
    {
      solver = onestepnspbxml->getSolver();

      string newSolvingMethod = onestepnspbxml->getSolverAlgorithmName();
      int MaxIter = onestepnspbxml->getSolverAlgorithmMaxIter();
      if (newSolvingMethod == "Lemke")
        fillSolvingMethod(newSolvingMethod, MaxIter);

      else if (newSolvingMethod == "Gsnl" || newSolvingMethod  == "Gcp")
      {
        double Tolerance = onestepnspbxml->getSolverAlgorithmTolerance();
        string NormType  = onestepnspbxml->getSolverAlgorithmNormType();
        fillSolvingMethod(newSolvingMethod, MaxIter, Tolerance, NormType);
      }

      else if (newSolvingMethod == "Latin")
      {
        double Tolerance = onestepnspbxml->getSolverAlgorithmTolerance();
        string NormType  = onestepnspbxml->getSolverAlgorithmNormType();
        double SearchDirection = onestepnspbxml->getSolverAlgorithmSearchDirection();
        fillSolvingMethod(newSolvingMethod, MaxIter, Tolerance, NormType, SearchDirection);
      }
      else RuntimeException::selfThrow("OneStepNSProblem::xml constructor, wrong solving method type");
    }
  }
  else RuntimeException::selfThrow("OneStepNSProblem::xml constructor, xml file=NULL");
  if (strategy != NULL)
  {
    interactionVector = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractions();
    ecVector = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getEqualityConstraints();
    // Default value for n  \warning: default size is the size of the first interaction in the vector
    n = interactionVector[0]->getNInteraction();
  }
  else cout << "OneStepNSPb xml-constructor - Warning: no strategy linked to OneStepPb" << endl;
}

// Constructor with given strategy and solving method parameters (optional)
OneStepNSProblem::OneStepNSProblem(Strategy * newStrat, const string& newSolver, const string& newSolvingMethod, const int& MaxIter,
                                   const double & Tolerance, const string & NormType, const double & SearchDirection):
  nspbType("none"), n(0), solver(""), strategy(newStrat), onestepnspbxml(NULL)
{
  if (strategy != NULL)
  {
    interactionVector = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractions();
    ecVector = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getEqualityConstraints();
    // Default value for n  \warning: default size is the size of the first interaction in the vector
    n = interactionVector[0]->getNInteraction();
    if (solver != "none")
    {
      solver = newSolver;
      if (newSolvingMethod == "Lemke")
        fillSolvingMethod(newSolvingMethod, MaxIter);

      else if (newSolvingMethod == "Gsnl" || newSolvingMethod  == "Gcp")
        fillSolvingMethod(newSolvingMethod, MaxIter, Tolerance, NormType);

      else if (newSolvingMethod == "Latin")
        fillSolvingMethod(newSolvingMethod, MaxIter, Tolerance, NormType, SearchDirection);

      else RuntimeException::selfThrow("OneStepNSProblem:: constructor from data, wrong solving method type");
    }
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem:: constructor from strategy, given strategy == NULL");
}

OneStepNSProblem::~OneStepNSProblem()
{}

Interaction* OneStepNSProblem::getInteractionPtr(const unsigned int& nb)
{
  if (nb >= interactionVector.size())
    RuntimeException::selfThrow("OneStepNSProblem::getInteractionPtr(const int& nb) - number greater than size of interaction vector");
  return interactionVector[nb];
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

  // if at least one relative degree is different from 0 or 1
  if (!topology->isTimeInvariant())
  {
    checkEffectiveOutput();
    // compute map of effective sizeOutput
    topology->computeEffectiveSizeOutput();
  }
  // else effectiveSizeOutput = SizeOutput
}

void OneStepNSProblem::checkEffectiveOutput()
{
  // get topology of the NonSmooth Dynamical System
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();
  // get time step
  double pasH = strategy->getTimeDiscretisationPtr()->getH();

  // loop over the interactions
  vector<Interaction*>::iterator it;
  for (it = interactionVector.begin(); it != interactionVector.end(); it++)
  {
    // get the output derivatives vector
    vector<SimpleVector *> yOld = (*it)->getYOld();

    // get relative degrees vector of this interaction
    vector<unsigned int> relativeDegree = topology->getRelativeDegrees(*it);

    unsigned int size = relativeDegree.size();

    // the prediction vector
    vector<SimpleVector *> yp;

    yp.resize(size, NULL);

    for (unsigned int i = 0; i < size ; i++)
      yp[i] = new SimpleVector(*yOld[i]);

    // \todo the way prediction is calculated should be defined by user elsewhere
    *(yp[0]) = *(yOld[0]) +  0.5 * pasH * *(yOld[1]) ;

    // loop from 0 to relative degree to find the first yp>0
    vector<unsigned int> indexMax;
    indexMax.resize(size, 0);
    for (unsigned int i = 0; i < size ; i++)
    {
      unsigned int j = 0;
      while ((*(yp[i]))(j) <= 0) j++;
      indexMax[i] = j;
    }
    topology->setIndexMax(*it, indexMax);
  }
}

void OneStepNSProblem::nextStep()
{
  vector<Interaction*>::iterator it;
  for (it = interactionVector.begin(); it != interactionVector.end(); it++)
    (*it)->swapInMemory();
  // get topology of the NonSmooth Dynamical System
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();
  // if relative degree different from 0 or 1
  if (!(topology->isTimeInvariant()))
  {
    checkEffectiveOutput();
    // compute map of effective sizeOutput
    topology->computeEffectiveSizeOutput();
  }
}

void OneStepNSProblem::updateInput()
{
  vector<Interaction*>::iterator it;
  double currentTime = strategy->getModelPtr()->getCurrentT();
  double pasH = strategy->getTimeDiscretisationPtr()->getH();
  for (it = interactionVector.begin(); it != interactionVector.end(); it++)
  {
    (*it)->update(currentTime, pasH);
    if (connectedInteractionMap.find(*it) != connectedInteractionMap.end())
      (*it)->getRelationPtr() -> computeInput(currentTime);
  }
}

void OneStepNSProblem::updateOutput()
{
  vector<Interaction*>::iterator it;
  double currentTime = strategy->getModelPtr()->getCurrentT();
  for (it = interactionVector.begin(); it != interactionVector.end(); it++)
    (*it)->getRelationPtr()->computeOutput(currentTime);
}


void OneStepNSProblem::checkInteraction()
{
  // --- check and update status of the interactions ---
  vector<Interaction*>::iterator it;
  double pasH = strategy->getTimeDiscretisationPtr()->getH();
  double time = strategy->getModelPtr()->getCurrentT();
  for (it = interactionVector.begin(); it != interactionVector.end(); it++)
    (*it)->check(time, pasH);
  updateConnectedInteractionMap();
}

void OneStepNSProblem::compute(const double& time)
{
  RuntimeException::selfThrow("OneStepNSProblem::compute - not yet implemented for problem type =" + getType());
}

void OneStepNSProblem::fillSolvingMethod(const string& newSolvingMethod, const int& MaxIter,
    const double & Tolerance, const string & NormType,
    const double & SearchDirection)
{
  if (solver ==  "LcpSolving")
    strcpy(solvingMethod.lcp.nom_method, newSolvingMethod.c_str());
  else if (solver == "RelayPrimalSolving")
    strcpy(solvingMethod.rp.nom_method, newSolvingMethod.c_str());
  else if (solver == "RelayDualSolving")
    strcpy(solvingMethod.rd.nom_method, newSolvingMethod.c_str());
  else if (solver == "ContactFrictionPrimalSolving")
    strcpy(solvingMethod.cfp.nom_method, newSolvingMethod.c_str());
  else if (solver == "ContactFrictionDualSolving")
    strcpy(solvingMethod.cfd.nom_method, newSolvingMethod.c_str());

  if (newSolvingMethod ==  "Lemke")
    setLemkeAlgorithm(solver, MaxIter);
  else if (newSolvingMethod == "Gsnl")
    setGsnlAlgorithm(solver, Tolerance, NormType, MaxIter);
  else if (newSolvingMethod == "Gcp")
    setGcpAlgorithm(solver, Tolerance, NormType, MaxIter);
  else if (newSolvingMethod == "Latin")
    setLatinAlgorithm(solver, Tolerance, NormType, MaxIter, SearchDirection);
  else
    RuntimeException::selfThrow("OneStepNSProblem::fillSolvingMethod, unknown method = " + newSolvingMethod);
}

void OneStepNSProblem::saveNSProblemToXML()
{
  IN("OneStepNSProblem::saveNSProblemToXML\n");
  if (onestepnspbxml != NULL)
  {
    onestepnspbxml->setN(n);
    vector<int> v;
    for (unsigned int i = 0; i < interactionVector.size(); i++)
      v.push_back(interactionVector[i]->getNumber());
    onestepnspbxml->setInteractionConcerned(v, allInteractionConcerned());

    /*
     * save of the solving method to XML
     */
    if (solver != "")
    {
      string methodName, normType;
      int maxIter;
      double tolerance, searchDirection;

      if (solver == OSNSP_LCPSOLVING)
      {
        methodName = solvingMethod.lcp.nom_method;
        normType = solvingMethod.lcp.normType;
        maxIter = solvingMethod.lcp.itermax;
        tolerance = solvingMethod.lcp.tol;
        searchDirection = solvingMethod.lcp.k_latin;
      }
      else if (solver == OSNSP_RPSOLVING)
      {
        methodName = solvingMethod.rp.nom_method;
        normType = solvingMethod.rp.normType;
        maxIter = solvingMethod.rp.itermax;
        tolerance = solvingMethod.rp.tol;
        searchDirection = solvingMethod.rp.k_latin;
      }
      else if (solver == OSNSP_RDSOLVING)
      {
        methodName = solvingMethod.rd.nom_method;
        normType = solvingMethod.rd.normType;
        maxIter = solvingMethod.rd.itermax;
        tolerance = solvingMethod.rd.tol;
        searchDirection = solvingMethod.rd.k_latin;
      }
      else if (solver == OSNSP_CFPSOLVING)
      {
        methodName = solvingMethod.cfp.nom_method;
        normType = solvingMethod.cfp.normType;
        maxIter = solvingMethod.cfp.itermax;
        tolerance = solvingMethod.cfp.tol;
        searchDirection = solvingMethod.cfp.k_latin;
      }
      else if (solver == OSNSP_CFDSOLVING)
      {
        methodName = solvingMethod.cfd.nom_method;
        normType = solvingMethod.cfd.normType;
        maxIter = solvingMethod.cfd.itermax;
        tolerance = solvingMethod.cfd.tol;
        searchDirection = solvingMethod.cfd.k_latin;
      }

      //      /*
      //       * according to the algorithm method used, only some attribute will be saved to XML
      //       * the other attributes will be put to the default value, so they won't be saved
      //       */
      //      if( methodName == OSNSP_LEMKE )
      //      {
      //        normType = DefaultAlgoNormType;
      //        maxIter = DefaultAlgoMaxIter;
      //        searchDirection = DefaultAlgoSearchDirection;
      //      }
      //      else if( methodName == OSNSP_GSNL )
      //      {
      //        searchDirection = DefaultAlgoSearchDirection;
      //      }
      //      else if( methodName == OSNSP_GCP )
      //      {
      //        searchDirection = DefaultAlgoSearchDirection;
      //      }
      //      else if( methodName == OSNSP_LATIN )
      //      {
      //        // all attributes are saved for Latin method
      //      }
      //      cout<<"OneStepNSproblem::saveNSProblemToXML => methodName == "<<methodName<<endl;

      onestepnspbxml->setSolver(solver,
                                methodName,
                                normType,
                                tolerance,
                                maxIter,
                                searchDirection);
    }
    else
    {
      cout << "# Warning : Can't save Solver tag, empty field" << endl;
    }
  }
  else RuntimeException::selfThrow("OneStepNSProblem::fillNSProblemWithNSProblemXML - OneStepNSProblemXML object not exists");
  OUT("OneStepNSProblem::saveNSProblemToXML\n");
}

bool OneStepNSProblem::allInteractionConcerned()
{
  bool res = false;

  if (interactionVector == getStrategyPtr()->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractions())
    res = true;

  return res;
}

void OneStepNSProblem::setLemkeAlgorithm(const string& meth,  const double& t)
{
  solver = meth;

  if (meth == OSNSP_LCPSOLVING)
  {
    strcpy(solvingMethod.lcp.nom_method, OSNSP_LEMKE.c_str());
    solvingMethod.lcp.tol = /*t*/ DefaultAlgoTolerance;
    strcpy(solvingMethod.lcp.normType, DefaultAlgoNormType.c_str());
    solvingMethod.lcp.itermax = /*DefaultAlgoMaxIter*/(int)t;
    solvingMethod.lcp.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_CFDSOLVING)
  {
    strcpy(solvingMethod.cfd.nom_method, OSNSP_LEMKE.c_str());
    solvingMethod.cfd.tol = /*t*/ DefaultAlgoTolerance;
    strcpy(solvingMethod.cfd.normType, DefaultAlgoNormType.c_str());
    solvingMethod.cfd.itermax = /*DefaultAlgoMaxIter*/(int)t;
    solvingMethod.cfd.k_latin = DefaultAlgoSearchDirection;
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem::setLemkeAlgorithm - solving method " + meth + " doesn't exists.");
}

void OneStepNSProblem::setGsnlAlgorithm(const string& meth,  const double& t,  const string& norm,  const int& iter)
{
  solver = meth;

  if (meth == OSNSP_LCPSOLVING)
  {
    strcpy(solvingMethod.lcp.nom_method, OSNSP_GSNL.c_str());
    solvingMethod.lcp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.lcp.normType, norm.c_str());
    solvingMethod.lcp.itermax = iter;
    solvingMethod.lcp.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_RPSOLVING)
  {
    strcpy(solvingMethod.rp.nom_method, OSNSP_GSNL.c_str());
    solvingMethod.rp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.rp.normType, norm.c_str());
    solvingMethod.rp.itermax = iter;
    solvingMethod.rp.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_RDSOLVING)
  {
    strcpy(solvingMethod.rd.nom_method, OSNSP_GSNL.c_str());
    solvingMethod.rd.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.rd.normType, norm.c_str());
    solvingMethod.rd.itermax = iter;
    solvingMethod.rd.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_CFPSOLVING)
  {
    strcpy(solvingMethod.cfp.nom_method, OSNSP_GSNL.c_str());
    solvingMethod.cfp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.cfp.normType, norm.c_str());
    solvingMethod.cfp.itermax = iter;
    solvingMethod.cfp.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_CFDSOLVING)
  {
    strcpy(solvingMethod.cfd.nom_method, OSNSP_GSNL.c_str());
    solvingMethod.cfd.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.cfd.normType, norm.c_str());
    solvingMethod.cfd.itermax = iter;
    solvingMethod.cfd.k_latin = DefaultAlgoSearchDirection;
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem::setGsnlAlgorithm - solving method " + meth + " doesn't exists.");
}

void OneStepNSProblem::setGcpAlgorithm(const string& meth,  const double& t,  const string& norm,  const int& iter)
{
  solver = meth;

  if (meth == OSNSP_LCPSOLVING)
  {
    strcpy(solvingMethod.lcp.nom_method, OSNSP_GCP.c_str());
    solvingMethod.lcp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.lcp.normType, norm.c_str());
    solvingMethod.lcp.itermax = iter;
    solvingMethod.lcp.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_RPSOLVING)
  {
    strcpy(solvingMethod.rp.nom_method, OSNSP_GCP.c_str());
    solvingMethod.rp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.rp.normType, norm.c_str());
    solvingMethod.rp.itermax = iter;
    solvingMethod.rp.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_RDSOLVING)
  {
    strcpy(solvingMethod.rd.nom_method, OSNSP_GCP.c_str());
    solvingMethod.rd.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.rd.normType, norm.c_str());
    solvingMethod.rd.itermax = iter;
    solvingMethod.rd.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_CFPSOLVING)
  {
    strcpy(solvingMethod.cfp.nom_method, OSNSP_GCP.c_str());
    solvingMethod.cfp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.cfp.normType, norm.c_str());
    solvingMethod.cfp.itermax = iter;
    solvingMethod.cfp.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_CFDSOLVING)
  {
    strcpy(solvingMethod.cfd.nom_method, OSNSP_GCP.c_str());
    solvingMethod.cfd.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.cfd.normType, norm.c_str());
    solvingMethod.cfd.itermax = iter;
    solvingMethod.cfd.k_latin = DefaultAlgoSearchDirection;
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem::setGcpAlgorithm - solving method " + meth + " doesn't exists.");
}

void OneStepNSProblem::setLatinAlgorithm(const string& meth, const double& t, const string& norm, const int& iter, const double& searchdirection)
{
  solver = meth;

  if (meth == OSNSP_LCPSOLVING)
  {
    strcpy(solvingMethod.lcp.nom_method, OSNSP_LATIN.c_str());
    solvingMethod.lcp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.lcp.normType, norm.c_str());
    solvingMethod.lcp.itermax = iter;
    solvingMethod.lcp.k_latin = searchdirection;
  }
  else if (meth == OSNSP_RPSOLVING)
  {
    strcpy(solvingMethod.rp.nom_method, OSNSP_LATIN.c_str());
    solvingMethod.rp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.rd.normType, norm.c_str());
    solvingMethod.rp.itermax = iter;
    solvingMethod.rp.k_latin = searchdirection;
  }
  else if (meth == OSNSP_RDSOLVING)
  {
    strcpy(solvingMethod.rd.nom_method, OSNSP_LATIN.c_str());
    solvingMethod.rd.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.rd.normType, norm.c_str());
    solvingMethod.rd.itermax = iter;
    solvingMethod.rd.k_latin = searchdirection;
  }
  else if (meth == OSNSP_CFPSOLVING)
  {
    strcpy(solvingMethod.cfp.nom_method, OSNSP_LATIN.c_str());
    solvingMethod.cfp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.cfp.normType, norm.c_str());
    solvingMethod.cfp.itermax = iter;
    solvingMethod.cfp.k_latin = searchdirection;
  }
  else if (meth == OSNSP_CFDSOLVING)
  {
    strcpy(solvingMethod.cfd.nom_method, OSNSP_LATIN.c_str());
    solvingMethod.cfd.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.cfd.normType, norm.c_str());
    solvingMethod.cfd.itermax = iter;
    solvingMethod.cfd.k_latin = searchdirection;
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem::setLatinAlgorithm - solving method " + meth + " doesn't exists.");
}

void OneStepNSProblem::updateConnectedInteractionMap()
{
  bool hasActiveConnection;
  vector<DynamicalSystem*> dsvect1, dsvect2;
  vector<int> status, status2;
  Connection *co;

  connectedInteractionMap.clear();

  // -- loop over the interactions --
  for (unsigned int i = 0; i < interactionVector.size(); i++)
  {
    // warning: only active interactions are put into the visibility table
    status = interactionVector[i]->getStatus();

    // -- loop over the status of the current interaction --
    for (unsigned int k = 0; k < status.size(); k++)
      if (status[k] == 1)
      {
        cout << "# interaction " << interactionVector[i]->getNumber() << " is active !" << endl;
        hasActiveConnection = false;
        for (unsigned int j = 0; j < interactionVector.size(); j++)
        {
          // we only put in the visibility table, the active interactions
          status2 = interactionVector[j]->getStatus();
          for (unsigned int l = 0; l < status2.size(); l++)
          {
            if (status2[l] == 1)
            {
              if (i != j)
              {
                hasActiveConnection = true;
                dsvect1 = interactionVector[i]->getDynamicalSystems();
                dsvect2 = interactionVector[j]->getDynamicalSystems();

                if ((dsvect1[0] == dsvect2[0]) || (dsvect1[0] == dsvect2[1]))
                {
                  // the interaction i and j are connected
                  co = new Connection();
                  co->status = 1;
                  co->originInteractionDSRank = 0;
                  if (dsvect1[0] == dsvect2[0]) co->connectedInteractionDSRank = 0;
                  else co->connectedInteractionDSRank = 1;
                  co->connected = interactionVector[j];

                  connectedInteractionMap[ interactionVector[i] ].push_back(co);
                }
                else if ((dsvect1[1] == dsvect2[0]) || (dsvect1[1] == dsvect2[1]))
                {
                  // the interaction i and j are connected
                  co = new Connection();
                  co->status = 1;
                  co->originInteractionDSRank = 1;
                  if (dsvect1[1] == dsvect2[0]) co->connectedInteractionDSRank = 0;
                  else co->connectedInteractionDSRank = 1;
                  co->connected = interactionVector[j];

                  connectedInteractionMap[ interactionVector[i] ].push_back(co);
                }
              }
            }
          }
        }
        if (!hasActiveConnection)
        {
          /*
           * if no active interaction is adjacent to the "interactionVector[i]"
           * active interaction, we add to the map "interactionVector[i]"
           * with an empty vector of Connection
           */
          //co = new Connection();
          connectedInteractionMap[ interactionVector[i] ].push_back(NULL);
        }
      }
  }
  //    displayConnectedInteractionMap();
  //    cout<<"#_#_  "<<connectedInteractionMap[ interactionVector[0] ].size()<<endl;
  //    cout<<"updateConnectedInteractionMap  <<press enter>>"<<endl;
  //    getchar();
}

void OneStepNSProblem::displayConnectedInteractionMap()
{
  map< Interaction*, vector<Connection*> >::iterator iter;

  //cout<<"#------OneStepNSProblem::displayConnectedInteractionMap-----------"<<endl;
  for (iter = connectedInteractionMap.begin(); iter != connectedInteractionMap.end(); iter++)
  {
    cout << "#-----------------" << endl;
    cout << "| Origin Interaction " << iter->first << endl;
    //cout<<"@ size of the second part of the map "<<iter->second.size()<<endl;
    if (iter->second[0] != NULL)
      for (unsigned int i = 0; i < iter->second.size(); i++)
      {
        cout << "|   + Connected Interaction " << iter->second[i]->connected << endl;
        cout << "|     status = " << iter->second[i]->status;
        cout << endl;
        cout << "|     originInteractionDSRank = " << iter->second[i]->originInteractionDSRank;
        cout << endl;
        cout << "|     connectedInteractionDSRank = " << iter->second[i]->connectedInteractionDSRank;
        cout << endl;
      }
    cout << "#-----------------" << endl;
  }
}

bool OneStepNSProblem::isOneStepNsProblemComplete()
{
  bool isComplete = true;

  if (nspbType != "LCP" || nspbType != "CFD" || nspbType != "QP" || nspbType != "Relay")
  {
    cout << "OneStepNSProblem is not complete: unknown problem type " << nspbType << endl;
    isComplete = false;
  }

  if (n == 0)
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
    vector< Interaction* >::iterator it;
    for (it = interactionVector.begin(); it != interactionVector.end(); it++)
      if (*it == NULL)
        cout << "OneStepNSProblem warning: an interaction points to NULL" << endl;
  }

  if (!(ecVector.size() > 0))
  {
    cout << "OneStepNSProblem warning: interaction vector is empty" << endl;
    isComplete = false;
  }
  else
  {
    vector< EqualityConstraint* >::iterator it;
    for (it = ecVector.begin(); it != ecVector.end(); it++)
      if (*it == NULL)
        cout << "OneStepNSProblem warning: an equalityConstraint of the problem points to NULL" << endl;
  }

  if (solver != "LcpSolving" || solver != "RelayPrimalSolving" || solver != "RelayDualSolving"
      || solver != "ContactFrictionPrimalSolving" || solver != "ContactFrictionDualSolving")
  {
    cout << "OneStepNSProblem is not complete: unknown solver type " << solver  << endl;
    isComplete = false;
  }

  if (!(connectedInteractionMap.size()) > 0)
  {
    cout << "OneStepNSProblem warning: connected interaction map is empty" << endl;
    isComplete = false;
  }
  else
  {
    map< Interaction*, std::vector<Connection*> >::iterator  mapit;
    mapit = connectedInteractionMap.find(NULL);
    if (mapit != connectedInteractionMap.end())
      cout << "OneStepNSProblem warning: a connected interaction points to NULL" << endl;
  }

  if (strategy == NULL)
  {
    cout << "OneStepNSProblem warning: no strategy linked with the problem" << endl;
    isComplete = false;
  }

  if (onestepnspbxml == NULL)
    cout << "OneStepNSProblem warning: xml linked-file == NULL" << endl;

  return isComplete;
}

