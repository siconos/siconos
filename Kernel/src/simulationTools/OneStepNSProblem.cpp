
#include "OneStepNSProblem.h"
using namespace std;

// --- CONSTRUCTORS/DESTRUCTOR ---

OneStepNSProblem::OneStepNSProblem():
  nspbType("none"), n(-1), solver(""), strategy(NULL), onestepnspbxml(NULL)
{}

OneStepNSProblem::OneStepNSProblem(OneStepNSProblemXML* osnspbxml, Strategy* newStrat):
  nspbType("none"), n(-1), solver(""), strategy(newStrat), onestepnspbxml(osnspbxml)
{
  if (onestepnspbxml != NULL)
  {
    if (onestepnspbxml->hasN()) n = onestepnspbxml->getN();
    if (onestepnspbxml->hasSolver())
    {
      solver = onestepnspbxml->getSolver();
      if (onestepnspbxml->hasSolverAlgorithm())
        fillSolvingMethod();
    }
  }
  else RuntimeException::selfThrow("OneStepNSProblem::xml constructor, xml file=NULL");
  if (strategy != NULL)
  {
    interactionVector = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractions();
    ecVector = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getEqualityConstraints();
    // Default value for n
    n = interactionVector[0]->getNInteraction();
  }
  else cout << "OneStepNSPb xml-constructor - Warning: no strategy linked to OneStepPb" << endl;
}

OneStepNSProblem::~OneStepNSProblem()
{}

Interaction* OneStepNSProblem::getInteractionPtr(const int& nb)
{
  if (nb < interactionVector.size())
    return interactionVector[nb];
  else RuntimeException::selfThrow("OneStepNSProblem::getInteraction(const int nb) - number greater than size of interaction vector");
}


void OneStepNSProblem::addInteraction(Interaction *interaction)
{
  interactionVector.push_back(interaction);
}

void OneStepNSProblem::updateState(void)
{
  vector<Interaction*>::iterator it;
  double currentTime = strategy->getModelPtr()->getCurrentT();
  double pasH = strategy->getTimeDiscretisationPtr()->getH();
  for (it = interactionVector.begin(); it != interactionVector.end(); ++it)
  {
    (*it)->update(currentTime, pasH);
    if (connectedInteractionMap.find(*it) != connectedInteractionMap.end())
      (*it)->getRelationPtr() -> computeInput(currentTime);
  }
}


void OneStepNSProblem::nextStep(void)
{
  vector<Interaction*>::iterator it;
  for (it = interactionVector.begin(); it != interactionVector.end(); ++it)
    (*it)->swapInMemory();
}

void OneStepNSProblem::checkInteraction()
{
  // --- check and update status of the interactions ---
  vector<Interaction*>::iterator it;
  double pasH = strategy->getTimeDiscretisationPtr()->getH();
  for (it = interactionVector.begin(); it != interactionVector.end(); ++it)
    (*it)->check(strategy->getModelPtr()->getCurrentT(), pasH);
  updateConnectedInteractionMap();
}

void OneStepNSProblem::formalize(const double& time)
{
  RuntimeException::selfThrow("OneStepNSProblem::formalize - not yet implemented for problem type =" + getType());
}


void OneStepNSProblem::compute(void)
{
  RuntimeException::selfThrow("OneStepNSProblem::compute - not yet implemented for problem type =" + getType());
}

void OneStepNSProblem::fillSolvingMethod()
{
  if (solver == OSNSP_LCPSOLVING)
    strcpy(solvingMethod.lcp.nom_method, onestepnspbxml->getSolverAlgorithmName().c_str());
  else if (solver == OSNSP_RPSOLVING)
    strcpy(solvingMethod.rp.nom_method, onestepnspbxml->getSolverAlgorithmName().c_str());
  else if (solver == OSNSP_RDSOLVING)
    strcpy(solvingMethod.rd.nom_method, onestepnspbxml->getSolverAlgorithmName().c_str());
  else if (solver == OSNSP_CFPSOLVING)
    strcpy(solvingMethod.cfp.nom_method, onestepnspbxml->getSolverAlgorithmName().c_str());
  else if (solver == OSNSP_CFDSOLVING)
    strcpy(solvingMethod.cfd.nom_method, onestepnspbxml->getSolverAlgorithmName().c_str());


  if (onestepnspbxml->getSolverAlgorithmName() == OSNSP_LEMKE)
  {
    setLemkeAlgorithm(solver,
                      /*onestepnspbxml->getSolverAlgorithmTolerance()*/
                      onestepnspbxml->getSolverAlgorithmMaxIter());
  }
  else if (onestepnspbxml->getSolverAlgorithmName() == OSNSP_GSNL)
  {
    setGsnlAlgorithm(solver,
                     onestepnspbxml->getSolverAlgorithmTolerance(),
                     onestepnspbxml->getSolverAlgorithmNormType(),
                     onestepnspbxml->getSolverAlgorithmMaxIter());
  }
  else if (onestepnspbxml->getSolverAlgorithmName() == OSNSP_GCP)
  {
    setGcpAlgorithm(solver,
                    onestepnspbxml->getSolverAlgorithmTolerance(),
                    onestepnspbxml->getSolverAlgorithmNormType(),
                    onestepnspbxml->getSolverAlgorithmMaxIter());
  }
  else if (onestepnspbxml->getSolverAlgorithmName() == OSNSP_LATIN)
  {
    setLatinAlgorithm(solver,
                      onestepnspbxml->getSolverAlgorithmTolerance(),
                      onestepnspbxml->getSolverAlgorithmNormType(),
                      onestepnspbxml->getSolverAlgorithmMaxIter(),
                      onestepnspbxml->getSolverAlgorithmSearchDirection());
  }
}

void OneStepNSProblem::saveNSProblemToXML()
{
  IN("OneStepNSProblem::saveNSProblemToXML\n");
  if (onestepnspbxml != NULL)
  {
    cout << " n===================== " << n << endl;
    onestepnspbxml->setN(n);

    vector<int> v;
    for (int i = 0; i < interactionVector.size(); i++)
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
    solvingMethod.lcp.itermax = /*DefaultAlgoMaxIter*/ t;
    solvingMethod.lcp.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_CFDSOLVING)
  {
    strcpy(solvingMethod.cfd.nom_method, OSNSP_LEMKE.c_str());
    solvingMethod.cfd.tol = /*t*/ DefaultAlgoTolerance;
    strcpy(solvingMethod.cfd.normType, DefaultAlgoNormType.c_str());
    solvingMethod.cfd.itermax = /*DefaultAlgoMaxIter*/ t;
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
    RuntimeException::selfThrow("OneStepNSProblem::setLemkeAlgorithm - solving method " + meth + " doesn't exists.");
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
    RuntimeException::selfThrow("OneStepNSProblem::setLemkeAlgorithm - solving method " + meth + " doesn't exists.");
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
    RuntimeException::selfThrow("OneStepNSProblem::setLemkeAlgorithm - solving method " + meth + " doesn't exists.");
}

void OneStepNSProblem::updateConnectedInteractionMap()
{
  bool hasActiveConnection;
  vector<DynamicalSystem*> dsvect1, dsvect2;
  vector<int> status, status2;
  Connection *co;

  connectedInteractionMap.clear();

  for (int i = 0; i < interactionVector.size(); i++)
  {
    //    for(int k=0; k<connectedInteractionMap[ interactionVector[i] ].size(); k++)
    //      delete connectedInteractionMap[ interactionVector[i] ][k];
    //    connectedInteractionMap[ interactionVector[i] ].clear();

    // we only put in the visibility table, the active interactions
    status = interactionVector[i]->getStatus();
    for (int k = 0; k < status.size(); k++)
      if (status[k] == 1)
      {
        cout << "# interaction " << i << " is active !" << endl;
        hasActiveConnection = false;
        for (int j = 0; j < interactionVector.size(); j++)
        {
          // we only put in the visibility table, the active interactions
          status2 = interactionVector[j]->getStatus();
          for (int l = 0; l < status2.size(); l++)
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
  for (iter = connectedInteractionMap.begin(); iter != connectedInteractionMap.end(); ++iter)
  {
    cout << "#-----------------" << endl;
    cout << "| Origin Interaction " << iter->first << endl;
    //cout<<"@ size of the second part of the map "<<iter->second.size()<<endl;
    if (iter->second[0] != NULL)
      for (int i = 0; i < iter->second.size(); i++)
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
