
#include "OneStepNSProblem.h"
#include "check.h"

OneStepNSProblem::OneStepNSProblem()
{
  this->onestepnspbxml = NULL;
  this->n = -1;
  this->solver = "";
}

OneStepNSProblem::OneStepNSProblem(OneStepNSProblemXML* osnspbxml)
{
  this->onestepnspbxml = osnspbxml;
  this->solver = "";
}

OneStepNSProblem::~OneStepNSProblem()
{}


Interaction* OneStepNSProblem::getInteraction(const int nb)
{
  if (nb < this->interactionVector.size())
    return this->interactionVector[nb];
  else RuntimeException::selfThrow("OneStepNSProblem::getInteraction(const int nb) - number greater than size of interaction vector");
}


void OneStepNSProblem::addInteraction(Interaction *interaction)
{
  this->interactionVector.push_back(interaction);
  //  cout<<"  - one of the Interactions Linked to the OneStepNSProblem ("<<this->nspbType<<") has the number "<<this->interactionVector[interactionVector.size()-1]->getNumber()<<" and his id is "<<this->interactionVector[interactionVector.size()-1]->getId()<<endl;
  //  cout<<" % the relation of the linked Interaction is : "<<this->interactionVector[interactionVector.size()-1]->getRelation()->getType()<<endl;
  //  cout<<" % this relation is linked to his XML object : "<<this->interactionVector[interactionVector.size()-1]->getRelation()->getRelationXML()<<endl;
}


//////////////////////////////

void OneStepNSProblem::initialize(void)
{
  // predict all the relations to see which ones have an effect
  Interaction *inter;
  vector<Interaction*>::iterator it;
  for (it = this->interactionVector.begin(); it != this->interactionVector.end(); ++it)
  {
    inter = *it;
    inter->initialize();
  }
}



void OneStepNSProblem::updateState(void)
{
  // predict all the relations to see which ones have an effect
  Interaction *inter;
  vector<Interaction*>::iterator it;
  double currentTime = this->strategy->getModel()->getCurrentT();

  for (it = this->interactionVector.begin(); it != this->interactionVector.end(); ++it)
  {
    inter = *it;
    inter->update(currentTime);

    if (this->connectedInteractionMap.find(inter) != this->connectedInteractionMap.end())
      inter->getRelation()->computeInput(currentTime);
  }
}


void OneStepNSProblem::nextStep(void)
{
  // predict all the relations to see which ones have an effect
  Interaction *inter;
  vector<Interaction*>::iterator it;
  for (it = this->interactionVector.begin(); it != this->interactionVector.end(); ++it)
  {
    inter = *it;
    inter->swapInMemory();
  }
}

void OneStepNSProblem::checkInteraction(void)
{
  // predict all the relations to see which ones have an effect
  Interaction *inter;
  vector<Interaction*>::iterator it;
  for (it = this->interactionVector.begin(); it != this->interactionVector.end(); ++it)
  {
    inter = *it;
    inter->check(this->strategy->getModel()->getCurrentT());
  }
}

void OneStepNSProblem::formalize(double time)
{
  //  Interaction *inter;
  //  vector<Interaction*>::iterator it;
  //  for( it = this->interactionVector.begin(); it != this->interactionVector.end(); ++it )
  //  {
  //    inter = *it;
  //    if( inter->isActive() )
  //    {
  //      /*
  //      inter->computeA();
  //      inter->computeB();
  //      */
  //    }
  //  }
}

void OneStepNSProblem::compute(void)
{
  // to do with NUMERICS
}

void OneStepNSProblem::fillNSProblemWithNSProblemXML()
{
  IN("OneStepNSProblem::fillNSProblemWithNSProblemXML\n");
  if (this->onestepnspbxml != NULL)
  {
    if (this->onestepnspbxml->hasN())
      this->n = this->onestepnspbxml->getN();
    else
    {
      //this->n =
      // \todo :  calculate the number "n" of the OneStepNSProblem
    }

    if (this->onestepnspbxml->hasSolver())
    {
      this->solver = this->onestepnspbxml->getSolver();
      if (this->onestepnspbxml->hasSolverAlgorithm())
      {
        // filling of the right method imported from Numerics
        this->fillSolvingMethod();
      }
    }
  }
  else RuntimeException::selfThrow("OneStepNSProblem::fillNSProblemWithNSProblemXML - OneStepNSProblemXML object not exists");
  OUT("OneStepNSProblem::fillNSProblemWithNSProblemXML\n");
}

void OneStepNSProblem::fillSolvingMethod()
{
  if (this->solver == OSNSP_LCPSOLVING)
    strcpy(this->solvingMethod.lcp.nom_method, this->onestepnspbxml->getSolverAlgorithmName().c_str());
  else if (this->solver == OSNSP_RPSOLVING)
    strcpy(this->solvingMethod.rp.nom_method, this->onestepnspbxml->getSolverAlgorithmName().c_str());
  else if (this->solver == OSNSP_RDSOLVING)
    strcpy(this->solvingMethod.rd.nom_method, this->onestepnspbxml->getSolverAlgorithmName().c_str());
  else if (this->solver == OSNSP_CFPSOLVING)
    strcpy(this->solvingMethod.cfp.nom_method, this->onestepnspbxml->getSolverAlgorithmName().c_str());
  else if (this->solver == OSNSP_CFDSOLVING)
    strcpy(this->solvingMethod.cfd.nom_method, this->onestepnspbxml->getSolverAlgorithmName().c_str());


  if (this->onestepnspbxml->getSolverAlgorithmName() == OSNSP_LEMKE)
  {
    this->setLemkeAlgorithm(this->solver,
                            /*this->onestepnspbxml->getSolverAlgorithmTolerance()*/
                            this->onestepnspbxml->getSolverAlgorithmMaxIter());
  }
  else if (this->onestepnspbxml->getSolverAlgorithmName() == OSNSP_GSNL)
  {
    this->setGsnlAlgorithm(this->solver,
                           this->onestepnspbxml->getSolverAlgorithmTolerance(),
                           this->onestepnspbxml->getSolverAlgorithmNormType(),
                           this->onestepnspbxml->getSolverAlgorithmMaxIter());
  }
  else if (this->onestepnspbxml->getSolverAlgorithmName() == OSNSP_GCP)
  {
    this->setGcpAlgorithm(this->solver,
                          this->onestepnspbxml->getSolverAlgorithmTolerance(),
                          this->onestepnspbxml->getSolverAlgorithmNormType(),
                          this->onestepnspbxml->getSolverAlgorithmMaxIter());
  }
  else if (this->onestepnspbxml->getSolverAlgorithmName() == OSNSP_LATIN)
  {
    this->setLatinAlgorithm(this->solver,
                            this->onestepnspbxml->getSolverAlgorithmTolerance(),
                            this->onestepnspbxml->getSolverAlgorithmNormType(),
                            this->onestepnspbxml->getSolverAlgorithmMaxIter(),
                            this->onestepnspbxml->getSolverAlgorithmSearchDirection());
  }
}

void OneStepNSProblem::saveNSProblemToXML()
{
  IN("OneStepNSProblem::saveNSProblemToXML\n");
  if (this->onestepnspbxml != NULL)
  {
    this->onestepnspbxml->setN(this->n);

    vector<int> v;
    for (int i = 0; i < this->interactionVector.size(); i++)
      v.push_back(this->interactionVector[i]->getNumber());
    this->onestepnspbxml->setInteractionConcerned(v, this->allInteractionConcerned());

    /*
     * save of the solving method to XML
     */
    if (this->solver != "")
    {
      string methodName, normType;
      int maxIter;
      double tolerance, searchDirection;

      if (this->solver == OSNSP_LCPSOLVING)
      {
        methodName = this->solvingMethod.lcp.nom_method;
        normType = this->solvingMethod.lcp.normType;
        maxIter = this->solvingMethod.lcp.itermax;
        tolerance = this->solvingMethod.lcp.tol;
        searchDirection = this->solvingMethod.lcp.k_latin;
      }
      else if (this->solver == OSNSP_RPSOLVING)
      {
        methodName = this->solvingMethod.rp.nom_method;
        normType = this->solvingMethod.rp.normType;
        maxIter = this->solvingMethod.rp.itermax;
        tolerance = this->solvingMethod.rp.tol;
        searchDirection = this->solvingMethod.rp.k_latin;
      }
      else if (this->solver == OSNSP_RDSOLVING)
      {
        methodName = this->solvingMethod.rd.nom_method;
        normType = this->solvingMethod.rd.normType;
        maxIter = this->solvingMethod.rd.itermax;
        tolerance = this->solvingMethod.rd.tol;
        searchDirection = this->solvingMethod.rd.k_latin;
      }
      else if (this->solver == OSNSP_CFPSOLVING)
      {
        methodName = this->solvingMethod.cfp.nom_method;
        normType = this->solvingMethod.cfp.normType;
        maxIter = this->solvingMethod.cfp.itermax;
        tolerance = this->solvingMethod.cfp.tol;
        searchDirection = this->solvingMethod.cfp.k_latin;
      }
      else if (this->solver == OSNSP_CFDSOLVING)
      {
        methodName = this->solvingMethod.cfd.nom_method;
        normType = this->solvingMethod.cfd.normType;
        maxIter = this->solvingMethod.cfd.itermax;
        tolerance = this->solvingMethod.cfd.tol;
        searchDirection = this->solvingMethod.cfd.k_latin;
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

      this->onestepnspbxml->setSolver(this->solver,
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


void OneStepNSProblem::fillInteractionVector()
{
  this->interactionVector = this->strategy->getModel()->getNonSmoothDynamicalSystem()->getInteractions();
  this->ecVector = this->strategy->getModel()->getNonSmoothDynamicalSystem()->getEqualityConstraints();
}

void OneStepNSProblem::init()
{
  if ((this->onestepnspbxml == NULL) || (! this->onestepnspbxml->hasN()))
  {
    /*
     * n is initialized with the value nInteraction of the first interaction of the NonSmoothDynamicalSystem
     */
    this->n = interactionVector[0]->getNInteraction();
  }
}

bool OneStepNSProblem::allInteractionConcerned()
{
  bool res = false;

  if (this->interactionVector == this->getStrategy()->getModel()->getNonSmoothDynamicalSystem()->getInteractions())
    res = true;

  return res;
}

void OneStepNSProblem::setLemkeAlgorithm(string meth, double t)
{
  this->solver = meth;

  if (meth == OSNSP_LCPSOLVING)
  {
    strcpy(this->solvingMethod.lcp.nom_method, OSNSP_LEMKE.c_str());
    this->solvingMethod.lcp.tol = /*t*/ DefaultAlgoTolerance;
    strcpy(this->solvingMethod.lcp.normType, DefaultAlgoNormType.c_str());
    this->solvingMethod.lcp.itermax = /*DefaultAlgoMaxIter*/ t;
    this->solvingMethod.lcp.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_CFDSOLVING)
  {
    strcpy(this->solvingMethod.cfd.nom_method, OSNSP_LEMKE.c_str());
    this->solvingMethod.cfd.tol = /*t*/ DefaultAlgoTolerance;
    strcpy(this->solvingMethod.cfd.normType, DefaultAlgoNormType.c_str());
    this->solvingMethod.cfd.itermax = /*DefaultAlgoMaxIter*/ t;
    this->solvingMethod.cfd.k_latin = DefaultAlgoSearchDirection;
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem::setLemkeAlgorithm - solving method " + meth + " doesn't exists.");
}

void OneStepNSProblem::setGsnlAlgorithm(string meth, double t, string norm, int iter)
{
  this->solver = meth;

  if (meth == OSNSP_LCPSOLVING)
  {
    strcpy(this->solvingMethod.lcp.nom_method, OSNSP_GSNL.c_str());
    this->solvingMethod.lcp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(this->solvingMethod.lcp.normType, norm.c_str());
    this->solvingMethod.lcp.itermax = iter;
    this->solvingMethod.lcp.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_RPSOLVING)
  {
    strcpy(this->solvingMethod.rp.nom_method, OSNSP_GSNL.c_str());
    this->solvingMethod.rp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(this->solvingMethod.rp.normType, norm.c_str());
    this->solvingMethod.rp.itermax = iter;
    this->solvingMethod.rp.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_RDSOLVING)
  {
    strcpy(this->solvingMethod.rd.nom_method, OSNSP_GSNL.c_str());
    this->solvingMethod.rd.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(this->solvingMethod.rd.normType, norm.c_str());
    this->solvingMethod.rd.itermax = iter;
    this->solvingMethod.rd.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_CFPSOLVING)
  {
    strcpy(this->solvingMethod.cfp.nom_method, OSNSP_GSNL.c_str());
    this->solvingMethod.cfp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(this->solvingMethod.cfp.normType, norm.c_str());
    this->solvingMethod.cfp.itermax = iter;
    this->solvingMethod.cfp.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_CFDSOLVING)
  {
    strcpy(this->solvingMethod.cfd.nom_method, OSNSP_GSNL.c_str());
    this->solvingMethod.cfd.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(this->solvingMethod.cfd.normType, norm.c_str());
    this->solvingMethod.cfd.itermax = iter;
    this->solvingMethod.cfd.k_latin = DefaultAlgoSearchDirection;
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem::setLemkeAlgorithm - solving method " + meth + " doesn't exists.");
}

void OneStepNSProblem::setGcpAlgorithm(string meth, double t, string norm, int iter)
{
  this->solver = meth;

  if (meth == OSNSP_LCPSOLVING)
  {
    strcpy(this->solvingMethod.lcp.nom_method, OSNSP_GCP.c_str());
    this->solvingMethod.lcp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(this->solvingMethod.lcp.normType, norm.c_str());
    this->solvingMethod.lcp.itermax = iter;
    this->solvingMethod.lcp.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_RPSOLVING)
  {
    strcpy(this->solvingMethod.rp.nom_method, OSNSP_GCP.c_str());
    this->solvingMethod.rp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(this->solvingMethod.rp.normType, norm.c_str());
    this->solvingMethod.rp.itermax = iter;
    this->solvingMethod.rp.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_RDSOLVING)
  {
    strcpy(this->solvingMethod.rd.nom_method, OSNSP_GCP.c_str());
    this->solvingMethod.rd.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(this->solvingMethod.rd.normType, norm.c_str());
    this->solvingMethod.rd.itermax = iter;
    this->solvingMethod.rd.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_CFPSOLVING)
  {
    strcpy(this->solvingMethod.cfp.nom_method, OSNSP_GCP.c_str());
    this->solvingMethod.cfp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(this->solvingMethod.cfp.normType, norm.c_str());
    this->solvingMethod.cfp.itermax = iter;
    this->solvingMethod.cfp.k_latin = DefaultAlgoSearchDirection;
  }
  else if (meth == OSNSP_CFDSOLVING)
  {
    strcpy(this->solvingMethod.cfd.nom_method, OSNSP_GCP.c_str());
    this->solvingMethod.cfd.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(this->solvingMethod.cfd.normType, norm.c_str());
    this->solvingMethod.cfd.itermax = iter;
    this->solvingMethod.cfd.k_latin = DefaultAlgoSearchDirection;
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem::setLemkeAlgorithm - solving method " + meth + " doesn't exists.");
}

void OneStepNSProblem::setLatinAlgorithm(string meth, double t, string norm, int iter, double searchdirection)
{
  this->solver = meth;

  if (meth == OSNSP_LCPSOLVING)
  {
    strcpy(this->solvingMethod.lcp.nom_method, OSNSP_LATIN.c_str());
    this->solvingMethod.lcp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(this->solvingMethod.lcp.normType, norm.c_str());
    this->solvingMethod.lcp.itermax = iter;
    this->solvingMethod.lcp.k_latin = searchdirection;
  }
  else if (meth == OSNSP_RPSOLVING)
  {
    strcpy(this->solvingMethod.rp.nom_method, OSNSP_LATIN.c_str());
    this->solvingMethod.rp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(this->solvingMethod.rd.normType, norm.c_str());
    this->solvingMethod.rp.itermax = iter;
    this->solvingMethod.rp.k_latin = searchdirection;
  }
  else if (meth == OSNSP_RDSOLVING)
  {
    strcpy(this->solvingMethod.rd.nom_method, OSNSP_LATIN.c_str());
    this->solvingMethod.rd.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(this->solvingMethod.rd.normType, norm.c_str());
    this->solvingMethod.rd.itermax = iter;
    this->solvingMethod.rd.k_latin = searchdirection;
  }
  else if (meth == OSNSP_CFPSOLVING)
  {
    strcpy(this->solvingMethod.cfp.nom_method, OSNSP_LATIN.c_str());
    this->solvingMethod.cfp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(this->solvingMethod.cfp.normType, norm.c_str());
    this->solvingMethod.cfp.itermax = iter;
    this->solvingMethod.cfp.k_latin = searchdirection;
  }
  else if (meth == OSNSP_CFDSOLVING)
  {
    strcpy(this->solvingMethod.cfd.nom_method, OSNSP_LATIN.c_str());
    this->solvingMethod.cfd.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(this->solvingMethod.cfd.normType, norm.c_str());
    this->solvingMethod.cfd.itermax = iter;
    this->solvingMethod.cfd.k_latin = searchdirection;
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


  this->connectedInteractionMap.clear();

  for (int i = 0; i < this->interactionVector.size(); i++)
  {
    //    for(int k=0; k<this->connectedInteractionMap[ interactionVector[i] ].size(); k++)
    //      delete this->connectedInteractionMap[ interactionVector[i] ][k];
    //    this->connectedInteractionMap[ interactionVector[i] ].clear();

    // we only put in the visibility table, the active interactions
    status = this->interactionVector[i]->getStatus();
    for (int k = 0; k < status.size(); k++)
      if (status[k] == 1)
      {
        cout << "# interaction " << i << " is active !" << endl;
        hasActiveConnection = false;
        for (int j = 0; j < this->interactionVector.size(); j++)
        {
          // we only put in the visibility table, the active interactions
          status2 = this->interactionVector[j]->getStatus();
          for (int l = 0; l < status2.size(); l++)
          {
            if (status2[l] == 1)
            {
              if (i != j)
              {
                hasActiveConnection = true;
                dsvect1 = this->interactionVector[i]->getDynamicalSystems();
                dsvect2 = this->interactionVector[j]->getDynamicalSystems();

                if ((dsvect1[0] == dsvect2[0]) || (dsvect1[0] == dsvect2[1]))
                {
                  // the interaction i and j are connected
                  co = new Connection();
                  co->status = 1;
                  co->originInteractionDSRank = 0;
                  if (dsvect1[0] == dsvect2[0]) co->connectedInteractionDSRank = 0;
                  else co->connectedInteractionDSRank = 1;
                  co->connected = interactionVector[j];

                  this->connectedInteractionMap[ interactionVector[i] ].push_back(co);
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

                  this->connectedInteractionMap[ interactionVector[i] ].push_back(co);
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
          this->connectedInteractionMap[ interactionVector[i] ].push_back(NULL);
        }
      }
  }
  //    this->displayConnectedInteractionMap();
  //    cout<<"#_#_  "<<this->connectedInteractionMap[ interactionVector[0] ].size()<<endl;
  //    cout<<"updateConnectedInteractionMap  <<press enter>>"<<endl;
  //    getchar();
}

void OneStepNSProblem::displayConnectedInteractionMap()
{
  map< Interaction*, vector<Connection*> >::iterator iter;

  //cout<<"#------OneStepNSProblem::displayConnectedInteractionMap-----------"<<endl;
  for (iter = this->connectedInteractionMap.begin(); iter != this->connectedInteractionMap.end(); ++iter)
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
