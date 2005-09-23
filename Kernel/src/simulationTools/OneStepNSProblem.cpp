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
      if (newSolvingMethod == "Lemke" || newSolvingMethod == "LexicoLemke")
        fillSolvingMethod(newSolvingMethod, MaxIter);

      else if (newSolvingMethod == "QP" || newSolvingMethod  == "NSQP")
      {
        double Tolerance = onestepnspbxml->getSolverAlgorithmTolerance();
        fillSolvingMethod(newSolvingMethod, 0, Tolerance);
      }
      else if (newSolvingMethod == "NLGS" || newSolvingMethod  == "CPG")
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
      if (newSolvingMethod == "Lemke" || newSolvingMethod == "LexicoLemke")
        fillSolvingMethod(newSolvingMethod, MaxIter);

      else if (newSolvingMethod == "QP" || newSolvingMethod  == "NSQP")
        fillSolvingMethod(newSolvingMethod, 0, Tolerance);

      else if (newSolvingMethod == "NLGS" || newSolvingMethod  == "CPG")
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
}

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
    unsigned int size = relativeDegree.size(); // this size corresponds to the interaction size, ie the number of relations

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
      indexMax.resize(size, 0);
      for (j = 0; j < size; j++)
      {
        for (i = 0; i < sizeYp; i++)
        {
          if ((*(yp[i]))(j) <= 0)
            indexMax[j]++;
          else
            break;
        }
      }
      topology->setIndexMax(*it, indexMax);

      for (i = 0; i < sizeYp ; i++)
        delete yp[i];

      // --- effective indexes ---

      // compute sizeOutput for the current interaction
      sizeOutput = topology->computeEffectiveSizeOutput(*it);

      vector<unsigned int> effectiveIndexes, indexMin;
      indexMin = topology->getIndexMin(*it);

      effectiveIndexes.resize(sizeOutput);

      for (k = 0; k < sizeOutput; k++)
      {
        for (j = 0; j < size; j++)
        {
          for (i = 0; i < (indexMax[j] - indexMin[j]); i++)
          {
            effectiveIndexes[k] = i + j * (relativeDegree[j] - indexMin[j]);
          }
        }
      }

      topology->setEffectiveIndexes(*it, effectiveIndexes);
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

void OneStepNSProblem::fillSolvingMethod(const string& newSolvingMethod, const int& MaxIter,
    const double & Tolerance, const string & NormType,
    const double & SearchDirection)
{
  if (solver ==  "LcpSolving")
    strcpy(solvingMethod.lcp.name, newSolvingMethod.c_str());
  else if (solver == "PrimalRelaySolving")
    strcpy(solvingMethod.pr.name, newSolvingMethod.c_str());
  else if (solver == "DualRelaySolving")
    strcpy(solvingMethod.dr.name, newSolvingMethod.c_str());
  else if (solver == "PrimalFrictionContact2DSolving")
    strcpy(solvingMethod.pfc_2D.name, newSolvingMethod.c_str());
  else if (solver == "DualFrictionContact2DSolving")
    strcpy(solvingMethod.dfc_2D.name, newSolvingMethod.c_str());

  if (newSolvingMethod ==  "Lemke")
    setLemkeAlgorithm(solver, MaxIter);
  else if (newSolvingMethod == "LexicoLemke")
    setLexicoLemkeAlgorithm(solver, MaxIter);
  else if (newSolvingMethod == "NLGS")
    setNLGSAlgorithm(solver, Tolerance, MaxIter);
  else if (newSolvingMethod == "QP")
    setQPAlgorithm(solver, Tolerance);
  else if (newSolvingMethod == "NSQP")
    setNSQPAlgorithm(solver, Tolerance);
  else if (newSolvingMethod == "CPG")
    setCPGAlgorithm(solver, Tolerance, MaxIter);
  else if (newSolvingMethod == "Latin")
    setLatinAlgorithm(solver, Tolerance, MaxIter, NormType, SearchDirection);
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
        methodName = solvingMethod.lcp.name;
        normType = solvingMethod.lcp.normType;
        maxIter = solvingMethod.lcp.itermax;
        tolerance = solvingMethod.lcp.tol;
        searchDirection = solvingMethod.lcp.k_latin;
      }
      else if (solver == OSNSP_PRSOLVING)
      {
        methodName = solvingMethod.pr.name;
        normType = solvingMethod.pr.normType;
        maxIter = solvingMethod.pr.itermax;
        tolerance = solvingMethod.pr.tol;
        searchDirection = solvingMethod.pr.k_latin;
      }
      else if (solver == OSNSP_DRSOLVING)
      {
        methodName = solvingMethod.dr.name;
        normType = solvingMethod.dr.normType;
        maxIter = solvingMethod.dr.itermax;
        tolerance = solvingMethod.dr.tol;
        searchDirection = solvingMethod.dr.k_latin;
      }
      else if (solver == OSNSP_PFC_2DSOLVING)
      {
        methodName = solvingMethod.pfc_2D.name;
        normType = solvingMethod.pfc_2D.normType;
        maxIter = solvingMethod.pfc_2D.itermax;
        tolerance = solvingMethod.pfc_2D.tol;
        searchDirection = solvingMethod.pfc_2D.k_latin;
      }
      else if (solver == OSNSP_DFC_2DSOLVING)
      {
        methodName = solvingMethod.dfc_2D.name;
        normType = solvingMethod.dfc_2D.normType;
        maxIter = solvingMethod.dfc_2D.itermax;
        tolerance = solvingMethod.dfc_2D.tol;
        searchDirection = solvingMethod.dfc_2D.k_latin;
      }

      onestepnspbxml->setSolver(solver, methodName, normType, tolerance, maxIter, searchDirection);
    }
    else
      cout << "# Warning : Can't save Solver tag, empty field" << endl;
  }
  else RuntimeException::selfThrow("OneStepNSProblem::saveNSProblemToXML - OneStepNSProblemXML object not exists");
  OUT("OneStepNSProblem::saveNSProblemToXML\n");
}

bool OneStepNSProblem::allInteractionConcerned()
{
  bool res = false;

  if (interactionVector == getStrategyPtr()->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractions())
    res = true;

  return res;
}

void OneStepNSProblem::setLemkeAlgorithm(const string& meth,  const unsigned int& iter)
{
  solver = meth;

  if (meth == OSNSP_LCPSOLVING)
  {
    strcpy(solvingMethod.lcp.name, OSNSP_LEMKE.c_str());
    solvingMethod.lcp.itermax = iter;
  }
  else if (meth == OSNSP_DFC_2DSOLVING)
  {
    strcpy(solvingMethod.dfc_2D.name, OSNSP_LEMKE.c_str());
    solvingMethod.dfc_2D.itermax = iter;
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem::setLemkeAlgorithm - solving method " + meth + " doesn't exists.");
}

void OneStepNSProblem::setLexicoLemkeAlgorithm(const string& meth,  const unsigned int& iter)
{
  solver = meth;

  if (meth == OSNSP_LCPSOLVING)
  {
    strcpy(solvingMethod.lcp.name, OSNSP_LEXICOLEMKE.c_str());
    solvingMethod.lcp.itermax = iter;
  }
  else if (meth == OSNSP_DFC_2DSOLVING)
  {
    strcpy(solvingMethod.dfc_2D.name, OSNSP_LEXICOLEMKE.c_str());
    solvingMethod.dfc_2D.itermax = iter;
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem::setLemkeAlgorithm - solving method " + meth + " doesn't exists.");
}

void OneStepNSProblem::setNLGSAlgorithm(const string& meth,  const double& tolerance,  const int& iter,  const string& norm)
{
  solver = meth;

  if (meth == OSNSP_LCPSOLVING)
  {
    strcpy(solvingMethod.lcp.name, OSNSP_NLGS.c_str());
    solvingMethod.lcp.tol = tolerance;
    /*##### normType is not yet implemented in Numerics  #####*/
    //strcpy( solvingMethod.lcp.normType, norm.c_str() );
    solvingMethod.lcp.itermax = iter;
  }
  else if (meth == OSNSP_PRSOLVING)
  {
    strcpy(solvingMethod.pr.name, OSNSP_NLGS.c_str());
    solvingMethod.pr.tol = tolerance;
    /*##### normType is not yet implemented in Numerics  #####*/
    //strcpy( solvingMethod.pr.normType, norm.c_str() );
    solvingMethod.pr.itermax = iter;
  }
  else if (meth == OSNSP_DRSOLVING)
  {
    strcpy(solvingMethod.dr.name, OSNSP_NLGS.c_str());
    solvingMethod.dr.tol = tolerance;
    /*##### normType is not yet implemented in Numerics  #####*/
    //strcpy( solvingMethod.dr.normType, norm.c_str() );
    solvingMethod.dr.itermax = iter;
  }
  else if (meth == OSNSP_PFC_2DSOLVING)
  {
    strcpy(solvingMethod.pfc_2D.name, OSNSP_NLGS.c_str());
    solvingMethod.pfc_2D.tol = tolerance;
    /*##### normType is not yet implemented in Numerics  #####*/
    //strcpy( solvingMethod.pfc_2D.normType, norm.c_str() );
    solvingMethod.pfc_2D.itermax = iter;
  }
  else if (meth == OSNSP_DFC_2DSOLVING)
  {
    strcpy(solvingMethod.dfc_2D.name, OSNSP_NLGS.c_str());
    solvingMethod.dfc_2D.tol = tolerance;
    /*##### normType is not yet implemented in Numerics  #####*/
    //strcpy( solvingMethod.dfc_2D.normType, norm.c_str() );
    solvingMethod.dfc_2D.itermax = iter;
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem::setNLGSAlgorithm - solving method " + meth + " doesn't exists.");
}

void OneStepNSProblem::setQPAlgorithm(const string& meth,  const double& tolerance)
{
  solver = meth;

  if (meth == OSNSP_LCPSOLVING)
  {
    strcpy(solvingMethod.lcp.name, OSNSP_QP.c_str());
    solvingMethod.lcp.tol = tolerance;
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem::setQPAlgorithm - solving method " + meth + " doesn't exists.");
}

void OneStepNSProblem::setNSQPAlgorithm(const string& meth,  const double& tolerance)
{
  solver = meth;

  if (meth == OSNSP_LCPSOLVING)
  {
    strcpy(solvingMethod.lcp.name, OSNSP_NSQP.c_str());
    solvingMethod.lcp.tol = tolerance;
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem::setNSQPAlgorithm - solving method " + meth + " doesn't exists.");
}

void OneStepNSProblem::setCPGAlgorithm(const string& meth,  const double& tolerance, const int& iter, const string& norm)
{
  solver = meth;

  if (meth == OSNSP_LCPSOLVING)
  {
    strcpy(solvingMethod.lcp.name, OSNSP_CPG.c_str());
    solvingMethod.lcp.tol = tolerance;
    /*##### normType is not yet implemented in Numerics  #####*/
    //strcpy( solvingMethod.lcp.normType, norm.c_str() );
    solvingMethod.lcp.itermax = iter;
  }
  else if (meth == OSNSP_PRSOLVING)
  {
    strcpy(solvingMethod.pr.name, OSNSP_CPG.c_str());
    solvingMethod.pr.tol = tolerance;
    /*##### normType is not yet implemented in Numerics  #####*/
    //strcpy( solvingMethod.pr.normType, norm.c_str() );
    solvingMethod.pr.itermax = iter;
  }
  else if (meth == OSNSP_DRSOLVING)
  {
    strcpy(solvingMethod.dr.name, OSNSP_CPG.c_str());
    solvingMethod.dr.tol = tolerance;
    /*##### normType is not yet implemented in Numerics  #####*/
    //strcpy( solvingMethod.dr.normType, norm.c_str() );
    solvingMethod.dr.itermax = iter;
  }
  else if (meth == OSNSP_PFC_2DSOLVING)
  {
    strcpy(solvingMethod.pfc_2D.name, OSNSP_CPG.c_str());
    solvingMethod.pfc_2D.tol = tolerance;
    /*##### normType is not yet implemented in Numerics  #####*/
    //strcpy( solvingMethod.pfc_2D.normType, norm.c_str() );
    solvingMethod.pfc_2D.itermax = iter;
  }
  else if (meth == OSNSP_DFC_2DSOLVING)
  {
    strcpy(solvingMethod.dfc_2D.name, OSNSP_CPG.c_str());
    solvingMethod.dfc_2D.tol = tolerance;
    /*##### normType is not yet implemented in Numerics  #####*/
    //strcpy( solvingMethod.dfc_2D.normType, norm.c_str() );
    solvingMethod.dfc_2D.itermax = iter;
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem::setCPGAlgorithm - solving method " + meth + " doesn't exists.");
}

void OneStepNSProblem::setLatinAlgorithm(const string& meth, const double& t, const int& iter, const string& norm, const double& searchdirection)
{
  solver = meth;

  if (meth == OSNSP_LCPSOLVING)
  {
    strcpy(solvingMethod.lcp.name, OSNSP_LATIN.c_str());
    solvingMethod.lcp.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.lcp.normType, norm.c_str());
    solvingMethod.lcp.itermax = iter;
    solvingMethod.lcp.k_latin = searchdirection;
  }
  else if (meth == OSNSP_PRSOLVING)
  {
    strcpy(solvingMethod.pr.name, OSNSP_LATIN.c_str());
    solvingMethod.pr.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.dr.normType, norm.c_str());
    solvingMethod.pr.itermax = iter;
    solvingMethod.pr.k_latin = searchdirection;
  }
  else if (meth == OSNSP_DRSOLVING)
  {
    strcpy(solvingMethod.dr.name, OSNSP_LATIN.c_str());
    solvingMethod.dr.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.dr.normType, norm.c_str());
    solvingMethod.dr.itermax = iter;
    solvingMethod.dr.k_latin = searchdirection;
  }
  else if (meth == OSNSP_PFC_2DSOLVING)
  {
    strcpy(solvingMethod.pfc_2D.name, OSNSP_LATIN.c_str());
    solvingMethod.pfc_2D.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.pfc_2D.normType, norm.c_str());
    solvingMethod.pfc_2D.itermax = iter;
    solvingMethod.pfc_2D.k_latin = searchdirection;
  }
  else if (meth == OSNSP_DFC_2DSOLVING)
  {
    strcpy(solvingMethod.dfc_2D.name, OSNSP_LATIN.c_str());
    solvingMethod.dfc_2D.tol = t;
    /*##### normType is not yet implemented in Numerics  #####*/
    strcpy(solvingMethod.dfc_2D.normType, norm.c_str());
    solvingMethod.dfc_2D.itermax = iter;
    solvingMethod.dfc_2D.k_latin = searchdirection;
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem::setLatinAlgorithm - solving method " + meth + " doesn't exists.");
}

bool OneStepNSProblem::isOneStepNsProblemComplete()
{
  bool isComplete = true;

  if (nspbType != "LCP" || nspbType != "DFC_2D" || nspbType != "QP" || nspbType != "Relay")
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

  if (solver != "LcpSolving" || solver != "PrimalRelaySolving" || solver != "DualRelaySolving"
      || solver != "PrimalFrictionContact2DSolving" || solver != "DualFrictionContact2DSolving")
  {
    cout << "OneStepNSProblem is not complete: unknown solver type " << solver  << endl;
    isComplete = false;
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

