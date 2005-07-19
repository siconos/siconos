#include "LCP.h"

// includes to be deleted thanks to factories
#include "Moreau.h"
#include "LagrangianLinearR.h"
#include "NewtonImpactLawNSL.h"
#include "LinearTIR.h"

using namespace std;


// Default constructor
LCP::LCP(): OneStepNSProblem(), nLcp(0), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  nspbType = LCP_OSNSP;
}

// xml constructor
LCP::LCP(OneStepNSProblemXML* osNsPbXml, Strategy* newStrat):
  OneStepNSProblem(osNsPbXml, newStrat), nLcp(0), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  nspbType = LCP_OSNSP;
  if (osNsPbXml != NULL)
  {
    LCPXML * xmllcp = (static_cast<LCPXML*>(osNsPbXml));
    // no getter-xml for nlcp ...
    if (xmllcp->hasM())
    {
      nLcp = (xmllcp->getM()).size(0);
      M = new SiconosMatrix(nLcp, nLcp);
      isMAllocatedIn = true;
      *M = xmllcp->getM();
    }
    if (xmllcp->hasQ())
    {
      nLcp = (xmllcp->getQ()).size(0);
      q = new SimpleVector(nLcp);
      isQAllocatedIn = true;
      *q = xmllcp->getQ();
    }
    if (xmllcp->hasQ() && xmllcp->hasM() && ((xmllcp->getM()).size(0) != (xmllcp->getQ()).size(0)))
      RuntimeException::selfThrow("LCP: xml constructor, inconsistent sizes between given q and M");
  }
  else RuntimeException::selfThrow("LCP: xml constructor, xml file=NULL");
}

// Constructor from a set of data
LCP::LCP(Strategy* newStrat, const string& newSolver, const string& newSolvingMethod,
         const int& MaxIter, const double & Tolerance, const string & NormType,
         const double & SearchDirection):
  OneStepNSProblem(newStrat, newSolver, newSolvingMethod, MaxIter, Tolerance, NormType, SearchDirection),
  nLcp(0), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  nspbType = LCP_OSNSP;
}

// destructor
LCP::~LCP()
{
  if (isWAllocatedIn)
  {
    delete w;
    w = NULL;
  }
  if (isZAllocatedIn)
  {
    delete z;
    z = NULL;
  }
  if (isMAllocatedIn)
  {
    delete M;
    M = NULL;
  }
  if (isQAllocatedIn)
  {
    delete q;
    q = NULL;
  }
}

void LCP::formalize(const double& time)
{
  IN("LCP::formalize(void)\n");
  // compute M and q operators for LCP problem
  computeM();
  computeQ(time);

  // check if size w and z are allright; if not, reallocate.
  if ((isWAllocatedIn) && (w->size() != q->size()))
  {
    delete w;
    if (isZAllocatedIn) delete z;
    w = new SimpleVector(nLcp);
    isWAllocatedIn = true;
    z = new SimpleVector(nLcp);
    isZAllocatedIn = true;
  }
  OUT("LCP::formalize(void)\n");
}

// Setters

void LCP::setW(const SimpleVector& newValue)
{
  if (nLcp != newValue.size())
    RuntimeException::selfThrow("LCP: setW, inconsistent size between given w size and problem size. You should set nLcp before");

  if (!isWAllocatedIn)
  {
    w = new SimpleVector(nLcp);
    isWAllocatedIn = true;
  }

  *w = newValue;
}

void LCP::setWPtr(SimpleVector* newPtr)
{
  if (isWAllocatedIn) delete w;
  w = newPtr;
  isWAllocatedIn = false;
}


void LCP::setZ(const SimpleVector& newValue)
{
  if (nLcp != newValue.size())
    RuntimeException::selfThrow("LCP: setZ, inconsistent size between given z size and problem size. You should set nLcp before");

  if (!isZAllocatedIn)
  {
    z = new SimpleVector(nLcp);
    isZAllocatedIn = true;
  }

  *z = newValue;
}

void LCP::setZPtr(SimpleVector* newPtr)
{
  if (isZAllocatedIn) delete z;
  z = newPtr;
  isZAllocatedIn = false;
}

void LCP::setM(const SiconosMatrix& newValue)
{
  if (nLcp != newValue.size(0) || nLcp != newValue.size(1))
    RuntimeException::selfThrow("LCP: setM, inconsistent size between given M size and problem size. You should set nLcp before");

  if (!isMAllocatedIn)
  {
    M = new SiconosMatrix(nLcp, nLcp);
    isMAllocatedIn = true;
  }

  *M = newValue;
}

void LCP::setMPtr(SiconosMatrix* newPtr)
{
  if (isMAllocatedIn) delete M;
  M = newPtr;
  isMAllocatedIn = false;
}


void LCP::setQ(const SimpleVector& newValue)
{
  if (nLcp != newValue.size())
    RuntimeException::selfThrow("LCP: setQ, inconsistent size between given q size and problem size. You should set nLcp before");

  if (!isQAllocatedIn)
  {
    q = new SimpleVector(nLcp);
    isQAllocatedIn = true;
  }

  *q = newValue;
}

void LCP::setQPtr(SimpleVector* newPtr)
{
  if (isQAllocatedIn) delete q;
  q = newPtr;
  isQAllocatedIn = false;
}

void LCP::compute()
{
  IN("LCP::compute(void)\n");
  int res;

  // --- Check that nLcp corresponds to the number of active interactions ---
  if (nLcp != connectedInteractionMap.size())
    RuntimeException::selfThrow("LCP: compute(), LCP problem size differs from active interaction number");

  if (!isWAllocatedIn)
  {
    w = new SimpleVector(nLcp);
    isWAllocatedIn = true;
  }
  if (!isZAllocatedIn)
  {
    z = new SimpleVector(nLcp);
    isZAllocatedIn = true ;
  }

  // Call Numerics solver
  if (nLcp != 0)
    res = solve_lcp(M->getArray(), q->getArray(), (int*)&nLcp, &solvingMethod, z->getArray(), w->getArray());

  // Update the relation
  SimpleVector *yDot, *lambda, *ynew;
  int activeInteraction = 0;

  if (nLcp != 0)
  {
    vector<Interaction*>::iterator it;
    for (it = interactionVector.begin(); it != interactionVector.end(); it++)
    {
      lambda = (*it) ->getLambdaPtr();
      lambda->zero();

      // if the current interaction is in the connected interaction map, set yDot and lambda to w and z
      if (connectedInteractionMap.find(*it) != connectedInteractionMap.end())
      {
        yDot = (*it)->getYDotPtr();
        ynew = (*it)-> getYPtr();

        Relation *R = (*it)->getRelationPtr();
        if (R->getType() == LAGRANGIANLINEARRELATION)
          (*yDot)(0) = (*w)(activeInteraction);
        else if (R->getType() == LINEARTIRELATION)
          (*ynew)(0) = (*w)(activeInteraction);
        else
          RuntimeException::selfThrow("LCP::compute not yet implemented for relation of type " + R->getType());
        (*lambda)(0) = (*z)(activeInteraction);
        activeInteraction++;
      }
    }
  }
  OUT("LCP::compute(void)\n");
}

void LCP::computeM()
{
  IN("LCP::computeM(void)\n");
  int number, orgDSRank, connectedDSRank;
  int currentActiveInteraction = 0;
  int interConnectedNumber = 0;
  vector<Connection*> vCo;
  SiconosMatrix *H, *WW ;
  bool isWWAllocatedIn = false;
  SiconosMatrix orgH, connectedH, wTmp, Mtmp;
  Relation *R, *RConnected;
  LagrangianLinearR *LLR ;
  unsigned int i;
  // --- Count the number of active interactions ---
  nLcp = connectedInteractionMap.size();

  if (isMAllocatedIn) delete M;
  M = new SiconosMatrix(nLcp, nLcp);
  isMAllocatedIn = true;
  M->zero();

  // \todo improve WW organization (WW block-structured)

  // --- For each interaction in the Map (ie active interaction) ... ---
  map<Interaction* , vector <Connection*> >::iterator iter;
  for (iter = connectedInteractionMap.begin(); iter != connectedInteractionMap.end(); ++iter)
  {
    Interaction *CurrentInteraction = iter->first;

    // --- Check if W matrix of Moreau's integrator is already inversed ---
    // At the time, W is inversed in initialize of Moreau.cpp
    // -> \todo improve this step using forward-backward to compute inverse of Mlcp blocks rather than W
    vector<DynamicalSystem*> vDS = CurrentInteraction ->getDynamicalSystems();
    unsigned int sizeDS = vDS.size();

    // Get W matrix of each DS concerned by the interaction and save in into v
    vector<SiconosMatrix*> v;
    vector<OneStepIntegrator*> OsiV;
    OsiV.reserve(sizeDS);
    for (i = 0; i < sizeDS; i++)
    {
      number = vDS[i]->getNumber();
      OsiV[i] = strategy->getIntegratorOfDSPtr(number);
      if (OsiV[i]->getType() == MOREAU_INTEGRATOR)
        v.push_back((static_cast<Moreau*>(OsiV[i]))->getWPtr());
      else
        RuntimeException::selfThrow("LCP::computeM not yet implemented for Integrator of type " + OsiV[i]->getType());
    }

    // Built block matrix WW with all the W previously saved in v..
    unsigned int size = v[0]->size(0);
    for (i = 1; i < sizeDS; i++)
      size += v[i]->size(0);
    WW = new SiconosMatrix(size, size);
    isWWAllocatedIn = true;
    *WW = BlockMatrixAssemble(v);

    // --- Get the relation parameters and compute M ---
    R = CurrentInteraction->getRelationPtr();
    if (R->getType() == LAGRANGIANLINEARRELATION)
    {
      //  compute H W Ht
      LLR = static_cast<LagrangianLinearR*>(R);
      H = LLR->getHPtr();
      Mtmp = *H * WW->multTranspose(*H);
      M->blockMatrixCopy(Mtmp, currentActiveInteraction, currentActiveInteraction);
    }
    else if (R->getType() == LINEARTIRELATION)
    {
      LinearTIR *LTIR = static_cast<LinearTIR*>(R);
      SiconosMatrix* Bloc = LTIR->getBPtr();
      SiconosMatrix* Cloc = LTIR->getCPtr();
      SiconosMatrix* Dloc = LTIR->getDPtr();
      double h = strategy->getTimeDiscretisationPtr()->getH(); // time step
      unsigned int sizeMtmp = Dloc->size(0);
      SiconosMatrix Mtmp(sizeMtmp, sizeMtmp);
      Mtmp = (h * (*Cloc * *WW * *Bloc)) + *Dloc;
      M->blockMatrixCopy(Mtmp, currentActiveInteraction, currentActiveInteraction);
    }
    else
      RuntimeException::selfThrow("LCP::computeM [level1] not yet implemented for relation of type " + R->getType());

    // --- Compute M for connected interactions ---
    interConnectedNumber = 0;
    if (iter ->second[0] != NULL)
    {
      // get from the map the connexion vector of the current interaction
      vCo = iter -> second ;
      for (unsigned int k = 0; k < vCo.size(); k++)
      {
        orgDSRank = vCo[k]->originInteractionDSRank;
        connectedDSRank = vCo[k]->connectedInteractionDSRank;

        // get W(Moreau) of common DS
        wTmp = *v[orgDSRank];

        // get H matrix of the common DS
        // /!\ we supose that all the DS have the same size !!!!
        if (R->getType() == LAGRANGIANLINEARRELATION)
        {
          LLR = static_cast<LagrangianLinearR*>(R);
          // /!\ copy of matrices !!!
          orgH = LLR->getHRelatingToDS(orgDSRank);
        }
        else
          RuntimeException::selfThrow("LCP::computeM [level2] not yet implemented for relation of type " + R->getType());

        // get H matrix of the connected DS
        RConnected = vCo[k]->connected->getRelationPtr();
        if (RConnected->getType() == LAGRANGIANLINEARRELATION)
        {
          LLR = static_cast<LagrangianLinearR*>(RConnected);
          // /!\ copy of matrices !!!
          connectedH = LLR->getHRelatingToDS(connectedDSRank);
        }
        else
          RuntimeException::selfThrow("LCP::computeM [level3] not yet implemented for relation of type " + RConnected->getType());

        Mtmp = orgH * wTmp.multTranspose(connectedH);
        //int interConnectedNumber = vCo[k]->connected->getNumber();

        // \todo : to be verified !!
        if (interConnectedNumber == currentActiveInteraction)
          interConnectedNumber++;
        M->blockMatrixCopy(Mtmp, currentActiveInteraction, interConnectedNumber);
        interConnectedNumber++;
      }
    }
    // incrementation of the number of active interaction
    currentActiveInteraction++;
  }
  if (isWWAllocatedIn) delete WW;
  OUT("LCP::computeM(void)\n");
}

void LCP::computeQ(const double& time)
{
  IN("LCP::computeQ(void)\n");
  Relation *R;
  NonSmoothLaw *nslaw;

  SimpleVector *qLCPtmp = new SimpleVector(nLcp);

  // --- Count the number of active interactions ---
  // (Although it is supposed to be already done in computeM() -> only to be cautious)
  nLcp = connectedInteractionMap.size();

  int qPos = 0;
  if (isQAllocatedIn) delete q;
  q = new SimpleVector(nLcp);
  isQAllocatedIn = true;
  q->zero();

  // --- For each interaction in the Map (ie active interaction) ... ---
  map<Interaction* , vector <Connection*> >::iterator iter;
  for (iter = connectedInteractionMap.begin(); iter != connectedInteractionMap.end(); ++iter)
  {
    // get current interaction ...
    Interaction *CurrentInteraction = iter->first;
    // get the corresponding relation and nslaw ...
    R = CurrentInteraction->getRelationPtr();
    nslaw = CurrentInteraction->getNonSmoothLawPtr();
    if (R->getType() == LAGRANGIANLINEARRELATION)
    {
      LagrangianLinearR *LLR = static_cast<LagrangianLinearR*>(R);
      if (nslaw->getType() == NEWTONIMPACTLAWNSLAW)
      {
        NewtonImpactLawNSL * newton = static_cast<NewtonImpactLawNSL*>(nslaw);
        double e = newton->getE();
        LLR->computeFreeOutput(time);
        *qLCPtmp = CurrentInteraction -> getYDot();
        *qLCPtmp += e * CurrentInteraction -> getYDotOld();
        // Assemble q
        //q = (-1)*qLCPtmp;

        /*
         * in a LCP problem the contribution of each interaction
         * to the vector 'q' is only a vector of dimension 1
         * so for the moment the assemblage of the q vector will be the copy
         * of 1 double value into 'q' for each active interaction
         */
        (*q)(qPos++) = -(*qLCPtmp)(0);
      }
      else
        RuntimeException::selfThrow("LCP::computeQ not yet implemented for NSlaw of type " + nslaw->getType() + "and for relation of type " + R->getType());
    }
    else if (R->getType() == LINEARTIRELATION)
    {
      LinearTIR *LTIR = static_cast<LinearTIR*>(R);
      if (nslaw->getType() == COMPLEMENTARITYCONDITIONNSLAW)
      {
        SiconosMatrix* Cloc = LTIR->getCPtr();
        //  intermediate version : only one dynamical system is handled
        SiconosVector* xfreeloc = (CurrentInteraction->getDynamicalSystems())[0]->getXFreePtr();
        (*q)(qPos++) = -(*Cloc * *xfreeloc)(0);
      }
      else
        RuntimeException::selfThrow("LCP::computeQ not yet implemented for NSlaw of type " + nslaw->getType() + "and for relation of type " + R->getType());
    }
    else
      RuntimeException::selfThrow("LCP::computeQ not yet implemented for relation of type " + R->getType());
  }
  delete qLCPtmp;
  OUT("LCP::computeQ(void)\n");
}

void LCP::display() const
{
  cout << "======= LCP display ======" << endl;
  cout << "____ data of the LCP " << endl;
  cout << "| nLcp : " << nLcp << endl;
  cout << "| LCP Matrix M  : " << endl;
  if (M != NULL) M->display();
  else cout << "-> NULL" << endl;
  cout << "| LCP vector q : " << endl;
  if (q != NULL) q->display();
  else cout << "-> NULL" << endl;
  cout << "==========================" << endl;

}

void LCP::saveNSProblemToXML()
{
  IN("LCP::saveNSProblemToXML\n");
  OneStepNSProblem::saveNSProblemToXML();
  if (onestepnspbxml != NULL)
  {
    (static_cast<LCPXML*>(onestepnspbxml))->setM(*M);
    (static_cast<LCPXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("LCP::saveNSProblemToXML - OneStepNSProblemXML object not exists");
  OUT("LCP::saveNSProblemToXML\n");
}

void LCP::saveMToXML()
{
  IN("LCP::saveMToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<LCPXML*>(onestepnspbxml))->setM(*M);
  }
  else RuntimeException::selfThrow("LCP::saveMToXML - OneStepNSProblemXML object not exists");
  OUT("LCP::saveMToXML\n");
}

void LCP::saveQToXML()
{
  IN("LCP::saveQToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<LCPXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("LCP::saveQToXML - OneStepNSProblemXML object not exists");
  OUT("LCP::saveQToXML\n");
}

LCP* LCP::convert(OneStepNSProblem* osnsp)
{
  cout << "LCP::convert (OneStepNSProblem* osnsp)" << endl;
  LCP* lcp = dynamic_cast<LCP*>(osnsp);
  return lcp;
}

