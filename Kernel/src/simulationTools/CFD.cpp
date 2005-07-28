#include "CFD.h"
using namespace std;

CFD::CFD(): OneStepNSProblem(), nCfd(0), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  nspbType = CFD_OSNSP;
}

CFD::CFD(OneStepNSProblemXML* osnspbxml, Strategy * newStrat):
  OneStepNSProblem(osnspbxml, newStrat), nCfd(0), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(true), isZAllocatedIn(true), isMAllocatedIn(true), isQAllocatedIn(true)
{
  nspbType = CFD_OSNSP;
  if (osnspbxml != NULL)
  {
    // no getter-xml for nCfd ...
    nCfd = ((static_cast<CFDXML*>(onestepnspbxml))->getQ()).size();
    n = nCfd;
    w = new SimpleVector(nCfd);
    z = new SimpleVector(nCfd);
    M = new SiconosMatrix(nCfd, nCfd);
    q = new SimpleVector(nCfd);
    *M = (static_cast<CFDXML*>(onestepnspbxml))->getM();
    *q = (static_cast<CFDXML*>(onestepnspbxml))->getQ();
  }
  else RuntimeException::selfThrow("CFD: xml constructor, xml file=NULL");
}

CFD::~CFD()
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

void CFD::preCFD(const double& time)
{
  IN("CFD::formalize(void)\n");

  double pasH = strategy->getTimeDiscretisationPtr()->getH();
  //cout<<"## CFD::formalize - interactionVector.size() == "<<this->interactionVector.size()<<endl;
  for (unsigned int i = 0; i < interactionVector.size(); i++)
  {
    interactionVector[i]->check(time, pasH);
  }
  updateConnectedInteractionMap();
  cout << " #~ CFD::updateConnectedInteractionMap done" << endl;
  cout << " #~ CFD::next : computeM, computeQ" << endl;
  computeM();
  computeQ(time);
  OUT("CFD::formalize(void)\n");

  // formalisations specific to CFD problem
  // ...
}


void CFD::compute(const double& time)
{
  IN("CFD::compute(void)\n");

  // pre-treatment for CFD
  preCFD(time);

  cout << " #%#% use of Numerics" << endl;
  /***********************************
   * integration of Siconos/Numerics *
   ***********************************/
  // now we'll use Numerics !!
  int res;
  int nn = nCfd;//n;

  if (nCfd == 0)
  {}
  else //if( nCfd == 1 )
  {
    res = solve_lcp(/*vec*/M->getArray(), /*q*/(q->getArray()), &nn, &(solvingMethod), z->getArray(), w->getArray());
  }

  // Update the relation
  SimpleVector *yDot, *lambda;
  int activeInteraction = 0;

  yDot = interactionVector[0]->getYPtr(1);
  lambda = interactionVector[0]->getLambdaPtr();

  if (nCfd == 0)
  {
    lambda->zero();
  }
  else
  {
    for (unsigned int i = 0; i < interactionVector.size(); i++)
    {
      lambda = interactionVector[i]->getLambdaPtr();
      lambda->zero();

      if (connectedInteractionMap.find(interactionVector[i]) != connectedInteractionMap.end())
      {
        yDot = interactionVector[i]->getYPtr(1);
        lambda = interactionVector[i]->getLambdaPtr();

        (*yDot)(0) = (*w)(activeInteraction);
        (*lambda)(0) = (*z)(activeInteraction);

        activeInteraction++;
      }
    }
  }
  OUT("CFD::compute(void)\n");
}


void CFD::computeM(void)
{
  IN("CFD::computeM(void)\n");
  unsigned int i, n;
  int orgDSRank, connectedDSRank;
  int currentActiveInteraction = 0;
  int interConnectedNumber = 0;
  vector<int> status;
  vector<SiconosMatrix*> v(2);
  vector<DynamicalSystem*> vDS;
  vector<Connection*> vCo;

  SiconosMatrix *W1, *W2;
  //SiconosMatrix W1, W2;

  SiconosMatrix WW;
  SiconosMatrix orgH, connectedH;
  SiconosMatrix wTmp, Mtmp;
  SiconosMatrix *H;

  Relation* R;
  Relation *RConnected;
  OneStepIntegrator *I, *I2;
  Moreau *M1, *M2;
  LagrangianLinearR *LLR;

  /* temporary creation of the matrix M, not good for other example than the BouncingBall */
  //M = SiconosMatrix::SiconosMatrix(0, 0);
  nCfd = 0;

  /* New initialisation that will be good */
  /*
   * we count the number of interaction that have a status == 1
   */
  int activeInteraction = 0;
  for (i = 0; i < interactionVector.size(); i++)
  {
    if (connectedInteractionMap.find(interactionVector[i]) != connectedInteractionMap.end())
      activeInteraction++;
  }
  // ?????? size of the M matrix in Contact Friction Dual ?????
  //  M = SiconosMatrix::SiconosMatrix(activeInteraction/**2*/, activeInteraction/**2*/);
  M->zero();

  // \todo using the visibility table instead of the interaction vector !!!!!
  for (i = 0; i < interactionVector.size(); i++)
  {
    if (connectedInteractionMap.find(interactionVector[i]) != connectedInteractionMap.end())
    {
      WW = SiconosMatrix::SiconosMatrix(0, 0);
      vDS.clear();
      vDS = interactionVector[i]->getDynamicalSystems();
      if (vDS.size() == 2)
      {
        n = vDS[0]->getNumber();
        I = strategy->getIntegratorOfDSPtr(n);
        if (I->getType() == MOREAU_INTEGRATOR)
        {
          M1 = static_cast<Moreau*>(I);
          W1 = M1->getWPtr();
          W1->PLUInverseInPlace();
        }
        else
          RuntimeException::selfThrow("CFD::computeA not yet implemented for Integrator of type " + I->getType());

        n = vDS[1]->getNumber();
        I2 = strategy->getIntegratorOfDSPtr(n);
        if (I2->getType() == MOREAU_INTEGRATOR)
        {
          M2 = static_cast<Moreau*>(I2);
          W2 = M2->getWPtr();
          W2->PLUInverseInPlace();
        }
        else
          RuntimeException::selfThrow("CFD::computeA not yet implemented for Integrator of type " + I->getType());

        // Assemble of W
        v[0] = W1;
        v[1] = W2;
        WW = BlockMatrixAssemble(v);

        W1->PLUInverseInPlace();
        W2->PLUInverseInPlace();
      }
      else
        RuntimeException::selfThrow("CFD::computeA not yet implemented for one DS in a Interaction ");

      R = interactionVector[i]->getRelationPtr();
      if (R->getType() == LAGRANGIANLINEARRELATION)
      {
        LLR = static_cast<LagrangianLinearR*>(R);
        H = LLR->getHPtr();
        Mtmp = *H * WW.multTranspose(*H);
        cout << "#_# " << currentActiveInteraction << " - " << currentActiveInteraction << endl;
        M->display();
        Mtmp.display();
        cout << "___________________" << endl;
        M->blockMatrixCopy(Mtmp, currentActiveInteraction, currentActiveInteraction);
      }
      else
        RuntimeException::selfThrow("CFD::computeM [level1] not yet implemented for relation of type " + R->getType());

      interConnectedNumber = 0;
      if (connectedInteractionMap[interactionVector[i]][0] != NULL)
      {
        for (unsigned int k = 0; k < connectedInteractionMap[interactionVector[i]].size(); k++)
        {
          vCo = connectedInteractionMap[interactionVector[i]];
          orgDSRank = vCo[k]->originInteractionDSRank;
          connectedDSRank = vCo[k]->connectedInteractionDSRank;

          // getting the right W
          wTmp = *v[orgDSRank];

          // getting H matrix of the common DS
          // /!\ we supose that all the DS have the same size !!!!
          if (R->getType() == LAGRANGIANLINEARRELATION)
          {
            LLR = static_cast<LagrangianLinearR*>(R);
            // /!\ copy of matrices !!!
            orgH = LLR->getHRelatingToDS(orgDSRank);
          }
          else
            RuntimeException::selfThrow("CFD::computeM [level2] not yet implemented for relation of type " + R->getType());

          RConnected = vCo[k]->connected->getRelationPtr();
          if (RConnected->getType() == LAGRANGIANLINEARRELATION)
          {
            LLR = static_cast<LagrangianLinearR*>(RConnected);
            // /!\ copy of matrices !!!
            connectedH = LLR->getHRelatingToDS(connectedDSRank);
          }
          else
            RuntimeException::selfThrow("CFD::computeM [level3] not yet implemented for relation of type " + R->getType());

          Mtmp = orgH * wTmp.multTranspose(connectedH);
          //int interConnectedNumber = vCo[k]->connected->getNumber();

          // \todo : to be verified !!
          if (interConnectedNumber == currentActiveInteraction)
            interConnectedNumber++;

          cout << "#_# " << currentActiveInteraction << " - " << interConnectedNumber << endl;
          M->blockMatrixCopy(Mtmp, currentActiveInteraction, interConnectedNumber);
          interConnectedNumber++;
        }
      }
      // incrementation of the number of active interaction for the M creation
      currentActiveInteraction++;
    } // test on status
  }

  nCfd = activeInteraction;
  OUT("CFD::computeM(void)\n");
}

void CFD::computeQ(const double& time)
{
  IN("CFD::computeQ(void)\n");
  Relation *R;
  LagrangianLinearR *LLR;
  NonSmoothLaw *nslaw;
  NewtonImpactLawNSL * newton;

  double e;
  SimpleVector qCFDtmp;

  /*
   * \warning : initialisation of "q" removed! It seems that that BouncingBall sample
   * is still running good ...
   */
  int qPos = 0;
  //q = SimpleVector::SimpleVector(qSize);
  q->zero();

  for (unsigned int i = 0; i < interactionVector.size(); i++)
  {
    if (connectedInteractionMap.find(interactionVector[i]) != connectedInteractionMap.end())
    {
      R = interactionVector[i]->getRelationPtr();
      nslaw = interactionVector[i]->getNonSmoothLawPtr();
      if (R->getType() == LAGRANGIANLINEARRELATION)
      {
        LLR = static_cast<LagrangianLinearR*>(R);

        if (nslaw->getType() == NEWTONIMPACTLAWNSLAW)
        {
          newton = static_cast<NewtonImpactLawNSL*>(nslaw);
          e = newton->getE();
          LLR->computeFreeOutput(time);
          qCFDtmp = interactionVector[i]->getY(1);

          qCFDtmp += e * interactionVector[i]->getYOld(1);

          // Assemble q
          //q = (-1)*qLCPtmp;

          /*
           * in a CFD problem the contribution of each interaction
           * to the vector 'q' is only a vector of dimension 1
           * so for the moment the assemblage of the q vector will be the copy
           * of 1 double value into 'q' for each active interaction
           */
          (*q)(qPos++) = -qCFDtmp(0);

          //            cout<<"### computeQ - CFD (Ufree, Uold) :"<<endl;
          //            interactionVector[i]->getYDot().display();
          //            interactionVector[i]->getYDotOld().display();
        }
        else
          RuntimeException::selfThrow("CFD::computeQ not yet implemented for NSlaw of type " + nslaw->getType() + "and for relation of type " + R->getType());
      }
      else
        RuntimeException::selfThrow("CFD::computeQ not yet implemented for relation of type " + R->getType());
    }
    //    }
  }
  OUT("CFD::computeQ(void)\n");
}

void CFD::display() const
{
  cout << "------------------------------------------------------" << endl;
  cout << "____ data of the CFD " << endl;
  cout << "| nCfd : " << nCfd << endl;
  cout << "| CFD Matrix M  : " << endl;
  if (M != NULL) M->display();
  else cout << "-> NULL" << endl;
  cout << "| CFD vector q : " << endl;
  if (q != NULL) q->display();
  else cout << "-> NULL" << endl;
  cout << "____________________________" << endl;
  cout << "------------------------------------------------------" << endl;
}

void CFD::saveNSProblemToXML()
{
  IN("CFD::saveNSProblemToXML\n");
  OneStepNSProblem::saveNSProblemToXML();
  if (onestepnspbxml != NULL)
  {
    (static_cast<CFDXML*>(onestepnspbxml))->setM(*M);
    (static_cast<CFDXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("CFD::saveNSProblemToXML - OneStepNSProblemXML object not exists");
  OUT("CFD::saveNSProblemToXML\n");
}

void CFD::saveMToXML()
{
  IN("CFD::saveMToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<CFDXML*>(onestepnspbxml))->setM(*M);
  }
  else RuntimeException::selfThrow("CFD::saveMToXML - OneStepNSProblemXML object not exists");
  OUT("CFD::saveMToXML\n");
}

void CFD::saveQToXML()
{
  IN("CFD::saveQToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<CFDXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("CFD::saveQToXML - OneStepNSProblemXML object not exists");
  OUT("CFD::saveQToXML\n");
}

CFD* CFD::convert(OneStepNSProblem* osnsp)
{
  cout << "CFD::convert (OneStepNSProblem* osnsp)" << endl;
  CFD* lcp = dynamic_cast<CFD*>(osnsp);
  return lcp;
}
