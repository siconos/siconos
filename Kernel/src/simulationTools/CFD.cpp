#include "CFD.h"
#include "CFDXML.h"
#include "check.h"


#include "Interaction.h"
// hazardous dependency
#include "Moreau.h"
#include "LagrangianLinearR.h"
#include "NewtonImpactLawNSL.h"
#include <stdio.h>


CFD::CFD(): OneStepNSProblem()
{
  this->nspbType = CFD_OSNSP;
}

CFD::CFD(OneStepNSProblemXML* osnspbxml): OneStepNSProblem(osnspbxml)
{
  this->nspbType = CFD_OSNSP;
}

CFD::~CFD()
{}



SiconosMatrix* CFD::getMPtr(void)
{
  return &this->M;
}

SimpleVector* CFD::getQPtr(void)
{
  return &this->q;
}


void CFD::formalize(double time)
{
  IN("CFD::formalize(void)\n");
  OneStepNSProblem::formalize(time);

  //cout<<"## CFD::formalize - interactionVector.size() == "<<this->interactionVector.size()<<endl;
  for (int i = 0; i < this->interactionVector.size(); i++)
  {
    this->interactionVector[i]->check(time);
  }
  this->updateConnectedInteractionMap();

  this->w = SimpleVector::SimpleVector(0);
  this->z = SimpleVector::SimpleVector(0);

  this->computeM();
  this->computeQ(time);

  this->w = SimpleVector::SimpleVector(this->nCfd);
  this->z = SimpleVector::SimpleVector(this->nCfd);

  OUT("CFD::formalize(void)\n");

  // formalisations specific to CFD problem
  // ...
}


void CFD::compute(void)
{
  IN("CFD::compute(void)\n");

  cout << " #%#% use of Numerics" << endl;
  /***********************************
   * integration of Siconos/Numerics *
   ***********************************/
  // now we'll use Numerics !!
  int res, i, j;
  int nn = this->nCfd;//this->n;

  if (this->nCfd == 0)
  {}
  else //if( this->nCfd == 1 )
  {
    res = solve_lcp(/*vec*/this->M.getArray(), /*q*/(this->q.getArray()), &nn, &(this->solvingMethod), this->z.getArray(), this->w.getArray());
  }

  // Update the relation
  SiconosVector *yDot, *lambda;
  int activeInteraction = 0;

  yDot = this->interactionVector[0]->getYDotPtr();
  lambda = this->interactionVector[0]->getLambdaPtr();

  if (this->nCfd == 0)
  {
    lambda->zero();
  }
  else
  {
    for (int i = 0; i < this->interactionVector.size(); i++)
    {
      lambda = this->interactionVector[i]->getLambdaPtr();
      lambda->zero();

      if (this->connectedInteractionMap.find(this->interactionVector[i]) != this->connectedInteractionMap.end())
      {
        yDot = this->interactionVector[i]->getYDotPtr();
        lambda = this->interactionVector[i]->getLambdaPtr();

        (*yDot)(0) = this->w(activeInteraction);
        (*lambda)(0) = this->z(activeInteraction);

        activeInteraction++;
      }
    }
  }
  OUT("CFD::compute(void)\n");
}


void CFD::computeM(void)
{
  IN("CFD::computeM(void)\n");
  int i, j, n;
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
  //this->M = SiconosMatrix::SiconosMatrix(0, 0);
  this->nCfd = 0;

  /* New initialisation that will be good */
  /*
   * we count the number of interaction that have a status == 1
   */
  int activeInteraction = 0;
  for (i = 0; i < this->interactionVector.size(); i++)
  {
    if (this->connectedInteractionMap.find(this->interactionVector[i]) != this->connectedInteractionMap.end())
      activeInteraction++;
  }
  this->M = SiconosMatrix::SiconosMatrix(activeInteraction, activeInteraction);
  this->M.zero();

  // \todo using the visibility table instead of the interaction vector !!!!!
  for (i = 0; i < this->interactionVector.size(); i++)
  {
    if (this->connectedInteractionMap.find(this->interactionVector[i]) != this->connectedInteractionMap.end())
    {
      WW = SiconosMatrix::SiconosMatrix(0, 0);
      vDS.clear();
      vDS = this->interactionVector[i]->getDynamicalSystems();
      if (vDS.size() == 2)
      {
        n = vDS[0]->getNumber();
        I = this->strategy->getIntegratorOfDS(n);
        if (I->getType() == MOREAU_INTEGRATOR)
        {
          M1 = static_cast<Moreau*>(I);
          W1 = M1->getWPtr();
          W1->PLUInverseInPlace();
        }
        else
          RuntimeException::selfThrow("CFD::computeA not yet implemented for Integrator of type " + I->getType());

        n = vDS[1]->getNumber();
        I2 = this->strategy->getIntegratorOfDS(n);
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

      R = this->interactionVector[i]->getRelation();
      if (R->getType() == LAGRANGIANLINEARRELATION)
      {
        LLR = static_cast<LagrangianLinearR*>(R);
        H = LLR->getHPtr();
        Mtmp = *H * WW.multTranspose(*H);
        cout << "#_# " << currentActiveInteraction << " - " << currentActiveInteraction << endl;
        this->M.blockMatrixCopy(Mtmp, currentActiveInteraction, currentActiveInteraction);
      }
      else
        RuntimeException::selfThrow("CFD::computeM [level1] not yet implemented for relation of type " + R->getType());

      interConnectedNumber = 0;
      if (connectedInteractionMap[this->interactionVector[i]][0] != NULL)
      {
        for (int k = 0; k < connectedInteractionMap[this->interactionVector[i]].size(); k++)
        {
          vCo = connectedInteractionMap[this->interactionVector[i]];
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

          RConnected = vCo[k]->connected->getRelation();
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
          this->M.blockMatrixCopy(Mtmp, currentActiveInteraction, interConnectedNumber);
          interConnectedNumber++;
        }
      }
      // incrementation of the number of active interaction for the M creation
      currentActiveInteraction++;
    } // test on status
  }

  this->nCfd = activeInteraction;
  OUT("CFD::computeM(void)\n");
}

void CFD::computeQ(double time)
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
  int qSize = this->nCfd;
  int qPos = 0;
  this->q = SimpleVector::SimpleVector(qSize);
  this->q.zero();

  for (int i = 0; i < this->interactionVector.size(); i++)
  {
    if (this->connectedInteractionMap.find(this->interactionVector[i]) != this->connectedInteractionMap.end())
    {
      R = this->interactionVector[i]->getRelation();
      nslaw = this->interactionVector[i]->getNonSmoothLaw();
      if (R->getType() == LAGRANGIANLINEARRELATION)
      {
        LLR = static_cast<LagrangianLinearR*>(R);

        if (nslaw->getType() == NEWTONIMPACTLAWNSLAW)
        {
          newton = static_cast<NewtonImpactLawNSL*>(nslaw);
          e = newton->getE();
          LLR->computeFreeOutput(time);
          qCFDtmp = this->interactionVector[i]->getYDot();

          qCFDtmp += e * this->interactionVector[i]->getYDotOld();

          // Assemble q
          //this->q = (-1)*qLCPtmp;

          /*
           * in a CFD problem the contribution of each interaction
           * to the vector 'q' is only a vector of dimension 1
           * so for the moment the assemblage of the q vector will be the copy
           * of 1 double value into 'q' for each active interaction
           */
          this->q(qPos++) = -qCFDtmp(0);

          //            cout<<"### computeQ - CFD (Ufree, Uold) :"<<endl;
          //            this->interactionVector[i]->getYDot().display();
          //            this->interactionVector[i]->getYDotOld().display();
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

void CFD::fillNSProblemWithNSProblemXML()
{
  OUT("CFD::fillNSProblemWithNSProblemXML\n");
  OneStepNSProblem::fillNSProblemWithNSProblemXML();
  if (this->onestepnspbxml != NULL)
  {
    this->M = (static_cast<CFDXML*>(this->onestepnspbxml))->getM();
    this->q = (static_cast<CFDXML*>(this->onestepnspbxml))->getQ();
  }
  else RuntimeException::selfThrow("CFD::fillNSProblemWithNSProblemXML - OneStepNSProblemXML object not exists");
}

void CFD::display() const
{
  cout << "------------------------------------------------------" << endl;
  cout << "____ data of the CFD " << endl;
  cout << "| nCfd : " << this->nCfd << endl;
  cout << "| CFD Matrix M  : " << endl;
  M.display();
  cout << "| CFD vector q : " << endl;
  q.display();
  cout << "____________________________" << endl;
  cout << "------------------------------------------------------" << endl;
}

void CFD::saveNSProblemToXML()
{
  IN("CFD::saveNSProblemToXML\n");
  OneStepNSProblem::saveNSProblemToXML();
  if (this->onestepnspbxml != NULL)
  {
    (static_cast<CFDXML*>(this->onestepnspbxml))->setM(&(this->M));
    (static_cast<CFDXML*>(this->onestepnspbxml))->setQ(&(this->q));
  }
  else RuntimeException::selfThrow("CFD::saveNSProblemToXML - OneStepNSProblemXML object not exists");
  OUT("CFD::saveNSProblemToXML\n");
}

void CFD::saveMToXML()
{
  IN("CFD::saveMToXML\n");
  if (this->onestepnspbxml != NULL)
  {
    (static_cast<CFDXML*>(this->onestepnspbxml))->setM(&(this->M));
  }
  else RuntimeException::selfThrow("CFD::saveMToXML - OneStepNSProblemXML object not exists");
  OUT("CFD::saveMToXML\n");
}

void CFD::saveQToXML()
{
  IN("CFD::saveQToXML\n");
  if (this->onestepnspbxml != NULL)
  {
    (static_cast<CFDXML*>(this->onestepnspbxml))->setQ(&(this->q));
  }
  else RuntimeException::selfThrow("CFD::saveQToXML - OneStepNSProblemXML object not exists");
  OUT("CFD::saveQToXML\n");
}

void CFD::createOneStepNSProblem(OneStepNSProblemXML * osnspbXML, Strategy * strategy)
{
  if (osnspbXML != NULL)
  {
    this->onestepnspbxml = osnspbXML;
    this->fillNSProblemWithNSProblemXML();
  }
  else
  {
    this->strategy = strategy;
    this->nspbType = CFD_OSNSP;
    this->fillInteractionVector();
  }
  this->init();
}


CFD* CFD::convert(OneStepNSProblem* osnsp)
{
  cout << "CFD::convert (OneStepNSProblem* osnsp)" << endl;
  CFD* lcp = dynamic_cast<CFD*>(osnsp);
  return lcp;
}
