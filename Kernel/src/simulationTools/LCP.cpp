#include "LCP.h"
#include "LCPXML.h"
#include "check.h"


#include "Interaction.h"
// hazardous dependency
#include "Moreau.h"
#include "LagrangianLinearR.h"
#include "NewtonImpactLawNSL.h"
#include <stdio.h>


LCP::LCP(): OneStepNSProblem()
{
  this->nspbType = LCP_OSNSP;
}

LCP::LCP(OneStepNSProblemXML* osnspbxml): OneStepNSProblem(osnspbxml)
{
  this->nspbType = LCP_OSNSP;
}

LCP::~LCP()
{}



SiconosMatrix* LCP::getMPtr(void)
{
  return &this->M;
}

SimpleVector* LCP::getQPtr(void)
{
  return &this->q;
}


void LCP::formalize(double time)
{
  IN("LCP::formalize(void)\n");
  OneStepNSProblem::formalize(time);

  //cout<<"## LCP::formalize - interactionVector.size() == "<<this->interactionVector.size()<<endl;
  for (int i = 0; i < this->interactionVector.size(); i++)
  {
    this->interactionVector[i]->check(time);
  }
  this->updateConnectedInteractionMap();

  this->w = SimpleVector::SimpleVector(0);
  this->z = SimpleVector::SimpleVector(0);

  this->computeM();
  this->computeQ(time);

  this->w = SimpleVector::SimpleVector(this->nLcp);
  this->z = SimpleVector::SimpleVector(this->nLcp);

  OUT("LCP::formalize(void)\n");

  // formalisations specific to LCP problem
  // ...
}


void LCP::compute(void)
{
  IN("LCP::compute(void)\n");

  //cout<<" #%#% use of Numerics"<<endl;
  /***********************************
   * integration of Siconos/Numerics *
   ***********************************/
  // now we'll use Numerics !!
  int res, i, j;
  int nn = this->nLcp;//this->n;

  if (this->nLcp == 0)
  {}
  else //if( this->nLcp == 1 )
  {
    res = solve_lcp(/*vec*/this->M.getArray(), /*q*/(this->q.getArray()), &nn, &(this->solvingMethod), this->z.getArray(), this->w.getArray());
    /**
     * Pretty hazardous copy !!!!!
     */

    //    res = solve_lcp(/*vec*/this->M.getArray(), /*q*/(this->q.getArray()), &n, &(this->solvingMethod), this->z.getArray(), this->w.getArray());
    //
    //
    //    for (i=0; i<n; i++)
    //    {
    //      this->w(i) = z[i];
    //      this->z(i) = w[i];
    //    }
  }
  //  else
  //  RuntimeException::selfThrow("LCP::compute not yet implemented for size > 1 ");

  // Update the relation
  SiconosVector *yDot, *lambda;
  int activeInteraction = 0;

  yDot = this->interactionVector[0]->getYDotPtr();
  lambda = this->interactionVector[0]->getLambdaPtr();

  if (this->nLcp == 0)
  {
    lambda->zero();
  }
  else //if (this->nLcp == 1)
  {
    //    cout<<"### M then q :"<<endl;
    //    this->M.display();
    //    this->q.display();
    //    cout<<"### yDot then lambda :"<<endl;
    //    this->w.display();
    //    this->z.display();
    //    cout<<"### "<<endl;
    //    this->displayConnectedInteractionMap();
    //    cout<<"### ---"<<endl;

    //    cout<<"                    <<Press Enter>>"<<endl;
    //    getchar();

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
  //cout<<"### quit LCP::compute"<<endl;
  OUT("LCP::compute(void)\n");
}


void LCP::computeM(void)
{
  IN("LCP::computeM(void)\n");
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
  this->nLcp = 0;

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
          if (!W1->isInversed())
          {
            W1->PLUInverseInPlace();
            if (W1->isInversed()) cout << "KAPOUeeeeeee ########################################################" << endl;
          }
          //  W1->display();


        }
        else
          RuntimeException::selfThrow("LCP::computeA not yet implemented for Integrator of type " + I->getType());

        n = vDS[1]->getNumber();
        I2 = this->strategy->getIntegratorOfDS(n);
        if (I2->getType() == MOREAU_INTEGRATOR)
        {
          M2 = static_cast<Moreau*>(I2);
          W2 = M2->getWPtr();
          if (!W2->isInversed()) W2->PLUInverseInPlace();
          // W2 ->display();


        }
        else
          RuntimeException::selfThrow("LCP::computeA not yet implemented for Integrator of type " + I->getType());

        // Assemble of W
        v[0] = W1;
        v[1] = W2;
        WW = BlockMatrixAssemble(v);

        W1->PLUInverseInPlace();
        W2->PLUInverseInPlace();
      }
      else
        RuntimeException::selfThrow("LCP::computeA not yet implemented for one DS in a Interaction ");

      R = this->interactionVector[i]->getRelation();
      if (R->getType() == LAGRANGIANLINEARRELATION)
      {
        LLR = static_cast<LagrangianLinearR*>(R);
        H = LLR->getHPtr();
        Mtmp = *H * WW.multTranspose(*H);
        //cout<<"#_# "<<currentActiveInteraction<<" - "<<currentActiveInteraction<<endl;
        this->M.blockMatrixCopy(Mtmp, currentActiveInteraction, currentActiveInteraction);
      }
      else
        RuntimeException::selfThrow("LCP::computeM [level1] not yet implemented for relation of type " + R->getType());

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
          //this->interactionVector[i].getDynamicalSystem( orgDSRank );
          if (R->getType() == LAGRANGIANLINEARRELATION)
          {
            LLR = static_cast<LagrangianLinearR*>(R);
            // /!\ copy of matrices !!!
            orgH = LLR->getHRelatingToDS(orgDSRank);
          }
          else
            RuntimeException::selfThrow("LCP::computeM [level2] not yet implemented for relation of type " + R->getType());

          RConnected = vCo[k]->connected->getRelation();
          if (RConnected->getType() == LAGRANGIANLINEARRELATION)
          {
            LLR = static_cast<LagrangianLinearR*>(RConnected);
            // /!\ copy of matrices !!!
            connectedH = LLR->getHRelatingToDS(connectedDSRank);
          }
          else
            RuntimeException::selfThrow("LCP::computeM [level3] not yet implemented for relation of type " + R->getType());

          Mtmp = orgH * wTmp.multTranspose(connectedH);
          //int interConnectedNumber = vCo[k]->connected->getNumber();

          // \todo : to be verified !!
          if (interConnectedNumber == currentActiveInteraction)
            interConnectedNumber++;

          //cout<<"#_# "<<currentActiveInteraction<<" - "<<interConnectedNumber<<endl;
          this->M.blockMatrixCopy(Mtmp, currentActiveInteraction, interConnectedNumber);
          interConnectedNumber++;
        }
      }
      // incrementation of the number of active interaction for the M creation
      currentActiveInteraction++;
    } // test on status
  }

  this->nLcp = activeInteraction;
  OUT("LCP::computeM(void)\n");
}

void LCP::computeQ(double time)
{
  IN("LCP::computeQ(void)\n");
  Relation *R;
  LagrangianLinearR *LLR;
  NonSmoothLaw *nslaw;
  NewtonImpactLawNSL * newton;

  double e;
  SimpleVector qLCPtmp;

  /*
   * \warning : initialisation of "q" removed! It seems that that BouncingBall sample
   * is still running good ...
   */
  int qSize = this->nLcp;
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
          qLCPtmp = this->interactionVector[i]->getYDot();

          qLCPtmp += e * this->interactionVector[i]->getYDotOld();

          // Assemble q
          //this->q = (-1)*qLCPtmp;

          /*
           * in a LCP problem the contribution of each interaction
           * to the vector 'q' is only a vector of dimension 1
           * so for the moment the assemblage of the q vector will be the copy
           * of 1 double value into 'q' for each active interaction
           */
          this->q(qPos++) = -qLCPtmp(0);

          //            cout<<"### computeQ - LCP (Ufree, Uold) :"<<endl;
          //            this->interactionVector[i]->getYDot().display();
          //            this->interactionVector[i]->getYDotOld().display();
        }
        else
          RuntimeException::selfThrow("LCP::computeQ not yet implemented for NSlaw of type " + nslaw->getType() + "and for relation of type " + R->getType());
      }
      else
        RuntimeException::selfThrow("LCP::computeQ not yet implemented for relation of type " + R->getType());
    }
    //    }
  }
  OUT("LCP::computeQ(void)\n");
}

void LCP::fillNSProblemWithNSProblemXML()
{
  OUT("LCP::fillNSProblemWithNSProblemXML\n");
  OneStepNSProblem::fillNSProblemWithNSProblemXML();
  if (this->onestepnspbxml != NULL)
  {
    this->M = (static_cast<LCPXML*>(this->onestepnspbxml))->getM();
    this->q = (static_cast<LCPXML*>(this->onestepnspbxml))->getQ();
  }
  else RuntimeException::selfThrow("LCP::fillNSProblemWithNSProblemXML - OneStepNSProblemXML object not exists");
}

void LCP::display() const
{
  cout << "------------------------------------------------------" << endl;
  cout << "____ data of the LCP " << endl;
  cout << "| nLcp : " << this->nLcp << endl;
  cout << "| LCP Matrix M  : " << endl;
  M.display();
  cout << "| LCP vector q : " << endl;
  q.display();
  cout << "____________________________" << endl;
  cout << "------------------------------------------------------" << endl;
}

void LCP::saveNSProblemToXML()
{
  IN("LCP::saveNSProblemToXML\n");
  OneStepNSProblem::saveNSProblemToXML();
  if (this->onestepnspbxml != NULL)
  {
    (static_cast<LCPXML*>(this->onestepnspbxml))->setM(&(this->M));
    (static_cast<LCPXML*>(this->onestepnspbxml))->setQ(&(this->q));
  }
  else RuntimeException::selfThrow("LCP::saveNSProblemToXML - OneStepNSProblemXML object not exists");
  OUT("LCP::saveNSProblemToXML\n");
}

void LCP::saveMToXML()
{
  IN("LCP::saveMToXML\n");
  if (this->onestepnspbxml != NULL)
  {
    (static_cast<LCPXML*>(this->onestepnspbxml))->setM(&(this->M));
  }
  else RuntimeException::selfThrow("LCP::saveMToXML - OneStepNSProblemXML object not exists");
  OUT("LCP::saveMToXML\n");
}

void LCP::saveQToXML()
{
  IN("LCP::saveQToXML\n");
  if (this->onestepnspbxml != NULL)
  {
    (static_cast<LCPXML*>(this->onestepnspbxml))->setQ(&(this->q));
  }
  else RuntimeException::selfThrow("LCP::saveQToXML - OneStepNSProblemXML object not exists");
  OUT("LCP::saveQToXML\n");
}

void LCP::createOneStepNSProblem(OneStepNSProblemXML * osnspbXML, Strategy * strategy)
{
  if (osnspbXML != NULL)
  {
    this->onestepnspbxml = osnspbXML;
    this->fillNSProblemWithNSProblemXML();
  }
  else
  {
    this->strategy = strategy;
    this->nspbType = LCP_OSNSP;
    this->fillInteractionVector();
  }
  this->init();
}


LCP* LCP::convert(OneStepNSProblem* osnsp)
{
  cout << "LCP::convert (OneStepNSProblem* osnsp)" << endl;
  LCP* lcp = dynamic_cast<LCP*>(osnsp);
  return lcp;
}
