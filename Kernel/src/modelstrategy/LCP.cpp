//$Id: LCP.cpp,v 1.67 2005/03/23 15:03:55 jbarbier Exp $
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

  cout << " #%#% use of Numerics" << endl;
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
          W1->PLUInverseInPlace();
        }
        else
          RuntimeException::selfThrow("LCP::computeA not yet implemented for Integrator of type " + I->getType());

        n = vDS[1]->getNumber();
        I2 = this->strategy->getIntegratorOfDS(n);
        if (I2->getType() == MOREAU_INTEGRATOR)
        {
          M2 = static_cast<Moreau*>(I2);
          W2 = M2->getWPtr();
          W2->PLUInverseInPlace();
        }
        else
          RuntimeException::selfThrow("LCP::computeA not yet implemented for Integrator of type " + I->getType());

        // Assemble of W
        v[0] = W1;
        v[1] = W2;
        WW = BlockMatrixAssemble(v);
      }
      else
        RuntimeException::selfThrow("LCP::computeA not yet implemented for one DS in a Interaction ");

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

          cout << "#_# " << currentActiveInteraction << " - " << interConnectedNumber << endl;
          this->M.blockMatrixCopy(Mtmp, currentActiveInteraction, interConnectedNumber);
          interConnectedNumber++;
        }
      }
      // incrementation of the number of active interaction for the M creation
      currentActiveInteraction++;
    } // test on status
  }

  //  if( activeInteraction > 0 )
  //  {
  //    this->M.display();
  //    cout<<"<<enter>>"<<endl;
  //    getchar();
  //  }
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

//$Log: LCP.cpp,v $
//Revision 1.67  2005/03/23 15:03:55  jbarbier
//- adaptation to the LMGC90 tags in non smooth dynamical system and strategy
//
//Revision 1.66  2005/03/22 15:55:05  jbarbier
//- class NewtonImpactFriction non smooth law added to the kernel
//
//- xml schema modified for this new class
//- xml schema modified to accept a "joker" for further use of a LMGC90 mechanical plugin
//
//- new test added for the loading/saving of a NewtonImpactFrictionNSL
//
//Revision 1.65  2005/03/17 10:15:15  jbarbier
//- python sample of the bouncingball without xml input file
//
//Revision 1.64  2005/03/08 14:23:43  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.63  2005/03/08 12:41:37  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.62  2005/03/07 13:17:20  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.61  2005/03/04 15:35:26  jbarbier
//- README files added for some samples
//
//- beginning of the refactoring of XML module constants
//
//Revision 1.60  2005/03/04 08:05:29  jbarbier
//- bug corrected when filling M matrix of the LCP
//
//Revision 1.59  2005/03/02 16:06:34  jbarbier
//- DoubleContact sample runnig successfully!
//
//- computeM and computeQ of LCP fixed
//
//Revision 1.58  2005/03/01 16:44:23  charlety
//
//_ configure : first version of library version checking ok
//
//Revision 1.57  2005/03/01 15:53:09  jbarbier
//- new sample in progress : 3 balls with 2 balls which are going to touch th third ball at the same time
//
//Revision 1.56  2005/03/01 10:38:20  jbarbier
//- RollingBalls sample is OK
//
//Revision 1.55  2005/02/28 16:22:33  jbarbier
//- rolling balls sample almost finished
//
//- in LCP, compute function now use all the interactions to make computations
//
//Revision 1.54  2005/02/24 15:50:19  jbarbier
//- LCP prepared to changes needed for several interactions
//
//- new function for the SiconosMatrices to copy a block matrix into another matrix
//
//- tag names of BoundaryConditionXML, DSInputOutputXML, DSXML, InteractionXML, LagrangianLinearRXML, LagrangianNLDSXML put in XMLTagNames.h
//
//Revision 1.53  2005/02/24 11:03:05  charlety
//
//_ New attempt with the autotools.
//
//_ /!\ THIS VERSION IS USABLE ONLY IF YOU HAVE INSTALLED THE EXTERNAL LIBRARIES (cppunit, libxml, lapack++, lapack, nana) in /usr/ or /usr/local.
//
//_ This version was only tested on Fedora core2 for the moment.
//
//Revision 1.52  2005/02/15 15:15:33  charlety
//
//_ modified some very slow functions to increase performance
//
//Revision 1.51  2005/02/14 09:52:22  charlety
//_ getters / setters put inline
//
//Revision 1.50  2005/02/10 10:35:19  jbarbier
//- new file regrouping all the const values of the model, modelingTools and numericalStrategy
//
//- new function in the LagrangianLinearR to get the H matrix corresponding to one of the 2 dynamical systems linked to the relation
//
//- new atribute of the OneStepNSProblem. A visibility table of the Interaction.
//
//Revision 1.49  2005/02/08 09:54:34  jbarbier
//- optimisation of the call of Numerics function in the LCP
//
//- new function in SiconosMatrix to get the double* corresponding to the array of double values of the matrix
//
//Revision 1.48  2005/02/04 14:52:44  jbarbier
//- Rolling balls in progress (contact is detected)
//
//- time data is given in parameter for computInput and Output in the Relation. Parameter is sent by methods of the OneStepNSProblem
//
//Revision 1.47  2005/02/04 07:46:21  jbarbier
//- last modification for RollingBalls
//
//Revision 1.46  2005/02/02 15:54:51  jbarbier
//- sample RollingBalls added
//
//- function getArray() added to SimpleVector to return the pointer on the array of double values
//
//Revision 1.45  2005/02/01 11:08:42  charlety
//
//_ some displays of values during computations suppressed.
//
//Revision 1.44  2005/02/01 10:43:48  jbarbier
//- BouncingBall sample taking account of the solver defined in the XML input file
//
//Revision 1.43  2005/01/31 16:26:25  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.42  2005/01/27 13:57:46  jbarbier
//- suppression of old LCP and QP structures
//
//Revision 1.41  2004/12/08 13:53:54  jbarbier
//- Numerics : name of the methode in structures of the SiconosNumerics.h changed
//from char* to char[64] to simplify the use
//
//- Kernel : t, T, and t0 initialised to pass through READ_UNINIT_MEM
//
//Revision 1.40  2004/12/08 12:49:38  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.39  2004/12/06 10:10:34  jbarbier
//- integration of Numerics and use of Numerics on the bouncing ball sample
//
//- Numerics is now needed to run the bouncing ball sample!
//
//Revision 1.38  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.37  2004/09/22 14:11:13  charlety
//
//  _ revision of Doxygen comments in modelstrategy
//
//Revision 1.36  2004/09/21 11:49:09  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.35  2004/09/16 11:35:25  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.34  2004/09/15 13:23:13  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.33  2004/09/14 13:49:54  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.32  2004/09/10 11:26:16  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.31  2004/09/09 08:57:44  jbarbier
//- functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//createTimeDiscretisation of the Strategy done.
//
//=> all functions to create manually the objects of the platform are done
//
//Revision 1.30  2004/08/18 14:37:19  jbarbier
//- creation of Model, NSDS, Strategy(TimeStepping and EventDriven) and
//DynamicalSystem available when the creation is in a command program
//
//Revision 1.29  2004/08/12 11:55:18  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.28  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.27  2004/07/28 07:43:45  jbarbier
//- all methods createObjectOfThePlatform(...) are now existing
//
//Revision 1.26  2004/07/06 14:54:49  acary
//Renaming NSLaw into NonSmoothLaw
//Renaming RelayNSLaw into RelayNSL
//Renaming CCNSLaw into ComplementarityConditionNSL
//Renaming NewtonImpactLaw into NewtonImpactLawNSL
//
//Revision 1.25  2004/07/06 08:09:09  acary
//Renaming Class LLRelation into LagrangianLinearR
//
//Revision 1.24  2004/07/05 12:38:31  charlety
//
//try of false plugin developed in LagrangianTIDS. The Moreau integrator considers it like a LagrangianNLDS, but this is not the plugin which is used to compute the external strength, but a function of the LagrangianTIDS.
//
//Revision 1.23  2004/07/02 14:50:31  acary
//Added correct PLU factorization and Inversion of the Matrix
//
//Revision 1.22  2004/06/30 13:35:56  acary
//Formalization and Computation of the LCP with the NewtonImpactLawNSL
//
//Revision 1.21  2004/06/29 15:12:02  acary
//Change in the naming comvention for the LCP
//The LCP Matrix is now denoted by M.
//The LCP Vector is now denoted by q.
//
