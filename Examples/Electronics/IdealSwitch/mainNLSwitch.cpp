
#include "SiconosKernel.hpp"
#include "circuit.h"
#include "elecRelation.h"
#include "myDS.h"
#include <stdio.h>
#include <stdlib.h>
#include "circuit.h"

using namespace std;

/************************************************************/
/************************************************************/
/************************************************************/
/*call back for the source*/
/*call back for the formulation with inversion*/
//void (bLDS) (double t, unsigned int N, double* b, unsigned int z, double*zz){
//}


/************************************************************/
/************************************************************/
/************************************************************/
/************************************************************/
/*main program*/


int main()
{
  string solverDirEnum = "DIRECT_ENUM" ;
  string solverDirPath = "DIRECT_PATH" ;
  string solverDirSimplex = "DIRECT_SIMPLEX" ;
  string solverEnum = "ENUM" ;
  string solverSimplex = "SIMPLEX" ;
  string solverPath = "PATH" ;
  string * solverName = 0;
  bool diodeIsOn;
  bool switchIsOn;
  bool stateChanged;
  // One Step non smooth problem

  double* floatWorkingMem = 0;
  int * intWorkingMem = 0;

  int freq = 1000;
  int Nfreq = 0;
  int cmp = 0;

  int NbDataMax = 10000;
  int NData = 0;

  /************************************************************/
  /************************************************************/
  /*Solver options*/
  solverName = &solverEnum;


  int dimX = 1;
  SimpleMatrix * M = 0;
  SimpleMatrix * A = 0;
  SimpleVector* X0 = 0;
  SimpleVector* As = 0;
  SimpleVector* mti = 0;
  SimpleVector* xti = 0;


  X0 = new SimpleVector(dimX);
  xti = new SimpleVector(dimX);
  xti->setValue(0, 0);

  int NBStep = (int) floor(sTf / sStep);

  //  NBStep = 3;
  //*****BUILD THE DYNAMIC SYSTEM
  SP::MyDS aDS ;
  aDS.reset(new MyDS(*xti));

  DynamicalSystemsSet  Inter_DS ;
  Inter_DS.insert(aDS);

  //******BUILD THE RELATION
  SimpleMatrix* C = 0;
  SimpleMatrix* D = 0;
  SimpleMatrix* B = 0;
  SP::elecRelation aR;
  aR.reset(new elecRelation());

  //*****BUILD THE NSLAW
  SP::NonSmoothLaw aNSL;
  aNSL.reset(new MixedComplementarityConditionNSL(sN, sM));
  /*
  if (ACE_FORMULATION==ACE_FORMULATION_SEMI_EXPLICT){
    NSLawSize=m+s;
    aNSL.reset(new MixedComplementarityConditionNSL(m,s));
  }else{
    NSLawSize=m;
    aNSL.reset(new ComplementarityConditionNSL(m));
  }
  */

  //****BUILD THE INTERACTION
  SP::Interaction aI(new Interaction("MLCP", Inter_DS, 1, sNSLawSize, aNSL, aR));
  //****BUILD THE SYSTEM
  SP::NonSmoothDynamicalSystem  aNSDS(new NonSmoothDynamicalSystem(aDS, aI));

  SP::Model  aM(new Model(0, sTf));
  aM->setNonSmoothDynamicalSystemPtr(aNSDS);
  SP::TimeDiscretisation  aTD(new TimeDiscretisation(0, sStep));
  SP::TimeStepping aS(new TimeStepping(aTD));
  aS->setComputeResiduY(true);
  aS->setUseRelativeConvergenceCriteron(false);
  //*****BUILD THE STEP INTEGRATOR
  SP::OneStepIntegrator  aMoreau ;
  aMoreau.reset(new Moreau(aDS, 0.5));
  aS->insertIntegrator(aMoreau);
  //  SP::NonSmoothSolver  mySolver( new NonSmoothSolver((*solverName),iparam,dparam,floatWorkingMem,intWorkingMem));

  //**** BUILD THE STEP NS PROBLEM
  SP::MLCP  aMLCP ;
  aMLCP.reset(new MLCP());
  aS->insertNonSmoothProblem(aMLCP);
  aM->initialize(aS);
  SP::SolverOptions numSolOptions = aMLCP->numericsSolverOptions();
  //  Alloc working mem for the solver
  int aux = mlcp_driver_get_iwork(aMLCP->getNumericsMLCP().get(), &*numSolOptions);
  intWorkingMem = (int*)malloc(aux * sizeof(int));
  numSolOptions->iWork = intWorkingMem;
  aux = mlcp_driver_get_dwork(aMLCP->getNumericsMLCP().get(), &*numSolOptions);
  floatWorkingMem = (double*)malloc(aux * sizeof(double));
  numSolOptions->dWork = floatWorkingMem;
  //  numSolOptions->dparams[0]=1e-12;

  mlcp_driver_init(aMLCP->getNumericsMLCP().get(), &*numSolOptions);
  //      setNumericsVerbose(1);


  SP::SiconosVector  x = aDS->x();
  SP::SiconosVector  y = aI->y(0);
  SP::SiconosVector  lambda = aI->lambda(0);
  ofstream * fout = new ofstream("simu.log");

  unsigned int count = 0; // events counter.
  // do simulation while events remains in the "future events" list of events manager.
  cout << " ==== Start of  simulation : " << NBStep << " steps====" << endl;
#ifdef CLSC_CIRCUIT
#else
  (*fout) << "C_charge " << "V1 " << "R(t)" << endl;
#endif

  for (int k = 0 ; k < NBStep ; k++)
  {
    //      if (cmp==150)
    //        setNumericsVerbose(1);
    //      else if (cmp==151)
    //        setNumericsVerbose(0);
    cout << "..." << cmp << endl;
    cmp++;
    // solve ...
    /*aS->computeOneStep();*/

    aS-> newtonSolve(1e-4, 20);
    aS->nextStep();
    x = aDS->x();
    lambda = aI->lambda(0);
#ifdef CLSC_CIRCUIT

    //    std::cout<<"x="<<x->getValue(0)<<" Is="<<lambda->getValue(0)<<" Id="<<lambda->getValue(1)<<" V3="<<lambda->getValue(2);
    //   std::cout<<" V4="<<lambda->getValue(3)<<" V5="<<lambda->getValue(4)<<" l6="<<lambda->getValue(5)<<" l7="<<lambda->getValue(6);
    //std::cout<<" l8="<<lambda->getValue(7)<<" l9="<<lambda->getValue(8)<<std::endl;
    stateChanged = false;
    if (lambda->getValue(6) > 1)
    {
      if (switchIsOn || k == 0)
      {
        switchIsOn = false;
        stateChanged = true;
      }
    }
    else
    {
      if (!switchIsOn || k == 0)
      {
        switchIsOn = true;
        stateChanged = true;
      }
    }
    if (lambda->getValue(8) > 1)
    {
      if (diodeIsOn || k == 0)
      {
        diodeIsOn = false;
        stateChanged = true;
      }
    }
    else
    {
      if (!diodeIsOn || k == 0)
      {
        diodeIsOn = true;
        stateChanged = true;
      }
    }
    if (stateChanged)
    {
      if (switchIsOn)
        std::cout << "SWTCH=ON";
      else
        std::cout << "SWTCH=OFF";
      if (diodeIsOn)
        std::cout << " DIODE=ON";
      else
        std::cout << " DIODE=OFF";
      std::cout << std::endl;
    }

    (*fout) << cmp << " " << x->getValue(0) << " " << lambda->getValue(0) << " " << lambda->getValue(1) << " " << lambda->getValue(2) << " " << lambda->getValue(3) << " " << lambda->getValue(4) << " " << lambda->getValue(5) << endl;
#else
    (*fout) << cmp << " " << x->getValue(0) << " " << lambda->getValue(0) << " " << lambda->getValue(3) + sR1 << endl;
#endif


  }
  aMLCP->reset();
  delete fout;
  cout << "===== End of simulation. ==== " << endl;
  return 0;

}
