
#include "SiconosKernel.hpp"
#include "adjointInput.h"
#include "myDS.h"
#include <stdio.h>
#include <stdlib.h>

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

static int sNSLawSize = 2;

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
  IntParameters iparam(10);
  DoubleParameters dparam(10);

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
  iparam[5] = 5; // nb config
  iparam[0] = 0; // verbose
  dparam[0] = 1e-11; // Tolerance
  solverName = &solverEnum;


  int dimX = 4;
  SimpleMatrix * M = 0;
  SimpleMatrix * A = 0;
  SimpleVector* X0 = 0;
  SimpleVector* As = 0;
  SimpleVector* mti = 0;
  SimpleVector* xti = 0;


  X0 = new SimpleVector(dimX);
  xti = new SimpleVector(dimX);
  xti->setValue(0, 0);

  double sT = 1e-2;
  double sStep = 0.1e-6;
  int NBStep = (int) floor(sT / sStep);

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
  SP::adjointInput aR;
  aR.reset(new adjointInput());

  int sN = 2;
  int sM = 0;


  //*****BUILD THE NSLAW
  SP::NonSmoothLaw aNSL;
  aNSL.reset(new MixedComplementarityConditionNSL(sN, sM));





  //****BUILD THE INTERACTION
  SP::Interaction aI(new Interaction("MLCP", Inter_DS, 1, sNSLawSize, aNSL, aR));
  //****BUILD THE SYSTEM
  SP::NonSmoothDynamicalSystem  aNSDS(new NonSmoothDynamicalSystem(aDS, aI));

  SP::Model  aM(new Model(0, sT));
  aM->setNonSmoothDynamicalSystemPtr(aNSDS);
  SP::TimeDiscretisation  aTD(new TimeDiscretisation(0, sStep));
  SP::TimeStepping aS(new TimeStepping(aTD));
  aS->setComputeResiduY(true);
  aS->setUseRelativeConvergenceCriteron(false);
  //*****BUILD THE STEP INTEGRATOR
  SP::OneStepIntegrator  aMoreau ;
  aMoreau.reset(new Moreau(aDS, 0.5));
  aS->insertIntegrator(aMoreau);
  SP::NonSmoothSolver  mySolver(new NonSmoothSolver((*solverName), iparam, dparam, floatWorkingMem, intWorkingMem));

  //**** BUILD THE STEP NS PROBLEM
  SP::MLCP  aMLCP ;
  aMLCP.reset(new MLCP(mySolver, "MLCP"));
  aS->insertNonSmoothProblem(aMLCP);
  aM->initialize(aS);
  //  Alloc working mem for the solver
  int aux = mlcp_driver_get_iwork(aMLCP->getNumericsMLCP().get(), mySolver->numericsSolverOptions().get());
  intWorkingMem = (int*)malloc(aux * sizeof(int));
  mySolver->numericsSolverOptions()->iWork = intWorkingMem;
  aux = mlcp_driver_get_dwork(aMLCP->getNumericsMLCP().get(), mySolver->numericsSolverOptions().get());
  floatWorkingMem = (double*)malloc(aux * sizeof(double));
  mySolver->numericsSolverOptions()->dWork = floatWorkingMem;

  mlcp_driver_init(aMLCP->getNumericsMLCP().get(), mySolver->numericsSolverOptions().get());
  //      setNumericsVerbose(1);


  SP::SiconosVector  x = aDS->x();
  SP::SiconosVector  y = aI->y(0);
  SP::SiconosVector  lambda = aI->lambda(0);
  ofstream * fout = new ofstream("simu.log");

  unsigned int count = 0; // events counter.
  // do simulation while events remains in the "future events" list of events manager.
  cout << " ==== Start of  simulation : " << NBStep << " steps====" << endl;

  for (int k = 0 ; k < NBStep ; k++)
  {

    cout << " ==== Step: " << endl;
    //      if (cmp==150)
    setNumericsVerbose(1);
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


  }
  aMLCP->reset();
  delete fout;
  cout << "===== End of simulation. ==== " << endl;
  return 0;

}
