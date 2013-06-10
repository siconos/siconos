
#include "SiconosKernel.hpp"
#include "adjointInput.hpp"
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
  string * solverName = 0;
  bool diodeIsOn;
  bool switchIsOn;
  bool stateChanged;

  double* floatWorkingMem = 0;
  int * intWorkingMem = 0;

  int freq = 1000;
  unsigned int Nfreq = 0;
  int cmp = 0;
  int cmp10 = 0;

  unsigned int NbDataMax = 10000;
  unsigned int NData = 0;

  /************************************************************/
  /************************************************************/


  int dimX = 4;
  SimpleMatrix * M = 0;
  SimpleMatrix * A = 0;
  SiconosVector* x0 = 0;
  SiconosVector* As = 0;
  SiconosVector* mti = 0;


  x0 = new SiconosVector(dimX);
  //Point de départ hors arc singulier
  x0->setValue(0, 3.4999939172);
  x0->setValue(1, -2.2788416237);
  x0->setValue(2, 1.1935988302);
  x0->setValue(3, -0.6365413023);



  double sT = 1;
  double sStep = 0.1e-4;
  unsigned int NBStep = floor(sT / sStep);
  //NBStep =2;

  //  NBStep = 3;
  //*****BUILD THE DYNAMIC SYSTEM
  SP::MyDS aDS ;
  aDS.reset(new MyDS(*x0));

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
  aNSL.reset(new ComplementarityConditionNSL(sN));





  //****BUILD THE INTERACTION
  SP::Interaction aI(new Interaction(sNSLawSize, aNSL, aR));
  //****BUILD THE SYSTEM
  SP::Model  aM(new Model(0, sT));
  aM->nonSmoothDynamicalSystem()->insertDynamicalSystem(aDS);
  aM->nonSmoothDynamicalSystem()->link(aI,aDS);
  SP::TimeDiscretisation  aTD(new TimeDiscretisation(0, sStep));
  SP::TimeStepping aS(new TimeStepping(aTD));
  aS->setComputeResiduY(true);
  aS->setUseRelativeConvergenceCriteron(false);
  //*****BUILD THE STEP INTEGRATOR
  SP::OneStepIntegrator  aMoreau ;
  aMoreau.reset(new Moreau(0.5));
  aMoreau->insertDynamicalSystem(aDS);
  aS->insertIntegrator(aMoreau);

  //**** BUILD THE STEP NS PROBLEM
  SP::LCP  aLCP ;


  aLCP.reset(new LCP(SICONOS_LCP_ENUM));

  aS->insertNonSmoothProblem(aLCP);
  aM->initialize(aS);

  setNumericsVerbose(0);


  SP::SiconosVector  x = aDS->x();
  SP::SiconosVector  y = aI->y(0);
  SP::SiconosVector  lambda = aI->lambda(0);
  ofstream * fout = new ofstream("simu.log");


  unsigned int outputSize = 10; // number of required data
  SimpleMatrix dataPlot(NBStep, outputSize);

  SP::SiconosVector z = aDS->x();
  SP::SiconosVector lambdaOut = aI->lambda(0);
  SP::SiconosVector yOut = aI->y(0);

  dataPlot(0, 0) = aM->t0(); // Initial time of the model
  dataPlot(0, 1) = (*z)(0);
  dataPlot(0, 2) = (*z)(1);
  dataPlot(0, 3) = (*z)(2);
  dataPlot(0, 4) = (*z)(3);
  dataPlot(0, 5) = (*lambdaOut)(0);
  dataPlot(0, 6) = (*lambdaOut)(1);
  dataPlot(0, 7) = (*yOut)(0);
  dataPlot(0, 8) = (*yOut)(1);




  unsigned int count = 0; // events counter.
  // do simulation while events remains in the "future events" list of events manager.
  cout << " ==== Start of  simulation : " << NBStep << " steps====" << endl;

  for (int k = 0 ; k < NBStep ; k++)
  {

    cout << " ==== Step: " << endl;
    //      if (cmp==150)
    // setNumericsVerbose(à);
    //      else if (cmp==151)
    setNumericsVerbose(0);

    cout << "..." << cmp << endl; // a garder


    cmp++;
    // solve ...
    /*aS->computeOneStep();*/

    aS-> newtonSolve(1e-4, 200);
    aS->nextStep();
    x = aDS->x();
    lambda = aI->lambda(0);
    dataPlot(k, 0) = aS->nextTime(); // Initial time of the model
    dataPlot(k, 1) = (*x)(0);
    dataPlot(k, 2) = (*x)(1);
    dataPlot(k, 3) = (*x)(2);
    dataPlot(k, 4) = (*x)(3);
    dataPlot(k, 5) = (*lambda)(0);
    dataPlot(k, 6) = (*lambda)(1);
    dataPlot(k, 7) = (*yOut)(0);
    dataPlot(k, 8) = (*yOut)(1);


  }


  // --- Output files ---
  cout << "====> Output file writing ..." << endl;
  ioMatrix::write("OptimalControl.dat", "ascii", dataPlot, "noDim");

  delete fout;
  cout << "===== End of simulation. ==== " << endl;
  return 0;

}
