
#include "NumericsVerbose.h"
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


int main()
{
  int cmp = 0;

  /************************************************************/
  /************************************************************/


  int dimX = 4;

  SP::SiconosVector x0(new SiconosVector(dimX));

  //Point de départ hors arc singulier
  x0->setValue(0, 3.4999939172);
  x0->setValue(1, -2.2788416237);
  x0->setValue(2, 1.1935988302);
  x0->setValue(3, -0.6365413023);



  double sT = 10;
  double sStep = 2e-3;
  unsigned int NBStep = floor(sT / sStep);
  //NBStep =2;

  //  NBStep = 3;
  //*****BUILD THE DYNAMIC SYSTEM
  SP::MyDS aDS ;
  aDS.reset(new MyDS(x0));

  //******BUILD THE RELATION
  SP::adjointInput aR;
  aR.reset(new adjointInput());

  int sN = 2;

  //*****BUILD THE NSLAW
  SP::NonSmoothLaw aNSL;
  aNSL.reset(new ComplementarityConditionNSL(sN));





  //****BUILD THE INTERACTION
  SP::Interaction aI(new Interaction(aNSL, aR));
  //****BUILD THE SYSTEM
  SP::NonSmoothDynamicalSystem  aM(new NonSmoothDynamicalSystem(0, sT));
  aM->insertDynamicalSystem(aDS);
  aM->link(aI,aDS);
  SP::TimeDiscretisation  aTD(new TimeDiscretisation(0, sStep));
  SP::TimeStepping aS(new TimeStepping(aM, aTD));
  aS->setComputeResiduY(true);
  aS->setComputeResiduR(true);
  aS->setUseRelativeConvergenceCriteron(false);
  aS->setNewtonTolerance(1.1e-11);
  aS->setNewtonMaxIteration(50);
  //*****BUILD THE STEP INTEGRATOR
  SP::OneStepIntegrator  aEulerMoreauOSI ;
  aEulerMoreauOSI.reset(new EulerMoreauOSI(0.5));
  aS->insertIntegrator(aEulerMoreauOSI);

  //**** BUILD THE STEP NS PROBLEM
  SP::LCP  aLCP ;


  aLCP.reset(new LCP(SICONOS_LCP_ENUM));
//  aLCP.reset(new LCP(SICONOS_LCP_NEWTONFB));

  aS->insertNonSmoothProblem(aLCP);

  numerics_set_verbose(0);


  SP::SiconosVector  x = aDS->x();
  SP::SiconosVector  y = aI->y(0);
  SP::SiconosVector  lambda = aI->lambda(0);


  unsigned int outputSize = 9; // number of required data
  SimpleMatrix dataPlot(NBStep+10, outputSize);

  SP::SiconosVector z = aDS->x();

  dataPlot(0, 0) = aM->t0(); // Initial time of the model
  dataPlot(0, 1) = (*z)(0);
  dataPlot(0, 2) = (*z)(1);
  dataPlot(0, 3) = (*z)(2);
  dataPlot(0, 4) = (*z)(3);
  dataPlot(0, 5) = (*lambda)(0);
  dataPlot(0, 6) = (*lambda)(1);
  dataPlot(0, 7) = (*y)(0);
  dataPlot(0, 8) = (*y)(1);

  // do simulation while events remains in the "future events" list of events manager.
  cout << " ==== Start of  simulation : " << NBStep << " steps====" << endl;
  boost::progress_display show_progress(NBStep);
  boost::timer time;
  time.restart();
  unsigned int k = 0;
  while (aS->hasNextEvent())
  {
    k++;
    //      if (cmp==150)
    // numerics_set_verbose(à);
    //      else if (cmp==151)
    numerics_set_verbose(0);
    ++show_progress;

    cmp++;

    // solve ...
    aS->computeOneStep();
    x = aDS->x();
    lambda = aI->lambda(0);
    dataPlot(k, 0) = aS->nextTime(); // Initial time of the model
    dataPlot(k, 1) = (*x)(0);
    dataPlot(k, 2) = (*x)(1);
    dataPlot(k, 3) = (*x)(2);
    dataPlot(k, 4) = (*x)(3);
    dataPlot(k, 5) = (*lambda)(0);
    dataPlot(k, 6) = (*lambda)(1);
    dataPlot(k, 7) = (*y)(0);
    dataPlot(k, 8) = (*y)(1);
    aS->nextStep();


  }

  cout << "===== End of simulation. ==== " << endl;
  dataPlot.resize(k+1, 9);

  // --- Output files ---
  cout << "====> Output file writing ..." << endl;
  ioMatrix::write("OptimalControl.dat", "ascii", dataPlot, "noDim");

  double error=0.0, eps=1e-08;
  if ((error=ioMatrix::compareRefFile(dataPlot, "OptimalControl.ref", eps)) >= 0.0
      && error > eps)
    return 1;

  return 0;

}
