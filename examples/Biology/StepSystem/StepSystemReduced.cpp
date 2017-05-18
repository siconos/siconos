#include "SiconosKernel.hpp"
#include "const.h"
#include "NonlinearRelationReduced2.h"
#include "myDS.h"
#include <stdio.h>
#include <stdlib.h>

//#define SICONOS_DEBUG
using namespace std;

/*main program*/

int main(int argc, char *argv[])
{
  //printf("argc %i\n", argc);
  int cmp=0;
  int dimX = 2;

  char  filename[50] = "simu.";

//***** Set the initial condition
  SP::SiconosVector xti(new SiconosVector(dimX));
  if (argc==1)
  {
    xti->setValue(0,1);
    xti->setValue(1,6);
    strncpy(&filename[5],"1.6.log",7);
  }
  else if (argc==3)
  {
    //printf("argv[0] %s\n", argv[0]);
    printf("xti(0) is set to %f\n", atof(argv[1]));
    printf("xti(1) is set to %f\n", atof(argv[2]));

    xti->setValue(0,atof(argv[1]));
    xti->setValue(1,atof(argv[2]));
    int sizeofargv1 = strlen(argv[1]);
    // printf("sizeofargv1 %i\n",sizeofargv1);
    strncpy(&filename[5],argv[1],sizeofargv1);
    int sizeofargv2 = strlen(argv[2]);
    //printf("sizeofargv2 %i\n",sizeofargv2);
    strncpy(&filename[5+sizeofargv1],".",1);

    strncpy(&filename[5+sizeofargv1+1],argv[2],sizeofargv2);
    strncpy(&filename[5+sizeofargv1+sizeofargv2+1],".log",4);

    printf("Output is written in filename %s\n",  filename);
  }
  else
  {
    cout << "wrong  number of arguments = " << argc << endl;
  }

  int NBStep = (int) floor(sTf/sStep);
  // NBStep =1;
//*****BUILD THE DYNAMICAL SYSTEM

  SP::MyDS aDS ;
  aDS.reset(new MyDS(xti));

//******BUILD THE RELATION
  SP::NonlinearRelationReduced2 aR;
  aR.reset(new NonlinearRelationReduced2());

//*****BUILD THE NSLAW
  SP::NonSmoothLaw aNSL;
  RelayNSL * rNSL = new RelayNSL(sNSLawSize);
  aNSL.reset(rNSL); //nr de lambda din pb LCP
  rNSL->setLb(-0.0);
  rNSL->setUb(1.0);
// unsigned int sNSLawSize = 4;

//****BUILD THE INTERACTION
  SP::Interaction aI(new Interaction(aNSL,aR));

//****BUILD THE SYSTEM
  SP::Model  aM(new Model(0,sTf));
  aM->nonSmoothDynamicalSystem()->insertDynamicalSystem(aDS);
  aM->nonSmoothDynamicalSystem()->link(aI,aDS);

// -- (1) OneStepIntegrators --
  SP::OneStepIntegrator  aEulerMoreauOSI ;
  aEulerMoreauOSI.reset(new EulerMoreauOSI(0.5));

// -- (2) Time discretisation --
  SP::TimeDiscretisation  aTD(new TimeDiscretisation(0,sStep));

// -- (3) Non smooth problem
  SP::Relay osnspb(new Relay(SICONOS_RELAY_LEMKE));
  osnspb->numericsSolverOptions()->dparam[0]=1e-08;
  osnspb->numericsSolverOptions()->iparam[0]=0;  // Multiple solutions 0 or 1
//  osnspb->numericsSolverOptions()->iparam[3]=48;

  osnspb->setNumericsVerboseMode(0);

// -- (4) Simulation setup with (1) (2) (3)
  SP::TimeStepping aS(new TimeStepping(aTD,aEulerMoreauOSI,osnspb));
  aS->setComputeResiduY(true);
  aS->setUseRelativeConvergenceCriteron(false);

// Initialization
  printf("-> Initialisation \n");
  aM->setSimulation(aS);
  aM->initialize();
  printf("-> End of initialization \n");

// BUILD THE STEP INTEGRATOR

  SP::SiconosVector  x = aDS->x();
  SP::SiconosVector  vectorfield = aDS->rhs();
  SP::SiconosVector  y = aI->y(0);
  SP::SiconosVector  lambda = aI->lambda(0);

  unsigned int outputSize = 9;
  SimpleMatrix dataPlot(NBStep+1,outputSize);





  cout << "=== Start of simulation: "<<NBStep<<" steps ===" << endl;

  printf("=== Start of simulation: %d steps ===  \n", NBStep);

  dataPlot(0, 0) = aM->t0();
  dataPlot(0,1) = x->getValue(0);
  dataPlot(0,2) = x->getValue(1);
  dataPlot(0, 3) = lambda->getValue(0);
  dataPlot(0, 7) = vectorfield->getValue(0);
  dataPlot(0, 8) = vectorfield->getValue(1);


  boost::progress_display show_progress(NBStep);

  boost::timer time;
  time.restart();

  for(int k = 0 ; k < NBStep ; k++)
  {
#ifdef SICONOS_DEBUG
    std::cout<<"-> Running step:"<<k<<std::endl;
#endif
    cmp++;

    aS->newtonSolve(1e-4, 200);
 
    dataPlot(cmp, 0) = aS->nextTime();
    dataPlot(cmp, 1) = x->getValue(0);
    dataPlot(cmp, 2) = x->getValue(1);
    dataPlot(cmp, 3) = lambda->getValue(0);

    aDS->computeRhs(aS->nextTime(),true);

    if (cmp==1) // tricks just for display to avoid the computation of the initial Rhs
    {
      dataPlot(cmp-1, 7) = vectorfield->getValue(0);
      dataPlot(cmp-1, 8) = vectorfield->getValue(1);
    }

    dataPlot(cmp, 7) = vectorfield->getValue(0);
    dataPlot(cmp, 8) = vectorfield->getValue(1);
    ++show_progress;
    aS->nextStep();
   
    // (*fout)<<cmp<<" "<<x->getValue(0)<<" "<<x->getValue(1)<<" "<<lambda->getValue(0)<<" "<<lambda->getValue(1)<<" "<<lambda->getValue(2)<<" "<<lambda->getValue(3)<<endl;
  }


  dataPlot.resize(cmp,outputSize);
  ioMatrix::write(filename, "ascii", dataPlot, "noDim");
  cout << "=== End of simulation. === " << endl;
  return 0;

}


