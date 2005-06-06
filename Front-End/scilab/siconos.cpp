#include <iostream>
#include "Model.h"
#include "check.h"

#include "LagrangianDS.h"
#include "LinearSystemDS.h"
#include "LagrangianLinearTIDS.h"
#include "LagrangianNonLinearR.h"

#include <libxml/parser.h>
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"
#include "SiconosDOMTreeTools.h"
#include <math.h>
#include <stdio.h>
#include "LCP.h"

#include <sys/time.h>



using namespace std;


void cartouche()
{
  cout << endl;
  cout << endl;
  cout << "+---------------------------------------------------------------+\n";
  cout << "+                                                               +\n";
  cout << "+                        SICONOS / WP2                          +\n";
  cout << "+                       INRIA - 2005 (c)                        +\n";
  cout << "+                                                               +\n";
  cout << "+---------------------------------------------------------------+\n";
  cout << endl;
  cout << endl;
}

//
// TODO: add a structure with array and manage index
//
Model *GLOB_MODEL;
Strategy *GLOB_STRATEGIE;

extern "C" void sicLoadModel(int *ret, char ModelXmlFile[])
{
  try
  {

    GLOB_MODEL = new Model(ModelXmlFile);
    *ret = 0;

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in sicLoadModel" << endl;
  }

}

extern "C" void sicInitStrategy()
{
  try
  {

    // TBD : verify ptr GLOB_STRATEGIE and GLOB_MODEL
    GLOB_STRATEGIE = GLOB_MODEL->getStrategy();
    GLOB_STRATEGIE->initialize();

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in sicInitStrategy" << endl;
  }
}

extern "C" void sicTimeGetH(double *H)
{
  *H = GLOB_STRATEGIE->getTimeDiscretisationPtr()->getH();
}

extern "C" void sicTimeGetN(int *N)
{
  *N = GLOB_STRATEGIE->getTimeDiscretisationPtr()->getN();
}

extern "C" void sicTimeGetK(int *K)
{
  *K = GLOB_STRATEGIE->getTimeDiscretisationPtr()->getK();
}

extern "C" void  sicSTNextStep(int *ret)
{
  GLOB_STRATEGIE->nextStep();

  *ret = 0;
}

extern "C" void sicSTComputeFreeState(int *ret)
{
  GLOB_STRATEGIE->computeFreeState();

  *ret = 0;
}

extern "C" void sicSTformalisePb(int *ret)
{

  try
  {

    GLOB_STRATEGIE->formaliseOneStepNSProblem();
    *ret = 0;

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in sicSTformalisePb" << endl;
  }


}

extern "C" void sicSTcomputePb(int *ret)
{
  GLOB_STRATEGIE->computeOneStepNSProblem();

  *ret = 0;
}

extern "C" void sicSTupdateState(int *ret)
{
  GLOB_STRATEGIE->updateState();

  *ret = 0;
}

extern "C" void sicModelgetQ(double *value, int *index)
{


  try
  {

    LagrangianDS* system = static_cast<LagrangianDS*>(GLOB_MODEL->getNonSmoothDynamicalSystem()->getDynamicalSystem(*index));

    if (system != NULL)
    {
      int size = (system->getQ()).size();
      for (int i = 0; i < size; i++)
        value[i] = system->getQ()(i);
    }

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in sicSTformalisePb" << endl;
  }

}

extern "C" void  simul()
{
  cartouche();

  try
  {

    // SCILAB
    // en global: le modele bouncingBall
    //
    Model bouncingBall("./BouncingBall_TIDS.xml");

    // SCILAB
    // InitStrategy()
    //

    cout << "\n *** BouncingBall.xml loaded ***" << endl;// getchar();
    Strategy* s = bouncingBall.getStrategy();

    cout << "the strategy will be initialized" << endl;
    s->initialize();
    cout << "\n **** the strategy is ready ****" << endl;

    // SCILAB
    //  GetStartTime(), GetEnd,BeginTime()
    //

    TimeDiscretisation* t = s->getTimeDiscretisationPtr();
    int k = t->getK();
    int N = t->getN();

    // Trace Values
    SiconosMatrix m(N + 1, 6);
    //time
    m(k, 0) = k * t->getH();
    // position
    LagrangianDS* ball = static_cast<LagrangianDS*>(bouncingBall.getNonSmoothDynamicalSystem()->getDynamicalSystem(0));
    m(k, 1) = (ball->getQ())(0);

    // position
    m(k, 2) = (ball->getVelocity())(0);
    LagrangianDS* ground = static_cast<LagrangianDS*>(bouncingBall.getNonSmoothDynamicalSystem()->getDynamicalSystem(1));
    m(k, 3) = (ground->getQ())(0);
    // position
    m(k, 4) = (ground->getVelocity())(0);
    /*SiconosVector*/
    SimpleVector vv;

    Interaction* I =  bouncingBall.getNonSmoothDynamicalSystem()->getInteraction(0);
    I->display();

    //m(k, 5) =   (bouncingBall.getNonSmoothDynamicalSystem()->getInteractionOnNumber(0)->getLambda())(0);
    m(k, 5) = 0.0;

    double t1, t2, elapsed;
    struct timeval tp;
    int rtn;

    clock_t start, end;
    double elapsed2;

    start = clock();
    rtn = gettimeofday(&tp, NULL);
    clock_t t11 = clock();
    t1 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;

    // SCILAB
    //  boucle scilab  NextStep()  ComputeFreeState(),...
    //

    while (k < N)
    {
      s->nextStep();
      cout << "NextStep done" << endl;
      k = t->getK();

      cout << "iteration : " << k << endl;
      s->computeFreeState();



      s->formaliseOneStepNSProblem();



      s->computeOneStepNSProblem();


      s->updateState();



      // Trace Values
      //time
      m(k, 0) = k * t->getH();
      // position
      LagrangianDS* ball = static_cast<LagrangianDS*>(bouncingBall.getNonSmoothDynamicalSystem()->getDynamicalSystem(0));
      m(k, 1) = (ball->getQ())(0);
      // position
      m(k, 2) = (ball->getVelocity())(0);

      m(k, 3) = (ground->getQ())(0);
      // position
      m(k, 4) = (ground->getVelocity())(0);

      m(k, 5) = (bouncingBall.getNonSmoothDynamicalSystem()->getInteraction(0)->getLambda())(0);


    }

    end = clock();
    rtn = gettimeofday(&tp, NULL);
    t2 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;
    elapsed = t2 - t1;
    elapsed2 = (end - start) / (double)CLOCKS_PER_SEC;

    cout << "time = " << elapsed << " --- cpu time " << elapsed2 << endl;

    cout << "iterations  done: " << k << endl;
    m.write("result.dat", "ascii");

    bouncingBall.saveToXMLFile("./BouncingBall_TIDS.xml.output");

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/BouncingBall\'" << endl;
  }
}
