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
  cout << "+                       INRIA - 2004 (c)                        +\n";
  cout << "+                                                               +\n";
  cout << "+---------------------------------------------------------------+\n";
  cout << endl;
  cout << endl;
}



int main(int argc, char* argv[])
{
  cartouche();

  try
  {

    Model sf("./SimpleFriction.xml");

    cout << "\n *** SimpleFriction.xml loaded ***" << endl;
    Strategy* s = sf.getStrategy();

    cout << "the strategy will be initialized" << endl;
    s->initialize();
    cout << "\n **** the strategy is ready ****" << endl;


    TimeDiscretisation* t = s->getTimeDiscretisation();
    int k = t->getK();
    int N = t->getN();

    // Trace Values
    SiconosMatrix m(N + 1, 6);
    //time
    m(k, 0) = k * t->getH();
    // position
    LagrangianDS* ball = static_cast<LagrangianDS*>(sf.getNonSmoothDynamicalSystem()->getDynamicalSystem(0));
    m(k, 1) = (ball->getQ())(0);
    // position
    m(k, 2) = (ball->getQ())(1);
    LagrangianDS* ground = static_cast<LagrangianDS*>(sf.getNonSmoothDynamicalSystem()->getDynamicalSystem(1));
    m(k, 3) = (ground->getQ())(0);
    // position
    m(k, 4) = (ground->getVelocity())(0);
    /*SiconosVector*/
    SimpleVector vv;

    Interaction* I =  sf.getNonSmoothDynamicalSystem()->getInteraction(0);
    I->display();

    //m(k, 5) =   (sf.getNonSmoothDynamicalSystem()->getInteractionOnNumber(0)->getLambda())(0);
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


    while (k < N)
    {
      s->nextStep();
      cout << "NextStep done" << endl;

      k = t->getK();
      cout << "iteration : " << k << endl;

      s->computeFreeState();
      cout << "computeFreeState done" << endl;


      s->formaliseOneStepNSProblem();
      cout << "formaliseOneStepNSProblem done" << endl;


      s->computeOneStepNSProblem();
      cout << "computeOneStepNSProblem done" << endl;

      s->updateState();
      cout << "updateState done" << endl;


      // Trace Values
      //time
      m(k, 0) = k * t->getH();
      // position
      LagrangianDS* ball = static_cast<LagrangianDS*>(sf.getNonSmoothDynamicalSystem()->getDynamicalSystem(0));
      m(k, 1) = (ball->getQ())(0);
      // position
      m(k, 2) = (ball->getQ())(1);

      m(k, 3) = (ground->getQ())(0);
      // position
      m(k, 4) = (ground->getVelocity())(0);

      m(k, 5) = (sf.getNonSmoothDynamicalSystem()->getInteraction(0)->getLambda())(0);


    }

    end = clock();
    rtn = gettimeofday(&tp, NULL);
    t2 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;
    elapsed = t2 - t1;
    elapsed2 = (end - start) / (double)CLOCKS_PER_SEC;

    cout << "time = " << elapsed << " --- cpu time " << elapsed2 << endl;

    cout << "iterations  done: " << k << endl;
    m.write("resultSF.dat", "ascii");

    sf.saveToXMLFile("./SimpleFriction.xml.output");

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/SimpleFriction\'" << endl;
  }
}
