#include <iostream>
#include "Model.h"
#include "check.h"

#include "LagrangianNLDS.h"
#include "LinearSystemDS.h"
#include "LagrangianTIDS.h"
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

    Model uBalls("./UltraBalls.xml");

    cout << "\n *** UltraBalls.xml loaded ***" << endl;
    //getchar();
    Strategy* s = uBalls.getStrategy();

    cout << "the strategy will be initialized" << endl;
    s->initialize();
    cout << "\n **** the strategy is ready ****" << endl;
    //cout<<"K "<<k<<endl;

    //        static_cast<LCP*>(s->getOneStepNSProblem())->display();
    //        cout<<" @ LCP = "<<s->getOneStepNSProblem()<<endl;
    //        cout<<"Press <<enter>> to continue"<<endl;
    //        getchar();


    TimeDiscretisation* t = s->getTimeDiscretisation();
    int k = t->getK();
    int N = t->getN();

    // Trace Values
    SiconosMatrix m(N + 1, 11);

    //time
    m(k, 0) = k * t->getH();

    // position
    LagrangianTIDS* ball = static_cast<LagrangianTIDS*>(uBalls.getNonSmoothDynamicalSystem()->getDynamicalSystem(0));
    m(k, 1) = (ball->getQ())(0);
    m(k, 2) = (ball->getVelocity())(0);

    // position
    LagrangianTIDS* ball2 = static_cast<LagrangianTIDS*>(uBalls.getNonSmoothDynamicalSystem()->getDynamicalSystem(1));
    m(k, 3) = (ball2->getQ())(0);
    m(k, 4) = (ball2->getVelocity())(0);

    // position
    LagrangianTIDS* ball3 = static_cast<LagrangianTIDS*>(uBalls.getNonSmoothDynamicalSystem()->getDynamicalSystem(2));
    m(k, 5) = (ball3->getQ())(0);
    m(k, 6) = (ball3->getVelocity())(0);

    // position
    LagrangianTIDS* ground = static_cast<LagrangianTIDS*>(uBalls.getNonSmoothDynamicalSystem()->getDynamicalSystem(3));
    m(k, 7) = (ground->getQ())(0);
    m(k, 8) = (ground->getVelocity())(0);

    // position
    LagrangianTIDS* ceiling = static_cast<LagrangianTIDS*>(uBalls.getNonSmoothDynamicalSystem()->getDynamicalSystem(4));
    m(k, 9) = (ceiling->getQ())(0);
    m(k, 10) = (ceiling->getVelocity())(0);

    /*SiconosVector*/
    SimpleVector vv;

    Interaction* I =  uBalls.getNonSmoothDynamicalSystem()->getInteraction(0);
    I->display();
    //cout << "vv " << vv << endl;

    //    m(k,7) =0.0;


    while (k < N)
    {
      s->nextStep();
      cout << "NextStep done" << endl;

      k = t->getK();
      cout << "iteration : " << k << endl;

      s->computeFreeState();

      //      cout << "\n **** going to formaliseOneStepNSProblem... ****" << endl;
      //      cout<<" @ LCP = "<<s->getOneStepNSProblem()<<endl;
      //      static_cast<LCP*>(s->getOneStepNSProblem())->getQPtr()->display();
      //      cout<<"Press <<enter>> to continue"<<endl;
      //      getchar();

      s->formaliseOneStepNSProblem();
      //      cout << "\n **** going to computeOneStepNSProblem... ****" << endl;
      //      cout<<" @ LCP = "<<s->getOneStepNSProblem()<<endl;
      //      static_cast<LCP*>(s->getOneStepNSProblem())->getQPtr()->display();
      //      cout<<"Press <<enter>> to continue"<<endl;
      //      getchar();

      //      cout << "\n **** going to computeOneStepNSProblem... ****" << endl;
      //      if( static_cast<LCP*>(s->getOneStepNSProblem())->getNLcp() == 1)
      //      {
      //        cout<<"[----------------------]"<<endl;
      //        cout<<" nLcp == "<<static_cast<LCP*>(s->getOneStepNSProblem())->getNLcp()<<endl;
      //        cout<<" -- q -- M -- w -- z --"<<endl;
      //        static_cast<LCP*>(s->getOneStepNSProblem())->getQ().display();
      //        static_cast<LCP*>(s->getOneStepNSProblem())->getM().display();
      //        static_cast<LCP*>(s->getOneStepNSProblem())->getW().display();
      //        static_cast<LCP*>(s->getOneStepNSProblem())->getZ().display();
      //        cout<<"[----------------------]"<<endl;
      //        //getchar();
      //      }

      s->computeOneStepNSProblem();
      //      cout << "\n **** ... end of computeOneStepNSProblem ****" << endl;
      //      if( static_cast<LCP*>(s->getOneStepNSProblem())->getNLcp() == 1)
      //      {
      //        cout<<"[----------------------]"<<endl;
      //        cout<<" -- q -- M -- w -- z --"<<endl;
      //        cout<<"_ "<<s<<endl;
      //        cout<<"_ "<<s->getOneStepNSProblem()<<endl;
      //        static_cast<LCP*>(s->getOneStepNSProblem())->getQ().display();
      //        static_cast<LCP*>(s->getOneStepNSProblem())->getM().display();
      //        static_cast<LCP*>(s->getOneStepNSProblem())->getW().display();
      //        static_cast<LCP*>(s->getOneStepNSProblem())->getZ().display();
      //        cout<<"[----------------------]"<<endl;
      //        //getchar();
      //      }

      s->updateState();


      // Trace Values
      //time
      m(k, 0) = k * t->getH();

      // position
      LagrangianNLDS* ball = static_cast<LagrangianNLDS*>(uBalls.getNonSmoothDynamicalSystem()->getDynamicalSystem(0));
      m(k, 1) = (ball->getQ())(0);
      m(k, 2) = (ball->getVelocity())(0);

      m(k, 3) = (ball2->getQ())(0);
      m(k, 4) = (ball2->getVelocity())(0);

      m(k, 5) = (ball3->getQ())(0);
      m(k, 6) = (ball3->getVelocity())(0);

      m(k, 7) = (ground->getQ())(0);
      m(k, 8) = (ground->getVelocity())(0);

      m(k, 9) = (ceiling->getQ())(0);
      m(k, 10) = (ceiling->getVelocity())(0);

      //  m(k, 5) = (uBalls.getNonSmoothDynamicalSystem()->getInteraction(0)->getLambda())(0);

      //ball->display();
      //ground->display();

      //      cout<<"Press a key to begin the computations   <<enter>>"<<endl;
      //getchar();

    }
    cout << "iterations  done: " << k << endl;

    m.write("resultUB.dat", "ascii");

    uBalls.saveToXMLFile("./UltraBalls.xml.output");

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/RollingBalls\'" << endl;
  }
}
