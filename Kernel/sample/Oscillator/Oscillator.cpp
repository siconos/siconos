#include <iostream>
#include "Model.h"
#include "check.h"

#include "DynamicalSystem.h"

#include <libxml/parser.h>
#include "NewSiconosVector.h"
#include "SiconosMatrix.h"
#include "SiconosDOMTreeTools.h"
#include <math.h>
#include <stdio.h>



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

    Model Oscillator("Oscillator.xml");
    //
    cout << "\n *** Oscillator.xml loaded ***" << endl;
    Strategy* s = Oscillator.getStrategy();



    cout << "the strategy will be initialized" << endl;
    s->initialize();
    cout << "\n **** the strategy is ready ****" << endl;


    cout << "Press <<enter>> to continue" << endl;
    getchar();


    TimeDiscretisation* t = s->getTimeDiscretisation();
    int k = t->getK();
    int N = t->getN();

    // Trace Values
    SiconosMatrix m(N + 1, 6);
    //time
    m(k, 0) = k * t->getH();
    // position
    DynamicalSystem* ds = Oscillator.getNonSmoothDynamicalSystem()->getDynamicalSystem(0);



    m(k, 5) = 0.0;
    ds->display();

    while (k < N)
    {
      s->nextStep();
      cout << "NextStep done" << endl;
      k = t->getK();

      cout << "iteration : " << k << endl;
      s->computeFreeState();
      //s->formaliseOneStepNSProblem();
      //s->computeOneStepNSProblem();
      s->updateState();



      // Trace Values
      //time
      m(k, 0) = k * t->getH();
      // position
      //      DynamicalSystem* ball = static_cast<LagrangianNLDS*> (Oscillator.getNonSmoothDynamicalSystem()->getDynamicalSystem(0));
      //      m(k, 1) = (ball->getQ())(0);
      //            // position
      //      m(k, 2) = (ball->getVelocity())(0);
      //
      //      m(k, 3) = (ground->getQ())(0);
      //            // position
      //      m(k, 4) = (ground->getVelocity())(0);
      //
      //      m(k, 5) =   (Oscillator.getNonSmoothDynamicalSystem()->getInteraction(0)->getLambda())(0);

      ds->display();

      //      cout<<"Press a key to begin the computations   <<enter>>"<<endl;
      //      getchar();

    }
    cout << "iterations  done: " << k << endl;

    m.write("result.dat", "ascii");

    Oscillator.saveToXMLFile("./Oscillator.xml.output");

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'main_siconos\'" << endl;
  }
}
