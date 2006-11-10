#include <iostream>
#include <math.h>
#include "Model.h"
#include "check.h"
using namespace std;

int main(int argc, char* argv[])
{
  try
  {
    // initialise the model with an xml file
    Model bouncingBall("./BouncingBall_TIDS.xml");
    cout << "\n *** BouncingBall.xml loaded ***" << endl;

    // get the strategy
    Strategy* strategy = bouncingBall.getStrategy();

    // compute number of time steps
    double time = bouncingBall.getT0();
    double timeStop = bouncingBall.getFinalT();
    int N = t->getN();
    if (N == 1)
    {
      N = floor((timeStop - time) / t->getH());
      t->setN(N);
    }
    // k is the current step
    int k = t->getK();

    strategy->initialize();
    cout << "\n **** the strategy is ready ****" << endl;

    // loop of simulation for one step
    while (k < N)
    {
      strategy->nextStep();
      k = t->getK();
      strategy->computeFreeState();
      strategy->formaliseOneStepNSProblem();
      strategy->computeOneStepNSProblem();
      strategy->updateState();
    }
    cout << "Simulation done." << endl;

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
