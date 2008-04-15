/ Header files
#include "Model.h"
#include "FirstOrderLinearTIDS.h"
#include "LinearTIR.h"
#include "ComplementarityConditionNSL.h"
#include "TimeStepping.h"
#include "ioMatrix.h"
#include <string>
#include <iostream>

using namespace std;

// main program
int main(int argc, char* argv[])
{
  // Exception handling
  try
  {
    // == User-defined parameters ==
    unsigned int ndof = 2;  // number of degrees of freedom of your system
    double t0 = 0.0;
    double T = 5e-3;        // Total simulation time
    double h = 1.0e-6;      // Time step
    string solverName = "Lemke"; // non smooth problem solver algo name.

    // ================= Creation of the model =======================

    // == Creation of the NonSmoothDynamicalSystem ==
    // DynamicalSystem(s)
    SiconosMatrix * A = new SimpleMatrix(ndof, ndof); // All components of A are automatically set to 0.
    (*A)(0, 0) = 1.0;
    (*A)(0, 1) = 1.0;
    (*A)(1, 0) = 3.0;
    (*A)(1, 1) = 1.0;
    SiconosVector * x0 = new SimpleVector(ndof);
    (*x0)(0) = Vinit;
    DynamicalSystem * process  = new FirstOrderLinearDS(1, *x0, *A);

    process->setComputeBFunction("ObserverPlugin.so", "computeU");

    DynamicalSystemsSet allDS;
    allDS.insert(process);



    // Relations
    unsigned int ninter = 1; // dimension of your Interaction = size of y and lambda vectors
    unsigned int noutput = 1;
    unsigned int ncontrol = 1;

    SiconosMatrix * B = new SimpleMatrix(ndof, ninter);
    (*B)(0, 0) = -1.0;
    (*B)(1, 0) = 1.0;
    SiconosMatrix * C = new SimpleMatrix(ninter, ndof);
    (*C)(0, 0) = -1.0;
    (*C)(0, 1) = 1.0;
    LinearTIR * myProcessRelation = new LinearTIR(*C, *B);

    myProcessrelation->setComputeEFunction("ObserverPlugin.so", "computeE");


    SiconosMatrix* D = new SimpleMatrix(ninter, ninter);
    (*D)(0, 0) = 1.0;
    myProcessRelation->setDPtr(D);
    SiconosMatrix * F = new SimpleMatrix(ninter, ncontrol);
    (*F)(0, 0) = 1.0;
    SiconosMatrix * E = new SimpleMatrix(ndof, ncontrol);
    (*E)(0, 0) = 1.0;
    (*E)(1, 0) = 2.0;
    SiconosMatrix * L = new SimpleMatrix(ndof, noutput);
    (*L)(0, 0) = 1.0;
    (*L)(1.0) = 1.0;
    SiconosMatrix * G = new SimpleMatrix(noutput, ndof);
    (*G)(0, 0) = 2.0;
    (*G)(0, 1) = 2.0;
    SiconosMatrix * hatA = new SimpleMatrix(ndof, ndof);

    (*hatA) = *A - *L *  *G;

    DynamicalSystem * Observer = new FirstOrderLinearDS(1, *x0, *hatA);


    LinearTIR * myObserverRelation = new LinearTIR(*C, *B);

    myObserverRelation->setComputeEFunction("ObserverPlugin.so", "computeE");

    myObserverRelation->setDPtr(D);



    // NonSmoothLaw
    unsigned int nslawSize = 1;
    NonSmoothLaw * myNslaw = new ComplementarityConditionNSL(nslawSize);



    // Choose a name and a number for your Interaction
    string nameInter = "processInteraction";
    unsigned int numInter = 1;
    Interaction* myProcessInteraction = new Interaction(nameInter, process, numInter, ninter, myNslaw, myProcessRelation);
    InteractionsSet allInteractions;
    allInteractions.insert(myProcessInteraction);


    string nameInter = "observerInteraction";
    unsigned int numInter = 2;
    Interaction* myObserverInteraction = new Interaction(nameInter, observer, numInter, ninter, myNslaw, myObserverRelation);

    allInteractions.insert(myObserverInteraction);

    // NonSmoothDynamicalSystem
    NonSmoothDynamicalSystem* myNSDS = new NonSmoothDynamicalSystem(allDS, allInteractions);



    // Model
    Model * ObserverLCS = new Model(t0, T);
    Observer->setNonSmoothDynamicalSystemPtr(myNSDS);
    // == Creation of the Simulation ==
    Simulation * s = new TimeStepping(ObserverLCS);
    // TimeDiscretisation
    TimeDiscretisation * td = new TimeDiscretisation(h, s);
    // OneStepIntegrator
    double theta = 0.5;
    // One Step Integrator
    Moreau* myIntegrator = new Moreau(allDS, theta, s);
    // OneStepNSProblem
    // One Step non smooth problem
    OneStepNSProblem* myLCP = new LCP(s, "LCP", solverName, 101, 0.0001, "max", 0.6);
    // ================================= Computation =================================

    // --- Initialisation of the simulation ---
    s->initialize();

    int k = td->getK(); // Current step
    int N = td->getNSteps(); // Number of time steps
    unsigned int outputSize = 7; // number of required data
    SimpleMatrix dataPlot(N, outputSize);
    SimpleVector * processLambda = myProcessInteraction->getLambdaPtr(0);
    SimpleVector * observerLambda = myObserverInteraction->getLambdaPtr(0);
    // We get values for the initial time step:
    // time
    dataPlot(k, 0) = k * h;
    dataPlot(k, 1) = (observer->getX())(0);
    dataPlot(k, 2) = (observer->getX())(1);
    dataPlot(k, 3) = (process->getX())(0);
    dataPlot(k, 4) = (process->getX())(1);
    dataPlot(k, 5) = processLambda(0);
    dataPlot(k, 6) = observerLambda(0);


    // Simulation loop
    while (k < N - 1)
    {
      s->nextStep();
      // get current time step
      k = td->getK();
      s->computeOneStep();
      dataPlot(k, 0) = k * h;
      dataPlot(k, 1) = (observer->getX())(0);
      dataPlot(k, 2) = (observer->getX())(1);
      dataPlot(k, 3) = (process->getX())(0);
      dataPlot(k, 4) = (process->getX())(1);
      dataPlot(k, 5) = processLambda(0);
      dataPlot(k, 6) = observerLambda(0);
    }

    // Write the results into the file "ObserverLCS.dat"
    ioMatrix io("ObserverLCS.dat", "ascii");
    io.write(dataPlot);

    // --- Time loop ---

  }






  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in ObserverLCS.cpp" << endl;
  }
}
