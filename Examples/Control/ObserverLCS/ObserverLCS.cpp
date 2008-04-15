
#include "SiconosKernel.h"

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
    double T = 1;        // Total simulation time
    double h = 1.0e-3;      // Time step
    double Vinit = 1.0;


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
    FirstOrderLinearDS * process  = new FirstOrderLinearDS(0, *x0, *A);

    process->setComputeBFunction("ObserverLCSPlugin.so", "computeU");


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
    FirstOrderLinearR * myProcessRelation = new FirstOrderLinearR(C, B);

    myProcessRelation->setComputeEFunction("ObserverLCSPlugin.so", "computeE");


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
    (*L)(1, 0) = 1.0;
    SiconosMatrix * G = new SimpleMatrix(noutput, ndof);
    (*G)(0, 0) = 2.0;
    (*G)(0, 1) = 2.0;
    SiconosMatrix * hatA = new SimpleMatrix(ndof, ndof);




    (*hatA) = (*A)     -   prod(*L, *G);

    FirstOrderLinearDS* observer = new FirstOrderLinearDS(1, *x0, *hatA);
    //observer->setComputeBFunction("ObserverLCSPlugin.so","computebobserver");
    SimpleVector * z = new SimpleVector(1);
    observer->setZPtr(z);


    allDS.insert(observer);


    FirstOrderLinearR * myObserverRelation = new FirstOrderLinearR(C, B);

    myObserverRelation->setComputeEFunction("ObserverLCSPlugin.so", "computeE");

    myObserverRelation->setDPtr(D);



    // NonSmoothLaw
    unsigned int nslawSize = 1;
    NonSmoothLaw * myNslaw = new ComplementarityConditionNSL(nslawSize);



    // Choose a name and a number for your Interaction
    string nameInter = "processInteraction";
    unsigned int numInter = 1;
    cout << "tutu " << endl;

    allDS.display();


    DynamicalSystemsSet dsObserverConcerned;
    dsObserverConcerned.insert(allDS.getPtr(0));
    cout << "tutu " << endl;

    Interaction* myProcessInteraction = new Interaction(nameInter, dsObserverConcerned, numInter, ninter, myNslaw, myProcessRelation);

    InteractionsSet allInteractions;
    allInteractions.insert(myProcessInteraction);
    cout << "tutu " << endl;


    string nameInter2 = "observerInteraction";
    unsigned int numInter2 = 2;

    DynamicalSystemsSet dsProcessConcerned;
    dsProcessConcerned.insert(allDS.getPtr(1));

    Interaction* myObserverInteraction = new Interaction(nameInter2, dsProcessConcerned, numInter2, ninter, myNslaw, myObserverRelation);

    //allInteractions.insert(myObserverInteraction);




    // NonSmoothDynamicalSystem
    NonSmoothDynamicalSystem* myNSDS = new NonSmoothDynamicalSystem(allDS, allInteractions);

    cout << "toto " << endl;


    // Model
    Model * ObserverLCS = new Model(t0, T);
    ObserverLCS->setNonSmoothDynamicalSystemPtr(myNSDS);
    // TimeDiscretisation
    TimeDiscretisation * td = new TimeDiscretisation(h, ObserverLCS);
    // == Creation of the Simulation ==
    TimeStepping * s = new TimeStepping(td);


    // OneStepIntegrator
    double theta = 0.5;
    // One Step Integrator
    Moreau* myIntegrator = new Moreau(allDS, theta, s);
    // OneStepNSProblem



    // One Step non smooth problem
    IntParameters iparam(5);
    iparam[0] = 1000; // Max number of iteration
    DoubleParameters dparam(5);
    dparam[0] = 1e-5; // Tolerance
    string solverName = "Lemke" ;
    NonSmoothSolver * mySolver = new NonSmoothSolver(solverName, iparam, dparam);
    LCP * osnspb = new LCP(s, mySolver);

    // ================================= Computation =================================

    // --- Initialisation of the simulation ---
    s->initialize();

    *z = *(myProcessInteraction->getYPtr(0)->getVectorPtr(0));
    int k = 0; // Current step
    int N = td->getNSteps(); // Number of time steps
    unsigned int outputSize = 7; // number of required data
    SimpleMatrix dataPlot(N, outputSize);
    SiconosVector * processLambda = myProcessInteraction->getLambdaPtr(0);
    SiconosVector * observerLambda = myObserverInteraction->getLambdaPtr(0);
    // We get values for the initial time step:
    // time
    dataPlot(k, 0) = s->getNextTime();
    dataPlot(k, 1) = (observer->getX())(0);
    dataPlot(k, 2) = (observer->getX())(1);
    dataPlot(k, 3) = (process->getX())(0);
    dataPlot(k, 4) = (process->getX())(1);
    dataPlot(k, 5) = (*processLambda)(0);
    dataPlot(k, 6) = (*observerLambda)(0);


    // Simulation loop
    while (k < N - 1)
    {
      k++;

      s->nextStep();
      // get current time step

      s->computeOneStep();
      *z = *(myProcessInteraction->getYPtr(0)->getVectorPtr(0));
      dataPlot(k, 0) = s->getNextTime();;
      dataPlot(k, 1) = (observer->getX())(0);
      dataPlot(k, 2) = (observer->getX())(1);
      dataPlot(k, 3) = (process->getX())(0);
      dataPlot(k, 4) = (process->getX())(1);
      dataPlot(k, 5) = (*processLambda)(0);
      dataPlot(k, 6) = (*observerLambda)(0);
    }

    // Write the results into the file "ObserverLCS.dat"
    ioMatrix io("ObserverLCS.dat", "ascii");
    io.write(dataPlot, "noDim");

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
