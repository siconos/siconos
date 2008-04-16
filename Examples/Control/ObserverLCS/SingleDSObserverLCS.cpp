
#include "SiconosKernel.h"

using namespace std;

// main program
int main(int argc, char* argv[])
{
  // Exception handling
  try
  {
    // == User-defined parameters ==
    unsigned int ndof = 4;  // number of degrees of freedom of your system
    double t0 = 0.0;
    double T = 2;        // Total simulation time
    double h = 1.0e-4;      // Time step
    double Vinit = 10.0;
    unsigned int noutput = 1;
    unsigned int ncontrol = 1;



    // ================= Creation of the model =======================

    // == Creation of the NonSmoothDynamicalSystem ==
    // DynamicalSystem(s)
    SiconosMatrix * A = new SimpleMatrix(2, 2); // All components of A are automatically set to 0.
    (*A)(0, 0) = 1.0;
    (*A)(0, 1) = 1.0;
    (*A)(1, 0) = 3.0;
    (*A)(1, 1) = 1.0;

    SiconosMatrix * TildeA = new SimpleMatrix(ndof, ndof); // All components of A are automatically set to 0.
    (*TildeA)(0, 0) = (*A)(0, 0);
    (*TildeA)(0, 1) = (*A)(0, 1);
    (*TildeA)(1, 0) = (*A)(1, 0);
    (*TildeA)(1, 1) = (*A)(1, 1) ;

    SiconosMatrix * L = new SimpleMatrix(2, noutput);
    (*L)(0, 0) = 1.0;
    (*L)(1, 0) = 1.0;
    SiconosMatrix * G = new SimpleMatrix(noutput, 2);
    (*G)(0, 0) = 2.0;
    (*G)(0, 1) = 2.0;

    SiconosMatrix * hatA = new SimpleMatrix(2, 2);
    (*hatA) = (*A)     -   prod(*L, *G);
    cout << "A-LG" << endl;
    hatA->display();
    (*TildeA)(2, 2) = (*hatA)(0, 0);
    (*TildeA)(2, 3) = (*hatA)(0, 1);
    (*TildeA)(3, 2) = (*hatA)(1, 0);
    (*TildeA)(3, 3) = (*hatA)(1, 1);

    SiconosMatrix * LG = new SimpleMatrix(2, 2);
    (*LG) =  prod(*L, *G);
    (*TildeA)(2, 0) = (*LG)(0, 0);
    (*TildeA)(3, 0) = (*LG)(1, 0);
    (*TildeA)(2, 1) = (*LG)(0, 1);
    (*TildeA)(3, 1) = (*LG)(1, 1);
    TildeA-> display();

    SiconosVector * x0 = new SimpleVector(ndof);
    (*x0)(0) = Vinit;
    FirstOrderLinearDS * processObserver  = new FirstOrderLinearDS(0, *x0, *TildeA);
    processObserver->setComputeBFunction("SingleDSObserverLCSPlugin.so", "computeU");
    DynamicalSystemsSet allDS;
    allDS.insert(processObserver);


    // Relations
    unsigned int ninter = 2; // dimension of your Interaction = size of y and lambda vectors
    SiconosMatrix * B = new SimpleMatrix(ndof, ninter);
    (*B)(0, 0) = -1.0;
    (*B)(1, 0) = 1.0;
    (*B)(2, 1) = -1.0;
    (*B)(3, 1) = 1.0;
    B->display();
    SiconosMatrix * C = new SimpleMatrix(ninter, ndof);
    (*C)(0, 0) = -1.0;
    (*C)(0, 1) = 1.0;
    (*C)(1, 2) = -1.0;
    (*C)(1, 3) = 1.0;
    C->display();

    FirstOrderLinearR * myProcessRelation = new FirstOrderLinearR(C, B);

    myProcessRelation->setComputeEFunction("SingleDSObserverLCSPlugin.so", "computeE");

    SiconosMatrix* D = new SimpleMatrix(ninter, ninter);
    (*D)(0, 0) = 1.0;
    (*D)(1, 1) = 1.0;
    myProcessRelation->setDPtr(D);
    D->display();
    //return 0;

    // NonSmoothLaw
    unsigned int nslawSize = 2;
    NonSmoothLaw * myNslaw = new ComplementarityConditionNSL(nslawSize);

    // Choose a name and a number for your Interaction
    string nameInter = "Interaction";
    unsigned int numInter = 1;
    cout << "tutu " << endl;

    allDS.display();


    DynamicalSystemsSet dsProcessConcerned;
    dsProcessConcerned.insert(allDS.getPtr(0));

    Interaction* myProcessInteraction = new Interaction(nameInter, dsProcessConcerned, numInter, ninter, myNslaw, myProcessRelation);

    InteractionsSet allInteractions;
    allInteractions.insert(myProcessInteraction);
    cout << "tutu " << endl;

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

    int k = 0; // Current step
    int N = td->getNSteps(); // Number of time steps
    unsigned int outputSize = 7; // number of required data
    SimpleMatrix dataPlot(N, outputSize);
    SiconosVector * processLambda = myProcessInteraction->getLambdaPtr(0);


    // We get values for the initial time step:
    // time
    dataPlot(k, 0) = s->getNextTime();;
    dataPlot(k, 1) = (processObserver->getX())(0);
    dataPlot(k, 2) = (processObserver->getX())(1);
    dataPlot(k, 3) = (processObserver->getX())(2);
    dataPlot(k, 4) = (processObserver->getX())(3);
    dataPlot(k, 5) = (*processLambda)(0);
    dataPlot(k, 6) = (*processLambda)(1);


    // Simulation loop
    while (k < N - 1)
    {
      k++;
      cout << "k = " << k << endl;
      s->nextStep();
      // get current time step

      s->computeOneStep();

      dataPlot(k, 0) = s->getNextTime();;
      dataPlot(k, 1) = (processObserver->getX())(0);
      dataPlot(k, 2) = (processObserver->getX())(1);
      dataPlot(k, 3) = (processObserver->getX())(2);
      dataPlot(k, 4) = (processObserver->getX())(3);
      dataPlot(k, 5) = (*processLambda)(0);
      dataPlot(k, 6) = (*processLambda)(1);
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
