
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
    double T = 2;        // Total simulation time
    double h = 1.0e-4;      // Time step
    double Vinit = 10.0;


    // ================= Creation of the model =======================
    // Steps:
    // - create some Dynamical Systems
    // - create some Interactions between those Dynamical Systems
    //   Interaction = some relations (constraints) and a NonSmoothLaw
    // - create a NonSmoothDynamicalSystem with the DynamicalSystems and the Interactions
    // - add this NonSmoothDynamicalSystem into a Model
    // - add a Simulation to the model
    //  Simulation = TimeDiscretisation + OneStepIntegrator and OneStepNSProblem

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // First System:
    // dx/dt = Ax + u(t) + r
    // x(0) = x0
    // Note: r = Blambda, B defines in relation below.

    SiconosMatrix * A = new SimpleMatrix(ndof, ndof);
    (*A)(0, 0) = 1.0;
    (*A)(0, 1) = 1.0;
    (*A)(1, 0) = 3.0;
    (*A)(1, 1) = 1.0;
    SiconosVector * x0 = new SimpleVector(ndof);
    (*x0)(0) = Vinit;
    FirstOrderLinearDS * process  = new FirstOrderLinearDS(0, *x0, *A);
    process->setComputeBFunction("ObserverLCSPlugin.so", "uProcess");

    // Second System, the observer:
    // dx/dt = A hatx + u(t) + L(y-haty)
    // hatx(0) = x0
    // y = Gx
    // u(t) + Ly  is computed with uObserver function

    unsigned int noutput = 1;
    SiconosMatrix * L = new SimpleMatrix(ndof, noutput);
    (*L)(0, 0) = 1.0;
    (*L)(1, 0) = 1.0;
    SiconosMatrix * G = new SimpleMatrix(noutput, ndof);
    (*G)(0, 0) = 2.0;
    (*G)(0, 1) = 2.0;

    // hatA is initialized with A
    SimpleMatrix* hatA = new SimpleMatrix(ndof, ndof);
    (*hatA)(0, 0) = -1.0;
    (*hatA)(0, 1) = -1.0;
    (*hatA)(1, 0) = 1.0;
    (*hatA)(1, 1) = -1.0;

    SiconosVector * obsX0 = new SimpleVector(ndof);
    FirstOrderLinearDS* observer = new FirstOrderLinearDS(1, *obsX0, *hatA);
    observer->setComputeBFunction("ObserverLCSPlugin.so", "uObserver");
    //    SimpleVector * z = new SimpleVector(1);
    observer->setZPtr(process->getXPtr());
    // The set of all DynamicalSystems
    DynamicalSystemsSet allDS;
    allDS.insert(process);
    allDS.insert(observer);

    // --------------------
    // --- Interactions ---
    // --------------------
    unsigned int ninter = 1; // dimension of your Interaction = size of y and lambda vectors
    unsigned int ncontrol = 1;

    // First relation, related to the process
    // y = Cx + Dlambda
    // r = Blambda
    SiconosMatrix * B = new SimpleMatrix(ndof, ninter);
    (*B)(0, 0) = -1.0;
    (*B)(1, 0) = 1.0;
    SiconosMatrix * C = new SimpleMatrix(ninter, ndof);
    (*C)(0, 0) = -1.0;
    (*C)(0, 1) = 1.0;
    FirstOrderLinearR * myProcessRelation = new FirstOrderLinearR(C, B);
    SiconosMatrix* D = new SimpleMatrix(ninter, ninter);
    (*D)(0, 0) = 1.0;

    myProcessRelation->setDPtr(D);
    myProcessRelation->setComputeEFunction("ObserverLCSPlugin.so", "computeE");

    // Second relation, related to the observer
    // haty = C hatX + D hatLambda + E
    // hatR = B hatLambda
    FirstOrderLinearR * myObserverRelation = new FirstOrderLinearR(C, B);
    myObserverRelation->setDPtr(D);
    myObserverRelation->setComputeEFunction("ObserverLCSPlugin.so", "computeE");

    // NonSmoothLaw
    unsigned int nslawSize = 1;
    NonSmoothLaw * myNslaw = new ComplementarityConditionNSL(nslawSize);

    // The Interaction which involves the first DS (the process)
    string nameInter = "processInteraction"; // Name
    unsigned int numInter = 1; // Dim of the interaction = dim of y and lambda vectors

    Interaction* myProcessInteraction = new Interaction(nameInter, process, numInter, ninter, myNslaw, myProcessRelation);

    // The Interaction which involves the second DS (the observer)
    string nameInter2 = "observerInteraction"; // Name
    unsigned int numInter2 = 2;

    Interaction* myObserverInteraction = new Interaction(nameInter2, observer, numInter2, ninter, myNslaw, myObserverRelation);

    // The set of all Interactions
    InteractionsSet allInteractions;
    allInteractions.insert(myProcessInteraction);
    allInteractions.insert(myObserverInteraction);

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    NonSmoothDynamicalSystem* myNSDS = new NonSmoothDynamicalSystem(allDS, allInteractions);

    // -------------
    // --- Model ---
    // -------------
    Model * ObserverLCS = new Model(t0, T);
    ObserverLCS->setNonSmoothDynamicalSystemPtr(myNSDS); // set NonSmoothDynamicalSystem of this model

    // ------------------
    // --- Simulation ---
    // ------------------
    // TimeDiscretisation
    TimeDiscretisation * td = new TimeDiscretisation(h, ObserverLCS);
    // == Creation of the Simulation ==
    TimeStepping * s = new TimeStepping(td);
    // -- OneStepIntegrators --
    double theta = 0.5;
    Moreau* myIntegrator = new Moreau(allDS, theta, s);

    // -- OneStepNsProblem --
    IntParameters iparam(5);
    iparam[0] = 1000; // Max number of iteration
    DoubleParameters dparam(5);
    dparam[0] = 1e-6; // Tolerance
    string solverName = "Lemke" ;
    NonSmoothSolver * mySolver = new NonSmoothSolver(solverName, iparam, dparam);
    LCP * osnspb = new LCP(s, mySolver);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Simulation initialisation ..." << endl << endl;

    s->initialize();

    // --- Get the values to be plotted ---
    unsigned int outputSize = 10; // number of required data
    int N = td->getNSteps(); // Number of time steps
    SimpleMatrix dataPlot(N, outputSize);

    SiconosVector * xProc = process->getXPtr();
    SiconosVector * xObs = observer->getXPtr();
    SiconosVector * lambdaProc = myProcessInteraction->getLambdaPtr(0);
    SiconosVector * lambdaObs = myObserverInteraction->getLambdaPtr(0);
    SiconosVector * yProc = myProcessInteraction->getYPtr(0);
    SiconosVector * yObs = myObserverInteraction->getYPtr(0);
    SiconosVector * z = observer->getZPtr();

    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = ObserverLCS->getT0(); // Initial time of the model
    dataPlot(0, 1) = (*xProc)(0);
    dataPlot(0, 2) = (*xProc)(1);
    dataPlot(0, 3) = (*xObs)(0);
    dataPlot(0, 4) = (*xObs)(1);
    dataPlot(0, 5) = (*lambdaProc)(0);
    dataPlot(0, 6) = (*lambdaObs)(0);
    dataPlot(0, 7) = (*yProc)(0);
    dataPlot(0, 8) = (*yObs)(0);
    dataPlot(0, 9) = (*z)(0);


    // ==== Simulation loop =====
    cout << "====> Start computation ... " << endl << endl;

    // *z = *(myProcessInteraction->getYPtr(0)->getVectorPtr(0));
    int k = 0; // Current step

    // Simulation loop
    boost::timer time;
    time.restart();
    while (k < N - 1)
    {
      k++;
      //  *z = *(myProcessInteraction->getYPtr(0)->getVectorPtr(0));
      s->computeOneStep();
      dataPlot(k, 0) = s->getNextTime();
      dataPlot(k, 1) = (*xProc)(0);
      dataPlot(k, 2) = (*xProc)(1);
      dataPlot(k, 3) = (*xObs)(0);
      dataPlot(k, 4) = (*xObs)(1);
      dataPlot(k, 5) = (*lambdaProc)(0);
      dataPlot(k, 6) = (*lambdaObs)(0);
      dataPlot(k, 7) = (*yProc)(0);
      dataPlot(k, 8) = (*yObs)(0);
      dataPlot(k, 9) = (*z)(0);
      s->nextStep();
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix io("ObserverLCS.dat", "ascii");
    io.write(dataPlot, "noDim");
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
