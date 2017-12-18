
#include "SiconosKernel.hpp"

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

    SP::SiconosMatrix A(new SimpleMatrix(ndof, ndof));
    (*A)(0, 0) = 1.0;
    (*A)(0, 1) = 1.0;
    (*A)(1, 0) = 3.0;
    (*A)(1, 1) = 1.0;
    SP::SiconosVector x0(new SiconosVector(ndof));
    (*x0)(0) = Vinit;
    SP::FirstOrderLinearDS process(new FirstOrderLinearDS(x0, A));
    process->setComputebFunction("ObserverLCSPlugin", "uProcess");

    // Second System, the observer:
    // dx/dt = A hatx + u(t) + L(y-haty)
    // hatx(0) = x0
    // y = Gx
    // u(t) + Ly  is computed with uObserver function

    unsigned int noutput = 1;
    SP::SiconosMatrix L(new SimpleMatrix(ndof, noutput));
    (*L)(0, 0) = 1.0;
    (*L)(1, 0) = 1.0;
    SP::SiconosMatrix G(new SimpleMatrix(noutput, ndof));
    (*G)(0, 0) = 2.0;
    (*G)(0, 1) = 2.0;

    // hatA is initialized with A
    SP::SimpleMatrix hatA(new SimpleMatrix(ndof, ndof));
    (*hatA)(0, 0) = -1.0;
    (*hatA)(0, 1) = -1.0;
    (*hatA)(1, 0) = 1.0;
    (*hatA)(1, 1) = -1.0;

    SP::SiconosVector obsX0(new SiconosVector(ndof));
    SP::FirstOrderLinearDS observer(new FirstOrderLinearDS(obsX0, hatA));
    observer->setComputebFunction("ObserverLCSPlugin", "uObserver");
    //    SiconosVector z(new SiconosVector(1);
    observer->setzPtr(process->x());
    // The set of all DynamicalSystems
    // --------------------
    // --- Interactions ---
    // --------------------
    unsigned int ninter = 1; // dimension of your Interaction = size of y and lambda vectors

    // First relation, related to the process
    // y = Cx + Dlambda
    // r = Blambda
    SP::SimpleMatrix B(new SimpleMatrix(ndof, ninter));
    (*B)(0, 0) = -1.0;
    (*B)(1, 0) = 1.0;
    SP::SimpleMatrix C(new SimpleMatrix(ninter, ndof));
    (*C)(0, 0) = -1.0;
    (*C)(0, 1) = 1.0;
    SP::FirstOrderLinearR myProcessRelation(new FirstOrderLinearR(C, B));
    SP::SimpleMatrix D(new SimpleMatrix(ninter, ninter));
    (*D)(0, 0) = 1.0;

    myProcessRelation->setDPtr(D);
    myProcessRelation->setComputeEFunction("ObserverLCSPlugin", "computeE");

    // Second relation, related to the observer
    // haty = C hatX + D hatLambda + E
    // hatR = B hatLambda
    SP::FirstOrderLinearR myObserverRelation(new FirstOrderLinearR(C, B));
    myObserverRelation->setDPtr(D);
    myObserverRelation->setComputeEFunction("ObserverLCSPlugin", "computeE");

    // NonSmoothLaw
    unsigned int nslawSize = 1;
    SP::NonSmoothLaw myNslaw(new ComplementarityConditionNSL(nslawSize));

    // The Interaction which involves the first DS (the process)
    SP::Interaction myProcessInteraction(new Interaction(myNslaw, myProcessRelation));

    // The Interaction which involves the second DS (the observer)
    SP::Interaction myObserverInteraction(new Interaction(myNslaw, myObserverRelation));

    // -------------
    // --- Model ---
    // -------------
    SP::Model ObserverLCS(new Model(t0, T));
    ObserverLCS->nonSmoothDynamicalSystem()->insertDynamicalSystem(process);
    ObserverLCS->nonSmoothDynamicalSystem()->insertDynamicalSystem(observer);
    ObserverLCS->nonSmoothDynamicalSystem()->link(myProcessInteraction, process);
    ObserverLCS->nonSmoothDynamicalSystem()->link(myObserverInteraction, observer);

    // ------------------
    // --- Simulation ---
    // ------------------
    // TimeDiscretisation
    SP::TimeDiscretisation td(new TimeDiscretisation(t0, h));
    // == Creation of the Simulation ==
    SP::TimeStepping s(new TimeStepping(td));
    // -- OneStepIntegrators --
    double theta = 0.5;
    SP::EulerMoreauOSI myIntegrator(new EulerMoreauOSI(theta));

    s->insertIntegrator(myIntegrator);

    // -- OneStepNsProblem --
    SP::LCP osnspb(new LCP());
    s->insertNonSmoothProblem(osnspb);
    ObserverLCS->setSimulation(s);
    
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Simulation initialisation ..." << endl << endl;

    ObserverLCS->initialize();

    // --- Get the values to be plotted ---
    unsigned int outputSize = 10; // number of required data
    unsigned int N = ceil((T - t0) / h); // Number of time steps

    SimpleMatrix dataPlot(N, outputSize);

    SP::SiconosVector xProc = process->x();
    SP::SiconosVector xObs = observer->x();
    SP::SiconosVector lambdaProc = myProcessInteraction->lambda(0);
    SP::SiconosVector lambdaObs = myObserverInteraction->lambda(0);
    SP::SiconosVector yProc = myProcessInteraction->y(0);
    SP::SiconosVector yObs = myObserverInteraction->y(0);
    SP::SiconosVector z = observer->z();

    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = ObserverLCS->t0(); // Initial time of the model
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

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    unsigned int k = 0; // Current step

    // Simulation loop
    boost::timer time;
    time.restart();
    while (k < N - 1)
    {
      k++;
      //  *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
      s->computeOneStep();
      dataPlot(k, 0) = s->nextTime();
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
    ioMatrix::write("ObserverLCS.dat", "ascii", dataPlot, "noDim");
    // Comparison with a reference file
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("ObserverLCS.ref", "ascii", dataPlotRef);

    std::cout << "error =" << (dataPlot - dataPlotRef).normInf() << std::endl;
    if ((dataPlot - dataPlotRef).normInf() > 1e-10)
    {

      std::cout << "Warning. The results is rather different from the reference file." << std::endl;
      return 1;
    }



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
