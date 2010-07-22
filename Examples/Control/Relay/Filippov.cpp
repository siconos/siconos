
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
    double T = 10;        // Total simulation times
    double h = 1.0e-2;      // Time step
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
    (*A)(0, 0) = 0.0;
    (*A)(0, 1) = 0.0;
    (*A)(1, 0) = 0.0;
    (*A)(1, 1) = 0.0;
    SP::SiconosVector x0(new SimpleVector(ndof));
    (*x0)(0) = Vinit;
    (*x0)(1) = Vinit;
    SP::FirstOrderLinearDS process(new FirstOrderLinearDS(x0, A));
    //    process->setComputebFunction("ObserverLCSPlugin.so","uProcess");

    // The set of all DynamicalSystems
    DynamicalSystemsSet allDS;
    allDS.insert(process);

    // --------------------
    // --- Interactions ---
    // --------------------
    unsigned int ninter = 2; // dimension of your Interaction = size of y and lambda vectors

    // First relation, related to the process
    // y = Cx + Dlambda
    // r = Blambda
    SP::SimpleMatrix B(new SimpleMatrix(ndof, ninter));
    (*B)(0, 0) = 1.0;
    (*B)(1, 0) = 2.0;
    (*B)(0, 1) = -2.0;
    (*B)(1, 1) = 1.0;
    *B = 2.0 * (*B);
    SP::SimpleMatrix C(new SimpleMatrix(ninter, ndof));
    (*C)(0, 0) = 1.0;
    (*C)(1, 0) = 0.0;
    (*C)(0, 1) = 0.0;
    (*C)(1, 1) = 1.0;
    SP::FirstOrderLinearR myProcessRelation(new FirstOrderLinearR(C, B));
    SP::SimpleMatrix D(new SimpleMatrix(ninter, ninter));
    (*D)(0, 0) = 0.0;
    (*D)(0, 1) = 0.0;
    (*D)(1, 0) = 0.0;
    (*D)(1, 1) = 0.0;

    myProcessRelation->setDPtr(D);
    //myProcessRelation->setComputeEFunction("ObserverLCSPlugin.so","computeE");

    // Second relation, related to the observer
    // haty = C hatX + D hatLambda + E
    // hatR = B hatLambda

    // NonSmoothLaw
    unsigned int nslawSize = 2;
    SP::NonSmoothLaw myNslaw(new RelayNSL(nslawSize));

    myNslaw->display();

    // The Interaction which involves the first DS (the process)
    string nameInter = "processInteraction"; // Name
    unsigned int numInter = 2; // Dim of the interaction = dim of y and lambda vectors

    SP::Interaction myProcessInteraction(new Interaction(nameInter, process, numInter, ninter, myNslaw, myProcessRelation));

    // The set of all Interactions
    InteractionsSet allInteractions;
    allInteractions.insert(myProcessInteraction);

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    SP::NonSmoothDynamicalSystem myNSDS(new NonSmoothDynamicalSystem(allDS, allInteractions));

    // -------------
    // --- Model ---
    // -------------
    SP::Model filippov(new Model(t0, T));
    filippov->setNonSmoothDynamicalSystemPtr(myNSDS); // set NonSmoothDynamicalSystem of this model

    // ------------------
    // --- Simulation ---
    // ------------------
    // TimeDiscretisation
    SP::TimeDiscretisation td(new TimeDiscretisation(t0, h));
    // == Creation of the Simulation ==
    SP::TimeStepping s(new TimeStepping(td));
    // -- OneStepIntegrators --
    double theta = 0.5;
    SP::Moreau myIntegrator(new Moreau(allDS, theta));
    s->insertIntegrator(myIntegrator);

    // -- OneStepNsProblem --

    SP::Relay osnspb(new Relay());


    osnspb->numericsSolverOptions()->solverId = SICONOS_RELAY_LEMKE;
    osnspb->numericsSolverOptions()->dparam[0] = 1e-08;


    s->insertNonSmoothProblem(osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Simulation initialisation ..." << endl << endl;

    filippov->initialize(s);


    //  (s->oneStepNSProblems)[0]->initialize(s);


    // --- Get the values to be plotted ---
    unsigned int outputSize = 7; // number of required data
    int N = (int)((T - t0) / h); // Number of time steps

    SimpleMatrix dataPlot(N, outputSize);

    SP::SiconosVector xProc = process->x();
    SP::SiconosVector lambdaProc = myProcessInteraction->lambda(0);
    SP::SiconosVector yProc = myProcessInteraction->y(0);

    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = filippov->t0(); // Initial time of the model
    dataPlot(0, 1) = (*xProc)(0);
    dataPlot(0, 2) = (*xProc)(1);
    dataPlot(0, 3) = -1.0;
    dataPlot(0, 4) = -1.0;
    dataPlot(0, 5) = (*xProc)(0);
    dataPlot(0, 6) = (*xProc)(1);



    // ==== Simulation loop =====
    cout << "====> Start computation ... " << endl << endl;

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    int k = 0; // Current step

    // Simulation loop
    boost::timer time;
    time.restart();
    while (k < N - 1)
    {
      k++;
      cout << "step --> " << k << endl;

      //  osnspb->setNumericsVerboseMode(1);

      //  *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
      s->computeOneStep();
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*xProc)(0);
      dataPlot(k, 2) = (*xProc)(1);
      dataPlot(k, 3) = (*lambdaProc)(0);
      dataPlot(k, 4) = (*lambdaProc)(1);
      dataPlot(k, 5) = (*yProc)(0);
      dataPlot(k, 6) = (*yProc)(1);
      s->nextStep();
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix io("Filippov.dat", "ascii");
    io.write(dataPlot, "noDim");
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in Fillipov.cpp" << endl;
  }
}
