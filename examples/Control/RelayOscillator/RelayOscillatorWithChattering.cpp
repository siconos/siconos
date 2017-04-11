#include "SiconosKernel.hpp"
#include <math.h>
using namespace std;

// main program
/* Example of a limit cycle in Relay system with sliding motion and multisliding motion
 * see
 * "SELF-OSCILLATIONS AND SLIDING IN RELAY FEEDBACK SYSTEMS: SYMMETRY AND BIFURCATIONS"
 * Di Bernardo, M. , Johansson, K.H. and Vasca, F.
 * International Journal of Bifurcation and Chaos, Vol. 11, No. 4 (2001) 1121–1140
 *
 */


int main(int argc, char* argv[])
{
  // Exception handling
  try
  {
    // == User-defined parameters ==
    unsigned int ndof = 4;  // number of degrees of freedom of your system
    double t0 = 0.0;
    double T = 5;        // Total simulation times
    double h = 1.0e-4;      // Time step

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

    (*A)(0, 0) = -4.0;
    (*A)(0, 1) = 1.0;
    (*A)(0, 2) = 0.0;
    (*A)(0, 3) = 0.0;

    (*A)(1, 0) = -6.0;
    (*A)(1, 1) = 0.0;
    (*A)(1, 2) = 1.0;
    (*A)(1, 3) = 0.0;

    (*A)(2, 0) = -4.0;
    (*A)(2, 1) = 0.0;
    (*A)(2, 2) = 0.0;
    (*A)(2, 3) = 1.0;

    (*A)(3, 0) = -1.0;
    (*A)(3, 1) = 0.0;
    (*A)(3, 2) = 0.0;
    (*A)(3, 3) = 0.0;

    A->display();
    SP::SiconosVector x0(new SiconosVector(ndof));
    (*x0)(0) = -0.01;
    (*x0)(1) = -0.02;
    (*x0)(2) = -1.0;
    (*x0)(3) = 0.5;

    SP::FirstOrderLinearDS process(new FirstOrderLinearDS(x0, A));
    //    process->setComputebFunction("ObserverLCSPlugin","uProcess");

    // --------------------
    // --- Interactions ---
    // --------------------
    unsigned int ninter = 1; // dimension of your Interaction = size of y and lambda vectors

    // First relation, related to the process
    // y = Cx + Dlambda
    // r = Blambda
    double xsi = 0.2;

    SP::SimpleMatrix B(new SimpleMatrix(ndof, ninter));
    (*B)(0, 0) = 0;
    (*B)(1, 0) = 1;
    (*B)(2, 0) = -2 * xsi;
    (*B)(3, 0) = xsi * xsi;

    SP::SimpleMatrix C(new SimpleMatrix(ninter, ndof));
    (*C)(0, 0) = 1.0;
    (*C)(0, 1) = 0.0;
    (*C)(0, 2) = 0.0;
    (*C)(0, 3) = 0.0;

    SP::FirstOrderLinearR myProcessRelation(new FirstOrderLinearR(C, B));
    SP::SimpleMatrix D(new SimpleMatrix(ninter, ninter));
    (*D)(0, 0) = 0.0;

    myProcessRelation->setDPtr(D);
    //myProcessRelation->setComputeEFunction("ObserverLCSPlugin","computeE");

    // Second relation, related to the observer
    // haty = C hatX + D hatLambda + E
    // hatR = B hatLambda

    // NonSmoothLaw
    unsigned int nslawSize = 1;
    SP::NonSmoothLaw myNslaw(new RelayNSL(nslawSize));

    myNslaw->display();

    // The Interaction which involves the first DS (the process)
    SP::Interaction myProcessInteraction(new Interaction(myNslaw, myProcessRelation));
    // -------------
    // --- Model ---
    // -------------
    SP::Model relayOscillatorWithChattering(new Model(t0, T));
    relayOscillatorWithChattering->nonSmoothDynamicalSystem()->insertDynamicalSystem(process);
    relayOscillatorWithChattering->nonSmoothDynamicalSystem()->link(myProcessInteraction, process);
 
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

    SP::Relay osnspb(new Relay());


    osnspb->setSolverId(SICONOS_RELAY_LEMKE);
    osnspb->numericsSolverOptions()->dparam[0] = 1e-08;


    s->insertNonSmoothProblem(osnspb);
    relayOscillatorWithChattering->setSimulation(s);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Simulation initialisation ..." << endl << endl;

    relayOscillatorWithChattering->initialize();


    //  (s->oneStepNSProblems)[0]->initialize();


    // --- Get the values to be plotted ---
    unsigned int outputSize = 8; // number of required data
    unsigned int N = ceil((T - t0) / h); // Number of time steps

    SimpleMatrix dataPlot(N, outputSize);

    SP::SiconosVector xProc = process->x();
    SP::SiconosVector lambdaProc = myProcessInteraction->lambda(0);
    SP::SiconosVector yProc = myProcessInteraction->y(0);
    unsigned int step = 0; // Current step


    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = relayOscillatorWithChattering->t0(); // Initial time of the model
    dataPlot(step, 1) = (*xProc)(0);
    dataPlot(step, 2) = (*xProc)(1);
    dataPlot(step, 3) = (*xProc)(2);
    dataPlot(step, 4) = (*xProc)(3);
    dataPlot(step, 5) = (*lambdaProc)(0);
    dataPlot(step, 6) = (*yProc)(0);


    // ==== Simulation loop =====
    cout << "====> Start computation ... " << endl << endl;

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    // Simulation loop
    boost::timer time;
    time.restart();
    while (step < N - 1)
    {
      step++;

      //  osnspb->setNumericsVerboseMode(1);

      //  *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
      s->computeOneStep();
      dataPlot(step, 0) = s->nextTime();
      dataPlot(step, 1) = (*xProc)(0);
      dataPlot(step, 2) = (*xProc)(1);
      dataPlot(step, 3) = (*xProc)(2);
      dataPlot(step, 4) = (*xProc)(3);
      dataPlot(step, 5) = (*lambdaProc)(0);
      dataPlot(step, 6) = (*yProc)(0);

      s->nextStep();
    }
    cout << endl << "End of computation - Number of iterations done: " << step - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix::write("RelayOscillatorWithChattering.dat", "ascii", dataPlot, "noDim");

    // // Comparison with a reference file
    // SimpleMatrix dataPlotRef(dataPlot);
    // dataPlotRef.zero();
    // ioMatrix::read("RelayOscillatorWithChattering.ref", "ascii", dataPlotRef);
    // //std::cout << (dataPlot-dataPlotRef).normInf() <<std::endl;
    // if ((dataPlot-dataPlotRef).normInf() > 1e-12)
    // {
    //   std::cout << "Warning. The results is rather different from the reference file."<< std::endl;
    //   return 1;
    // }


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
