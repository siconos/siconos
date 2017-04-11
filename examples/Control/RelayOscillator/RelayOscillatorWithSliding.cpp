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
    unsigned int ndof = 3;  // number of degrees of freedom of your system
    double t0 = 0.0;
    double T = 100;        // Total simulation times
    double h = 1.0e-2;      // Time step

    // Parameters of the Dynamics with default values Limit cycles with sliding
    double omega = 1.0, xsi = 1.0, lambda = 1.0;
    double k = 1.0, sigma = -1.0, rho = 1.0;

    // Parameters of the Dynamics for multi-sliding periods.

    rho = 1.0;
    xsi = 0.05;
    omega = 8.0;

    // Parameters of the Dynamics for an a periodic attractor.
    rho = 1.0;
    xsi = -0.07;
    omega = 10.0;
    lambda = 0.05;


    double xinit = 0.5;

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
    (*A)(0, 0) = -(2 * xsi * omega + lambda);
    (*A)(0, 1) = 1.0;
    (*A)(0, 2) = 0.0;
    (*A)(1, 0) = -(2 * xsi * omega * lambda + omega * omega);
    (*A)(1, 1) = 0.0;
    (*A)(1, 2) = 1.0;
    (*A)(2, 0) = -lambda * omega * omega;
    (*A)(2, 1) = 0.0;
    (*A)(2, 2) = 0.0;
    SP::SiconosVector x0(new SiconosVector(ndof));
    (*x0)(0) = 0.0;
    (*x0)(1) = xinit;
    (*x0)(2) = 0.0;

    SP::FirstOrderLinearDS process(new FirstOrderLinearDS(x0, A));
    //    process->setComputebFunction("ObserverLCSPlugin","uProcess");

    // --------------------
    // --- Interactions ---
    // --------------------
    unsigned int ninter = 1; // dimension of your Interaction = size of y and lambda vectors

    // First relation, related to the process
    // y = Cx + Dlambda
    // r = Blambda
    SP::SimpleMatrix B(new SimpleMatrix(ndof, ninter));
    (*B)(0, 0) = k;
    (*B)(1, 0) = 2 * k * sigma * rho;
    (*B)(2, 0) = k * rho * rho;

    SP::SimpleMatrix C(new SimpleMatrix(ninter, ndof));
    (*C)(0, 0) = 1.0;
    (*C)(0, 1) = 0.0;
    (*C)(0, 2) = 0.0;

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
    SP::Model relayOscillatorWithSliding(new Model(t0, T));
    relayOscillatorWithSliding->nonSmoothDynamicalSystem()->insertDynamicalSystem(process);
    relayOscillatorWithSliding->nonSmoothDynamicalSystem()->link(myProcessInteraction, process);

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
    relayOscillatorWithSliding->setSimulation(s);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Simulation initialisation ..." << endl << endl;

    relayOscillatorWithSliding->initialize();


    //  (s->oneStepNSProblems)[0]->initialize();


    // --- Get the values to be plotted ---
    unsigned int outputSize = 7; // number of required data
    unsigned int N = ceil((T - t0) / h); // Number of time steps

    SimpleMatrix dataPlot(N, outputSize);

    SP::SiconosVector xProc = process->x();
    SP::SiconosVector lambdaProc = myProcessInteraction->lambda(0);
    SP::SiconosVector yProc = myProcessInteraction->y(0);
    unsigned int step = 0; // Current step


    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = relayOscillatorWithSliding->t0(); // Initial time of the model
    dataPlot(step, 1) = (*xProc)(0);
    dataPlot(step, 2) = (*xProc)(1);
    dataPlot(step, 3) = (*xProc)(2);
    dataPlot(step, 4) = (*lambdaProc)(0);
    dataPlot(step, 5) = (*yProc)(0);


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
      dataPlot(step, 4) = (*lambdaProc)(0);
      dataPlot(step, 5) = (*yProc)(0);
      s->nextStep();
    }
    cout << endl << "End of computation - Number of iterations done: " << step - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix::write("RelayOscillatorWithSliding.dat", "ascii", dataPlot, "noDim");

    // // Comparison with a reference file
    // SimpleMatrix dataPlotRef(dataPlot);
    // dataPlotRef.zero();
    // ioMatrix::read("RelayOscillatorWithSliding.ref", "ascii", dataPlotRef);
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
