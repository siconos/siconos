#include "SiconosKernel.hpp"
#include <math.h>
using namespace std;

// main program
/* Example of a limit cycle in Relay system predicted by the describing function approach
 * see Example 7.9 in "Nonlinear Systems" H. Khalil. Third Edition. Prentice Hall, 2002 ISBN 0-13-067389-7
 *
 * The transfer function of the linear part is $G(s) 1/(s(s+1)(s+2))$ and the describing function is $\Psi(a)=4/(\pi a)$
 * The predicted limit cycle is for the output $y(t) = a \sin (\omega t)$ with $a = 2/(3\pi)$ and $\omega = \sqrt(2)$.
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
    double xinit = 3.0 * sqrt(2.0) / (2.0 * M_PI);


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
    (*A)(0, 1) = 1.0;
    (*A)(0, 2) = 0.0;
    (*A)(1, 0) = 0.0;
    (*A)(1, 1) = 0.0;
    (*A)(1, 2) = 1.0;
    (*A)(2, 0) = 0.0;
    (*A)(2, 1) = -3.0;
    (*A)(2, 2) = -2.0;
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
    (*B)(0, 0) = 0.0;
    (*B)(1, 0) = 0.0;
    (*B)(2, 0) = 1.0;

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
    SP::Model relayOscillator(new Model(t0, T));
    relayOscillator->nonSmoothDynamicalSystem()->insertDynamicalSystem(process);
    relayOscillator->nonSmoothDynamicalSystem()->link(myProcessInteraction, process);

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

    relayOscillator->setSimulation(s);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Simulation initialisation ..." << endl << endl;

    relayOscillator->initialize();


    //  (s->oneStepNSProblems)[0]->initialize();


    // --- Get the values to be plotted ---
    unsigned int outputSize = 7; // number of required data
    unsigned int N = ceil((T - t0) / h); // Number of time steps

    SimpleMatrix dataPlot(N, outputSize);

    SP::SiconosVector xProc = process->x();
    SP::SiconosVector lambdaProc = myProcessInteraction->lambda(0);
    SP::SiconosVector yProc = myProcessInteraction->y(0);
    unsigned int k = 0; // Current step


    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = relayOscillator->t0(); // Initial time of the model
    dataPlot(k, 1) = (*xProc)(0);
    dataPlot(k, 2) = (*xProc)(1);
    dataPlot(k, 3) = (*xProc)(2);
    dataPlot(k, 4) = (*lambdaProc)(0);
    dataPlot(k, 5) = (*yProc)(0);


    // ==== Simulation loop =====
    cout << "====> Start computation ... " << endl << endl;

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    // Simulation loop
    boost::timer time;
    time.restart();
    while (k < N - 1)
    {
      k++;

      //  osnspb->setNumericsVerboseMode(1);

      //  *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
      s->computeOneStep();
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*xProc)(0);
      dataPlot(k, 2) = (*xProc)(1);
      dataPlot(k, 3) = (*xProc)(2);
      dataPlot(k, 4) = (*lambdaProc)(0);
      dataPlot(k, 5) = (*yProc)(0);
      s->nextStep();
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix::write("RelayOscillator.dat", "ascii", dataPlot, "noDim");

    // Comparison with a reference file
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("RelayOscillator.ref", "ascii", dataPlotRef);
    std::cout << (dataPlot-dataPlotRef).normInf() <<std::endl;
    if ((dataPlot - dataPlotRef).normInf() > 1e-12)
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
    cout << "Exception caught in Fillipov.cpp" << endl;
  }
}
