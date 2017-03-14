
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
    double T = 20.0;        // Total simulation times
    double h = 1.0e-1;      // Time step
    double Vinit = 10.0;

    double G = 10.0;
    double beta = .3;

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
    (*A)(0, 1) = 1.0;
    SP::SiconosVector x0(new SiconosVector(ndof));
    (*x0)(0) = Vinit;
    (*x0)(1) = Vinit;
    SP::FirstOrderLinearTIDS doubleIntegrator(new FirstOrderLinearTIDS(x0, A));

    // --------------------
    // --- Interactions ---
    // --------------------
    unsigned int ninter = 2; // dimension of your Interaction = size of y and lambda vectors

    // First relation, related to the doubleIntegrator
    // y = Cx + Dlambda
    // r = Blambda
    SP::SimpleMatrix B(new SimpleMatrix(ndof, ninter));
    (*B)(1, 0) = G;
    (*B)(1, 1) = G*beta;
    SP::SimpleMatrix C(new SimpleMatrix(ninter, ndof));
    C->eye();
    SP::FirstOrderLinearTIR twistingRelation(new FirstOrderLinearTIR(C, B));

    // NonSmoothLaw
    unsigned int nslawSize = 2;
    SP::SimpleMatrix H(new SimpleMatrix(4, 2));
    (*H)(0, 0) = 1.0;
    (*H)(1, 0) = -h/2.0;
    (*H)(2, 0) = -1.0;
    (*H)(3, 0) = h/2.0;
    (*H)(1, 1) = 1.0;
    (*H)(3, 1) = -1.0;

    SP::SiconosVector K(new SiconosVector(4));
    (*K)(0) = -1.0;
    (*K)(1) = -1.0;
    (*K)(2) = -1.0;
    (*K)(3) = -1.0;
    SP::NonSmoothLaw nslaw(new NormalConeNSL(nslawSize, H, K));

    SP::Interaction twistingInteraction(new Interaction(nslaw, twistingRelation));

    // -------------
    // --- Model ---
    // -------------
    SP::Model itw(new Model(t0, T));
    itw->nonSmoothDynamicalSystem()->insertDynamicalSystem(doubleIntegrator);
    itw->nonSmoothDynamicalSystem()->link(twistingInteraction,doubleIntegrator);

    // ------------------
    // --- Simulation ---
    // ------------------
    // TimeDiscretisation
    SP::TimeDiscretisation td(new TimeDiscretisation(t0, h));
    // == Creation of the Simulation ==
    SP::TimeStepping s(new TimeStepping(td));
    // -- OneStepIntegrators --
    double theta = 0.5;
    SP::EulerMoreauOSI integrator(new EulerMoreauOSI(theta));
    s->insertIntegrator(integrator);
    itw->setSimulation(s);
    // -- OneStepNsProblem --

    SP::AVI osnspb(new AVI());

    s->insertNonSmoothProblem(osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Simulation initialisation ..." << endl << endl;

    itw->initialize();


    //  (s->oneStepNSProblems)[0]->initialize();


    // --- Get the values to be plotted ---
    unsigned outputSize = 5; // number of required data
    unsigned N = ceil((T - t0) / h) + 10; // Number of time steps

    SP::SiconosMatrix dataPlot(new SimpleMatrix(N, outputSize));

    SiconosVector& xProc = *doubleIntegrator->x();
    SiconosVector& lambdaProc = *twistingInteraction->lambda(0);

    // -> saved in a matrix dataPlot
    (*dataPlot)(0, 0) = itw->t0(); // Initial time of the model
    (*dataPlot)(0, 1) = xProc(0);
    (*dataPlot)(0, 2) = xProc(1);
    (*dataPlot)(0, 3) = -1.0;


    (*dataPlot)(0, 4) = -1.0;



    // ==== Simulation loop =====
    cout << "====> Start computation ... " << endl << endl;

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    unsigned int k = 0; // Current step

    // Simulation loop
    boost::timer time;
    time.restart();
    while (s->hasNextEvent())
    {
      k++;
      //  osnspb->setNumericsVerboseMode(1);

      //  *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
      s->computeOneStep();
      (*dataPlot)(k, 0) = s->nextTime();
      (*dataPlot)(k, 1) = xProc(0);
      (*dataPlot)(k, 2) = xProc(1);
      (*dataPlot)(k, 3) = lambdaProc(0);
      (*dataPlot)(k, 4) = lambdaProc(1);
      s->nextStep();
    }
    dataPlot->resize(k, dataPlot->size(1));

    cout << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix::write("Twisting.dat", "ascii", *dataPlot, "noDim");

    // Comparison with a reference file
    SP::SiconosMatrix dataPlotRef( new SimpleMatrix(*dataPlot));
    dataPlotRef->zero();
    ioMatrix::read("Twisting.ref", "ascii", *dataPlotRef);
    // We do not compare the Lagrange multiplier that are very sensitive to numerical approximations
    Index * coord= new Index(4);
    (*coord)[0]= 0;
    (*coord)[1]= k;
    (*coord)[2]= 0;
    (*coord)[3]= 3;
    SP::SiconosMatrix  stateDataPlot(new SimpleMatrix(*dataPlot, *coord));
    SP::SiconosMatrix  stateDataPlotRef(new SimpleMatrix(*dataPlotRef, *coord));

    double error =  (*stateDataPlotRef-*stateDataPlot).normInf() ;
    std::cout << "Error =" << error <<std::endl;
    if (error > 1e-12)
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
