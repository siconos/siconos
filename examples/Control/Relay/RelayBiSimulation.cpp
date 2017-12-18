
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
    double T = 1;        // Total simulation time
    double h = 1.0e-4;      // Time step
    double hcontroller = 1.0e-2;      // Time step
    double Vinit = 1.0;

    if (h > hcontroller)
    {
      RuntimeException::selfThrow("hcontroller must be larger than h");
    }


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

    SP::SiconosVector x0(new SiconosVector(ndof));
    (*x0)(0) = Vinit;
    (*x0)(1) = -Vinit;

    SP::FirstOrderLinearDS processDS(new FirstOrderLinearDS(x0, A));
    processDS->setComputebFunction("plugins", "computeB");

    SP::FirstOrderLinearDS controllerDS(new FirstOrderLinearDS(x0, A));

    // --------------------
    // --- Interactions ---
    // --------------------
    unsigned int ninter = 2; // dimension of your Interaction = size of y and lambda vectors

    // First relation, related to the process
    // y = Cx + Dlambda
    // r = Blambda
    SP::SimpleMatrix B(new SimpleMatrix(ndof, ninter));
    (*B)(0, 0) = 1.0;
    (*B)(1, 0) = 0.0;
    (*B)(0, 1) = 0.0;
    (*B)(1, 1) = 1.0;
    *B = 2.0 * (*B);
    SP::SimpleMatrix C(new SimpleMatrix(ninter, ndof));
    (*C)(0, 0) = 1.0;
    (*C)(1, 0) = 0.0;
    (*C)(0, 1) = 0.0;
    (*C)(1, 1) = 1.0;

    SP::SimpleMatrix D(new SimpleMatrix(ninter, ninter));
    (*D)(0, 0) = 0.0;
    (*D)(0, 1) = 0.0;
    (*D)(1, 0) = 0.0;
    (*D)(1, 1) = 0.0;


    //     SP::FirstOrderLinearR myProcessRelation(new FirstOrderLinearR(C,B));
    //     myProcessRelation->setDPtr(D);
    //myProcessRelation->setComputeEFunction("ObserverLCSPlugin","computeE");

    SP::FirstOrderLinearR myControllerRelation(new FirstOrderLinearR(C, B));
    myControllerRelation->setDPtr(D);

    // NonSmoothLaw
    unsigned int nslawSize = 2;
    SP::NonSmoothLaw myNslaw(new RelayNSL(nslawSize));
    myNslaw->display();

    // The Interaction which involves the first DS (the process)
    string nameInter = "processInteraction"; // Name

    //    SP::Interaction myProcessInteraction(new Interaction(myNslaw, myProcessRelation));

    SP::Interaction myControllerInteraction(new Interaction(myNslaw, myControllerRelation));



    // -------------
    // --- Model process ---
    // -------------
    SP::Model process(new Model(t0, T));
    process->nonSmoothDynamicalSystem()->insertDynamicalSystem(processDS);
    //     process->nonSmoothDynamicalSystem()->link(myProcessInteraction,processDS);

    // ------------------
    // --- Simulation ---
    // ------------------
    // TimeDiscretisation
    SP::TimeDiscretisation processTD(new TimeDiscretisation(t0, h));
    // == Creation of the Simulation ==
    SP::TimeStepping processSimulation(new TimeStepping(processTD, 0));
    // -- OneStepIntegrators --
    double theta = 0.5;
    SP::EulerMoreauOSI processIntegrator(new EulerMoreauOSI(theta));


    // -------------
    // --- Model controller ---
    // -------------
    SP::Model controller(new Model(t0, T));
    controller->nonSmoothDynamicalSystem()->insertDynamicalSystem(controllerDS);
    controller->nonSmoothDynamicalSystem()->link(myControllerInteraction, controllerDS);

    // ------------------
    // --- Simulation ---
    // ------------------
    // TimeDiscretisation
    SP::TimeDiscretisation controllerTD(new TimeDiscretisation(t0, hcontroller));
    // == Creation of the Simulation ==
    SP::TimeStepping controllerSimulation(new TimeStepping(controllerTD));
    // -- OneStepIntegrators --
    double controllertheta = 0.5;
    SP::EulerMoreauOSI controllerIntegrator(new EulerMoreauOSI(controllertheta));


    // -- OneStepNsProblem --
    SP::LCP controllerLCP(new LCP());

    SP::Relay controllerOSNSPB(new Relay(SICONOS_RELAY_PGS));
    controllerSimulation->insertNonSmoothProblem(controllerOSNSPB);


    // coupling the simulation

    SP::SiconosVector sampledControl(new SiconosVector(2));
    processDS->setzPtr(sampledControl);




    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Simulation initialisation ..." << endl << endl;
    process->setSimulation(processSimulation);
    process->initialize();
    controller->setSimulation(controllerSimulation);
    controller->initialize();

    // --- Get the values to be plotted ---
    unsigned int outputSize = 10; // number of required data
    unsigned int N = ceil((T - t0) / h) + 10; // Number of time steps
    SimpleMatrix dataPlot(N, outputSize);

    SP::SiconosVector xProc = processDS->x();

    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = process->t0(); // Initial time of the model
    dataPlot(0, 1) = (*xProc)(0);

    unsigned int Ncontroller = ceil((T - t0) / hcontroller) + 10; // Number of time steps
    SimpleMatrix dataPlotController(Ncontroller, outputSize);

    SP::SiconosVector xController = controllerDS->x();
    SP::SiconosVector lambda = myControllerInteraction->lambda(0);
    SP::SiconosVector y = myControllerInteraction->y(0);

    // -> saved in a matrix dataPlot
    dataPlotController(0, 0) = controller->t0(); // Initial time of the model
    dataPlotController(0, 1) = (*xController)(0);
    dataPlotController(0, 2) = (*xController)(1);
    dataPlotController(0, 5) = (*lambda)(0);
    dataPlotController(0, 6) = (*lambda)(1);
    dataPlotController(0, 7) = (*y)(0);
    dataPlotController(0, 8) = (*y)(1);



    // ==== Simulation loop =====
    cout << "====> Start computation ... " << endl << endl;

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    int k = 0; // Current step
    int kcontroller = 0;
    // Simulation loop
    boost::timer time;
    time.restart();

    while (controllerSimulation->hasNextEvent())
    {
      kcontroller ++ ;
//      cout << "step controller--> " << kcontroller << " at time t =" << controllerSimulation->nextTime() << endl;

      // Computation of the controller over the sampling time
      controllerSimulation->computeOneStep();

      //  input of the controller in the process thanks to z and sampledControl
      prod(1.0, *B, *lambda, *sampledControl, true);

      while (processSimulation->hasNextEvent() && processSimulation->nextTime() < controllerSimulation->nextTime())
      {
        k++;
//        cout << "         step --> " << k  << " at time t =" << processSimulation->nextTime() << endl;

        processSimulation->computeOneStep();
        dataPlot(k, 0) = processSimulation->nextTime();
        dataPlot(k, 1) = (*xProc)(0);
        dataPlot(k, 2) = (*xProc)(1);
        processSimulation->nextStep();
      }
      dataPlotController(kcontroller, 0) = controllerSimulation->nextTime() ; // Initial time of the model
      dataPlotController(kcontroller, 1) = (*xController)(0);
      dataPlotController(kcontroller, 2) = (*xController)(1);
      dataPlotController(kcontroller, 5) = (*lambda)(0);
      dataPlotController(kcontroller, 6) = (*lambda)(1);
      dataPlotController(kcontroller, 7) = (*y)(0);
      dataPlotController(kcontroller, 8) = (*y)(1);

      // feedback output of the measures for the controller
      *(xController) = *(xProc);

      controllerSimulation->nextStep();
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    dataPlotController.resize(kcontroller, outputSize);
    ioMatrix::write("RelayBiSimulation-Controller.dat", "ascii", dataPlotController, "noDim");


    // Comparison with a reference file
    //    SimpleMatrix dataPlotRef(dataPlot);
    //     dataPlotRef.zero();
    //     ioMatrix::read("RelayBiSimulation.ref", "ascii", dataPlotRef);
    //     //std::cout << (dataPlot-dataPlotRef).normInf() <<std::endl;

    //     if ((dataPlot-dataPlotRef).normInf() > 1e-12)
    //     {
    //       std::cout << "Warning. The results is rather different from the reference file."<< std::endl;
    //       return 1;
    //     }

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in RelayBiSimulation.cpp" << endl;
  }
}
