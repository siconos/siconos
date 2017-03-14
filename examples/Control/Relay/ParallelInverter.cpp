#include "SiconosKernel.hpp"
#include <math.h>
using namespace std;
static unsigned int ndof = 4;
/* 02/2010 --> 08/2010*/
/*Author: Yamen ABDENNADHER */
/*Exemple from : Rafael Ramos, Dominigo Biel, Enric Fossas and Franisco Guinjoan. Interleaving Quasi-Sliding-Mode Control of Parallel-Connected Buck-based Inverters. IEEE vol 55 N11 Novembre 2008. For more detail, please see Yamen's report.*/
// main program
int main(int argc, char* argv[])
{
  // Exception handling
  try
  {
    // == User-defined parameters ==
    //unsigned int ndof = 4;  // number of degrees of freedom of your system
    double t0 = 0.0;
    double T = 0.02;       // Total simulation time
    double h = 1.0e-6;      // Time step
    //double Vinit= 1.0;
    double Cp = 60e-6;
    double leq = 0.385e-3;
    double Rl = 5;
    double m = 3;

    double L1 = 1.5e-3;
    double L2 = 1.22e-3;
    double L3 = 0.9e-3;
    //leq=1/((1/L1)+(1/L2)+(1/L3));

    double rl1 = 94e-3;
    double rl2 = 116e-3;
    double rl3 = 100e-3;

    double E1 = 70;
    double E2 = 70;
    double E3 = 70;

    double k1 = 1;
    double k2 = 6e-5;


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
    (*A)(0, 0) = 0;
    (*A)(0, 1) = 1;
    (*A)(0, 2) = 0;
    (*A)(0, 3) = 0;
    (*A)(1, 0) = (-1 / (Cp * leq)) - (rl1 / (Cp * L1 * Rl));
    (*A)(1, 1) = (-rl1 / L1) - (1 / (Cp * Rl));
    (*A)(1, 2) = (rl1 / (Cp * L1)) - (rl2 / (Cp * L2));
    (*A)(1, 3) = (rl1 / (Cp * L1)) - (rl3 / (Cp * L3));
    (*A)(2, 0) = -1 / L2;
    (*A)(2, 1) = 0;
    (*A)(2, 2) = -rl2 / L2;
    (*A)(2, 3) = 0;
    (*A)(3, 0) = -1 / L3;
    (*A)(3, 1) = 0;
    (*A)(3, 2) = 0;
    (*A)(3, 3) = -rl3 / L3;

    A->display();


    SP::SiconosVector x0(new SiconosVector(ndof));
    (*x0)(0) = 50;
    (*x0)(1) = 7;
    (*x0)(2) = 4;
    (*x0)(3) = 4;

    SP::FirstOrderLinearDS process(new FirstOrderLinearDS(x0, A));
    SP::SiconosVector zProc(new SiconosVector(1, 0));
    process->setzPtr(zProc);

    // --------------------
    // --- Interactions ---
    // --------------------
    unsigned int ninter = 3; // dimension of your Interaction = size of y and lambda vectors

    // First relation, related to the process
    // y = Cx + Dlambda +eDLS
    // r = Blambda
    SP::SimpleMatrix B(new SimpleMatrix(ndof, ninter));
    (*B)(0, 0) = 0;
    (*B)(0, 1) = 0;
    (*B)(0, 2) = 0;
    (*B)(1, 0) = E1 / (Cp * L1);
    (*B)(1, 1) = E2 / (Cp * L2);
    (*B)(1, 2) = E3 / (Cp * L3);
    (*B)(2, 0) = 0;
    (*B)(2, 1) = E2 / L2;
    (*B)(2, 2) = 0;
    (*B)(3, 0) = 0;
    (*B)(3, 1) = 0;
    (*B)(3, 2) = E3 / L3;

    B->display();

    *B = 1 * (*B);
    SP::SimpleMatrix C(new SimpleMatrix(ninter, ndof));
    (*C)(0, 0) = (Cp * k1 / k2) + ((m - 1) / Rl);
    (*C)(0, 1) = m * Cp;
    (*C)(0, 2) = -m;
    (*C)(0, 3) = -m;
    (*C)(1, 0) = (Cp * k1 / k2) - (1 / Rl);
    (*C)(1, 1) = 0;
    (*C)(1, 2) = m;
    (*C)(1, 3) = 0;
    (*C)(2, 0) = (Cp * k1 / k2) - (1 / Rl);
    (*C)(2, 1) = 0;
    (*C)(2, 2) = 0;
    (*C)(2, 3) = m;
    C->display();
    //((*C)*(*B))->display();
    SP::FirstOrderLinearR myProcessRelation(new FirstOrderLinearR(C, B));
    SP::SimpleMatrix D(new SimpleMatrix(ninter, ninter));
    (*D)(0, 0) = 0.0;
    (*D)(0, 1) = 0.0;
    (*D)(0, 2) = 0.0;
    (*D)(1, 0) = 0.0;
    (*D)(1, 1) = 0.0;
    (*D)(1, 2) = 0.0;
    (*D)(2, 0) = 0.0;
    (*D)(2, 1) = 0.0;
    (*D)(2, 2) = 0.0;
    myProcessRelation->setComputeEFunction("plugins", "eLDS");

    myProcessRelation->setDPtr(D);
    //myProcessRelation->setComputeEFunction("ObserverLCSPlugin","computeE");

    // Second relation, related to the observer
    // haty = C hatX + D hatLambda + E
    // hatR = B hatLambda

    // NonSmoothLaw
    unsigned int nslawSize = 3;
    SP::NonSmoothLaw myNslaw(new RelayNSL(nslawSize));

    myNslaw->display();

    // The Interaction which involves the first DS (the process)
    SP::Interaction myProcessInteraction(new Interaction(myNslaw, myProcessRelation));

    // -------------
    // --- Model ---
    // -------------
    SP::Model simpleExampleRelay(new Model(t0, T));
    simpleExampleRelay->nonSmoothDynamicalSystem()->insertDynamicalSystem(process);
    simpleExampleRelay->nonSmoothDynamicalSystem()->link(myProcessInteraction, process);
    
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
    // -- OneStepNsProblem --

    SP::Relay osnspb(new Relay(SICONOS_RELAY_ENUM));


    //osnspb->setNumericsSolverName("Lemke");
    //osnspb->numericsSolverOptions()->dparam[0]=1e-08;


    s->insertNonSmoothProblem(osnspb);


    //SP::LCP osnspb1(new LCP());

    //SP::Relay osnspb(new Relay("Lemke"));
    //s->insertNonSmoothProblem(osnspb);
    simpleExampleRelay->setSimulation(s);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Simulation initialisation ..." << endl << endl;

    simpleExampleRelay->initialize();


    //  (s->oneStepNSProblems)[0]->initialize();


    // --- Get the values to be plotted ---
    unsigned int outputSize = 17; // number of required data
    unsigned int N = ceil((T - t0) / h); // Number of time steps

    SimpleMatrix dataPlot(N, outputSize);

    SP::SiconosVector xProc = process->x();
    SP::SiconosVector lambdaProc = myProcessInteraction->lambda(0);
    SP::SiconosVector yProc = myProcessInteraction->y(0);

    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = simpleExampleRelay->t0(); // Initial time of the model
    dataPlot(0, 1) = (*xProc)(0);
    dataPlot(0, 2) = (*xProc)(1);
    dataPlot(0, 3) = (*xProc)(2);
    dataPlot(0, 4) = (*xProc)(3);
    dataPlot(0, 5) = (*lambdaProc)(0);
    dataPlot(0, 6) = (*lambdaProc)(1);
    dataPlot(0, 7) = (*lambdaProc)(2);
    dataPlot(0, 8) = (*yProc)(0);
    dataPlot(0, 9) = (*yProc)(1);
    dataPlot(0, 10) = (*yProc)(2);
    dataPlot(0, 11) = (*zProc)(0);            //v_ref
    dataPlot(0, 12) = (*zProc)(0) / Rl;       //i_ref
    dataPlot(0, 13) = 0;           //voltage_error
    dataPlot(0, 14) = 0;           //current_error
    dataPlot(0, 15) = 0;           //il1
    dataPlot(0, 16) = (*xProc)(1);



    // ==== Simulation loop =====
    cout << "====> Start computation ... " << endl << endl;

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    unsigned int k = 0; // Current step

    // Simulation loop
    boost::timer time;
    time.restart();

    unsigned int i = 0;
    int j = 0;
    SP::SiconosVector err(new SiconosVector(2));
    SP::SiconosVector il1(new SiconosVector(1));
    SP::SiconosVector temps(new SiconosVector(N + 1));


    while (k < N - 1)
    {
      k++;
      //    osnspb->setNumericsVerboseMode(1);
      //  *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
      s->computeOneStep();
      (*err)(0) = abs((*zProc)(0) - (*xProc)(0)); //voltage error
      (*err)(1) = abs(((*zProc)(0) / Rl) - ((*xProc)(0) / Rl)); //current error
      (*il1)(0) = ((*xProc)(0) / Rl) - ((*xProc)(2) + (*xProc)(3)) + Cp * (*xProc)(1); //current through L_1


      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*xProc)(0);          //output voltage v0
      dataPlot(k, 2) = (*xProc)(0) / Rl;     //output current i0
      dataPlot(k, 3) = (*xProc)(2);          //il2
      dataPlot(k, 4) = (*xProc)(3);          //il3
      dataPlot(k, 5) = (*lambdaProc)(0);     //u1
      dataPlot(k, 6) = (*lambdaProc)(1);     //u2
      dataPlot(k, 7) = (*lambdaProc)(2);     //u3
      dataPlot(k, 8) = (*yProc)(0);          //y1
      dataPlot(k, 9) = (*yProc)(1);          //y2
      dataPlot(k, 10) = (*yProc)(2);         //y3
      dataPlot(k, 11) = (*zProc)(0);            //v_ref
      dataPlot(k, 12) = (*zProc)(0) / Rl;       //i_ref
      dataPlot(k, 13) = (*err)(0);           //voltage_error
      dataPlot(k, 14) = (*err)(1);           //current_error
      dataPlot(k, 15) = (*il1)(0);           //il1
      dataPlot(k, 16) = (*xProc)(1);           //v'0
      s->nextStep();

      //////////////////////////////////////////////////////////////////////////////////
      if (((*yProc)(0) > 1e-8) || ((*yProc)(0) < -1e-8))
      {
        (*temps)(k) = (*yProc)(0);
      }
      ///////////////////////////////////////////////////////////////////////////////////
    }

    while (i < N - 1)
    {
      if ((*temps)(i) == 0)
        j = j + 1;

      i++;
    }
    cout << "The sliding mode appears at the step number : \n" << N - 1 - j << endl;

    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix::write("ParallelInverter.dat", "ascii", dataPlot, "noDim");
    std::cout << "Comparison with a reference file" << std::endl;
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("ParallelInverter.ref", "ascii", dataPlotRef);
    SP::SiconosVector errSim = compareMatrices(dataPlot, dataPlotRef);
    std::cout << "errSim display :" << std::endl;
    errSim->display();
    double error = ((dataPlot - dataPlotRef).normInf())/(dataPlotRef.normInf());
    std::cout << "error =" << error << std::endl;
    if (error > 1e-12) // some data are > 1e4
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
    cout << "Exception caught in SimpleExampleRelay.cpp" << endl;
  }
}
