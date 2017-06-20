
#include "SiconosKernel.hpp"

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
    double T = 25;        // Total simulation time
    double h = 1.0e-3;      // Time step
    double Vinit = 10.0;
    unsigned int noutput = 1;

    // ================= Creation of the model =======================

    // == Creation of the NonSmoothDynamicalSystem ==
    // DynamicalSystem(s)
    SimpleMatrix A(2, 2); // All components of A are automatically set to 0.
    A(0, 0) = 1.0;
    A(0, 1) = 1.0;
    A(1, 0) = 3.0;
    A(1, 1) = 1.0;
    A = 0.1 * A;
    SimpleMatrix TildeA(ndof, ndof); // All components of A are automatically set to 0.
    TildeA(0, 0) =  A(0, 0);
    TildeA(0, 1) =  A(0, 1);
    TildeA(1, 0) =  A(1, 0);
    TildeA(1, 1) =  A(1, 1) ;

    SimpleMatrix L(2, noutput);
    L(0, 0) = 1.0;
    L(1, 0) = 1.0;
    L = 0.1 * L;
    SimpleMatrix G(noutput, 2);
    G(0, 0) = 2.0;
    G(0, 1) = 2.0;

    SimpleMatrix hatA(2, 2);
    hatA = A     -   prod(L, G);
    TildeA(2, 2) = hatA(0, 0);
    TildeA(2, 3) = hatA(0, 1);
    TildeA(3, 2) = hatA(1, 0);
    TildeA(3, 3) = hatA(1, 1);

    SimpleMatrix LG(2, 2);
    LG =  prod(L, G);
    TildeA(2, 0) = LG(0, 0);
    TildeA(3, 0) = LG(1, 0);
    TildeA(2, 1) = LG(0, 1);
    TildeA(3, 1) = LG(1, 1);

    SP::SiconosVector x0(new SiconosVector(ndof));
    (*x0)(0) = Vinit;
    SP::FirstOrderLinearDS processObserver(new FirstOrderLinearDS(x0, createSPtrSimpleMatrix(TildeA)));
    processObserver->setComputebFunction("SingleDSObserverLCSPlugin", "computeU");

    // Relations
    unsigned int ninter = 2; // dimension of your Interaction = size of y and lambda vectors
    SimpleMatrix B(ndof, ninter);
    B(0, 0) = -1.0;
    B(1, 0) = 1.0;
    B(2, 1) = -1.0;
    B(3, 1) = 1.0;
    SimpleMatrix C(ninter, ndof);
    C(0, 0) = -1.0;
    C(0, 1) = 1.0;
    C(1, 2) = -1.0;
    C(1, 3) = 1.0;

    SP::FirstOrderLinearR myProcessRelation(new FirstOrderLinearR(createSPtrSimpleMatrix(C), createSPtrSimpleMatrix(B)));

    myProcessRelation->setComputeEFunction("SingleDSObserverLCSPlugin", "computeE");

    SimpleMatrix D(ninter, ninter);
    D(0, 0) = 1.0;
    D(1, 1) = 1.0;
    //myProcessRelation->setD(D);
    //return 0;

    // NonSmoothLaw
    unsigned int nslawSize = 2;
    SP::NonSmoothLaw myNslaw(new ComplementarityConditionNSL(nslawSize));

    SP::Interaction myProcessInteraction(new Interaction(myNslaw, myProcessRelation));

    // Model
    SP::Model ObserverLCS(new Model(t0, T));
    ObserverLCS->nonSmoothDynamicalSystem()->insertDynamicalSystem(processObserver);
    ObserverLCS->nonSmoothDynamicalSystem()->link(myProcessInteraction, processObserver);
    // TimeDiscretisation
    SP::TimeDiscretisation td(new TimeDiscretisation(t0, h));
    // == Creation of the Simulation ==
    SP::TimeStepping s(new TimeStepping(td));


    // OneStepIntegrator
    double theta = 0.5;
    // One Step Integrator
    SP::EulerMoreauOSI myIntegrator(new EulerMoreauOSI(theta));
    s->insertIntegrator(myIntegrator);
    // OneStepNSProblem



    // One Step non smooth problem

    SP::LCP osnspb(new LCP());
    s->insertNonSmoothProblem(osnspb);

    ObserverLCS->setSimulation(s);
    // ================================= Computation =================================

    // --- Initialisation of the simulation ---
    ObserverLCS->initialize();

    int k = 0; // Current step
    unsigned int N = ceil((T - t0) / h) + 1; // Number of time steps
    unsigned int outputSize = 10; // number of required data
    SimpleMatrix dataPlot(N, outputSize);
    SP::SiconosVector processLambda = myProcessInteraction->lambda(0);


    // We get values for the initial time step:
    // time
    dataPlot(k, 0) = s->nextTime();;
    dataPlot(k, 1) = (*processObserver->x())(0); // Observer x(1)
    dataPlot(k, 2) = (*processObserver->x())(1); //Observer x(2)
    dataPlot(k, 3) = (*processObserver->x())(2);// Process x(1)
    dataPlot(k, 4) = (*processObserver->x())(3);// Process x(2)
    dataPlot(k, 5) = (*processLambda)(0);
    dataPlot(k, 6) = (*processLambda)(1);
    dataPlot(k, 7) = (*processObserver->b())(0);
    dataPlot(k, 8) = abs((*processObserver->x())(0) - (*processObserver->x())(2)) ;
    dataPlot(k, 9) = abs((*processObserver->x())(1) - (*processObserver->x())(3))  ;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();
    // Simulation loop
    while (s->hasNextEvent())
    {
      k++;

      // get current time step

      s->computeOneStep();

      dataPlot(k, 0) = s->nextTime();;
      dataPlot(k, 1) = (*processObserver->x())(0);
      dataPlot(k, 2) = (*processObserver->x())(1);
      dataPlot(k, 3) = (*processObserver->x())(2);
      dataPlot(k, 4) = (*processObserver->x())(3);
      dataPlot(k, 5) = (*processLambda)(0);
      dataPlot(k, 6) = (*processLambda)(1);
      dataPlot(k, 7) = (*processObserver->b())(0);
      dataPlot(k, 8) = abs((*processObserver->x())(0) - (*processObserver->x())(2)) ;
      dataPlot(k, 9) = abs((*processObserver->x())(1) - (*processObserver->x())(3))  ;
      ++show_progress;

      s->nextStep();
    }

    // Write the results into the file "ObserverLCS.dat"
    ioMatrix::write("SingleDSObserverLCS.dat", "ascii", dataPlot, "noDim");

    // --- Time loop ---
    // Comparison with a reference file
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("SingleDSObserverLCS.ref", "ascii", dataPlotRef);

    std::cout << "Error =" << (dataPlot - dataPlotRef).normInf() << std::endl;
    if ((dataPlot - dataPlotRef).normInf() > 1e-09)
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
