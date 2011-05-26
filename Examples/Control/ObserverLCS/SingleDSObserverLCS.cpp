
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
    unsigned int ncontrol = 1;



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
    cout << "A-LG" << endl;
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

    SP::SimpleVector x0(new SimpleVector(ndof));
    (*x0)(0) = Vinit;
    SP::FirstOrderLinearDS processObserver(new FirstOrderLinearDS(x0, createSPtrSimpleMatrix(TildeA)));
    processObserver->setComputebFunction("SingleDSObserverLCSPlugin.so", "computeU");
    DynamicalSystemsSet allDS;
    allDS.insert(processObserver);


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

    myProcessRelation->setComputeEFunction("SingleDSObserverLCSPlugin.so", "computeE");

    SimpleMatrix D(ninter, ninter);
    D(0, 0) = 1.0;
    D(1, 1) = 1.0;
    //myProcessRelation->setD(D);
    //return 0;

    // NonSmoothLaw
    unsigned int nslawSize = 2;
    SP::NonSmoothLaw myNslaw(new ComplementarityConditionNSL(nslawSize));

    // Choose a name and a number for your Interaction
    string nameInter = "Interaction";
    unsigned int numInter = 1;


    DynamicalSystemsSet dsProcessConcerned;
    dsProcessConcerned.insert(allDS.getPtr(0));

    SP::Interaction myProcessInteraction(new Interaction(nameInter, dsProcessConcerned, numInter, ninter, myNslaw, myProcessRelation));

    InteractionsSet allInteractions;
    allInteractions.insert(myProcessInteraction);
    cout << "tutu " << endl;

    // NonSmoothDynamicalSystem
    SP::NonSmoothDynamicalSystem myNSDS(new NonSmoothDynamicalSystem(allDS, allInteractions));


    // Model
    SP::Model ObserverLCS(new Model(t0, T));
    ObserverLCS->setNonSmoothDynamicalSystemPtr(myNSDS);
    // TimeDiscretisation
    SP::TimeDiscretisation td(new TimeDiscretisation(t0, h));
    // == Creation of the Simulation ==
    SP::TimeStepping s(new TimeStepping(td));


    // OneStepIntegrator
    double theta = 0.5;
    // One Step Integrator
    SP::Moreau myIntegrator(new Moreau(allDS, theta));
    s->insertIntegrator(myIntegrator);
    // OneStepNSProblem



    // One Step non smooth problem

    SP::LCP osnspb(new LCP());
    s->insertNonSmoothProblem(osnspb);
    // ================================= Computation =================================

    // --- Initialisation of the simulation ---
    ObserverLCS->initialize(s);

    int k = 0; // Current step
    int N = (int)((T - t0) / h); // Number of time steps
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
    while (k < N - 1)
    {
      k++;

      s->nextStep();
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

    }

    // Write the results into the file "ObserverLCS.dat"
    ioMatrix io("SingleDSObserverLCS.dat", "ascii");
    io.write(dataPlot, "noDim");

    // --- Time loop ---

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
