
// =============================== Robot sample (Bip) ===============================
//
// Keywords: LagrangianDS, LagrangianLinear relation, Moreau TimeStepping, LCP.
//
// =============================================================================================

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 21;           // degrees of freedom for robot
    double t0 = 0;                   // initial computation time
    double T = 1.0;//0.005;                   // final computation time
    double h = 0.001;                // time step
    double criterion = 0.5;
    unsigned int maxIter = 5000;
    double e = 0.0;                  // nslaw

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // unsigned int i;
    DynamicalSystemsSet allDS; // the list of DS

    // --- DS: robot bip ---

    // The dof are angles between differents parts of the robot.

    // Initial position (angles in radian)
    SimpleVector q0(nDof), v0(nDof);
    q0(1) = -0.1;
    q0(2) = 0.2;
    q0(3) = -0.1;
    q0(5) = -0.1;
    q0(6) = 0.2;
    q0(7) = -0.1;

    LagrangianDS * bip = new LagrangianDS(1, q0, v0);

    // external plug-in
    bip->setComputeMassFunction("RobotPlugin.so", "mass");
    bip->setComputeNNLFunction("RobotPlugin.so", "NNL");
    bip->setComputeJacobianNNLFunction(0, "RobotPlugin.so", "jacobianQNNL");
    bip->setComputeJacobianNNLFunction(1, "RobotPlugin.so", "jacobianVNNL");

    allDS.insert(bip);

    // -------------------
    // --- Interactions---
    // -------------------

    // Two interactions:
    //  - one with Lagrangian non linear relation to define contact with ground
    //  - the other to define angles limitations (articular stops), with lagrangian linear relation
    //  Both with newton impact nslaw.

    InteractionsSet allInteractions;

    // -- relations --

    NonSmoothLaw * nslaw = new NewtonImpactNSL(e);
    string G = "RobotPlugin:G0";
    Relation * relation = new LagrangianScleronomousR("RobotPlugin:h0", G);
    Interaction * inter = new Interaction("floor-feet", allDS, 0, 23, nslaw, relation);

    //The linear contraint corresponding to joints limits (hq+b>0)
    SimpleMatrix H(30, 21);
    SimpleVector b(30);
    H.zero();
    H(0, 0) = -1;
    H(1, 0) = 1;
    b(0) = 0.21;
    b(1) = 0.21;
    H(2, 1) = -1;
    H(3, 1) = 1;
    b(2) = 0.64;
    b(3) = 0.64;
    H(4, 2) = -1;
    H(5, 2) = 1;
    b(4) = 1.41;
    b(5) = -0.11;
    H(6, 3) = -1;
    H(7, 3) = 1;
    b(6) = 0.2;
    b(7) = 1.17;
    H(8, 4) = -1;
    H(9, 4) = 1;
    b(8) = 0.21;
    b(9) = 0.21;
    H(10, 5) = -1;
    H(11, 5) = 1;
    b(10) = 0.64;
    b(11) = 0.64;
    H(12, 6) = -1;
    H(13, 6) = 1;
    b(12) = 1.41;
    b(13) = -0.11;
    H(14, 7) = -1;
    H(15, 7) = 1;
    b(14) = 0.2;
    b(15) = 1.17;
    H(16, 8) = -1;
    H(17, 8) = 1;
    b(16) = 0.17;
    b(17) = 0.17;
    H(18, 9) = -1;
    H(19, 9) = 1;
    b(18) = 0.14;
    b(19) = 0.14;
    H(20, 10) = -1;
    H(21, 10) = 1;
    b(20) = 0.17;
    b(21) = 0.17;
    H(22, 11) = -1;
    H(23, 11) = 1;
    b(22) = 0.14;
    b(23) = 0.14;
    H(24, 12) = -1;
    H(25, 12) = 1;
    b(24) = 0.17;
    b(25) = 0.17;
    H(26, 13) = -1;
    H(27, 13) = 1;
    b(26) = 0.08;
    b(27) = 0.08;
    H(28, 14) = -1;
    H(29, 14) = 1;
    b(28) = 0.08;
    b(29) = 0.21;


    Relation * relation2 = new LagrangianLinearR(H, b);
    Interaction * inter2 =  new Interaction("joint-limits", allDS, 1, 30, nslaw, relation2);


    allInteractions.insert(inter);
    allInteractions.insert(inter2);

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------

    bool isBVP = 0;
    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions, isBVP);

    // -------------
    // --- Model ---
    // -------------

    Model * Robot = new Model(t0, T);
    Robot->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    TimeDiscretisation * t = new TimeDiscretisation(h, Robot);

    TimeStepping* s = new TimeStepping(t);


    // -- OneStepIntegrators --
    OneStepIntegrator * OSI =  new Moreau(bip, 0.500001, s);

    // -- OneStepNsProblem --
    //OneStepNSProblem * osnspb = new LCP(s,"name","Lemke",2001, 0.0005);
    LCP * osnspb = new LCP(s, "name", "PGS", 2001, 0.0005);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition ===========================


    // ================================= Computation =================================

    // --- Simulation initialization ---
    s->initialize();
    cout << "End of simulation initialisation" << endl;

    int k = 0; // Current step
    int N = t->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 22;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time
    dataPlot(k, 0) = Robot->getT0();

    for (unsigned int i = 1; i < 22; i++)
      dataPlot(k, i) = bip->getQ()(i - 1);

    // --- Compute elapsed time ---
    double t1, t2, elapsed;
    struct timeval tp;
    int rtn;
    clock_t start, end;
    double elapsed2;
    start = clock();
    rtn = gettimeofday(&tp, NULL);
    t1 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;

    EventsManager * eventsManager = s->getEventsManagerPtr();


    // --- Time loop ---
    cout << "Start computation ... " << endl;
    while (eventsManager->getNextEventPtr() != NULL)
    {
      // get current time step
      k++;
      cout << k << endl;
      s->newtonSolve(criterion, maxIter);
      s->nextStep();
      dataPlot(k, 0) = s->getStartingTime();

      for (unsigned int i = 1; i < 22; i++)
        dataPlot(k, i) = bip->getQ()(i - 1);
    }
    end = clock();
    rtn = gettimeofday(&tp, NULL);
    t2 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;
    elapsed = t2 - t1;
    elapsed2 = (end - start) / (double)CLOCKS_PER_SEC;
    cout << "time = " << elapsed << " --- cpu time " << elapsed2 << endl;
    cout << "End of computation - Number of iterations done: " << k << endl;

    // --- Output files ---
    ioMatrix out("result.dat", "ascii");
    out.write(dataPlot, "noDim");

    // --- Free memory ---
    delete osnspb;
    delete t;
    delete OSI;
    delete s;
    delete Robot;
    delete nsds;
    delete inter;
    delete inter2;
    delete relation;
    delete nslaw;
    delete relation2;
    delete bip;
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/MultiBeadsColumn\'" << endl;
  }
}
