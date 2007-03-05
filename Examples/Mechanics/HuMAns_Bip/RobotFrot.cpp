
// =============================== Robot sample (Bip) ===============================
//
// Keywords: LagrangianDS, LagrangianLinear relation, Moreau TimeStepping, LCP.
//
// =============================================================================================

#include "SiconosKernel.h"

using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 21;           // degrees of freedom for robot
    double t0 = 0;                   // initial computation time
    double T = 1.3;//0.005;                   // final computation time
    double h = 0.001;                // time step
    double criterion = 0.5;
    unsigned int maxIter = 5000;
    double en = 0.0, et = 0.0, mu = 0.1;      // nslaw

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
    bip->setComputeMassFunction("RobotFrotPlugin.so", "mass");
    bip->setComputeNNLFunction("RobotFrotPlugin.so", "NNL");
    bip->setComputeJacobianNNLFunction(0, "RobotFrotPlugin.so", "jacobianQNNL");
    bip->setComputeJacobianNNLFunction(1, "RobotFrotPlugin.so", "jacobianVNNL");

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

    NonSmoothLaw * nslaw = new NewtonImpactFrictionNSL(en, et, mu, 3);
    string G = "RobotFrotPlugin:G0";
    Relation * relation = new LagrangianScleronomousR("RobotFrotPlugin:h0", G);
    Interaction * inter = new Interaction("floor-feet", allDS, 0, 69, nslaw, relation);

    //The linear contraint corresponding to joints limits (hq+b>0)
    NonSmoothLaw * nslaw2 = new NewtonImpactFrictionNSL(en, et, 0, 3);
    SimpleMatrix H(90, 21);
    SimpleVector b(90);
    H.zero();
    b.zero();
    for (unsigned int i = 0; i < 15; ++i)
    {
      for (unsigned int j = 0; j < 3; ++j)
      {
        H(6 * i + j, i) = -1;
        H(6 * i + j + 3, i) = 1;
      }
    }

    unsigned int i = 0;
    b(3 * i + 0) = 0.21;
    b(3 * i + 1) = 0.21;
    i++;
    b(3 * i + 2) = 0.64;
    b(3 * i + 3) = 0.64;
    i++;
    b(3 * i + 4) = 1.41;
    b(3 * i + 5) = -0.11;
    i++;
    b(3 * i + 6) = 0.20;
    b(3 * i + 7) = 1.17;
    i++;
    b(3 * i + 8) = 0.21;
    b(3 * i + 9) = 0.21;
    i++;
    b(3 * i + 10) = 0.64;
    b(3 * i + 11) = 0.64;
    i++;
    b(3 * i + 12) = 1.41;
    b(3 * i + 13) = -0.11;
    i++;
    b(3 * i + 14) = 0.2;
    b(3 * i + 15) = 1.17;
    i++;
    b(3 * i + 16) = 0.17;
    b(3 * i + 17) = 0.17;
    i++;
    b(3 * i + 18) = 0.14;
    b(3 * i + 19) = 0.14;
    i++;
    b(3 * i + 20) = 0.17;
    b(3 * i + 21) = 0.17;
    i++;
    b(3 * i + 22) = 0.14;
    b(3 * i + 23) = 0.14;
    i++;
    b(3 * i + 24) = 0.17;
    b(3 * i + 25) = 0.17;
    i++;
    b(3 * i + 26) = 0.08;
    b(3 * i + 27) = 0.08;
    i++;
    b(3 * i + 28) = 0.08;
    b(3 * i + 29) = 0.21;


    Relation * relation2 = new LagrangianLinearR(H, b);
    Interaction * inter2 =  new Interaction("joint-limits", allDS, 1, 90, nslaw2, relation2);


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
    OneStepNSProblem * osnspb = new FrictionContact3D(s, "name", "NLGS", 2000, 0.0005);

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
    dataPlot(k, 0) = Robot->getCurrentT();

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

    //   em->display();

    // --- Time loop ---
    cout << "Start computation ... " << endl;
    while (eventsManager->getNextEventPtr() != NULL)
    {
      // get current time step
      k++;
      cout << k << endl;
      s->newtonSolve(criterion, maxIter);
      dataPlot(k, 0) = Robot->getCurrentT();

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
